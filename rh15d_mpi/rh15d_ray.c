/* $Id$ */
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"
#include "io.h"

#define COMMENT_CHAR    "#"
#define RAY_INPUT_FILE  "ray.input"
#define SVN_ID "$Id$"

/* --- Function prototypes --                          -------------- */
void init_ncdf_ray(void);
void init_ncdf_ray_new(void);
void init_ncdf_ray_old(void);
void writeRay(void);
void close_ncdf_ray(void);
void calculate_ray(void);
void writeAtmos_p(void);

/* --- Global variables --                             -------------- */

enum Topology topology = ONE_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
NCDF_Atmos_file infile;
CommandLine commandline;
char messageStr[MAX_MESSAGE_LENGTH];
BackgroundData bgdat;
MPI_data  mpi;
IO_data   io;
IO_buffer iobuf;

/* ------- begin -------------------------- rhf1d.c ----------------- */

int main(int argc, char *argv[])
{
  bool_t write_analyze_output, equilibria_only, run_ray, writej, exit_on_EOF;
  int    niter, nact, i, Nspect, Ntest, k, Nread, Nrequired, ierror,
         checkPoint, save_Nrays, *wave_index = NULL, conv_iter;
  double muz, save_muz, save_mux, save_muy, save_wmu;
  Atom *atom;
  Molecule *molecule;
  AtomicLine *line;
  ActiveSet *as;
  FILE  *fp_ray;

  time_t     curtime;
  struct tm *loctime;
  char timestr[ARR_STRLEN];
  char  inputLine[MAX_LINE_SIZE];



  /* --- Set up MPI ----------------------             -------------- */
  initParallel(&argc, &argv, run_ray=FALSE);

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  /* Direct log stream into MPI log files */
  mpi.main_logfile     = commandline.logfile;
  commandline.logfile  = mpi.logfile;
  strcpy(mpi.svn_id, SVN_ID); /* save SVN version */

  /* --- Read input data and initialize --             -------------- */
  readInput();
  spectrum.updateJ = TRUE;

  getCPU(1, TIME_START, NULL);
  init_ncdf_atmos(&atmos, &geometry, &infile);

  /* Find out the work load for each process */
  distribute_jobs();
  
  
  /* --- Read ray.input --                            --------------- */
  /* --- Read direction cosine for ray --              -------------- */
  if ((fp_ray = fopen(RAY_INPUT_FILE, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", RAY_INPUT_FILE);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }
  
  getLine(fp_ray, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%lf", &muz);
  checkNread(Nread, Nrequired=1, argv[0], checkPoint=1);

  if (muz <= 0.0  ||  muz > 1.0) {
    sprintf(messageStr,
	    "Value of muz = %f does not lie in interval <0.0, 1.0]\n", muz);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }

  /* --- read how many points to write detailed S, chi, eta, etc ---- */
  Nread = fscanf(fp_ray, "%d", &Nspect);
  checkNread(Nread, 1, argv[0], checkPoint=2);
  io.ray_nwave_sel = Nspect;

   /* --- Read wavelength indices for which chi and S are to be
       written out for the specified direction --    -------------- */

  if (Nspect > 0) {
    io.ray_wave_idx = (int *) malloc(Nspect * sizeof(int));
    Nread = 0;
    while (fscanf(fp_ray, "%d", &io.ray_wave_idx[Nread]) != EOF) Nread++;
    checkNread(Nread, Nspect, argv[0], checkPoint=3);
    fclose(fp_ray);

    wave_index = io.ray_wave_idx;
  }
  
  /* --- Save geometry values to change back after --    ------------ */
  save_Nrays = atmos.Nrays;   save_wmu = geometry.wmu[0];
  save_muz = geometry.muz[0]; save_mux = geometry.mux[0]; save_muy = geometry.muy[0];


  /* Main loop over tasks */
  for (mpi.task = 0; mpi.task < mpi.Ntasks; mpi.task++) {
    
      mpi.isfirst = (mpi.task == 0); 
      
      if (mpi.stop) mpi.stop = FALSE;
      

      /* Indices of x and y */
      mpi.ix = mpi.taskmap[mpi.task + mpi.my_start][0];
      mpi.iy = mpi.taskmap[mpi.task + mpi.my_start][1];
    
      /* Printout some info */
      sprintf(messageStr,
        "Process %3d: --- START task %3ld [of %ld], (xi,yi) = (%3d,%3d)\n",
         mpi.rank, mpi.task+1, mpi.Ntasks, mpi.xnum[mpi.ix], mpi.ynum[mpi.iy]);
      fprintf(mpi.main_logfile, messageStr);
      Error(MESSAGE, "main", messageStr);
      
            
      /* Read atmosphere column */
      readAtmos_ncdf(mpi.xnum[mpi.ix],mpi.ynum[mpi.iy], &atmos, &geometry, &infile);
      
      if (atmos.Stokes) Bproject();
      
      
      /* --- Run only once --                                  --------- */
      if (mpi.isfirst) {
        readAtomicModels();   
        readMolecularModels();
	
        SortLambda();
	
	/* Check if wavelength indices make sense */
	for (i = 0;  i < Nspect;  i++) {
	  if (wave_index[i] < 0  ||  wave_index[i] >= spectrum.Nspect) {
  	      sprintf(messageStr, "Illegal value of wave_index[n]: %4d\n"
  	      "Value has to be between 0 and %4d\n", 
  	      wave_index[i], spectrum.Nspect);
	      Error(ERROR_LEVEL_2, argv[0], messageStr);
	      continue;
	  }
	}
	
        initParallelIO(run_ray=FALSE, writej=FALSE);
        init_ncdf_ray();
      
      } else {
        /* Update quantities that depend on atmosphere and initialise others */
        UpdateAtmosDep();
      }
      
      /* --- Calculate background opacities --             ------------- */
      Background_p(write_analyze_output=TRUE, equilibria_only=FALSE);
      
      getProfiles();
      initSolution_p();
      initScatter();
      
      getCPU(1, TIME_POLL, "Total Initialize");
      
      /* --- Solve radiative transfer for active ingredients -- --------- */
      Iterate_p(input.NmaxIter, input.iterLimit);
      
      
      /* Treat odd cases as a crash */
      if (isnan(mpi.dpopsmax[mpi.task]) || isinf(mpi.dpopsmax[mpi.task]) || 
	  (mpi.dpopsmax[mpi.task] < 0) || ((mpi.dpopsmax[mpi.task] == 0) && (input.NmaxIter > 0)))
        mpi.stop = TRUE;
      

    /* In case of crash, write dummy data and proceed to next task */
    if (mpi.stop) {
      sprintf(messageStr,
	      "Process %3d: *** SKIP  task %3ld (crashed after %d iterations)\n",
	      mpi.rank, mpi.task+1, mpi.niter[mpi.task]);
      fprintf(mpi.main_logfile, messageStr);
      Error(MESSAGE, "main", messageStr);

      mpi.ncrash++;
      mpi.stop = FALSE;
      mpi.dpopsmax[mpi.task] = 0.0;
      mpi.convergence[mpi.task] = -1;

      continue;
    }



    /* Printout some info, finished iter */
    if (mpi.convergence[mpi.task]) {
      sprintf(messageStr,
       "Process %3d: *** END   task %3ld iter, iterations = %3d, CONVERGED\n",
       mpi.rank, mpi.task+1, mpi.niter[mpi.task]);
      mpi.nconv++;
    } else {
      sprintf(messageStr,
       "Process %3d: *** END   task %3ld iter, iterations = %3d, NO convergence\n",
       mpi.rank, mpi.task+1, mpi.niter[mpi.task]);
      mpi.nnoconv++;
    }

    fprintf(mpi.main_logfile, messageStr);
    Error(MESSAGE, "main", messageStr);
      
    /* Lambda iterate mean radiation field */
    adjustStokesMode();
    niter = 0;
    while (niter < input.NmaxScatter) {
      if (solveSpectrum(FALSE, FALSE) <= input.iterLimit) break;
      niter++;
    }
      
    copyBufVars(writej=FALSE);
      
 
    
    if (mpi.convergence[mpi.task]) {
      /* Redefine geometry just for this ray */
      atmos.Nrays     = 1;
      geometry.Nrays  = 1;
      geometry.muz[0] = muz;
      geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
      geometry.muy[0] = 0.0;
      geometry.wmu[0] = 1.0;
      spectrum.updateJ = FALSE;
      
      calculate_ray();
      writeRay();

      
      /* Put back previous values for geometry  */
      atmos.Nrays     = geometry.Nrays = save_Nrays;
      geometry.muz[0] = save_muz;
      geometry.mux[0] = save_mux;
      geometry.muy[0] = save_muy;
      geometry.wmu[0] = save_wmu;
      spectrum.updateJ = TRUE;
      
    }

    /* --- Write output files --                         -------------- */
    getCPU(1, TIME_START, NULL);
    writeAtmos_p();


    if (input.p15d_wspec) 
      writeSpectrum_p(); /* replaces writeSpectrum, writeFlux */

    getCPU(1, TIME_POLL, "Write output");

  } /* End of main task loop */


  writeOutput(writej=FALSE);

  closeParallelIO(run_ray=FALSE, writej=FALSE);
  close_ncdf_ray();

  /* Frees from memory stuff used for job control */
  finish_jobs();

  sprintf(messageStr,
   "*** Job ending. Total %ld 1-D columns: %ld converged, %ld not converged, %ld crashed.\n%s",
	  mpi.Ntasks, mpi.nconv, mpi.nnoconv, mpi.ncrash,
	  "*** RH finished gracefully.\n");	  
  if (mpi.rank == 0) fprintf(mpi.main_logfile, messageStr);
  Error(MESSAGE,"main",messageStr);

  printTotalCPU();
  MPI_Finalize();

  return 0;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
