#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>

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

#ifndef REV_ID
#define REV_ID "UNKNOWN"
#endif

/* --- Function prototypes --                          -------------- */
void init_ncdf_ray(void);
void writeRay(void);
void close_ncdf_ray(void);

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
  bool_t write_analyze_output, equilibria_only, run_ray, exit_on_EOF;
  int    niter, nact, i, Ntest, k, Nspect, Nread, Nrequired,
         checkPoint, *wave_index = NULL;
  double muz;
  Atom *atom;
  Molecule *molecule;
  AtomicLine *line;
  FILE   *fp_out, *fp_ray, *fp_stokes;
  ActiveSet *as;

  time_t     curtime;
  struct tm *loctime;
  char timestr[ARR_STRLEN], rayFileName[14], inputLine[MAX_LINE_SIZE];


  /* --- Set up MPI ----------------------             -------------- */
  initParallel(&argc, &argv, run_ray=FALSE);

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  /* Direct log stream into MPI log files */
  mpi.main_logfile     = commandline.logfile;
  commandline.logfile  = mpi.logfile;
  strcpy(mpi.rev_id, REV_ID); /* save revision */

  /* --- Read input data and initialize --             -------------- */
  readInput();
  spectrum.updateJ = FALSE;

  getCPU(1, TIME_START, NULL);
  init_ncdf_atmos(&atmos, &geometry, &infile);

  /* Find out the work load for each process */
  distribute_jobs();



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

  if (input.StokesMode == FIELD_FREE ||
      input.StokesMode == POLARIZATION_FREE) {
    input.StokesMode = FULL_STOKES;
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

  /* --- redefine geometry for just this one ray --    -------------- */

  atmos.Nrays = geometry.Nrays = 1;
  geometry.muz[0] = muz;
  geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
  geometry.muy[0] = 0.0;
  geometry.wmu[0] = 1.0;


  /* --- Force LTE populations for all active atoms -- -------------- */
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    atom->initial_solution = LTE_POPULATIONS;
  }


  /* Main loop over tasks */
  for (mpi.task = 0; mpi.task < mpi.Ntasks; mpi.task++) {

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
    if (mpi.task == 0) {
      readAtomicModels();   
      readMolecularModels();

      SortLambda();
      
      bgdat.write_BRS = FALSE;
      
      /* --- START stuff from initParallelIO, just getting the needed parts --- */
      init_Background();
      
      mpi.StokesMode_save = input.StokesMode;
      
      /* Get file position of atom files (to re-read collisions) */
      io.atom_file_pos = (long *) malloc(atmos.Nactiveatom * sizeof(long));
      mpi.zcut_hist    = (int *)    calloc(mpi.Ntasks , sizeof(int));

      for (nact = 0; nact < atmos.Nactiveatom; nact++) {
	atom = atmos.activeatoms[nact];
	io.atom_file_pos[nact] = ftell(atom->fp_input);
      }

      mpi.single_log       = FALSE;  
    
      /* --- END of stuff from initParalleIO ---*/
      
      init_ncdf_ray();

    } else {
      /* Update quantities that depend on atmosphere and initialise others */
      UpdateAtmosDep();
    }

    /* --- Calculate background opacities --             ------------- */
    Background_p(write_analyze_output=FALSE, equilibria_only=FALSE);

    getProfiles();
    initSolution_p();
    initScatter();

    getCPU(1, TIME_POLL, "Total Initialize");

    
    /* --- Solve radiative transfer equations --         -------------- */
    solveSpectrum(FALSE, FALSE);

    /* --- Write emergent spectrum to output file --     -------------- */
    writeRay();

   
    sprintf(messageStr,
      "Process %3d: *** END   task %3ld\n",
	    mpi.rank, mpi.task+1);
    fprintf(mpi.main_logfile, messageStr);
    Error(MESSAGE, "main", messageStr);	

    mpi.nconv++;

  } /* End of main task loop */


  /* --- Stuff that was on closeParallelIO --- */
  close_ncdf_atmos(&atmos, &geometry, &infile);
  close_Background();
  free(io.atom_file_pos);
  /* --- END of stuff from closeParallelIO ---*/
  
  
  close_ncdf_ray();
  
  /* Frees from memory stuff used for job control */
  finish_jobs();

  sprintf(messageStr,
	  "*** Job ending. Total %ld 1-D columns: %ld computed, %ld skipped.\n%s",
	  mpi.Ntasks, mpi.nconv, mpi.ncrash,
	  "*** RH lte ray finished gracefully.\n");	  
  if (mpi.rank == 0) fprintf(mpi.main_logfile, messageStr);
  Error(MESSAGE,"main",messageStr);

  printTotalCPU();
  MPI_Finalize();

  return 0;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
