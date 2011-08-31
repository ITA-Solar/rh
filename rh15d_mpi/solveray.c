/* ------- file: -------------------------- solveray.c --------------

       Version:       rh2.0, 1.5-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Fri Jan 14 16:00:00 2011 --

       --------------------------                      ----------RH-- */

/* --- Solves radiative transfer for given atmosphere and model atom
       along a ray with arbitrary \mu_z, assuming the atom's population
       numbers and angle-averaged radiation field is given.


       Expects input file ``ray.input'' containing two lines of the form

         muz
         Nspect  wave_index1  ....   wave_indexNspect
       --                                              -------------- */

#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "inputs.h"
#include "error.h"
#include "parallel.h"
#include "io.h"

#define COMMENT_CHAR    "#"
#define RAY_INPUT_FILE  "ray.input"


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
char messageStr[MAX_LINE_SIZE];
BackgroundData bgdat;
MPI_data mpi;
IO_data io;
IO_buffer iobuf;

/* ------- begin -------------------------- solveray.c -------------- */

int main(int argc, char *argv[])
{
  register int n, k;

  char    rayFileName[14], inputLine[MAX_LINE_SIZE];
  bool_t  result, exit_on_EOF, to_obs, initialize, crosscoupling,
          analyze_output, equilibria_only, run_ray, writej;
  int     Nspect, Nread, Nrequired, checkPoint, *wave_index = NULL;
  double  muz, *S, *chi, *J;
  FILE   *fp_out, *fp_ray, *fp_stokes;
  ActiveSet *as;

  /* --- Set up MPI ----------------------             -------------- */
  initParallel(&argc, &argv, run_ray=TRUE);

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  /* Direct log stream into MPI log files */
  mpi.main_logfile     = commandline.logfile;
  commandline.logfile  = mpi.logfile;

  /* --- Read input data and initialize --             -------------- */

  readInput();
  spectrum.updateJ = FALSE;

  /* --- Read input data for atmosphere --             -------------- */

  getCPU(1, TIME_START, NULL);
  init_ncdf_atmos(&atmos, &geometry, &infile);

  /* --- Find out the work load for each process --    -------------- */
  distribute_jobs();

  /* --- Find out which columns converged in the RH run ------------- */
  readConvergence();


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

  input.startJ = OLD_J;

  /* Main loop over tasks */
  for (mpi.task = 0; mpi.task < mpi.Ntasks; mpi.task++) {

    /* Indices of x and y */
    mpi.ix = mpi.taskmap[mpi.task + mpi.my_start][0];
    mpi.iy = mpi.taskmap[mpi.task + mpi.my_start][1];

    /* Printout some info */
    sprintf(messageStr,
      "Process %3d: --- START task %ld [of %ld], (xi,yi) = (%d,%d)\n",
	    mpi.rank, mpi.task+1, mpi.Ntasks, mpi.xnum[mpi.ix], mpi.ynum[mpi.iy]);
    fprintf(mpi.main_logfile, messageStr);
    Error(MESSAGE, "main", messageStr);


    /* Read atmosphere column */
    readAtmos_ncdf(mpi.xnum[mpi.ix],mpi.ynum[mpi.iy], &atmos, &geometry, &infile);

    
    if (atmos.Stokes) Bproject();

    if (mpi.task == 0) {
      readAtomicModels();
      readMolecularModels();

      SortLambda();

      bgdat.write_BRS = FALSE;

      /* writej=TRUE here is only for reading old J, not writing anything */
      initParallelIO(run_ray=TRUE, writej=TRUE);
      init_ncdf_ray();

    } else {
      /* Update quantities that depend on atmosphere and initialise others */
      UpdateAtmosDep();
    }


    /* --- Calculate background opacities --             ------------- */
    Background_p(analyze_output=FALSE, equilibria_only=FALSE); 

    getProfiles();
    initSolution_p();
    initScatter();

    /* If RH did not converge, do not calculate for this column */
    if (mpi.rh_converged[mpi.ix][mpi.iy] < 1) {
      sprintf(messageStr,
	      "Process %3d: *** SKIP  task %ld (RH did not converge)\n",
	      mpi.rank,mpi.task+1);
      fprintf(mpi.main_logfile, messageStr);
      Error(MESSAGE, "main", messageStr);	
      mpi.ncrash++;
      continue;
    }

    getCPU(1, TIME_POLL, "Total initialize");

    /* --- Solve radiative transfer equations --         -------------- */
    solveSpectrum(FALSE, FALSE);

    /* --- Write emergent spectrum to output file --     -------------- */
    writeRay();
   
    sprintf(messageStr,
      "Process %3d: *** END   task %ld\n",
	    mpi.rank, mpi.task+1, mpi.Ntasks);
    fprintf(mpi.main_logfile, messageStr);
    Error(MESSAGE, "main", messageStr);	


    mpi.nconv++;
  }  /* End of main task loop */

  closeParallelIO(run_ray=TRUE, writej=TRUE);
  close_ncdf_ray();

  finish_jobs();
  freeMatrix((void **) mpi.rh_converged);

  sprintf(messageStr,
	  "*** Job ending. Total %ld 1-D columns: %ld computed, %ld skipped.\n%s",
	  mpi.Ntasks, mpi.nconv, mpi.ncrash,
	  "*** Solveray finished gracefully.\n");	  
  if (mpi.rank == 0) fprintf(mpi.main_logfile, messageStr);
  Error(MESSAGE,"main",messageStr);


  printTotalCPU();
  MPI_Finalize();

  return 0;
}
/* ------- end ---------------------------- solveray.c -------------- */


