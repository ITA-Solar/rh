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

#ifndef REV_ID
#define REV_ID "UNKNOWN"
#endif

/* --- Global variables --                             -------------- */
enum Topology topology = ONE_D_PLANE;
Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
Input_Atmos_file infile;
CommandLine commandline;
char messageStr[MAX_MESSAGE_LENGTH];
BackgroundData bgdat;
MPI_data  mpi;
IO_data   io;
IO_buffer iobuf;
const float FILLVALUE = FILL;  /* Default fill value for HDF5 */

int main(int argc, char *argv[])
{
  bool_t write_analyze_output, equilibria_only, run_ray, writej;
  int    niter;

  /* --- Set up MPI ----------------------             -------------- */
  /* Initialize main_logfile to stderr as safe default during MPI setup */
  mpi.main_logfile = stderr;
  initParallel(&argc, &argv, run_ray=FALSE);
  memset(&spectrum,0,sizeof(spectrum));

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();
  /* Direct log stream into MPI log files */
  mpi.main_logfile     = commandline.logfile;
  commandline.logfile  = mpi.logfile;
  strcpy(mpi.rev_id, REV_ID);  /* save revision */
  /* --- Read input data and initialize --             -------------- */
  readInput(NULL);
  if (input.p15d_rerun) readSavedKeywords();
  spectrum.updateJ = TRUE;
  getCPU(1, TIME_START, NULL);
  init_atmos(&atmos, &geometry, &infile);
  distribute_jobs();  /* Find out the work load for each process */
  /* Saved input overrides any current options */
  if (input.p15d_rerun) {
      readSavedInput();
  } else {
      readRayInput();
  }
  atmos.moving = TRUE;  /* To prevent moving change from column [0, 0] */
   /* Read first atmosphere column just to get dimensions */
  readAtmos(0, 0, &atmos, &geometry, &infile);
  if (atmos.Stokes) Bproject();
  readAtomicModels();
  readMolecularModels();
  SortLambda();
  checkValuesRayInput();
  initParallelIO(run_ray=FALSE, writej=FALSE);

  /* Main loop over tasks */
  for (mpi.task = 0; mpi.task < mpi.Ntasks; mpi.task++) {
    mpi.isfirst = (mpi.task == 0);
    if (mpi.stop) mpi.stop = FALSE;
    /* Indices of x and y */
    mpi.ix = mpi.taskmap[mpi.task + mpi.my_start][0];
    mpi.iy = mpi.taskmap[mpi.task + mpi.my_start][1];
    /* Printout some info */
    sprintf(messageStr,
            "Process %4d: --- START task %3ld [of %ld], (xi,yi) = (%3d,%3d)\n",
            mpi.rank, mpi.task+1, mpi.Ntasks, mpi.xnum[mpi.ix], mpi.ynum[mpi.iy]);
    fprintf(mpi.main_logfile, "%s", messageStr);
    Error(MESSAGE, "main", messageStr);
    /* Read atmosphere column */
    readAtmos(mpi.xnum[mpi.ix],mpi.ynum[mpi.iy], &atmos, &geometry, &infile);
    /* Update quantities that depend on atmosphere and initialise others */
    UpdateAtmosDep();
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
      	(mpi.dpopsmax[mpi.task] < 0) || ((mpi.dpopsmax[mpi.task] == 0) &&
        (input.NmaxIter > 0))) mpi.stop = TRUE;
    /* In case of crash, write dummy data and proceed to next task */
    if (mpi.stop) {
      sprintf(messageStr,
	          "Process %4d: *** SKIP  task %3ld (crashed after %d iterations)\n",
	          mpi.rank, mpi.task+1, mpi.niter[mpi.task]);
      fprintf(mpi.main_logfile, "%s", messageStr);
      Error(MESSAGE, "main", messageStr);
      close_Background();  /* To avoid many open files */
      mpi.ncrash++;
      mpi.stop = FALSE;
      mpi.dpopsmax[mpi.task] = 0.0;
      mpi.convergence[mpi.task] = -1;
      continue;
    }
    /* Printout some info, finished iter */
    if (mpi.convergence[mpi.task]) {
      sprintf(messageStr,
              "Process %4d: *** END   task %3ld iter, iterations = %3d,"
              " CONVERGED\n",  mpi.rank, mpi.task+1, mpi.niter[mpi.task]);
      mpi.nconv++;
    } else {
      sprintf(messageStr,
              "Process %4d: *** END   task %3ld iter, iterations = %3d,"
              " NO convergence\n", mpi.rank, mpi.task+1, mpi.niter[mpi.task]);
      mpi.nnoconv++;
    }
    fprintf(mpi.main_logfile, "%s", messageStr);
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
        geometry.muz[0] = io.ray_muz;
        geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
        geometry.muy[0] = 0.0;
        geometry.wmu[0] = 1.0;
        spectrum.updateJ = FALSE;
        calculate_ray();
        writeRay();
        /* Put back previous values for geometry  */
        atmos.Nrays     = geometry.Nrays = geometry.save_Nrays;
        geometry.muz[0] = geometry.save_muz;
        geometry.mux[0] = geometry.save_mux;
        geometry.muy[0] = geometry.save_muy;
        geometry.wmu[0] = geometry.save_wmu;
        spectrum.updateJ = TRUE;
    }

    /* --- Write output files --                         -------------- */
    getCPU(1, TIME_START, NULL);
    writeAtmos_p();
    getCPU(1, TIME_POLL, "Write output");
  } /* End of main task loop */
  writeOutput(writej=FALSE);
  closeParallelIO(run_ray=FALSE, writej=FALSE);
  /* Frees from memory stuff used for job control */
  finish_jobs();
  sprintf(messageStr,
          "*** Job ending. Total %ld 1-D columns: %ld converged, "
          "%ld did not converge, %ld crashed.\n%s", mpi.Ntasks, mpi.nconv,
          mpi.nnoconv, mpi.ncrash, "*** RH finished gracefully.\n");
  if (mpi.rank == 0) fprintf(mpi.main_logfile, "%s", messageStr);
  Error(MESSAGE, "main", messageStr);
  printTotalCPU();
  MPI_Finalize();
  return 0;
}
