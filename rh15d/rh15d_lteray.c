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
/* Default fill value for HDF5 */
const float FILLVALUE = FILL;

/* ------- begin -------------------------- rhf1d.c ----------------- */

int main(int argc, char *argv[])
{
  bool_t write_analyze_output, equilibria_only, run_ray;
  int    nact;
  Atom *atom;

  /* --- Set up MPI ----------------------             -------------- */
  /* Initialize main_logfile to stderr as safe default during MPI setup */
  mpi.main_logfile = stderr;
  initParallel(&argc, &argv, run_ray=FALSE);

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  /* Direct log stream into MPI log files */
  mpi.main_logfile     = commandline.logfile;
  commandline.logfile  = mpi.logfile;
  strcpy(mpi.rev_id, REV_ID); /* save revision */

  /* --- Read input data and initialize --             -------------- */
  readInput(NULL);
  spectrum.updateJ = FALSE;
  getCPU(1, TIME_START, NULL);
  init_atmos(&atmos, &geometry, &infile);
  distribute_jobs();  /* Find out the work load for each process */
  readRayInput();

  if (input.StokesMode == FIELD_FREE ||
      input.StokesMode == POLARIZATION_FREE) {
    input.StokesMode = FULL_STOKES;
  }

  /* --- redefine geometry for just this one ray --    -------------- */
  atmos.Nrays = geometry.Nrays = 1;
  geometry.muz[0] = io.ray_muz;
  geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
  geometry.muy[0] = 0.0;
  geometry.wmu[0] = 1.0;

  atmos.moving = TRUE;  /* To prevent moving change from column [0, 0] */
   /* Read first atmosphere column just to get dimensions */
  readAtmos(0, 0, &atmos, &geometry, &infile);
  readAtomicModels();
  readMolecularModels();
  SortLambda();
  checkValuesRayInput();
  /* --- Force LTE populations for all active atoms -- ------------ */
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    atom->initial_solution = LTE_POPULATIONS;
  }

  /* --- START stuff from initParallelIO, just getting the needed parts --- */
  init_Background();
  mpi.StokesMode_save = input.StokesMode;
  /* Get file position of atom files (to re-read collisions) */
  io.atom_file_pos = (long *) malloc(atmos.Nactiveatom * sizeof(long));
  mpi.zcut_hist    = (int *)    calloc(mpi.Ntasks , sizeof(int));
  mpi.single_log       = FALSE;
  init_hdf5_ray();
  /* --- END of stuff from initParalleIO ---*/

  /* Main loop over tasks */
  for (mpi.task = 0; mpi.task < mpi.Ntasks; mpi.task++) {

    /* Indices of x and y */
    mpi.ix = mpi.taskmap[mpi.task + mpi.my_start][0];
    mpi.iy = mpi.taskmap[mpi.task + mpi.my_start][1];

    /* Printout some info */
    sprintf(messageStr,  "Process %4d: --- START task %3ld [of %ld], "
            "(xi,yi) = (%3d,%3d)\n", mpi.rank, mpi.task+1, mpi.Ntasks,
            mpi.xnum[mpi.ix], mpi.ynum[mpi.iy]);
    fprintf(mpi.main_logfile, "%s", messageStr);
    Error(MESSAGE, "main", messageStr);

    /* Read atmosphere column */
    readAtmos(mpi.xnum[mpi.ix],mpi.ynum[mpi.iy], &atmos, &geometry, &infile);

    /* Update quantities that depend on atmosphere and initialise others */
    UpdateAtmosDep();

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
      "Process %4d: *** END   task %3ld\n",
	    mpi.rank, mpi.task+1);
    fprintf(mpi.main_logfile, "%s", messageStr);
    Error(MESSAGE, "main", messageStr);

    mpi.nconv++;
  } /* End of main task loop */

  if (mpi.Ntasks == 0) {
    sprintf(messageStr, "Process %4d: *** NO WORK (more processes than"
            "tasks!)\n", mpi.rank);
    fprintf(mpi.main_logfile, "%s", messageStr);
    Error(MESSAGE, "main", messageStr);
  }
  /* --- Stuff that was on closeParallelIO --- */
  close_atmos(&atmos, &geometry, &infile);
  free(io.atom_file_pos);
  close_hdf5_ray();
  /* --- END of stuff from closeParallelIO ---*/

  /* Frees from memory stuff used for job control */
  finish_jobs();
  sprintf(messageStr, "*** Job ending. Total %ld 1-D columns: %ld computed, "
          "%ld skipped.\n%s", mpi.Ntasks, mpi.nconv, mpi.ncrash,
	      "*** RH lte ray finished gracefully.\n");
  if (mpi.rank == 0) fprintf(mpi.main_logfile, "%s", messageStr);
  Error(MESSAGE, "main", messageStr);

  printTotalCPU();
  MPI_Finalize();

  return 0;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
