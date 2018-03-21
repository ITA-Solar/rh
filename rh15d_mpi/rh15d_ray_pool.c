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

#define WORKTAG 1
#define DIETAG 2

#ifndef REV_ID
#define REV_ID "UNKNOWN"
#endif

/* --- Function prototypes --                          -------------- */
void overlord(void);
void drone(void);


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

/* ------- begin -------------------------- rhf1d.c ----------------- */

int main(int argc, char *argv[])
{
  bool_t run_ray, writej;

  /* --- Set up MPI ----------------------             -------------- */
  initParallel(&argc, &argv, run_ray=FALSE);
  mpi.size -= 1;  /* Remove overlord from count, as it is not doing work */
  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();
  /* Direct log stream into MPI log files */
  mpi.main_logfile     = commandline.logfile;
  commandline.logfile  = mpi.logfile;
  if (mpi.size == 0) {
    sprintf(messageStr,
            "Must run rh15d_ray_pool with more than one process. Aborting.");
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }
  strcpy(mpi.rev_id, REV_ID); /* save revision */
  /* --- Read input data and initialize --             -------------- */
  readInput();
  spectrum.updateJ = TRUE;
  getCPU(1, TIME_START, NULL);
  init_atmos(&atmos, &geometry, &infile);
  /* Find out the work load for each process, put only one task for pool */
  distribute_jobs();
  /* Saved input overrides any current options */
  if (input.p15d_rerun) {
      readSavedInput();
  } else {
      readRayInput();
  }
  mpi.Ntasks = 1;
  atmos.moving = TRUE;  /* To prevent moving change from column [0, 0] */
  /* Read first atmosphere column just to get dimensions */
  readAtmos(0, 0, &atmos, &geometry, &infile);
  if (atmos.Stokes) Bproject();
  readAtomicModels();
  readMolecularModels();
  SortLambda();
  checkValuesRayInput();
  initParallelIO(run_ray=FALSE, writej=FALSE);
  /*//////////////////////
  ////////////////////////
  //////////////////////*/
  if (mpi.rank == 0) {

    overlord();

  } else {

    drone();

  }
  /*//////////////////////
  ////////////////////////
  //////////////////////*/
  /* writeOutput(writej=FALSE); */
  closeParallelIO(run_ray=FALSE, writej=FALSE);
  /* Frees from memory stuff used for job control */
  finish_jobs();
  sprintf(messageStr, "*** Job ending. Total %ld 1-D columns: %ld converged, "
          " %ld did not converge, %ld crashed.\n%s", mpi.total_tasks, mpi.nconv,
          mpi.nnoconv, mpi.ncrash, "*** RH finished gracefully.\n");
  if (mpi.rank == 0) fprintf(mpi.main_logfile, "%s", messageStr);
  Error(MESSAGE,"main",messageStr);
  printTotalCPU();
  MPI_Finalize();
  return 0;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */


/* ------- start ---------------------------- overlord.c ------------ */
void overlord(void) {
  MPI_Status status;
  int result;
  long rank, current_task=0;

  if (mpi.total_tasks == 1) {
      /* Send task to first process */
      MPI_Sendrecv(&current_task, 1, MPI_LONG, 1, WORKTAG,
		   &result, 1, MPI_INT, 1, MPI_ANY_TAG,
		   MPI_COMM_WORLD, &status);
      /* Tell all the drones to exit by sending an empty message with the DIETAG. */
     for (rank = 1; rank <= mpi.size; ++rank) {
       MPI_Ssend(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
     }
  } else {
    /* Seed the drones; send one unit of work to each drone. */
    for (rank = 1; rank <= mpi.size; ++rank) {
      if (rank > mpi.total_tasks)
	break;

      /* Send it to each rank */
      MPI_Ssend(&current_task, 1, MPI_LONG, rank, WORKTAG, MPI_COMM_WORLD);
      ++current_task;
    }

    /* Loop over getting new work requests until there is no more work to be done */
    while (current_task < mpi.total_tasks) {
      /* Receive results from a drone */
      MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
	       MPI_COMM_WORLD, &status);
      /* Send the drone a new work unit */
      MPI_Ssend(&current_task, 1, MPI_LONG, status.MPI_SOURCE, WORKTAG,
	        MPI_COMM_WORLD);
      /* Get the next unit of work to be done */
      ++current_task;
    }

    /* There's no more work to be done, so receive all the outstanding results
       from the drones. */
    for (rank = 1; rank <= mpi.size; ++rank) {
      if (rank > mpi.total_tasks)
	break;
      MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE,
	       MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    /* Tell all the drones to exit by sending an empty message with the DIETAG. */
    for (rank = 1; rank <= mpi.size; ++rank) {
      MPI_Ssend(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    }
  }
}
/* ------- end   ---------------------------- overlord.c ------------ */

/* ------- start ---------------------------- drone.c --------------- */
void drone(void) {
    MPI_Status status;
    bool_t write_analyze_output, equilibria_only;
    int niter, result=1;
    long task=1;

    mpi.isfirst = TRUE;
    /* Main loop over tasks */
    while (1) {
        if (mpi.stop) mpi.stop = FALSE;
        /* Receive a message from the overlord */
        MPI_Recv(&mpi.task, 1, MPI_LONG, 0, MPI_ANY_TAG,
                 MPI_COMM_WORLD, &status);
        /* Check the tag of the received message. */
        if (status.MPI_TAG == DIETAG) return;

        ++task;

        /* Do the work */
        /* Indices of x and y */
        mpi.ix = mpi.taskmap[mpi.task][0];
        mpi.iy = mpi.taskmap[mpi.task][1];
        /* To use only first element of Ntasks arrays, set mpi.task to zero */
        mpi.task = 0;
        /* Printout some info */
        sprintf(messageStr,
                "Process %4d: --- START task %3ld, (xi,yi) = (%3d,%3d)\n",
                mpi.rank, task-1, mpi.xnum[mpi.ix], mpi.ynum[mpi.iy]);
        fprintf(mpi.main_logfile, "%s", messageStr);
        Error(MESSAGE, "main", messageStr);

        /* Read atmosphere column */
        readAtmos(mpi.xnum[mpi.ix],mpi.ynum[mpi.iy], &atmos, &geometry,
                  &infile);
        /* Update quantities that depend on atmosphere and initialise others */
        UpdateAtmosDep();
        /* --- Calculate background opacities --             ------------- */
        Background_p(write_analyze_output=TRUE, equilibria_only=FALSE);
        getProfiles();
        initSolution_p();
        initScatter();
        mpi.isfirst = FALSE; /* Put down here because initSolution_p uses it */
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
                  "Process %4d: *** SKIP  task %3ld (crashed after %d "
                  "iterations)\n", mpi.rank, task-1, mpi.niter[mpi.task]);
          fprintf(mpi.main_logfile, "%s", messageStr);
          Error(MESSAGE, "main", messageStr);
          close_Background();  /* To avoid many open files */
          mpi.ncrash++;
          mpi.stop = FALSE;
          mpi.dpopsmax[mpi.task] = 0.0;
          mpi.convergence[mpi.task] = -1;
          /* Write MPI output and send result to overlord*/
          writeMPI_p(task);
          MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
          continue;
        }
        /* Printout some info, finished iter */
        if (mpi.convergence[mpi.task]) {
            sprintf(messageStr,
                    "Process %4d: *** END   task %3ld iter, iterations = %3d,"
                    " CONVERGED\n", mpi.rank, task-1, mpi.niter[mpi.task]);
            mpi.nconv++;
        } else {
          sprintf(messageStr,
                  "Process %4d: *** END   task %3ld iter, iterations = %3d,"
                  " NO convergence\n", mpi.rank, task-1, mpi.niter[mpi.task]);
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
        if (mpi.convergence[mpi.task]) {
            /* Make sure aux written before ray redefined */
            writeAux_p();
            writeAtmos_p();
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
        /* --- Write output MPI group, send result to overlord ---------- */
        writeMPI_p(task);
        MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    } /* End of main task loop */
}
/* ------- end   ---------------------------- drone.c ------------ */
