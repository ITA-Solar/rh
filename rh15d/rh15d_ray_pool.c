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

#define WORKTAG  1
#define DIETAG   2
#define BATCHTAG 3

#ifndef REV_ID
#define REV_ID "UNKNOWN"
#endif

/* --- Function prototypes --                          -------------- */
void overlord_batch_nodeaware(long *node_cursor, const long *node_end,
                              long flush_interval, long *dispatched_out);
int  drone_batch(PoolOutputBuf *poolbuf);
static void compute_node_task_ranges(long *node_start, long *node_end);


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
  bool_t run_ray = FALSE, writej = FALSE;

  /* --- Set up MPI ----------------------             -------------- */
  /* Initialize main_logfile to stderr as safe default during MPI setup */
  mpi.main_logfile = stderr;
  initParallel(&argc, &argv, run_ray=FALSE);
  memset(&spectrum,0,sizeof(spectrum));

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
  readInput(NULL);
  if (input.p15d_rerun) readSavedKeywords();
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

  /* --- Populate node-shared atmosphere cache.  Must come AFTER
     distribute_jobs (sets node_ix0/ix1, mpi.ny) and BEFORE the first
     readAtmos so that even the dimension-discovery call is served
     from the cache.                                                 */
  init_atmos_node_cache(&atmos, &infile);

  /* Read first atmosphere column just to get dimensions.  Use a
     node-local column so every rank's first readAtmos hits the
     cache instead of falling back to file I/O.                      */
  readAtmos(mpi.xnum[mpi.node_ix0], mpi.ynum[0], &atmos, &geometry, &infile);
  if (atmos.Stokes) Bproject();
  readAtomicModels();
  readMolecularModels();
  SortLambda();
  checkValuesRayInput();
  initParallelIO(run_ray=FALSE, writej=FALSE);

  /* --- Node-aware batched pool ------------------------------------
     Every drone is pinned to its node's row range.  The overlord keeps
     a per-node cursor into the global taskmap; when a drone reports
     done, the overlord serves the next task from THAT drone's node's
     cursor.  Drones on drained nodes sit out the rest of the run at
     BATCHTAG but still participate in the collective flush with an
     empty poolbuf.  Because every assigned (ix, iy) lies in the
     drone's own node row range, readAtmos_hdf5 always hits the cache
     and never touches Lustre during the run loop.                   */

  int flush_interval = input.p15d_flush_interval;
  if (flush_interval < 1) flush_interval = 1;
  if ((long)flush_interval * mpi.size > mpi.total_tasks)
    flush_interval = (int)(mpi.total_tasks / mpi.size) + 1;
  if (flush_interval < 1) flush_interval = 1;

  /* Per-node cursor and end over the global taskmap.  Only rank 0
     uses these; drones don't need to know the task layout. */
  long *node_cursor = NULL;
  long *node_end    = NULL;
  if (mpi.rank == 0) {
    node_cursor = (long *) calloc((size_t) mpi.n_nodes, sizeof(long));
    node_end    = (long *) calloc((size_t) mpi.n_nodes, sizeof(long));
    compute_node_task_ranges(node_cursor, node_end);
  }

  PoolOutputBuf poolbuf;
  poolbuf_init(&poolbuf, flush_interval + 4);

  if (mpi.rank == 0) {
    sprintf(messageStr,
            "Pool mode (node-aware): %ld total tasks, %d nodes, "
            "flush every %d columns/rank\n",
            mpi.total_tasks, mpi.n_nodes, flush_interval);
    fprintf(mpi.main_logfile, "%s", messageStr);
    Error(MESSAGE, "main", messageStr);
    /* Log per-node task counts for diagnostics */
    for (int k = 0; k < mpi.n_nodes; k++) {
      sprintf(messageStr,
              "  node %2d: %ld tasks (taskmap [%ld, %ld))\n",
              k, node_end[k] - node_cursor[k], node_cursor[k], node_end[k]);
      fprintf(mpi.main_logfile, "%s", messageStr);
    }
  }

  int any_work_left = 1;
  while (any_work_left) {
    if (mpi.rank == 0) {
      long dispatched = 0;
      overlord_batch_nodeaware(node_cursor, node_end,
                               (long) flush_interval, &dispatched);
      /* Any node still has work? */
      any_work_left = 0;
      for (int k = 0; k < mpi.n_nodes; k++) {
        if (node_cursor[k] < node_end[k]) { any_work_left = 1; break; }
      }
    } else {
      drone_batch(&poolbuf);
    }

    /* All ranks collectively write buffered columns to disk and flush */
    writeCollective_pool(&poolbuf, TRUE);
    poolbuf_reset(&poolbuf);

    /* Broadcast termination decision */
    MPI_Bcast(&any_work_left, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  if (mpi.rank == 0) {
    free(node_cursor);
    free(node_end);
  }

  /* Tell all drones to exit */
  if (mpi.rank == 0) {
    long rank;
    for (rank = 1; rank <= mpi.size; ++rank)
      MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  } else {
    MPI_Status status;
    long dummy;
    MPI_Recv(&dummy, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }

  poolbuf_free(&poolbuf);
  closeParallelIO(run_ray=FALSE, writej=FALSE);
  /* Frees from memory stuff used for job control */
  finish_jobs();
  sprintf(messageStr, "*** Job ending. Total %ld 1-D columns: %ld converged, "
          " %ld did not converge, %ld crashed.\n%s", mpi.total_tasks, mpi.nconv,
          mpi.nnoconv, mpi.ncrash, "*** RH finished gracefully.\n");
  if (mpi.rank == 0) fprintf(mpi.main_logfile, "%s", messageStr);
  Error(MESSAGE,"main",messageStr);
  printTotalCPU();
  MPI_Info_free(&mpi.info);
  MPI_Finalize();
  return 0;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */


/* ------- start --------- compute_node_task_ranges --------------- */
static void compute_node_task_ranges(long *node_start, long *node_end) {
/* Compute per-node [start, end) ranges in the global taskmap.  Relies
   on get_taskmap's row-major iteration order: rows with the same ix
   are contiguous, so all rows owned by node k form a contiguous slice.
   Rank 0 calls this once at startup.

   node_ix1[k] is derived from the same deterministic base/remainder
   split used in distribute_jobs, so rank 0 does not need to
   communicate with other ranks to learn the boundaries.             */
  int  k;
  long t;
  int *node_ix1_arr = (int *) malloc((size_t) mpi.n_nodes * sizeof(int));
  int  base = mpi.nx / mpi.n_nodes;
  int  rem  = mpi.nx % mpi.n_nodes;
  for (k = 0; k < mpi.n_nodes; k++) {
    int ix0 = k * base + ((k < rem) ? k : rem);
    int ix1 = ix0 + base + ((k < rem) ? 1 : 0);
    node_ix1_arr[k] = ix1;
    node_start[k]   = 0;
    node_end[k]     = 0;
  }
  /* Single linear scan of the taskmap.  current_k tracks the node
     whose slice we're accumulating; we advance it whenever we see an
     ix that has crossed into the next node's row range. */
  if (mpi.taskmap == NULL || mpi.total_tasks == 0) {
    free(node_ix1_arr);
    return;
  }
  k = 0;
  node_start[0] = 0;
  for (t = 0; t < mpi.total_tasks; t++) {
    int ix = (int) mpi.taskmap[t][0];
    while (k < mpi.n_nodes && ix >= node_ix1_arr[k]) {
      node_end[k] = t;
      k++;
      if (k < mpi.n_nodes) node_start[k] = t;
    }
    if (k >= mpi.n_nodes) break;
  }
  if (k < mpi.n_nodes) node_end[k] = mpi.total_tasks;
  /* Any trailing nodes with no tasks get [end, end) — zero-length slice. */
  for (k = k + 1; k < mpi.n_nodes; k++) {
    node_start[k] = mpi.total_tasks;
    node_end[k]   = mpi.total_tasks;
  }
  free(node_ix1_arr);
}
/* ------- end   --------- compute_node_task_ranges --------------- */

/* ------- start ---------------------------- overlord.c ------------ */
void overlord_batch_nodeaware(long *node_cursor, const long *node_end,
                              long flush_interval, long *dispatched_out)
{
/* Node-aware dynamic dispatch for one batch.  Each drone is served at
   most `flush_interval` tasks from its own node's cursor.  Drones
   whose node is out of tasks (or whose node has nothing left for this
   batch) are sent BATCHTAG immediately and sit out until the next
   collective flush.

   State: node_cursor[k] advances monotonically through [start, end[k])
   across batches, so work already dispatched in previous batches is
   never re-sent.                                                     */
  MPI_Status status;
  int  result;
  long rank;
  long dispatched = 0;
  int  active = 0;

  /* drone_got[r] counts how many tasks rank r has received this batch.
     Used to enforce the per-drone flush_interval budget. */
  int *drone_got = (int *) calloc((size_t) mpi.size + 1, sizeof(int));

  /* --- Seed: send every drone its first task (if available) --- */
  for (rank = 1; rank <= mpi.size; rank++) {
    int nk = mpi.rank_node[rank];
    if (nk < mpi.n_nodes && node_cursor[nk] < node_end[nk]) {
      long t = node_cursor[nk]++;
      MPI_Send(&t, 1, MPI_LONG, rank, WORKTAG, MPI_COMM_WORLD);
      drone_got[rank] = 1;
      dispatched++;
      active++;
    } else {
      /* Drone's node has no work left — excuse it from this batch. */
      MPI_Send(0, 0, MPI_INT, rank, BATCHTAG, MPI_COMM_WORLD);
    }
  }

  /* --- Main dispatch loop: serve more work as drones report done --- */
  while (active > 0) {
    MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);
    int src = status.MPI_SOURCE;
    int nk  = mpi.rank_node[src];
    if (drone_got[src] < (int) flush_interval &&
        nk < mpi.n_nodes && node_cursor[nk] < node_end[nk]) {
      long t = node_cursor[nk]++;
      MPI_Send(&t, 1, MPI_LONG, src, WORKTAG, MPI_COMM_WORLD);
      drone_got[src]++;
      dispatched++;
    } else {
      /* Drone is done for this batch (either hit flush budget or its
         node is drained).  Send BATCHTAG so it returns for flush. */
      MPI_Send(0, 0, MPI_INT, src, BATCHTAG, MPI_COMM_WORLD);
      active--;
    }
  }

  free(drone_got);
  if (dispatched_out != NULL) *dispatched_out = dispatched;
}
/* ------- end   ---------------------------- overlord.c ------------ */

/* ------- start ---------------------------- drone.c --------------- */
int drone_batch(PoolOutputBuf *poolbuf) {
/* Processes columns until receiving a BATCHTAG (end of batch) or
   DIETAG (shutdown).  Returns the tag that ended the loop. */

    MPI_Status status;
    bool_t write_analyze_output, equilibria_only;
    int niter, result=1;
    static long task = 0;
    static bool_t first_call = TRUE;

    if (first_call) {
      mpi.isfirst = TRUE;
      first_call = FALSE;
    }

    /* Main loop over tasks in this batch */
    while (1) {
        if (mpi.stop) mpi.stop = FALSE;
        /* Receive a message from the overlord */
        MPI_Recv(&mpi.task, 1, MPI_LONG, 0, MPI_ANY_TAG,
                 MPI_COMM_WORLD, &status);
        /* Check the tag of the received message. */
        if (status.MPI_TAG == BATCHTAG) return BATCHTAG;
        if (status.MPI_TAG == DIETAG)   return DIETAG;

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
                mpi.rank, task, mpi.xnum[mpi.ix], mpi.ynum[mpi.iy]);
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
        /* In case of crash, buffer metadata and proceed to next task */
        if (mpi.stop) {
          sprintf(messageStr,
                  "Process %4d: *** SKIP  task %3ld (crashed after %d "
                  "iterations)\n", mpi.rank, task, mpi.niter[mpi.task]);
          fprintf(mpi.main_logfile, "%s", messageStr);
          Error(MESSAGE, "main", messageStr);
          close_Background();  /* To avoid many open files */
          mpi.ncrash++;
          mpi.stop = FALSE;
          mpi.dpopsmax[mpi.task] = 0.0;
          mpi.convergence[mpi.task] = -1;
          /* Buffer MPI metadata only (no ray/aux for crashed columns) */
          poolbuf_store_mpi(poolbuf, task);
          MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
          continue;
        }
        /* Printout some info, finished iter */
        if (mpi.convergence[mpi.task]) {
            sprintf(messageStr,
                    "Process %4d: *** END   task %3ld iter, iterations = %3d,"
                    " CONVERGED\n", mpi.rank, task, mpi.niter[mpi.task]);
            mpi.nconv++;
        } else {
          sprintf(messageStr,
                  "Process %4d: *** END   task %3ld iter, iterations = %3d,"
                  " NO convergence\n", mpi.rank, task, mpi.niter[mpi.task]);
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
            /* Buffer aux + atmos data BEFORE geometry is redefined */
            poolbuf_store_mpi(poolbuf, task);
            poolbuf_store_aux_atmos(poolbuf);
            /* Redefine geometry just for this ray */
            atmos.Nrays     = 1;
            geometry.Nrays  = 1;
            geometry.muz[0] = io.ray_muz;
            geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
            geometry.muy[0] = 0.0;
            geometry.wmu[0] = 1.0;
            spectrum.updateJ = FALSE;
            calculate_ray();
            /* Buffer ray data after calculate_ray() populates spectrum.I */
            poolbuf_store_ray(poolbuf);
            /* Put back previous values for geometry  */
            atmos.Nrays     = geometry.Nrays = geometry.save_Nrays;
            geometry.muz[0] = geometry.save_muz;
            geometry.mux[0] = geometry.save_mux;
            geometry.muy[0] = geometry.save_muy;
            geometry.wmu[0] = geometry.save_wmu;
            spectrum.updateJ = TRUE;
        } else {
            /* Non-converged: buffer MPI metadata only */
            poolbuf_store_mpi(poolbuf, task);
        }
        /* --- Send result to overlord --- */
        MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    } /* End of batch task loop */
}
/* ------- end   ---------------------------- drone.c ------------ */
