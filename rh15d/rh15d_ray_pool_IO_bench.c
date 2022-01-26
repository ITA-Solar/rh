#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

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
void get_stats(double,const char*);


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

double tic,toc,dt_tictoc;
double global_drone_min_stats[5];
double global_drone_max_stats[5];


/* ------- begin -------------------------- rhf1d.c ----------------- */

int main(int argc, char *argv[]) {
  bool_t run_ray = FALSE, writej = FALSE;

  for (size_t i = 0; i < 5; ++i) {
    global_drone_max_stats[i] = 0.0;
    global_drone_min_stats[i] = 1e6;
  }

  /* --- Set up MPI ----------------------             -------------- */
  initParallel(&argc, &argv, run_ray=FALSE);

  tic = MPI_Wtime();

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

  toc = MPI_Wtime();

  get_stats(toc - tic, "Initialization");

  /* --- Read input data and initialize --             -------------- */
  
  tic = MPI_Wtime();

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

  toc = MPI_Wtime();
  
  get_stats(toc - tic, "Read RayInput");


  tic = MPI_Wtime();

  /* Read first atmosphere column just to get dimensions */
  readAtmos(0, 0, &atmos, &geometry, &infile);
  if (atmos.Stokes) Bproject();
  readAtomicModels();
  readMolecularModels();

  toc = MPI_Wtime();

  get_stats(toc - tic, "read[Atomic/Molecular]Models");

  SortLambda();
  checkValuesRayInput();

  tic = MPI_Wtime();
  initParallelIO(run_ray=FALSE, writej=FALSE);
  toc = MPI_Wtime();
  
  get_stats(toc - tic, "initParallelIO");

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

void get_stats(double dt_local, const char* label) {
  double dt_min = 0.0;
  double dt_max = 0.0;

  MPI_Reduce(&dt_local,&dt_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&dt_local,&dt_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  if (mpi.rank == 0) {
    printf("[bench] %64s min/max: %0.3f [s] / %0.3f [s]\n", label, dt_min, dt_max);
  }
}


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

    double dt_readAtoms_min = 1e6;
    double dt_readAtoms_max = 0.0;

    double dt_writeAux_min  = 1e6;
    double dt_writeAux_max  = 0.0;

    double dt_writeRay_min  = 1e6;
    double dt_writeRay_max  = 0.0;

    double dt_writeAtmos_min = 1e6;
    double dt_writeAtmos_max = 0.0;

    double dt_writeMPI_min = 1e6;
    double dt_writeMPI_max = 0.0;

    /* Main loop over tasks */
    while (1) {
        if (mpi.stop) mpi.stop = FALSE;
        /* Receive a message from the overlord */
        MPI_Recv(&mpi.task, 1, MPI_LONG, 0, MPI_ANY_TAG,
                 MPI_COMM_WORLD, &status);
        /* Check the tag of the received message. */
        if (status.MPI_TAG == DIETAG) {
          printf("[bench][%04d] %50s min/max: %0.3f [s] / %0.3f [s]\n", mpi.rank, "readAtoms", dt_readAtoms_min, dt_readAtoms_max);
          printf("[bench][%04d] %50s min/max: %0.3f [s] / %0.3f [s]\n", mpi.rank, "writeAux", dt_writeAux_min, dt_writeAux_max);
          printf("[bench][%04d] %50s min/max: %0.3f [s] / %0.3f [s]\n", mpi.rank, "writeRay", dt_writeRay_min, dt_writeRay_max);
          printf("[bench][%04d] %50s min/max: %0.3f [s] / %0.3f [s]\n", mpi.rank, "writeAtmos", dt_writeAtmos_min, dt_writeAtmos_max);
          printf("[bench][%04d] %50s min/max: %0.3f [s] / %0.3f [s]\n", mpi.rank, "writeMPI", dt_writeMPI_min, dt_writeMPI_max);
         return; 
        }

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
        tic = MPI_Wtime(); 
        readAtmos(mpi.xnum[mpi.ix],mpi.ynum[mpi.iy], &atmos, &geometry,
                  &infile);
        toc = MPI_Wtime();

        dt_tictoc = toc - tic;
        if (dt_tictoc <= dt_readAtoms_min) dt_readAtoms_min = dt_tictoc;
        if (dt_tictoc >= dt_readAtoms_max) dt_readAtoms_max = dt_tictoc;

         // Update quantities that depend on atmosphere and initialise others 
        UpdateAtmosDep();
        /* --- Calculate background opacities --             ------------- */
        Background_p(write_analyze_output=TRUE, equilibria_only=FALSE);
        getProfiles();
        initSolution_p();
        initScatter();
        mpi.isfirst = FALSE; /* Put down here because initSolution_p uses it */
        getCPU(1, TIME_POLL, "Total Initialize");
        /* Skip calulcation and just say it worked */
        mpi.convergence[mpi.task] = true;

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
        // adjustStokesMode();
        niter = 1;
        // Dp not calculate just wait random n seconds [ 0, 1 or 3]
        int wait = rand() % 3;
        int milli_seconds = 1000 * wait;
        clock_t start_time = clock();
        while (clock() < start_time + milli_seconds) 
          ; 
        // while (niter < input.NmaxScatter) {
            // if (solveSpectrum(FALSE, FALSE) <= input.iterLimit) break;
            // niter++;
        // }
        if (mpi.convergence[mpi.task]) {
            tic = MPI_Wtime();
            /* Make sure aux written before ray redefined */
            writeAux_p();
            toc = MPI_Wtime();

            dt_tictoc = toc - tic;
            if (dt_tictoc <= dt_writeAux_min) dt_writeAux_min = dt_tictoc;
            if (dt_tictoc >= dt_writeAux_max) dt_writeAux_max = dt_tictoc;

            tic = MPI_Wtime();
            writeAtmos_p();
            toc = MPI_Wtime();

            dt_tictoc = toc - tic;
            if (dt_tictoc <= dt_writeAtmos_min) dt_writeAtmos_min = dt_tictoc;
            if (dt_tictoc >= dt_writeAtmos_max) dt_writeAtmos_max = dt_tictoc;


            // get_stats(toc-tic,"Write Atmos");
            /* Redefine geometry just for this ray */
            atmos.Nrays     = 1;
            geometry.Nrays  = 1;
            geometry.muz[0] = io.ray_muz;
            geometry.mux[0] = sqrt(1.0 - SQ(geometry.muz[0]));
            geometry.muy[0] = 0.0;
            geometry.wmu[0] = 1.0;
            spectrum.updateJ = FALSE;
            calculate_ray();
            tic = MPI_Wtime();
            writeRay();
            toc = MPI_Wtime();

            dt_tictoc = toc - tic;
            if (dt_tictoc <= dt_writeRay_min) dt_writeRay_min = dt_tictoc;
            if (dt_tictoc >= dt_writeRay_max) dt_writeRay_max = dt_tictoc;
            // printf("[id=%03d][task=%04ld] writeRay elapsed time: %0.2f [s]\n",mpi.rank, mpi.task, toc - tic);

            /* Put back previous values for geometry  */
            atmos.Nrays     = geometry.Nrays = geometry.save_Nrays;
            geometry.muz[0] = geometry.save_muz;
            geometry.mux[0] = geometry.save_mux;
            geometry.muy[0] = geometry.save_muy;
            geometry.wmu[0] = geometry.save_wmu;
            spectrum.updateJ = TRUE;
        }
        /* --- Write output MPI group, send result to overlord ---------- */
        tic = MPI_Wtime();
        writeMPI_p(task);
        toc = MPI_Wtime();

        dt_tictoc = toc - tic;
        if (dt_tictoc <= dt_writeMPI_min) dt_writeMPI_min = dt_tictoc;
        if (dt_tictoc >= dt_writeMPI_max) dt_writeMPI_max = dt_tictoc;

        MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    } /* End of main task loop */

}
/* ------- end   ---------------------------- drone.c ------------ */
