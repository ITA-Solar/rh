/* $Id: rh15d_ray.c 42 2011-12-02 02:13:29Z tiago $ */
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
#define SVN_ID "$Id: rh15d_ray.c 42 2011-12-02 02:13:29Z tiago $"
#define WORKTAG 1
#define DIETAG 2

/* --- Function prototypes --                          -------------- */
void init_ncdf_ray(void);
void init_ncdf_ray_new(void);
void init_ncdf_ray_old(void);
void writeRay(void);
void close_ncdf_ray(void);
void calculate_ray(void);
void overlord(void);
void drone(void);

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
int Nspect, *wave_index = NULL, save_Nrays;
double muz, save_muz, save_mux, save_muy, save_wmu;


/* ------- begin -------------------------- rhf1d.c ----------------- */

int main(int argc, char *argv[])
{
  bool_t analyze_output, equilibria_only, run_ray, writej, exit_on_EOF;
  int    nact, i, Ntest, k, Nread, Nrequired, ierror, checkPoint, conv_iter;
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
  mpi.size -= 1; /* Remove overlord from count, as it is not doing work */

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

  /* Find out the work load for each process, put only one task for pool */
  distribute_jobs();
  mpi.Ntasks = 1; 
  
         
  /* Read first atmosphere column just to get dimensions */
  readAtmos_ncdf(0, 0, &atmos, &geometry, &infile);
      
  if (atmos.Stokes) Bproject();
  
  readAtomicModels();   
  readMolecularModels();
	
  SortLambda();
	
  /* Check if wavelength indices make sense */
    for (i = 0;  i < Nspect;  i++) {
      if (wave_index[i] < 0  ||  wave_index[i] >= spectrum.Nspect) {
	  sprintf(messageStr, "Illegal value of wave_index[n]: %4d\n"
	  "Value has to be between 0 and %4d\n", 
	  wave_index[i], spectrum.Nspect);
	  Error(ERROR_LEVEL_2, "main", messageStr);
      }
    }
	
  initParallelIO(run_ray=FALSE, writej=FALSE);
  init_ncdf_ray();

  ////////////////////////
  ////////////////////////
  ////////////////////////
  if (mpi.rank == 0) {

    overlord();
    
  } else {

    drone();
    
  }
  ////////////////////////
  ////////////////////////
  ////////////////////////

  //writeOutput(writej=FALSE);

  closeParallelIO(run_ray=FALSE, writej=FALSE);
  close_ncdf_ray();

  /* Frees from memory stuff used for job control */
  finish_jobs();

  sprintf(messageStr,
   "*** Job ending. Total %ld 1-D columns: %ld converged, %ld not converged, %ld crashed.\n%s",
	  mpi.total_tasks, mpi.nconv, mpi.nnoconv, mpi.ncrash,
	  "*** RH finished gracefully.\n");	  
  if (mpi.rank == 0) fprintf(mpi.main_logfile, messageStr);
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
  long rank, proc_tasks, current_task=0;
  
  proc_tasks = MIN(mpi.total_tasks, mpi.size+1);

  
  /* Seed the drones; send one unit of work to each drone. */
  for (rank = 1; rank < proc_tasks; ++rank) {

    /* Send it to each rank */
    MPI_Send(&current_task,     /* message buffer */
             1,                 /* one data item */
             MPI_LONG,          /* data item is an integer */
             rank,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
    
    ++current_task;
  }

  /* Loop over getting new work requests until there is no more work to be done */
  while (current_task < mpi.total_tasks) {

    
    /* Receive results from a drone */
    MPI_Recv(&result,           /* message buffer */
             1,                 /* one data item */
             MPI_INT,           /* of type double real */
             MPI_ANY_SOURCE,    /* receive from any sender */
             MPI_ANY_TAG,       /* any type of message */
             MPI_COMM_WORLD,    /* default communicator */
             &status);          /* info about the received message */

    /* Send the drone a new work unit */
    MPI_Send(&current_task,     /* message buffer */
             1,                 /* one data item */
             MPI_LONG,          /* data item is an integer */
             status.MPI_SOURCE, /* to who we just received from */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */

    /* Get the next unit of work to be done */
    ++current_task;
  }
  
  /* There's no more work to be done, so receive all the outstanding results
     from the drones. */
  for (rank = 1; rank < proc_tasks; ++rank) {
    MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE,
             MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  }

  /* Tell all the drones to exit by sending an empty message with the DIETAG. */
  for (rank = 1; rank < proc_tasks; ++rank) {
    MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  }
}
/* ------- end   ---------------------------- overlord.c ------------ */

/* ------- start ---------------------------- drone.c --------------- */
void drone(void) {
  MPI_Status status;
  bool_t analyze_output, equilibria_only, run_ray, writej, exit_on_EOF;
  int i, niter, result=1;
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
        "Process %3d: --- START task %3ld, (xi,yi) = (%3d,%3d)\n",
         mpi.rank, task-1, mpi.xnum[mpi.ix], mpi.ynum[mpi.iy]);
      fprintf(mpi.main_logfile, messageStr);
      Error(MESSAGE, "main", messageStr);

            
      /* Read atmosphere column */
      readAtmos_ncdf(mpi.xnum[mpi.ix],mpi.ynum[mpi.iy], &atmos, &geometry, &infile);
      
      if (atmos.Stokes) Bproject();

      /* Update quantities that depend on atmosphere and initialise others */
      UpdateAtmosDep();
      
      /* --- Calculate background opacities --             ------------- */
      Background_p(analyze_output=TRUE, equilibria_only=FALSE);
      
      getProfiles();
      initSolution_p();
      initScatter();
      
      mpi.isfirst = FALSE; /* Put down here because initSolution_p uses it */
      
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
	      mpi.rank, task-1, mpi.niter[mpi.task]);
      fprintf(mpi.main_logfile, messageStr);
      Error(MESSAGE, "main", messageStr);
      
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
       "Process %3d: *** END   task %3ld iter, iterations = %3d, CONVERGED\n",
       mpi.rank, task-1, mpi.niter[mpi.task]);
      mpi.nconv++;
    } else {
      sprintf(messageStr,
       "Process %3d: *** END   task %3ld iter, iterations = %3d, NO convergence\n",
       mpi.rank, task-1, mpi.niter[mpi.task]);
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
      
    //copyBufVars(writej=FALSE);
    
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

    /* --- Write output files, send result to overlord ---------- */
    writeMPI_p(task);
    writeAtmos_p();
    MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

  } /* End of main task loop */

  
  
}
/* ------- end   ---------------------------- drone.c ------------ */

