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

/* --- Function prototypes --                          -------------- */


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
MPI_data mpi;
IO_data io;

/* ------- begin -------------------------- rhf1d.c ----------------- */

int main(int argc, char *argv[])
{
  bool_t analyze_output, equilibria_only, run_ray;
  int    niter, nact, i, Ntest, k;

  Atom *atom;
  Molecule *molecule;
  AtomicLine *line;


  /* --- Set up MPI ----------------------             -------------- */
  initParallel(&argc, &argv, run_ray=FALSE);

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  /* Direct log stream into MPI log files */
  mpi.main_logfile     = commandline.logfile;
  commandline.logfile  = mpi.logfile;

  /* --- Read input data and initialize --             -------------- */
  readInput();
  spectrum.updateJ = TRUE;

  getCPU(1, TIME_START, NULL);
  init_ncdf_atmos(&atmos, &geometry, &infile);

  /* Find out the work load for each process */
  distribute_jobs();

  // Temporary
  //mpi.Ntasks = 2;

  /* Main loop over tasks */
  for (mpi.task = 0; mpi.task < mpi.Ntasks; mpi.task++) {

    /* Indices of x and y */
    mpi.ix = mpi.taskmap[mpi.task + mpi.my_start][0];
    mpi.iy = mpi.taskmap[mpi.task + mpi.my_start][1];
   
    /* Printout some info */
    sprintf(messageStr,
      "Process %d: --- START task %ld [of %ld], (xi,yi) = (%3d,%3d)\n",
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
      initParallelIO(run_ray=FALSE);

    } else {
      /* Update quantities that depend on atmosphere and initialise others */
      UpdateAtmosDep();
    }
    
    /* --- Calculate background opacities --             ------------- */
    Background_p(analyze_output=TRUE, equilibria_only=FALSE);

    getProfiles();
    initSolution_p();
    initScatter();

    getCPU(1, TIME_POLL, "Total Initialize");

    /* --- Solve radiative transfer for active ingredients -- --------- */
    Iterate_p(input.NmaxIter, input.iterLimit);

    /* Treat odd cases as a crash */
    if (isnan(mpi.dpopsmax) || isinf(mpi.dpopsmax) || (mpi.dpopsmax <= 0))
      mpi.stop = TRUE;

    /* In case of crash, write dummy data and proceed to next task */
    if (mpi.stop) {
      sprintf(messageStr,
	      "Process %d: *** SKIP  task %ld (crashed after %d iterations)\n",
	      mpi.rank, mpi.task+1, mpi.niter);
      fprintf(mpi.main_logfile, messageStr);
      Error(MESSAGE, "main", messageStr);

      mpi.ncrash++;
      mpi.stop = FALSE;
      mpi.dpopsmax = 0.0;
      mpi.convergence = -1;

      writeAtmos_p();    
      writeMPI_p();

      continue;
    }



    /* Printout some info, finished iter */
    if (mpi.convergence) {
      sprintf(messageStr,
       "Process %d: *** END   task %ld iter, iterations = %3d, CONVERGED\n",
       mpi.rank, mpi.task+1, mpi.niter);
      mpi.nconv++;
    } else {
      sprintf(messageStr,
       "Process %d: *** END   task %ld iter, iterations = %3d, NO convergence\n",
       mpi.rank, mpi.task+1, mpi.niter);
      mpi.nnoconv++;
    }

    fprintf(mpi.main_logfile, messageStr);
    Error(MESSAGE, "main", messageStr);


    adjustStokesMode();
    niter = 0;
    while (niter < input.NmaxScatter) {
      if (solveSpectrum(FALSE, FALSE) <= input.iterLimit) break;
      niter++;
    }

    /* --- Write output files --                         -------------- */
    getCPU(1, TIME_START, NULL);

    if (input.p15d_wspec) 
      writeSpectrum_p(); /* replaces writeSpectrum, writeFlux */
    writeAtmos_p();      /* replaces writeInput, writeAtmos, writeGeometry */
    writeJ_p();
    writeMPI_p();
    writeAux_p();        /* replaces writeAtom, writePopulations, writeRadrate, 
                            writeCollisionRate, writeDamping, and writeOpacity */

    getCPU(1, TIME_POLL, "Write output");

  } /* End of main task loop */



  closeParallelIO(run_ray = FALSE);

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
