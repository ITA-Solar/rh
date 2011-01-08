#include <string.h>
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
#include "xdr.h"
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
  bool_t analyze_output, equilibria_only;
  int    niter, nact, i, Ntest;

  Atom *atom;
  Molecule *molecule;
  AtomicLine *line;


  /* --- Set up MPI ----------------------             -------------- */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi.size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
  MPI_Get_processor_name(mpi.name, &mpi.namelen);

  mpi.comm = MPI_COMM_WORLD;
  mpi.info = MPI_INFO_NULL;

  /* --- Read input data and initialize --             -------------- */

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  readInput();
  spectrum.updateJ = TRUE;

  getCPU(1, TIME_START, NULL);
  init_ncdf_atmos(&atmos, &geometry, &infile);

  /* Find out the work load for each process */
  distribute_jobs();

  // Temporary
  mpi.Ntasks = 2;

  /* Main loop over tasks */
  for (mpi.task = 0; mpi.task < mpi.Ntasks; mpi.task++) {
    printf("#############################\n");
    printf("######    TASK  %d     #######\n",mpi.task);
    printf("#############################\n");

    /* Indices of x and y */
    mpi.ix = mpi.taskmap[mpi.task + mpi.my_start][0];
    mpi.iy = mpi.taskmap[mpi.task + mpi.my_start][1];

    /* Read atmosphere column */
    //readAtmos_ncdf(mpi.xnum[mpi.ix],mpi.ynum[mpi.iy], &atmos, &geometry, &infile);

    if (mpi.rank == 0)
      readAtmos_ncdf(214,220, &atmos, &geometry, &infile);
    if (mpi.rank == 1)
      readAtmos_ncdf(280,100, &atmos, &geometry, &infile);

    if (atmos.Stokes) Bproject();

  
    /* --- Run only once --                                  --------- */
    if (mpi.task == 0) {
      readAtomicModels();   
      readMolecularModels();

      SortLambda();
      initParallelIO();
    } else {
      /* Update quantities that depend on atmosphere and initialise others */
      UpdateAtmosDep();
    }
    

    Background_p(analyze_output=TRUE, equilibria_only=FALSE);

    getProfiles();
    initSolution_p();


   

    initScatter();

    /*
    printf("    n[0][i]        n[1][i]        n[2][i]         n[3][i]          \n");
    atom = atmos.activeatoms[0];
    for(i=0; i<175; i+=5){
      printf(" %10.4e  %10.4e  %10.4e  %10.4e\n",
	     atom->n[0][i],
	     atom->n[1][i],
	     atom->n[2][i],
	     atom->n[3][i]);
    }
    */

  
    // TIAGO: problem is NOT in collisional rates, but in stuff added further to Gamma
    //        (radiation?). Differences occur in solvespectrum. Must look in formal.c 
    //        and iterate.c. For non-PRD atom there doesn't seem to be an issue!!!
    //        Check source function, chi, eta, etc. in formal.c

  // end testing

    getCPU(1, TIME_POLL, "Total Initialize");

    /* --- Solve radiative transfer for active ingredients -- --------- */
    Iterate(input.NmaxIter, input.iterLimit);

 

    adjustStokesMode();
    niter = 0;
    while (niter < input.NmaxScatter) {
      if (solveSpectrum(FALSE, FALSE) <= input.iterLimit) break;
      niter++;
    }
    /* --- Write output files --                         -------------- */
    getCPU(1, TIME_START, NULL);

    writeSpectrum_p(); /* replaces writeSpectrum, writeFlux */
    writeAtmos_p();    /* replaces writeInput, writeAtmos, writeGeometry */
    writeMPI_p();
    writeAux_p();      /* replaces writeAtom, writePopulations, writeRadrate, 
                          writeCollisionRate, writeDamping, and writeOpacity */

    //freeAtom(&atmos.atoms[0]);
    //free(spectrum.lambda);
   
    // This will not work because active atoms cannot be freed with freeAtom,
    // as the line->lambda pointer has be repointed to spectrum.lambda. Must
    // create freeActiveaAtom. Use atom->active for this. 
 
    /*
    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      atom = atmos.activeatoms[nact];
      printf("### Freeing atom %s\n",atom->ID);
      freeAtom(atom);
    }
    */
    

    // TODO: must write stuff to put the molecules
   
    /*
    for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
      molecule = atmos.activemols[nact];
      writeMolPops(molecule);
    }
    */
   

    getCPU(1, TIME_POLL, "Write output");

    // TODO: here must free active atoms/molecules so that there is no problem
    //       in reallocating them in readAtomicModels

  } /* End of main task loop */


  close_ncdf_atmos(&atmos, &geometry, &infile);
  closeParallelIO();

  /* Frees from memory stuff used for job control */
  finish_jobs();
    
  printTotalCPU();
  MPI_Finalize();

  return 0;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
