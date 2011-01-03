#include <string.h>

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


  /* --- Set up MPI ----------------------             -------------- */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi.size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
  MPI_Get_processor_name(mpi.name, &mpi.namelen);

  mpi.comm = MPI_COMM_WORLD;
  mpi.info = MPI_INFO_NULL;


  // Temporary: 
  mpi.Ntasks = 1;
  mpi.task   = 0; 
  mpi.backgrrecno = 0; 

  /* --- Read input data and initialize --             -------------- */

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  readInput();
  spectrum.updateJ = TRUE;

  getCPU(1, TIME_START, NULL);
  init_ncdf_atmos(&atmos, &geometry, &infile);
  // init_io should go in here, and should be the following order:
  /* init_ncdf_atmos
   * init_ncdf_J
   * init_Background
   * init_ncdf_input
   * init_ncdf_aux 
   * init_ncdf_spec
  */

  /* Find out the work load for each process */
  distribute_jobs();

  // Temporary
  mpi.Ntasks = 1;

  printf("MPI nx = %d\n",mpi.nx);
  printf("MPI ny = %d\n",mpi.ny); 
  printf("Atmos moving = %d\n",atmos.moving);
  printf("Atmos ID = %s\n",atmos.ID);
  printf("Atmos ID len = %d\n",strlen(atmos.ID)); 


  /* Main loop over tasks */
  for (mpi.task = 0; mpi.task < mpi.Ntasks; mpi.task++) {
    /* Indices of x and y */
    mpi.ix = mpi.taskmap[mpi.task + mpi.my_start][0];
    mpi.iy = mpi.taskmap[mpi.task + mpi.my_start][1];

    /* Read atmosphere column */
    //readAtmos_ncdf(mpi.xnum[mpi.ix],mpi.ynum[mpi.iy], &atmos, &geometry, &infile);

    if (mpi.rank == 0)
      readAtmos_ncdf(214,220, &atmos, &geometry, &infile);
    if (mpi.rank == 1)
      readAtmos_ncdf(280,100, &atmos, &geometry, &infile);

    /*
      printf("    height      temp        ne           vz         vturb\n");
      //printf("    nh(0)       nh(1)       nh(2)        nh(3)      nh(4)        nh(5)       nhtot\n");
      for(i=0; i<175; i++){
      //printf(" %10.4e  %10.4e  %10.4e  %10.4e %10.4e  %10.4e  %10.4e\n",
      //   atmos.nH[0][i],atmos.nH[1][i],atmos.nH[2][i],atmos.nH[3][i],atmos.nH[4][i],
      //   atmos.nH[5][i],atmos.nHtot[i]);
      
      printf(" %10.4e  %10.4e  %10.4e  %10.4e %10.4e\n",
      geometry.height[i],atmos.T[i],atmos.ne[i],geometry.vel[i],atmos.vturb[i]);
      }
      exit(0);
    */

    if (atmos.Stokes) Bproject();


    readAtomicModels();
    readMolecularModels();
    SortLambda();
  
    // Maybe we should have an init_io where all these files are initialised?
    init_ncdf_J();
    init_Background();
    init_ncdf_spec();
    Background_p(analyze_output=TRUE, equilibria_only=FALSE);

    init_ncdf_indata();

    getProfiles();
    initSolution_p();
    initScatter();

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

    writeAtmos_p();
    /* These will be superseed by writeInputData_p() 
       writeInput();
       writeAtmos(&atmos);
       writeGeometry(&geometry); 
    */
    writeSpectrum_p();
    /* These two were superseeded by writeSpectrum_p()
       writeSpectrum(&spectrum);
       writeFlux(FLUX_DOT_OUT);
    */

    // Tiago: commenting out the writing of output files for now
    /*
      for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      atom = atmos.activeatoms[nact];

      writeAtom(atom);
      writePopulations(atom);
      writeRadRate(atom);
      writeCollisionRate(atom);
      writeDamping(atom);
      } 
      for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
      molecule = atmos.activemols[nact];
      writeMolPops(molecule);
      }

      writeOpacity();
      
      getCPU(1, TIME_POLL, "Write output");
    */

  } /* End of main task loop */

  close_ncdf_J();
  close_Background();
  close_ncdf_atmos(&atmos, &geometry, &infile);
  close_ncdf_spec();
  close_ncdf_indata();


  /* Frees from memory stuff used for job control */
  finish_jobs();
    
  printTotalCPU();
  MPI_Finalize();

  return 0;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
