#include <string.h>
#include <mpi.h>

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
MPI_Comm mpi_comm;
MPI_Info mpi_info;

/* ------- begin -------------------------- rhf1d.c ----------------- */

int main(int argc, char *argv[])
{
  bool_t analyze_output, equilibria_only;
  int    niter, nact, i;

  Atom *atom;
  Molecule *molecule;

  /* --- Set up MPI ----------------------             -------------- */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Get_processor_name(mpi_name, &mpi_namelen);

  mpi_comm = MPI_COMM_WORLD;
  mpi_info = MPI_INFO_NULL;

  /* --- Read input data and initialize --             -------------- */

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  readInput();
  spectrum.updateJ = TRUE;

  getCPU(1, TIME_START, NULL);
  init_ncdf(&atmos, &geometry, &infile);
  readAtmos_ncdf(214,220, &atmos, &geometry, &infile);

  printf("Atmos moving = %d\n",atmos.moving);
  printf("Atmos ID = %s\n",atmos.ID);

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

  
  //MULTIatmos(&atmos, &geometry);
  if (atmos.Stokes) Bproject();


  readAtomicModels();
  readMolecularModels();
  SortLambda();
  
  //getBoundary(&geometry);
  
  Background(analyze_output=TRUE, equilibria_only=FALSE);
  //convertScales(&atmos, &geometry);

  getProfiles();
  initSolution();
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

  writeInput();
  writeAtmos(&atmos);
  writeGeometry(&geometry);
  writeSpectrum(&spectrum);
  writeFlux(FLUX_DOT_OUT);

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


  close_ncdf(&atmos, &geometry, &infile);

  printTotalCPU();
  MPI_Finalize();

  return 0;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
