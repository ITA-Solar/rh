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

/* ------- begin -------------------------- rhf1d.c ----------------- */

int main(int argc, char *argv[])
{
  bool_t analyze_output, equilibria_only;
  int    niter, nact, i;

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
  init_ncdf(&atmos, &geometry, &infile);
  readAtmos_ncdf(214,220, &atmos, &geometry, &infile);

  // Temporary, to be put in parallel init file
  if ((input.p15d_x1 <= 0) || (input.p15d_x1 > infile.nx))
    input.p15d_x1 = infile.nx;
  if ((input.p15d_y1 <= 0) || (input.p15d_y1 > infile.ny))
    input.p15d_y1 = infile.ny;

  mpi.nx = (input.p15d_x1 - input.p15d_x0) / input.p15d_xst + 
           (input.p15d_x1 - input.p15d_x0) % input.p15d_xst;

  mpi.ny = (input.p15d_y1 - input.p15d_y0) / input.p15d_yst + 
           (input.p15d_y1 - input.p15d_y0) % input.p15d_yst;


  printf("MPI nx = %d\n",mpi.nx);
  printf("MPI ny = %d\n",mpi.ny); 
  printf("Atmos moving = %d\n",atmos.moving);
  printf("Atmos ID = %s\n",atmos.ID);
  printf("Atmos ID len = %d\n",strlen(atmos.ID)); 

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
  Background_p(analyze_output=TRUE, equilibria_only=FALSE);

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


  close_ncdf_J();
  close_Background();
  close_atmos_ncdf(&atmos, &geometry, &infile);

  printTotalCPU();
  MPI_Finalize();

  return 0;
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
