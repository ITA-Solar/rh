/* ------- file: -------------------------- backgrcontr.c -----------

       Version:       rh2.0, 2-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Apr 22 09:48:50 2009 --

       --------------------------                      ----------RH-- */

#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "background.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"

#define CONTR_INPUT_FILE "contribute.input"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

enum Topology topology = TWO_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
InputData input;
ProgramStats stats;
CommandLine commandline;
char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- backgrcontr.c ----------- */

int main(int argc, char *argv[])
{
  register int n;

  int     Nlambda;
  double *lambda;
  FILE   *fp;
  Atom *atom;
  Molecule *molecule;

  /* --- Read input data and initialize --             -------------- */

  setOptions(argc, argv);
  SetFPEtraps();

  readInput();
  readAtmos(&atmos, &geometry);
  readAtomicModels();
  readMolecularModels();

  if ((fp = fopen(CONTR_INPUT_FILE, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", CONTR_INPUT_FILE);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }

  fscanf(fp, "%d", &Nlambda);
  lambda = (double *) malloc(Nlambda * sizeof(double));
  for (n = 0;  n < Nlambda;  n++) fscanf(fp, "%lf", &lambda[n]);
  fclose(fp);

  backgrOpac(Nlambda, lambda);

  free(lambda);
}
/* ------- end ---------------------------- backgrcontr.c ----------- */
