/* ------- file: -------------------------- readb_xdr.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Jan 26 11:32:34 2000 --

       --------------------------                      ----------RH-- */

/* --- Reads the atmospheric magnetic field.
       XDR (external data representation) version. --  -------------- */


/* --- Read the absolute value of the magnetic field, and the angles
       gamma (the angle of the field with the normal direction to the
       atmosphere) and chi (the angle of the field with the positive
       x-axis, measured counter-clockwise).

 Note: Magnetic field strength should be given in Tesla (1 T = 10^4 Gauss) 
       Angles should be given in radians.

       --                                              -------------- */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- readB.c ----------------- */

bool_t readB(Atmosphere *atmos)
{
  const char routineName[] = "readB";

  bool_t  result = TRUE;
  FILE   *fp_stokes;
  XDR     xdrs;

  if (!strcmp(input.Stokes_input, "none")) return FALSE;

  if ((fp_stokes = fopen(input.Stokes_input, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", input.Stokes_input);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  xdrstdio_create(&xdrs, fp_stokes, XDR_DECODE);

  atmos->B = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->gamma_B = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->chi_B   = (double *) malloc(atmos->Nspace * sizeof(double));

  result &= xdr_vector(&xdrs, (char *) atmos->B, atmos->Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result &= xdr_vector(&xdrs, (char *) atmos->gamma_B, atmos->Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result &= xdr_vector(&xdrs, (char *) atmos->chi_B, atmos->Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);

  if (!result) {
    sprintf(messageStr, "Unable to read from input file %s",
	    input.Stokes_input);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp_stokes);

  return result;
}
/* ------- end ---------------------------- readB.c ----------------- */
