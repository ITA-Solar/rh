/* ------- file: -------------------------- writespect_p.c --------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Jun 19 17:02:19 2006 --

       --------------------------                      ----------RH-- */

/* --- Writes spectroscopic data to output file.
       XDR (external data representation) version. --  -------------- */


#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "parallel.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern enum Topology topology;

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- writeSpectrum.c --------- */

void writeSpectrum(Spectrum *spectrum)
{
  const char routineName[] = "writeSpectrum";
  register int nspect;

  bool_t  result = TRUE;
  int     Nintensity;
  double *lambda_air, vacuum_to_air_limit = VACUUM_TO_AIR_LIMIT;
  FILE   *fp_spectrum;
  XDR     xdrs;

  if (!strcmp(input.spectrum_output, "none")) return;

  if ((fp_spectrum = fopen(input.spectrum_output, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s",
	    input.spectrum_output);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_spectrum, XDR_ENCODE);

  result &= xdr_int(&xdrs, &spectrum->Nspect);

  if (spectrum->vacuum_to_air) {
    lambda_air = (double *) malloc(spectrum->Nspect * sizeof(double));
    vacuum_to_air(spectrum->Nspect, spectrum->lambda, lambda_air);
    result &= xdr_vector(&xdrs, (char *) lambda_air, spectrum->Nspect,
			 sizeof(double), (xdrproc_t) xdr_double);
    free(lambda_air);
  } else
    result &= xdr_vector(&xdrs, (char *) spectrum->lambda, spectrum->Nspect,
			 sizeof(double), (xdrproc_t) xdr_double);

  switch (topology) {
  case ONE_D_PLANE:
    Nintensity = atmos.Nrays * spectrum->Nspect;
    break;
  case TWO_D_PLANE:
    Nintensity = atmos.N[0] * atmos.Nrays * spectrum->Nspect;
    break;
  case THREE_D_PLANE:
    Nintensity = atmos.N[0]*atmos.N[1] * atmos.Nrays * spectrum->Nspect;
    break;
  case SPHERICAL_SYMMETRIC:
    Nintensity = atmos.Nrays * spectrum->Nspect;
    break;
  default:
    Nintensity = 0;
    break;
  }
  result &= xdr_vector(&xdrs, (char *) spectrum->I[0], Nintensity, 
		       sizeof(double), (xdrproc_t) xdr_double);

  result &= xdr_bool(&xdrs, &spectrum->vacuum_to_air);
  result &= xdr_double(&xdrs, &vacuum_to_air_limit);

  if (atmos.Stokes || input.backgr_pol) {
    result &= xdr_vector(&xdrs, (char *) spectrum->Stokes_Q[0], Nintensity, 
			 sizeof(double), (xdrproc_t) xdr_double);
    result &= xdr_vector(&xdrs, (char *) spectrum->Stokes_U[0], Nintensity, 
			 sizeof(double), (xdrproc_t) xdr_double);
    result &= xdr_vector(&xdrs, (char *) spectrum->Stokes_V[0], Nintensity, 
			 sizeof(double), (xdrproc_t) xdr_double);
  }

  if (!result) {
    sprintf(messageStr, "Unable to write proper amount to output file %s",
	    input.spectrum_output);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }
  xdr_destroy(&xdrs);
  fclose(fp_spectrum);

  /* --- Write angle-averaged mean intensity to file -- ------------- */

  if (spectrum->updateJ && !input.limit_memory) {
    for (nspect = 0;  nspect < spectrum->Nspect;  nspect++)
      writeJlambda_ncdf(nspect, spectrum->J[nspect]);

    /* --- Write the anisotropy J^2_0 in the z-direction -- --------- */
    if (input.backgr_pol) {
      for (nspect = 0;  nspect < spectrum->Nspect;  nspect++)
	writeJ20_ncdf(nspect, spectrum->J20[nspect]);
    }
  }
}
/* ------- end ---------------------------- writeSpectrum.c --------- */
