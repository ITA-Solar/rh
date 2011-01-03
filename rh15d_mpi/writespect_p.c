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
#include "geometry.h"
#include "spectrum.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "parallel.h"
#include "io.h"

#define SPEC_FILE "output_spectrum.ncdf"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern enum Topology topology;

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];
extern MPI_data mpi;
extern IO_data io; 

/* ------- begin -------------------------- init_ncdf_spec.c -------- */
void init_ncdf_spec(void)
/* Creates the netCDF file for the spectrum */
{
  const char routineName[] = "init_ncdf_spec";
  int     ierror, ncid, nx_id, ny_id, nspect_id, nrays_id, wave_var,
          intensity_var, flux_var, stokes_u_var, stokes_q_var, 
          stokes_v_var, dimids[4];
  FILE   *test;

  /* Check if we can open the file */
  if ((test = fopen(SPEC_FILE, "w")) == NULL) {
    sprintf(messageStr, "Unable to open spectrum output file %s", SPEC_FILE);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(test);
  }

  /* Create the file  */
  if ((ierror = nc_create_par(SPEC_FILE, NC_NETCDF4 | NC_CLOBBER | NC_MPIIO, 
			      mpi.comm, mpi.info, &ncid))) ERR(ierror,routineName);
  
  /* Write atmos.ID as global attribute */
  if ((ierror = nc_put_att_text(ncid, NC_GLOBAL, "atmosID", strlen(atmos.ID),
				atmos.ID ))) ERR(ierror,routineName);

  // More stuff to add as global attributes:
  // * date of creation (master node only?)
  // * hostname (master node only?)
  // * version of the code

  /* Create dimensions */ 
  if ((ierror = nc_def_dim( ncid, "nx",    mpi.nx,          &nx_id     ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim( ncid, "ny",    mpi.ny,          &ny_id     ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim( ncid, "nrays", atmos.Nrays,     &nrays_id  ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim( ncid, "nwave", spectrum.Nspect, &nspect_id ))) 
    ERR(ierror,routineName);
  
  /* Create variables */
  dimids[0] = nx_id;
  dimids[1] = ny_id;
  dimids[2] = nrays_id;
  dimids[3] = nspect_id;

  /* Intensity */
  if ((ierror = nc_def_var( ncid,  "Intensity", NC_FLOAT, 4, dimids, &intensity_var)))
    ERR(ierror,routineName);

  /* Other Stokes parameters, if available */
  if (atmos.Stokes || input.backgr_pol) {
    if ((ierror = nc_def_var( ncid, "Stokes_Q", NC_FLOAT, 4, dimids, &stokes_q_var)))
      ERR(ierror,routineName); 
    if ((ierror = nc_def_var( ncid, "Stokes_U", NC_FLOAT, 4, dimids, &stokes_u_var)))
      ERR(ierror,routineName); 
    if ((ierror = nc_def_var( ncid, "Stokes_V", NC_FLOAT, 4, dimids, &stokes_v_var)))
      ERR(ierror,routineName); 
  }

  /* Flux in z-direction */
  dimids[0] = nx_id;
  dimids[1] = ny_id;
  dimids[2] = nspect_id;

  if ((ierror = nc_def_var( ncid,  "Flux", NC_FLOAT, 3, dimids, &flux_var)))
    ERR(ierror,routineName);

  // A variable for the ray intensity should be added in solveray... 

  // TODO: Put unit attributes in variables!!!

  /* Array with wavelengths */
  dimids[0] = nspect_id;

  if ((ierror = nc_def_var( ncid,  "Wavelength", NC_DOUBLE, 1, dimids, &wave_var)))
    ERR(ierror,routineName);

  /* End define mode */
  if ((ierror = nc_enddef(ncid))) ERR(ierror,routineName);

  
  /* Copy stuff to the IO data struct */
  io.spec_ncid         = ncid;
  io.spec_wave_var     = wave_var;
  io.spec_int_var      = intensity_var;
  io.spec_flux_var     = flux_var;
  io.spec_stokes_q_var = stokes_q_var;
  io.spec_stokes_u_var = stokes_u_var;
  io.spec_stokes_v_var = stokes_v_var;


  return; 

}

/* ------- end   -------------------------- init_ncdf_spec.c -------- */

/* ------- begin -------------------------- close_ncdf_spec.c -------- */
void close_ncdf_spec(void)
/* Closes the spec netCDF file */ 
{
  const char routineName[] = "close_ncdf_spec";
  int        ierror;

  if ((ierror = nc_close(io.spec_ncid))) ERR(ierror,routineName);
  return; 
}
/* ------- end   -------------------------- close_ncdf_spec.c -------- */


/* ------- begin -------------------------- writeSpectrum_p.c --------- */
void writeSpectrum_p(void)
/* Writes spectrum to netCDF file. NOTE: also writes flux, avoiding the need for flux.out */
{
  const char routineName[] = "writeSpectrum_p";
  register int mu, nspect;
  int     ierror;
  double *lambda_air, *flux, *wmuz, vacuum_to_air_limit = VACUUM_TO_AIR_LIMIT;
  size_t  start[] = {0, 0, 0, 0};
  size_t  count[] = {1, 1, 1, 1};

  start[0] = mpi.ix; 
  start[1] = mpi.iy;

  /* Write wavelength */
  if (spectrum.vacuum_to_air) {
    lambda_air = (double *) malloc(spectrum.Nspect * sizeof(double));
    vacuum_to_air(spectrum.Nspect, spectrum.lambda, lambda_air);
    if ((ierror = nc_put_var_double(io.spec_ncid, io.spec_wave_var, lambda_air )))
      ERR(ierror,routineName);
    free(lambda_air);
  } else
    if ((ierror = nc_put_var_double(io.spec_ncid, io.spec_wave_var, spectrum.lambda )))
      ERR(ierror,routineName);

  /* Write intensity */
  count[2] = atmos.Nrays;
  for(nspect = 0; nspect < spectrum.Nspect; nspect++){
    start[3] = nspect;
    if ((ierror = nc_put_vara_double(io.spec_ncid, io.spec_int_var, start, count,
				     spectrum.I[nspect] ))) ERR(ierror,routineName);
    /* If necessary, write rest of Stokes vector */
    if (atmos.Stokes || input.backgr_pol) {
      if ((ierror = nc_put_vara_double(io.spec_ncid, io.spec_stokes_q_var, start, count,
		 	       spectrum.Stokes_Q[nspect] ))) ERR(ierror,routineName);
      if ((ierror = nc_put_vara_double(io.spec_ncid, io.spec_stokes_u_var, start, count,
			       spectrum.Stokes_U[nspect] ))) ERR(ierror,routineName);  
      if ((ierror = nc_put_vara_double(io.spec_ncid, io.spec_stokes_v_var, start, count,
			       spectrum.Stokes_V[nspect] ))) ERR(ierror,routineName);  
    }
  }

  /* Calculate flux */
  wmuz = (double *) malloc(atmos.Nrays * sizeof(double));
  for (mu = 0;  mu < atmos.Nrays;  mu++)
    wmuz[mu] = geometry.muz[mu] * geometry.wmu[mu];

  flux = (double *) calloc(spectrum.Nspect, sizeof(double));
  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
    for (mu = 0;  mu < atmos.Nrays;  mu++)
      flux[nspect] += spectrum.I[nspect][mu] * wmuz[mu];
    flux[nspect] *= 2.0 * PI;
  }

  /* Write flux */
  count[2] = spectrum.Nspect;
  if ((ierror = nc_put_vara_double(io.spec_ncid, io.spec_flux_var, start, count,
				   flux ))) ERR(ierror,routineName);  
  free(wmuz); free(flux);


  /* --- Write angle-averaged mean intensity to file -- ------------- */

  if (spectrum.updateJ && !input.limit_memory) {
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
      writeJlambda_ncdf(nspect, spectrum.J[nspect]);

    /* --- Write the anisotropy J^2_0 in the z-direction -- --------- */
    if (input.backgr_pol) {
      for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
	writeJ20_ncdf(nspect, spectrum.J20[nspect]);
    }
  }


  return;
}
/* ------- end   -------------------------- writeSpectrum_p.c --------- */


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
