/* ------- file: -------------------------- writespect_p.c --------

       Version:       rh2.0, 1.5-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Tue Jan 18 17:28:25 2011 --

       --------------------------                      ----------RH-- */

/* --- Writes spectroscopic data to output file. NetCDF version. ---- */


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
#include "parallel.h"
#include "io.h"


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
  if ((ierror = nc_def_var( ncid,  INT_NAME, NC_FLOAT, 4, dimids, &intensity_var)))
    ERR(ierror,routineName);

  /* Other Stokes parameters, if available */
  if (atmos.Stokes || input.backgr_pol) {
    if ((ierror = nc_def_var( ncid, STOKES_Q, NC_FLOAT, 4, dimids, &stokes_q_var)))
      ERR(ierror,routineName); 
    if ((ierror = nc_def_var( ncid, STOKES_U, NC_FLOAT, 4, dimids, &stokes_u_var)))
      ERR(ierror,routineName); 
    if ((ierror = nc_def_var( ncid, STOKES_V, NC_FLOAT, 4, dimids, &stokes_v_var)))
      ERR(ierror,routineName); 
  }

  /* Flux in z-direction */
  dimids[0] = nx_id;
  dimids[1] = ny_id;
  dimids[2] = nspect_id;

  if ((ierror = nc_def_var( ncid,  FLUX_NAME, NC_FLOAT, 3, dimids, &flux_var)))
    ERR(ierror,routineName);

  /* Array with wavelengths */
  dimids[0] = nspect_id;

  if ((ierror = nc_def_var( ncid, WAVE_NAME, NC_DOUBLE, 1, dimids, &wave_var)))
    ERR(ierror,routineName);

  /* Write units */
  if ((ierror = nc_put_att_text( ncid, wave_var,       "units",  2,
                         "nm" )))                      ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid, intensity_var,  "units",  23,
                         "J s^-1 m^-2 Hz^-1 sr^-1" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid, flux_var,       "units",  17,
                         "J s^-1 m^-2 Hz^-1" )))       ERR(ierror,routineName);
  
  if (atmos.Stokes || input.backgr_pol) {
    if ((ierror = nc_put_att_text( ncid, stokes_q_var, "units",  23,
                         "J s^-1 m^-2 Hz^-1 sr^-1" ))) ERR(ierror,routineName);
    if ((ierror = nc_put_att_text( ncid, stokes_u_var, "units",  23,
                         "J s^-1 m^-2 Hz^-1 sr^-1" ))) ERR(ierror,routineName);
    if ((ierror = nc_put_att_text( ncid, stokes_v_var, "units",  23,
                         "J s^-1 m^-2 Hz^-1 sr^-1" ))) ERR(ierror,routineName);
  }

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


  return;
}
/* ------- end   -------------------------- writeSpectrum_p.c --------- */
