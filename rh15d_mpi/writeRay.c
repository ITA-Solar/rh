/* ------- file: -------------------------- writeray.c ----- --------

       Version:       rh2.0
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Thu Jan 13 16:35:11 2011 --

       --------------------------                      ----------RH-- */

/* --- Writes ray data into output file  --            -------------- */

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "constant.h"
#include "background.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"
#include "io.h"


/* --- Function prototypes --                          -------------- */
void init_hdf5_ray_new(void);
void init_hdf5_ray_existing(void);
void loadBackground(int la, int mu, bool_t to_obs);

/* --- Global variables --                             -------------- */

extern enum Topology topology;

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern Geometry geometry;
extern char messageStr[];
extern Input_Atmos_file infile;
extern MPI_data mpi;
extern IO_data io;


/* ------- begin --------------------------   init_hdf5_ray.c  ------ */
void init_hdf5_ray(void) {
  /* Wrapper to find out if we should use old file or create new one. */

  if (input.p15d_rerun) init_hdf5_ray_existing(); else init_hdf5_ray_new();

  return;
}
/* ------- end   --------------------------   init_hdf5_ray.c  ------ */

/* ------- begin -------------------------- init_hdf5_ray_new.c ----- */
void init_hdf5_ray_new(void)
/* Creates the file for the ray */
{
  const char routineName[] = "init_hdf5_ray_new";
  int     k;
  hid_t   plist, ncid, file_dspace, id_x, id_y, id_z, id_wave, id_wave_sel;
  hsize_t dims[4];
  bool_t  write_xtra;
  double *lambda_air, *wave_selected, *tmp;
  char    timestr[ARR_STRLEN];
  time_t  curtime;
  struct tm *loctime;

  write_xtra = (io.ray_nwave_sel > 0);

  /* Create the file with parallel MPI-IO access */
  if (( plist = H5Pcreate(H5P_FILE_ACCESS )) < 0) HERR(routineName);
  if (( H5Pset_fapl_mpio(plist, mpi.comm, mpi.info) ) < 0) HERR(routineName);
  if (( ncid = H5Fcreate(RAY_FILE, H5F_ACC_TRUNC, H5P_DEFAULT, plist) ) < 0)
    HERR(routineName);
  if (( H5Pclose(plist) ) < 0) HERR(routineName);

  /* Global attributes */
  if (( H5LTset_attribute_string(ncid, "/", "atmosID", atmos.ID)) < 0)
    HERR(routineName);
  if (( H5LTset_attribute_ushort(ncid, "/", "snapshot_number",
			 (const unsigned short *) &mpi.snap_number, 1 ) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid, "/", "rev_id", mpi.rev_id) ) < 0)
    HERR(routineName);
  /* Write old netcdf dimensions as global attributes */
  if (( H5LTset_attribute_int(ncid, "/", "nx", &mpi.nx, 1) ) < 0)
      HERR(routineName);
  if (( H5LTset_attribute_int(ncid, "/", "ny", &mpi.ny, 1) ) < 0)
      HERR(routineName);
  if (( H5LTset_attribute_int(ncid, "/", "nz", &infile.nz, 1 )) < 0)
      HERR(routineName);
  if (( H5LTset_attribute_int(ncid, "/", "nwave", &spectrum.Nspect, 1)) < 0)
      HERR(routineName);
  if (write_xtra) {
    if (( H5LTset_attribute_int(ncid, "/", WAVE_SEL,
				&io.ray_nwave_sel, 1) ) < 0) HERR(routineName);
  }

  /* --- Write dimensions: x, y, z, wavelength, wavelength selected --- */
  dims[0] = mpi.nx;
  if (( H5LTmake_dataset(ncid, X_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                         geometry.xscale) ) < 0)  HERR(routineName);
  if (( id_x = H5Dopen2(ncid, X_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
  dims[0] = mpi.ny;
  if (( H5LTmake_dataset(ncid, Y_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                         geometry.yscale) ) < 0)  HERR(routineName);
  if (( id_y = H5Dopen2(ncid, Y_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
  dims[0] = infile.nz;
  tmp = (double *) calloc(infile.nz , sizeof(double));
  if (( H5LTmake_dataset(ncid, ZOUT_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                         tmp) ) < 0)  HERR(routineName);
  /* Define as dimension only (no coordinates) in netCDF standard */
  if (( H5LTset_attribute_string(ncid, ZOUT_NAME, "NAME",
                                 NETCDF_COMPAT) ) < 0) HERR(routineName);
  if (( id_z = H5Dopen2(ncid, ZOUT_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
  free(tmp);

  wave_selected = (double *) malloc(io.ray_nwave_sel * sizeof(double));
  dims[0] = spectrum.Nspect;
  if (spectrum.vacuum_to_air) {
    lambda_air = (double *) malloc(spectrum.Nspect * sizeof(double));
    vacuum_to_air(spectrum.Nspect, spectrum.lambda, lambda_air);
    if (( H5LTmake_dataset(ncid, WAVE_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                           lambda_air) ) < 0) HERR(routineName);
    for (k = 0;  k < io.ray_nwave_sel;  k++)
        wave_selected[k] = lambda_air[io.ray_wave_idx[k]];
    free(lambda_air);
  } else {
    if (( H5LTmake_dataset(ncid, WAVE_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                           spectrum.lambda) ) < 0) HERR(routineName);
    for (k = 0;  k < io.ray_nwave_sel;  k++)
        wave_selected[k] = spectrum.lambda[io.ray_wave_idx[k]];
  }
  if (( id_wave = H5Dopen2(ncid, WAVE_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
  if (write_xtra) {
    dims[0] = io.ray_nwave_sel;
    if (( H5LTmake_dataset(ncid, WAVE_SEL_IDX, 1, dims, H5T_NATIVE_INT,
                           io.ray_wave_idx) ) < 0)  HERR(routineName);
    if (( H5LTmake_dataset(ncid, WAVE_SEL_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                           wave_selected) ) < 0)  HERR(routineName);
    if (( id_wave_sel = H5Dopen2(ncid, WAVE_SEL_NAME,
                                 H5P_DEFAULT)) < 0) HERR(routineName);
  }
  free(wave_selected);

  /* Create dataspace with same fill value as netcdf */
  dims[0] = mpi.nx;
  dims[1] = mpi.ny;
  dims[2] = spectrum.Nspect;
  if (( file_dspace = H5Screate_simple(3, dims, NULL) ) < 0) HERR(routineName);
  if (( plist = H5Pcreate(H5P_DATASET_CREATE) ) < 0) HERR(routineName);
  if (( H5Pset_fill_value(plist, H5T_NATIVE_FLOAT, &FILLVALUE) ) < 0)
    HERR(routineName);
  if (( H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY) ) < 0) HERR(routineName);
  if (( H5Pset_fill_time(plist, H5D_FILL_TIME_ALLOC) ) < 0) HERR(routineName);
  /* Intensity */
  if (( io.ray_int_var = H5Dcreate(ncid, INT_NAME, H5T_NATIVE_FLOAT,
				    file_dspace, H5P_DEFAULT, plist,
            H5P_DEFAULT)) < 0) HERR(routineName);
  /* Attach dimension scales */
  if (( H5DSattach_scale(io.ray_int_var, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.ray_int_var, id_y, 1)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.ray_int_var, id_wave, 2)) < 0) HERR(routineName);


  /* Other Stokes parameters, if available */
  if (atmos.Stokes || input.backgr_pol) {
    if (( io.ray_stokes_q_var = H5Dcreate(ncid, STOKES_Q, H5T_NATIVE_FLOAT,
				          file_dspace, H5P_DEFAULT, plist,
				          H5P_DEFAULT)) < 0) HERR(routineName);
    /* Attach dimension scales */
    if (( H5DSattach_scale(io.ray_stokes_q_var, id_x, 0)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_stokes_q_var, id_y, 1)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_stokes_q_var, id_wave, 2)) < 0) HERR(routineName);

    if (( io.ray_stokes_u_var = H5Dcreate(ncid, STOKES_U, H5T_NATIVE_FLOAT,
				          file_dspace, H5P_DEFAULT, plist,
				          H5P_DEFAULT)) < 0) HERR(routineName);
    /* Attach dimension scales */
    if (( H5DSattach_scale(io.ray_stokes_u_var, id_x, 0)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_stokes_u_var, id_y, 1)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_stokes_u_var, id_wave, 2)) < 0) HERR(routineName);

    if (( io.ray_stokes_v_var = H5Dcreate(ncid, STOKES_V, H5T_NATIVE_FLOAT,
				          file_dspace, H5P_DEFAULT, plist,
				          H5P_DEFAULT)) < 0) HERR(routineName);
    /* Attach dimension scales */
    if (( H5DSattach_scale(io.ray_stokes_v_var, id_x, 0)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_stokes_v_var, id_y, 1)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_stokes_v_var, id_wave, 2)) < 0) HERR(routineName);
  }
  if (input.p15d_wtau) {
    if (( io.ray_tau1_var = H5Dcreate(ncid, TAU1_NAME, H5T_NATIVE_FLOAT,
				                      file_dspace, H5P_DEFAULT, plist,
				                      H5P_DEFAULT)) < 0) HERR(routineName);
    /* Attach dimension scales */
    if (( H5DSattach_scale(io.ray_tau1_var, id_x, 0)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_tau1_var, id_y, 1)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_tau1_var, id_wave, 2)) < 0) HERR(routineName);
  }
  if (write_xtra) {
    dims[0] = mpi.nx;
    dims[1] = mpi.ny;
    dims[2] = infile.nz;
    dims[3] = io.ray_nwave_sel;
    if (( file_dspace = H5Screate_simple(4, dims, NULL) ) < 0) HERR(routineName);
    /* Source function, opacity and emissivity, line and continuum */
    if (( io.ray_chi_var = H5Dcreate(ncid, CHI_NAME, H5T_NATIVE_FLOAT,
				                     file_dspace, H5P_DEFAULT, plist,
				                     H5P_DEFAULT)) < 0) HERR(routineName);
    /* Attach dimension scales */
    if (( H5DSattach_scale(io.ray_chi_var, id_x, 0)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_chi_var, id_y, 1)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_chi_var, id_z, 2)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_chi_var, id_wave_sel, 3)) < 0) HERR(routineName);
    if (( io.ray_S_var = H5Dcreate(ncid, S_NAME, H5T_NATIVE_FLOAT,
				                   file_dspace, H5P_DEFAULT, plist,
				                   H5P_DEFAULT)) < 0) HERR(routineName);
    /* Attach dimension scales */
    if (( H5DSattach_scale(io.ray_S_var, id_x, 0)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_S_var, id_y, 1)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_S_var, id_z, 2)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_S_var, id_wave_sel, 3)) < 0) HERR(routineName);
    if (( io.ray_j_var = H5Dcreate(ncid, "Jlambda", H5T_NATIVE_FLOAT,
				                   file_dspace, H5P_DEFAULT, plist,
				                   H5P_DEFAULT)) < 0) HERR(routineName);
    /* Attach dimension scales */
    if (( H5DSattach_scale(io.ray_j_var, id_x, 0)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_j_var, id_y, 1)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_j_var, id_z, 2)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_j_var, id_wave_sel, 3)) < 0) HERR(routineName);
    if (( io.ray_sca_c_var = H5Dcreate(ncid, SCA_C_NAME, H5T_NATIVE_FLOAT,
				                       file_dspace, H5P_DEFAULT, plist,
				                       H5P_DEFAULT)) < 0) HERR(routineName);
    /* Attach dimension scales */
    if (( H5DSattach_scale(io.ray_sca_c_var, id_x, 0)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_sca_c_var, id_y, 1)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_sca_c_var, id_z, 2)) < 0) HERR(routineName);
    if (( H5DSattach_scale(io.ray_sca_c_var, id_wave_sel, 3)) < 0) HERR(routineName);
  }

  /* --- Write attributes --- */
  /* Time of creation in ISO 8601 */
  curtime = time(NULL);
  loctime = localtime(&curtime);
  strftime(timestr, ARR_STRLEN, "%Y-%m-%dT%H:%M:%S%z", loctime);
  if (( H5LTset_attribute_string(ncid, "/", "creation_time",
                            (const char *) &timestr) ) < 0) HERR(routineName);

  /*  units and fill value  */
  if (( H5LTset_attribute_string(ncid, X_NAME, "units",
                                 "m") ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid, Y_NAME, "units",
                                 "m") ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid, ZOUT_NAME, "units",
                                 "m") ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid, WAVE_NAME, "units",
                                 "nm") ) < 0) HERR(routineName);
  if (write_xtra) if (( H5LTset_attribute_string(ncid, WAVE_SEL_NAME, "units",
                                                 "nm") ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid, INT_NAME, "units",
                           "J s^-1 m^-2 Hz^-1 sr^-1") ) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid, INT_NAME, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);


  if (atmos.Stokes || input.backgr_pol) {
    if (( H5LTset_attribute_string(ncid, STOKES_Q, "units",
                           "J s^-1 m^-2 Hz^-1 sr^-1") ) < 0) HERR(routineName);
    if (( H5LTset_attribute_float(ncid, STOKES_Q, "_FillValue",
                                  &FILLVALUE, 1) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid, STOKES_U, "units",
                           "J s^-1 m^-2 Hz^-1 sr^-1") ) < 0) HERR(routineName);
    if (( H5LTset_attribute_float(ncid, STOKES_U, "_FillValue",
                                  &FILLVALUE, 1) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid, STOKES_V, "units",
                           "J s^-1 m^-2 Hz^-1 sr^-1") ) < 0) HERR(routineName);
    if (( H5LTset_attribute_float(ncid, STOKES_V, "_FillValue",
                                  &FILLVALUE, 1) ) < 0) HERR(routineName);
  }

  if (input.p15d_wtau) {
    if (( H5LTset_attribute_string(ncid, TAU1_NAME, "units",
                                   "m") ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid, TAU1_NAME, DESC_NAME,
                     "Height of optical depth unity") ) < 0) HERR(routineName);
  }

  if (write_xtra) {
    if (( H5LTset_attribute_string(ncid, S_NAME, "units",
                          "J s^-1 m^-2 Hz^-1 sr^-1" ) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid, S_NAME, DESC_NAME,
         "Total source function (line + continuum)" ) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_float(ncid, S_NAME, "_FillValue",
                                  &FILLVALUE, 1) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid, CHI_NAME, "units",
                                   "m^-1" ) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid, CHI_NAME, DESC_NAME,
              "Total absorption (line + continuum)" ) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_float(ncid, CHI_NAME, "_FillValue",
                                  &FILLVALUE, 1) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid, "/Jlambda", "units",
                          "J s^-1 m^-2 Hz^-1 sr^-1" ) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid, "/Jlambda", DESC_NAME,
                             "Mean radiation field" ) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_float(ncid, "/Jlambda", "_FillValue",
                                  &FILLVALUE, 1) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid, SCA_C_NAME, DESC_NAME,
            "Scattering term multiplied by Jlambda" ) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_float(ncid, SCA_C_NAME, "_FillValue",
                                  &FILLVALUE, 1) ) < 0) HERR(routineName);
  }
  if (( H5Pclose(plist) ) < 0 ) HERR(routineName);
  if (( H5Dclose(id_x) ) < 0) HERR(routineName);
  if (( H5Dclose(id_y) ) < 0) HERR(routineName);
  if (( H5Dclose(id_z) ) < 0) HERR(routineName);
  if (( H5Dclose(id_wave) ) < 0) HERR(routineName);
  if (write_xtra) if (( H5Dclose(id_wave_sel) ) < 0) HERR(routineName);

  /* Flush ensures file is created in case of crash */
  if (( H5Fflush(ncid, H5F_SCOPE_LOCAL) ) < 0) HERR(routineName);
  /* Copy stuff to the IO data struct */
  io.ray_ncid = ncid;
  return;
}

/* ------- end   -------------------------- init_hdf5_ray_new.c ------- */


/* ------- begin   -------------------- init_hdf5_ray_existing.c ------ */
void init_hdf5_ray_existing(void)
/* Opens an existing ray file */
{
  const char routineName[] = "init_hdf5_ray_existing";
  int     ncid;
  bool_t  write_xtra;
  size_t  attr_size;
  hid_t   plist;
  char   *atmosID;
  H5T_class_t type_class;

  write_xtra = (io.ray_nwave_sel > 0);

  /* Open the file with parallel MPI-IO access */
  if (( plist = H5Pcreate(H5P_FILE_ACCESS )) < 0) HERR(routineName);
  if (( H5Pset_fapl_mpio(plist, mpi.comm, mpi.info) ) < 0) HERR(routineName);
  if (( ncid = H5Fopen(RAY_FILE, H5F_ACC_RDWR, plist) ) < 0) HERR(routineName);
  if (( H5Pclose(plist) ) < 0) HERR(routineName);

  io.ray_ncid = ncid;
  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if (( H5LTget_attribute_info(ncid, "/", "atmosID", NULL, &type_class,
                               &attr_size) ) < 0) HERR(routineName);
  atmosID = (char *) malloc(attr_size + 1);
  if (( H5LTget_attribute_string(ncid, "/", "atmosID", atmosID) ) < 0)
    HERR(routineName);
  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,
         "Ray file was calculated for different atmosphere (%s) than current",
	    atmosID);
    Error(WARNING, routineName, messageStr);
    }
  free(atmosID);
  /* --- Open datasets collectively ---*/
  if (( io.ray_int_var = H5Dopen(ncid, INT_NAME, H5P_DEFAULT) ) < 0)
    HERR(routineName);
  if (input.p15d_wtau) {
    if (( io.ray_tau1_var = H5Dopen(ncid, TAU1_NAME, H5P_DEFAULT) ) < 0)
      HERR(routineName);
  }
  if (atmos.Stokes || input.backgr_pol) {
    if (( io.ray_stokes_q_var = H5Dopen(ncid, STOKES_Q, H5P_DEFAULT) ) < 0)
      HERR(routineName);
    if (( io.ray_stokes_u_var = H5Dopen(ncid, STOKES_U, H5P_DEFAULT) ) < 0)
      HERR(routineName);
    if (( io.ray_stokes_v_var = H5Dopen(ncid, STOKES_V, H5P_DEFAULT) ) < 0)
      HERR(routineName);
  }
  if (write_xtra) {
    if (( io.ray_chi_var = H5Dopen(ncid, CHI_NAME, H5P_DEFAULT) ) < 0)
      HERR(routineName);
    if (( io.ray_S_var = H5Dopen(ncid, S_NAME, H5P_DEFAULT) ) < 0)
      HERR(routineName);
    if (( io.ray_j_var = H5Dopen(ncid, "Jlambda", H5P_DEFAULT) ) < 0)
      HERR(routineName);
    if (( io.ray_sca_c_var = H5Dopen(ncid, SCA_C_NAME, H5P_DEFAULT) ) < 0)
      HERR(routineName);
  }
  return;
}

/* ------- end   --------------------- init_hdf5_ray_existing.c ------ */


/* ------- begin -------------------------- close_hdf5_ray.c --------- */
void close_hdf5_ray(void) {
  /* Closes the ray file */
  const char routineName[] = "close_hdf5_ray";
  bool_t     write_xtra;

  write_xtra = (io.ray_nwave_sel > 0);
  /* Close all datasets */
  if (( H5Dclose(io.ray_int_var) ) < 0) HERR(routineName);
  if (atmos.Stokes || input.backgr_pol) {
    if (( H5Dclose(io.ray_stokes_q_var) ) < 0) HERR(routineName);
    if (( H5Dclose(io.ray_stokes_u_var) ) < 0) HERR(routineName);
    if (( H5Dclose(io.ray_stokes_v_var) ) < 0) HERR(routineName);
  }
  if (input.p15d_wtau) {
    if (( H5Dclose(io.ray_tau1_var) ) < 0) HERR(routineName);
  }
  if (write_xtra) {
    if (( H5Dclose(io.ray_chi_var ) ) < 0) HERR(routineName);
    if (( H5Dclose(io.ray_S_var ) ) < 0) HERR(routineName);
    if (( H5Dclose(io.ray_j_var ) ) < 0) HERR(routineName);
    if (( H5Dclose(io.ray_sca_c_var ) ) < 0) HERR(routineName);
  }
  if (( H5Fclose(io.ray_ncid) ) < 0) HERR(routineName);
  return;
}
/* ------- end   -------------------------- close_hdf5_ray.c --------- */


/* ------- begin -------------------------- writeRay.c --------------- */
void writeRay(void) {
  /* Writes ray data to file. */
  const char routineName[] = "writeRay";
  int        idx, ncid, k, l, nspect;
  double    *J;
  float    **chi, **S, **sca, *tau_one, tau_cur, tau_prev, tmp, *chi_tmp;
  float    **Jnu;
  hsize_t    offset[] = {0, 0, 0, 0};
  hsize_t    count[] = {1, 1, 1, 1};
  hsize_t    dims[4];
  bool_t     write_xtra, crosscoupling, to_obs, initialize,prdh_limit_mem_save;
  ActiveSet *as;
  hid_t      file_dataspace, mem_dataspace;

  write_xtra = (io.ray_nwave_sel > 0);
  ncid = io.ray_ncid;

  /* Memory dataspace */
  dims[0] = spectrum.Nspect;
  if (( mem_dataspace = H5Screate_simple(1, dims, NULL) ) < 0)
    HERR(routineName);
  /* File dataspace */
  offset[0] = mpi.ix;
  offset[1] = mpi.iy;
  count[2] = spectrum.Nspect;
  if (( file_dataspace = H5Dget_space(io.ray_int_var) ) < 0) HERR(routineName);
  if (( H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset,
                            NULL, count, NULL) ) < 0) HERR(routineName);

  /* Write intensity */
  if (( H5Dwrite(io.ray_int_var, H5T_NATIVE_DOUBLE, mem_dataspace,
          file_dataspace, H5P_DEFAULT, spectrum.I[0]) ) < 0) HERR(routineName);

  /* Calculate height of tau=1 and write to file*/
  if (input.p15d_wtau) {
    tau_one = (float *) calloc(spectrum.Nspect, sizeof(float));

    /* set switch so that shift of rho_prd is done with a fresh interpolation */
    prdh_limit_mem_save = FALSE;
    if (input.prdh_limit_mem) prdh_limit_mem_save = TRUE;
    input.prdh_limit_mem = TRUE;

    for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
      as  = &spectrum.as[nspect];
      alloc_as(nspect, crosscoupling=FALSE);
      Opacity(nspect, 0, to_obs=TRUE, initialize=TRUE);
      chi_tmp = (float *) calloc(infile.nz, sizeof(float));

      if (input.backgr_in_mem) {
      	loadBackground(nspect, 0, to_obs=TRUE);
      } else {
      	readBackground(nspect, 0, to_obs=TRUE);
      }

      tau_prev = 0.0;
      tau_cur  = 0.0;

      for (k = 0;  k < atmos.Nspace;  k++) {
        l = k + mpi.zcut;
        chi_tmp[l] = (float) (as->chi[k] + as->chi_c[k]);

        /* Calculate tau=1 depth, manual linear interpolation */
        if (k > 0) {
          /* Manually integrate tau, trapezoidal rule*/
          tau_prev = tau_cur;
          tmp = 0.5 * (geometry.height[k - 1] - geometry.height[k]);
          tau_cur  = tau_prev + tmp * (chi_tmp[l] + chi_tmp[l - 1]);

          if (((tau_cur > 1.) && (k == 1)) || ((tau_cur < 1.) &&
                                               (k == atmos.Nspace - 1))) {
            tau_one[nspect] = geometry.height[k];
          }
          else if ((tau_cur > 1.) && (tau_prev < 1.)) {
            /* Manual linear interpolation for tau=1 */
            tmp = (1. - tau_prev) / (tau_cur - tau_prev);
            tau_one[nspect] = geometry.height[k - 1] + tmp *
              (geometry.height[k] - geometry.height[k - 1]);
          }
        }
      }
      free_as(nspect, crosscoupling=FALSE);
      free(chi_tmp);
    }
    /* Write to file */
    if (( H5Dwrite(io.ray_tau1_var, H5T_NATIVE_FLOAT, mem_dataspace,
                 file_dataspace, H5P_DEFAULT, tau_one) ) < 0) HERR(routineName);
    /* set back PRD input option */
    if (input.PRD_angle_dep == PRD_ANGLE_APPROX && atmos.NPRDactive > 0)
      input.prdh_limit_mem = prdh_limit_mem_save ;
    free(tau_one);
  }

  if (atmos.Stokes || input.backgr_pol) { /* Write rest of Stokes vector */
    if (( H5Dwrite(io.ray_stokes_q_var, H5T_NATIVE_DOUBLE, mem_dataspace,
                   file_dataspace, H5P_DEFAULT,
                   spectrum.Stokes_Q[0]) ) < 0) HERR(routineName);
    if (( H5Dwrite(io.ray_stokes_u_var, H5T_NATIVE_DOUBLE, mem_dataspace,
                   file_dataspace, H5P_DEFAULT,
                   spectrum.Stokes_U[0]) ) < 0) HERR(routineName);
    if (( H5Dwrite(io.ray_stokes_v_var, H5T_NATIVE_DOUBLE, mem_dataspace,
                   file_dataspace, H5P_DEFAULT,
                   spectrum.Stokes_V[0]) ) < 0) HERR(routineName);
  }
  /* release dataspace resources */
  if (( H5Sclose(mem_dataspace) ) < 0) HERR(routineName);
  if (( H5Sclose(file_dataspace) ) < 0) HERR(routineName);

  if (write_xtra) {
    /* Write opacity and emissivity for line and continuum */
    chi = matrix_float(infile.nz, io.ray_nwave_sel);
    S   = matrix_float(infile.nz, io.ray_nwave_sel);
    sca = matrix_float(infile.nz, io.ray_nwave_sel);
    Jnu = matrix_float(infile.nz, io.ray_nwave_sel);

    if (input.limit_memory) J = (double *) malloc(atmos.Nspace *
                                                  sizeof(double));
    /* set switch so that shift of rho_prd is done with a fresh
       interpolation */
    prdh_limit_mem_save = FALSE;
    if (input.prdh_limit_mem) prdh_limit_mem_save = TRUE;
    input.prdh_limit_mem = TRUE;

    for (nspect = 0; nspect < io.ray_nwave_sel; nspect++) {
      idx = io.ray_wave_idx[nspect];
      as  = &spectrum.as[idx];

      alloc_as(idx, crosscoupling=FALSE);
      Opacity(idx, 0, to_obs=TRUE, initialize=TRUE);

      if (input.backgr_in_mem) {
        loadBackground(idx, 0, to_obs=TRUE);
      } else {
        readBackground(idx, 0, to_obs=TRUE);
      }

      if (input.limit_memory) {
        //readJlambda_single(idx, J);
      } else
        J = spectrum.J[idx];

      /* Zero S and chi  */
      for (k = 0; k < infile.nz; k++) {
        chi[k][nspect] = FILLVALUE;
        S[k][nspect]   = FILLVALUE;
        sca[k][nspect] = FILLVALUE;
        Jnu[k][nspect] = FILLVALUE;
      }

      for (k = 0;  k < atmos.Nspace;  k++) {
        l = k + mpi.zcut;
        chi[l][nspect] = (float) (as->chi[k] + as->chi_c[k]);
      	sca[l][nspect] = (float) (as->sca_c[k] * J[k]);
      	S[l][nspect]   = (float) ((as->eta[k] + as->eta_c[k] +
                           as->sca_c[k] * J[k]) / (as->chi[k] + as->chi_c[k]));
        Jnu[l][nspect] = (float) J[k];
      }
      free_as(idx, crosscoupling=FALSE);
    }

    /* set back PRD input option */
    if (input.PRD_angle_dep == PRD_ANGLE_APPROX && atmos.NPRDactive > 0)
      input.prdh_limit_mem = prdh_limit_mem_save;

    /* Write variables */
    /* Memory dataspace */
    dims[0] = infile.nz;
    dims[1] = io.ray_nwave_sel;
    if (( mem_dataspace = H5Screate_simple(2, dims, NULL) ) < 0)
      HERR(routineName);
    /* File dataspace */
    offset[0] = mpi.ix;  count[0] = 1;
    offset[1] = mpi.iy;  count[1] = 1;
    offset[2] = 0;       count[2] = infile.nz;
    offset[3] = 0;       count[3] = io.ray_nwave_sel;
    if (( file_dataspace = H5Dget_space(io.ray_S_var) ) < 0) HERR(routineName);
    if (( H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, offset,
                              NULL, count, NULL) ) < 0) HERR(routineName);
    if (( H5Dwrite(io.ray_chi_var, H5T_NATIVE_FLOAT, mem_dataspace,
                 file_dataspace, H5P_DEFAULT, chi[0]) ) < 0) HERR(routineName);
    if (( H5Dwrite(io.ray_S_var, H5T_NATIVE_FLOAT, mem_dataspace,
                 file_dataspace, H5P_DEFAULT, S[0]) ) < 0) HERR(routineName);
    if (( H5Dwrite(io.ray_j_var, H5T_NATIVE_FLOAT, mem_dataspace,
                 file_dataspace, H5P_DEFAULT, Jnu[0]) ) < 0) HERR(routineName);
    if (( H5Dwrite(io.ray_sca_c_var, H5T_NATIVE_FLOAT, mem_dataspace,
                 file_dataspace, H5P_DEFAULT, sca[0]) ) < 0) HERR(routineName);
    freeMatrix((void **) chi);
    freeMatrix((void **) S);
    freeMatrix((void **) Jnu);
    freeMatrix((void **) sca);
    /* release dataspace resources */
    if (( H5Sclose(mem_dataspace) ) < 0) HERR(routineName);
    if (( H5Sclose(file_dataspace) ) < 0) HERR(routineName);
    if (input.limit_memory) free(J);
  }
  close_Background();  /* To avoid many open files */
  return;
}
/* ------- end   -------------------------- writeRay.c -------------- */


/* ------- begin -------------------------- calculate_ray.c --------- */
void calculate_ray(void) {
  /* performs necessary reallocations and inits, and solves for ray */
  int i, nact, ierror, mu, k;
      bool_t analyze_output, equilibria_only, prdh_limit_mem_save;
      Atom *atom;
      AtomicLine *line;

  close_Background();  /* Was opened previously, will be opened here again */

  /* Recalculate magnetic field projections */
  if (atmos.Stokes) {
    if (atmos.cos_gamma != NULL) {
      freeMatrix((void **) atmos.cos_gamma);
      atmos.cos_gamma = NULL;
    }
    if (atmos.cos_2chi != NULL) {
      freeMatrix((void **) atmos.cos_2chi);
      atmos.cos_2chi = NULL;
    }
    if (atmos.sin_2chi != NULL) {
      freeMatrix((void **) atmos.sin_2chi);
      atmos.sin_2chi = NULL;
    }
    Bproject();
  }

  /* Must calculate background opacity for new ray, need some prep first */
  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];

    // TIAGO: commented below
    /* Rewind atom files to point just before collisional data
    if ((ierror = fseek(atom->fp_input, io.atom_file_pos[nact], SEEK_SET))) {
      sprintf(messageStr, "Unable to rewind atom file for %s", atom->ID);
      Error(ERROR_LEVEL_2, "rh15d_ray", messageStr);
    }
    */

    /* Free collisions matrix, going to be re-allocated in background */
    if (atom->C != NULL) {
      freeMatrix((void **) atom->C);
      atom->C = NULL;
    }

    /* Recalculate line profiles for new angle */
    for (i = 0; i < atom->Nline; i++) {
      line = &atom->line[i];

      if (line->phi  != NULL) {
        freeMatrix((void **) line->phi);
        line->phi = NULL;
      }
      if (line->wphi != NULL) {
        free(line->wphi);
        line->wphi = NULL;
      }

      if (atmos.moving && line->polarizable && (input.StokesMode>FIELD_FREE)) {
        if (line->phi_Q != NULL) {
          freeMatrix((void **) line->phi_Q);
          line->phi_Q = NULL;
        }
        if (line->phi_U != NULL) {
          freeMatrix((void **) line->phi_U);
          line->phi_U = NULL;
        }
        if (line->phi_V != NULL) {
          freeMatrix((void **) line->phi_V);
          line->phi_V = NULL;
        }
        if (input.magneto_optical) {
          if (line->psi_Q != NULL) {
            freeMatrix((void **) line->psi_Q);
            line->psi_Q = NULL;
          }
          if (line->psi_U != NULL) {
            freeMatrix((void **) line->psi_U);
            line->psi_U = NULL;
          }
          if (line->psi_V != NULL) {
            freeMatrix((void **) line->psi_V);
            line->psi_V = NULL;
          }
        }
      }
      Profile(line);
    }
  }

  /* reallocate intensities for correct number of rays */
  if (spectrum.I != NULL) freeMatrix((void **) spectrum.I);
  spectrum.I = matrix_double(spectrum.Nspect, atmos.Nrays);
  if (atmos.Stokes || input.backgr_pol) {
    if (spectrum.Stokes_Q != NULL) freeMatrix((void **) spectrum.Stokes_Q);
    spectrum.Stokes_Q = matrix_double(spectrum.Nspect, atmos.Nrays);
    if (spectrum.Stokes_U != NULL) freeMatrix((void **) spectrum.Stokes_U);
    spectrum.Stokes_U = matrix_double(spectrum.Nspect, atmos.Nrays);
    if (spectrum.Stokes_V != NULL) freeMatrix((void **) spectrum.Stokes_V);
    spectrum.Stokes_V = matrix_double(spectrum.Nspect, atmos.Nrays);
  }

  if (input.PRD_angle_dep == PRD_ANGLE_APPROX &&  atmos.NPRDactive > 0) {
    // recalculate line-of-sight velocity
    if (spectrum.v_los != NULL) freeMatrix((void **) spectrum.v_los);
    spectrum.v_los = matrix_double( atmos.Nrays, atmos.Nspace);
    for (mu = 0;  mu < atmos.Nrays;  mu++) {
      for (k = 0;  k < atmos.Nspace;  k++) {
        spectrum.v_los[mu][k] = vproject(k, mu);
      }
    }

    /* set switch so that shift of rho_prd is done with a fresh
       interpolation */
    prdh_limit_mem_save = FALSE;
    if (input.prdh_limit_mem) prdh_limit_mem_save = TRUE;
    input.prdh_limit_mem = TRUE;
  }

  Background_p(analyze_output=FALSE, equilibria_only=FALSE);

  /* --- Solve radiative transfer for ray --           -------------- */
  solveSpectrum(FALSE, FALSE);

  // set back PRD input option
  if (input.PRD_angle_dep == PRD_ANGLE_APPROX && atmos.NPRDactive > 0)
    input.prdh_limit_mem = prdh_limit_mem_save ;

}

/* ------- end   -------------------------- calculate_ray.c --------- */
