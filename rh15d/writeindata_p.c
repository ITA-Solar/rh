/* ------- file: -------------------------- writeindata_p.c ------------

       Version:       rh2.0, 1.5-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Thu Dec 30 22:58:15 2010 --

       --------------------------                      -----------RH-- */

/* --- Writes input data (including input, atmos, geom) to output file */
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"
#include "io.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern InputData input;
extern char messageStr[];
extern Input_Atmos_file infile;
extern MPI_data mpi;
extern IO_data io;


/* ------- begin --------------------------   init_hdf5_indata.c  --- */
void init_hdf5_indata(void) {
  /* Wrapper to find out if we should use old file or create new one */

  if (input.p15d_rerun) init_hdf5_indata_existing(); else init_hdf5_indata_new();

  return;
}
/* ------- end   --------------------------   init_hdf5_indata.c  --- */


/* ------- begin --------------------------   init_hdf5_indata.c  --- */
void init_hdf5_indata_new(void)
/* Creates the file for the input data */
{
  const char routineName[] = "init_hdf5_indata_new";
  unsigned int *tmp;
  int     i, PRD_angle_dep;
  double  *eweight, *eabund, *tmp_double;;
  /* This value is harcoded for efficiency.
     Maximum number of iterations ever needed */
  int     NMaxIter = 1500;
  int     atom_name_size = 8;
  hid_t   plist, ncid, file_dspace, ncid_input, ncid_atmos, ncid_mpi, ncid_tmp;
  hid_t   id_x, id_y, id_z, id_n, id_tmp;
  hsize_t dims[4];
  bool_t  XRD;
  char    startJ[MAX_LINE_SIZE], StokesMode[MAX_LINE_SIZE];
  char    angleSet[MAX_LINE_SIZE], group_name[ARR_STRLEN], **atom_names;
  Atom   *atom;

  if (input.NmaxIter > NMaxIter) {  /* Adjust NmaxIter to avoid hdf5 writing error */
      sprintf(messageStr,  
              "Increasing internal NMaxIter (%d) to N_MAX_ITER (%d).\n", 
              NMaxIter, input.NmaxIter);
      Error(WARNING, routineName, messageStr);
      NMaxIter = input.NmaxIter;
  }

  /* Create the file  */
  if (( plist = H5Pcreate(H5P_FILE_ACCESS )) < 0) HERR(routineName);
  if (( H5Pset_fapl_mpio(plist, mpi.comm, mpi.info) ) < 0) HERR(routineName);
  if (( ncid = H5Fcreate(INPUTDATA_FILE, H5F_ACC_TRUNC, H5P_DEFAULT,
                         plist) ) < 0) HERR(routineName);
  if (( H5Pclose(plist) ) < 0) HERR(routineName);

  /* Create groups */
  if (( ncid_input = H5Gcreate(ncid, "/input", H5P_DEFAULT, H5P_DEFAULT,
                               H5P_DEFAULT) ) < 0) HERR(routineName);
  if (( ncid_atmos = H5Gcreate(ncid, "/atmos", H5P_DEFAULT, H5P_DEFAULT,
                               H5P_DEFAULT) ) < 0) HERR(routineName);
  if (( ncid_mpi = H5Gcreate(ncid, "/mpi", H5P_DEFAULT, H5P_DEFAULT,
                               H5P_DEFAULT) ) < 0) HERR(routineName);

  /* --- Definitions for the root group --- */
  /* dimensions as attributes */
  if (( H5LTset_attribute_int(ncid, "/", "nx", &mpi.nx, 1) ) < 0)
      HERR(routineName);
  if (( H5LTset_attribute_int(ncid, "/", "ny", &mpi.ny, 1) ) < 0)
      HERR(routineName);
  if (( H5LTset_attribute_int(ncid, "/", "nz", (int *) &infile.nz, 1 )) < 0)
      HERR(routineName);
  /* attributes */
  if (( H5LTset_attribute_string(ncid, "/", "atmosID", atmos.ID)) < 0)
    HERR(routineName);
  if (( H5LTset_attribute_string(ncid, "/", "rev_id", mpi.rev_id) ) < 0)
    HERR(routineName);
  /* --- dimension datasets, all values set to zero --- */
  dims[0] = infile.nz;
  tmp_double = (double *) calloc(infile.nz , sizeof(double));
  if (( H5LTmake_dataset(ncid, ZOUT_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                         tmp_double) ) < 0)  HERR(routineName);
  free(tmp_double);
  if (( id_z = H5Dopen2(ncid, ZOUT_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
  /* For compatibility with netCDF readers, only use dataset as dimension */
  if (( H5LTset_attribute_string(ncid, ZOUT_NAME, "NAME",
                                 NETCDF_COMPAT) ) < 0) HERR(routineName);

  /* --- Definitions for the INPUT group --- */
  /* attributes */
  if ( atmos.NPRDactive > 0)
    PRD_angle_dep = input.PRD_angle_dep;
  else
    PRD_angle_dep=0;

  XRD = (input.XRD  &&  atmos.NPRDactive > 0);

  if (( H5LTset_attribute_uchar(ncid_input, ".", "Magneto_optical",
          (unsigned char *) &input.magneto_optical, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_uchar(ncid_input, ".", "PRD_angle_dep",
          (unsigned char *) &PRD_angle_dep, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_uchar(ncid_input, ".", "XRD",
          (unsigned char *) &XRD, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_uchar(ncid_input, ".", "Background_polarization",
          (unsigned char *) &input.backgr_pol, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_input, ".", "natoms", &atmos.Natom, 1) ) < 0)
      HERR(routineName);

  switch (input.startJ) {
  case UNKNOWN:
    strcpy(startJ, "Unknown");
    break;
  case LTE_POPULATIONS:
    strcpy(startJ, "LTE_POPULATIONS");
    break;
  case ZERO_RADIATION:
    strcpy(startJ, "ZERO_RADIATION");
    break;
  case OLD_POPULATIONS:
    strcpy(startJ, "OLD_POPULATIONS");
    break;
  case ESCAPE_PROBABILITY:
    strcpy(startJ, "ESCAPE_PROBABILITY");
    break;
  case NEW_J:
    strcpy(startJ, "NEW_J");
    break;
  case OLD_J:
    strcpy(startJ, "OLD_J");
    break;
  }
  if (( H5LTset_attribute_string(ncid_input, ".", "Start_J", startJ)) < 0)
    HERR(routineName);

  switch (input.StokesMode) {
  case NO_STOKES:
    strcpy(StokesMode, "NO_STOKES");
    break;
  case FIELD_FREE:
    strcpy(StokesMode, "FIELD_FREE");
    break;
  case POLARIZATION_FREE:
    strcpy(StokesMode, "POLARIZATION_FREE");
    break;
  case FULL_STOKES:
    strcpy(StokesMode, "FULL_STOKES");
    break;
  }
  if (( H5LTset_attribute_string(ncid_input, ".", "Stokes_mode",
                                 StokesMode) ) < 0) HERR(routineName);

  switch (atmos.angleSet.set) {
  case SET_VERTICAL:
    strcpy(angleSet, "SET_VERTICAL");
    break;
  case SET_GL:
    strcpy(angleSet, "SET_GL");
    break;
  case SET_A2:
    strcpy(angleSet, "SET_A2");
    break;
  case SET_A4:
    strcpy(angleSet, "SET_A4");
    break;
  case SET_A6:
    strcpy(angleSet, "SET_A6");
    break;
  case SET_A8:
    strcpy(angleSet, "SET_A8");
    break;
  case SET_B4:
    strcpy(angleSet, "SET_B4");
    break;
  case SET_B6:
    strcpy(angleSet, "SET_B6");
    break;
  case SET_B8:
    strcpy(angleSet, "SET_B8");
    break;
  case NO_SET:
    strcpy(angleSet, "NO_SET");
    break;
  }
  if (( H5LTset_attribute_string(ncid_input, ".", "Angle_set",
                                 angleSet) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid_input, ".", "Atmos_file",
                                 input.atmos_input) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid_input, ".", "Abundances_file",
                                 input.abund_input) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid_input, ".", "Kurucz_PF_data",
                                 input.pfData) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_double(ncid_input, ".", "Iteration_limit",
                                 &input.iterLimit, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_double(ncid_input, ".", "PRD_Iteration_limit",
                              &input.PRDiterLimit, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_input, ".", "N_max_iter",
                              &input.NmaxIter, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_input, ".", "Ng_delay",
                              &input.Ngdelay, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_input, ".", "Ng_order",
                              &input.Ngorder, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_input, ".", "Ng_period",
                              &input.Ngperiod, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_input, ".", "PRD_N_max_iter",
                              &input.PRD_NmaxIter, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_input, ".", "PRD_Ng_delay",
                              &input.PRD_Ngdelay, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_input, ".", "PRD_Ng_order",
                              &input.PRD_Ngorder, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_input, ".", "PRD_Ng_period",
                              &input.PRD_Ngperiod, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_double(ncid_input, ".", "Metallicity",
                               &input.metallicity, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_double(ncid_input, ".", "Lambda_reference",
                                &atmos.lambda_ref, 1) ) < 0) HERR(routineName);
  /* contents of input files as variables */
  if (( H5LTmake_dataset_string(ncid_input, "keyword_file_contents",
                          input.keyword_file_contents) ) < 0)  HERR(routineName);
  if (( H5LTmake_dataset_string(ncid_input, "atoms_file_contents",
                          input.atoms_file_contents) ) < 0)  HERR(routineName);
  free(input.atoms_file_contents);
  /* Make sub groups for storing contents of atom files */
  atom_names = matrix_char(atmos.Natom, atom_name_size);
  for (i=0; i < atmos.Natom; i++) {
      atom = &atmos.atoms[i];
      sprintf(group_name, (atom->ID[1] == ' ') ? "atom_%.1s" : "atom_%.2s",
              atom->ID);
      strncpy(atom_names[i], group_name, atom_name_size - 1);
      if (( ncid_tmp = H5Gcreate(ncid_input, group_name, H5P_DEFAULT,
                            H5P_DEFAULT, H5P_DEFAULT) ) < 0) HERR(routineName);
      if (( H5LTset_attribute_string(ncid_tmp, ".", "file_name",
                                     atom->atom_file)) < 0) HERR(routineName);
      if (( H5LTmake_dataset_string(ncid_tmp, "file_contents",
                                    atom->fp_input) ) < 0)  HERR(routineName);
      if (( H5Gclose(ncid_tmp) ) < 0) HERR(routineName);
  }
  dims[0] = atmos.Natom;
  dims[1] = atom_name_size;
  if (( H5LTmake_dataset(ncid_input, "atom_groups", 2, dims,
                         H5T_C_S1, atom_names[0]) ) < 0)  HERR(routineName);
  freeMatrix((void **) atom_names);
  /* Information from ray.input */
  if (( H5LTset_attribute_double(ncid_input, ".", INPUT_MU,
                                 &io.ray_muz, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_uint(ncid_input, ".", WAVE_SEL,
                (unsigned int *) &io.ray_nwave_sel, 1) ) < 0) HERR(routineName);
  if (io.ray_nwave_sel > 0) {
      dims[0] = io.ray_nwave_sel;
      if (( H5LTmake_dataset(ncid_input, WAVE_SEL_IDX, 1, dims,
                   H5T_NATIVE_UINT, io.ray_wave_idx) ) < 0)  HERR(routineName);
  }
  /* Information from wavetable */
  if (input.wavetable != NULL) {
      dims[0] = input.Nxwave;
      if (( H5LTmake_dataset_double(ncid_input, WAVETABLE, 1, dims,
                                    input.wavetable) ) < 0)  HERR(routineName);
      if (( H5LTset_attribute_uint(ncid_input, ".", NXWAVE,
                                   &input.Nxwave, 1) ) < 0) HERR(routineName);
      free(input.wavetable);
  }
  /* Information from Kurucz tables */
  if (( H5LTset_attribute_uint(ncid_input, ".", NKURUCZ,
            (unsigned int *) &input.Nkurucz_files, 1) ) < 0) HERR(routineName);
  if (input.Nkurucz_files > 0) {
      if (( H5LTmake_dataset_string(ncid_input, "kurucz_file_contents",
              input.kurucz_file_contents)) < 0)  HERR(routineName);
      for (i=0; i < input.Nkurucz_files; i++) {
          sprintf(group_name, KURUCZ_LINE_FILE, i);
          if (( ncid_tmp = H5Gcreate(ncid_input, group_name, H5P_DEFAULT,
                            H5P_DEFAULT, H5P_DEFAULT) ) < 0) HERR(routineName);
          if (( H5LTset_attribute_string(ncid_tmp, ".", "file_name",
                       input.kurucz_line_file_name[i])) < 0) HERR(routineName);
          if (( H5LTmake_dataset_string(ncid_tmp, "file_contents",
                  input.kurucz_line_file_contents[i])) < 0)  HERR(routineName);
          if (( H5Gclose(ncid_tmp) ) < 0) HERR(routineName);
      }
  }

  /* --- Definitions for the ATMOS group --- */
  /* dimensions */
  if (( H5LTset_attribute_int(ncid_atmos, ".", "nelements",
                              &atmos.Nelem, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_atmos, ".", "nrays",
                              &geometry.Nrays, 1) ) < 0) HERR(routineName);
  /* --- dimension datasets --- */
  dims[0] = mpi.nx;
  if (( H5LTmake_dataset(ncid_atmos, X_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                         geometry.xscale) ) < 0)  HERR(routineName);
  if (( id_x = H5Dopen2(ncid_atmos, X_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
  dims[0] = mpi.ny;
  if (( H5LTmake_dataset(ncid_atmos, Y_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                         geometry.yscale) ) < 0)  HERR(routineName);
  if (( id_y = H5Dopen2(ncid_atmos, Y_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
  dims[0] = atmos.Nelem;
  tmp = (unsigned int *) calloc(atmos.Nelem , sizeof(unsigned int));
  if (( H5LTmake_dataset(ncid_atmos, ELEM_NAME, 1, dims, H5T_NATIVE_UINT,
                         tmp) ) < 0)  HERR(routineName);
  free(tmp);
  dims[0] = geometry.Nrays;
  tmp = (unsigned int *) calloc(geometry.Nrays , sizeof(unsigned int));
  if (( H5LTmake_dataset(ncid_atmos, RAY_NAME, 1, dims, H5T_NATIVE_UINT,
                         tmp) ) < 0)  HERR(routineName);
  free(tmp);
  /* For compatibility with netCDF readers, only use dataset as dimension */
  if (( H5LTset_attribute_string(ncid_atmos, ELEM_NAME, "NAME",
                                 NETCDF_COMPAT) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid_atmos, RAY_NAME, "NAME",
                                 NETCDF_COMPAT) ) < 0) HERR(routineName);

  /* variables*/
  dims[0] = mpi.nx;
  dims[1] = mpi.ny;
  dims[2] = infile.nz;
  if (( file_dspace = H5Screate_simple(3, dims, NULL) ) < 0) HERR(routineName);
  if (( plist = H5Pcreate(H5P_DATASET_CREATE) ) < 0) HERR(routineName);
  if (( H5Pset_fill_value(plist, H5T_NATIVE_FLOAT, &FILLVALUE) ) < 0)
    HERR(routineName);
  if (( H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY) ) < 0) HERR(routineName);
  if (( H5Pset_fill_time(plist, H5D_FILL_TIME_ALLOC) ) < 0) HERR(routineName);
  if (( io.in_atmos_T = H5Dcreate(ncid_atmos, TEMP_NAME, H5T_NATIVE_FLOAT,
         file_dspace, H5P_DEFAULT, plist, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_T, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_T, id_y, 1)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_T, id_z, 2)) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid_atmos, TEMP_NAME, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);
  if (( io.in_atmos_vz = H5Dcreate(ncid_atmos, VZ_NAME, H5T_NATIVE_FLOAT,
         file_dspace, H5P_DEFAULT, plist, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_vz, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_vz, id_y, 1)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_vz, id_z, 2)) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid_atmos, VZ_NAME, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);
  if (( io.in_atmos_z = H5Dcreate(ncid_atmos, ZH_NAME, H5T_NATIVE_FLOAT,
         file_dspace, H5P_DEFAULT, plist, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_z, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_z, id_y, 1)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_z, id_z, 2)) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid_atmos, ZH_NAME, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);
  if (( io.in_atmos_ne = H5Dcreate(ncid_atmos, NE_NAME, H5T_NATIVE_DOUBLE,
         file_dspace, H5P_DEFAULT, plist, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_ne, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_ne, id_y, 1)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_atmos_ne, id_z, 2)) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid_atmos, NE_NAME, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);

  if (( H5Pclose(plist) ) < 0) HERR(routineName);
  if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
  /* --- Write some data that does not depend on xi, yi, ATMOS group --- */
  /* arrays of number of elements */
  eweight = (double *) malloc(atmos.Nelem * sizeof(double));
  eabund = (double *) malloc(atmos.Nelem * sizeof(double));
  for (i=0; i < atmos.Nelem; i++) {
    eweight[i] = atmos.elements[i].weight;
    eabund[i] = atmos.elements[i].abund;
  }
  dims[0] = atmos.Nelem;
  if (( id_n = H5Dopen2(ncid_atmos, ELEM_NAME,
                        H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5LTmake_dataset(ncid_atmos, "element_weight", 1, dims,
                H5T_NATIVE_DOUBLE, eweight) ) < 0) HERR(routineName);
  if (( id_tmp = H5Dopen2(ncid_atmos, "element_weight",
                          H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
  if (( H5Dclose(id_tmp) ) < 0) HERR(routineName);
  if (( H5LTmake_dataset(ncid_atmos, "element_abundance", 1, dims,
                H5T_NATIVE_DOUBLE, eabund) ) < 0) HERR(routineName);
  if (( id_tmp = H5Dopen2(ncid_atmos, "element_abundance",
                          H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
  if (( H5Dclose(id_tmp) ) < 0) HERR(routineName);
  if (( H5Dclose(id_n) ) < 0) HERR(routineName);
  /* Not writing element_id for now
  dims[1] = strlen;
  if (( H5LTmake_dataset(ncid_atmos, "element_id", 2, dims,
                H5T_C_S1, eID) ) < 0) HERR(routineName);
  */
  free(eweight);
  free(eabund);

  dims[0] = geometry.Nrays;
  if (( id_n = H5Dopen2(ncid_atmos, RAY_NAME,
                        H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5LTmake_dataset(ncid_atmos, "muz", 1, dims,
              H5T_NATIVE_DOUBLE, geometry.muz) ) < 0) HERR(routineName);
  if (( id_tmp = H5Dopen2(ncid_atmos, "muz", H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
  if (( H5Dclose(id_tmp) ) < 0) HERR(routineName);
  if (( H5LTmake_dataset(ncid_atmos, "wmu", 1, dims,
              H5T_NATIVE_DOUBLE, geometry.wmu) ) < 0) HERR(routineName);
  if (( id_tmp = H5Dopen2(ncid_atmos, "wmu", H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
  if (( H5Dclose(id_tmp) ) < 0) HERR(routineName);
  if (( H5Dclose(id_n) ) < 0) HERR(routineName);

  /* attributes */
  if (( H5LTset_attribute_uchar(ncid_atmos, ".", "moving",
                   (unsigned char *) &atmos.moving, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_uchar(ncid_atmos, ".", "stokes",
                   (unsigned char *) &atmos.Stokes, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid_atmos, TEMP_NAME, "units",
                                 UNIT_TEMP) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid_atmos, VZ_NAME, "units",
                                 UNIT_VELOCITY) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid_atmos, ZH_NAME, "units",
                                 UNIT_LENGTH) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid_atmos, NE_NAME, "units",
                                 UNIT_PER_VOLUME) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_string(ncid_atmos,  "element_weight", "units",
                                 UNIT_AMU) ) < 0) HERR(routineName);
  if (( H5Dclose(id_x) ) < 0) HERR(routineName);
  if (( H5Dclose(id_y) ) < 0) HERR(routineName);

  /* --- Definitions for the MPI group --- */
  /* dimensions */
  if (( H5LTset_attribute_int(ncid_mpi, ".", "nprocesses",
                              &mpi.size, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_mpi, ".", "niterations",
                              &NMaxIter, 1) ) < 0) HERR(routineName);
  /* --- dimension datasets --- */
  dims[0] = mpi.nx;
  if (( H5LTmake_dataset(ncid_mpi, X_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                         geometry.xscale) ) < 0)  HERR(routineName);
  if (( id_x = H5Dopen2(ncid_mpi, X_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
  dims[0] = mpi.ny;
  if (( H5LTmake_dataset(ncid_mpi, Y_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                         geometry.yscale) ) < 0)  HERR(routineName);
  if (( id_y = H5Dopen2(ncid_mpi, Y_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
  dims[0] = NMaxIter;
  tmp = (unsigned int *) calloc(NMaxIter , sizeof(unsigned int));
  if (( H5LTmake_dataset(ncid_mpi, IT_NAME, 1, dims, H5T_NATIVE_UINT,
                         tmp) ) < 0)  HERR(routineName);
  free(tmp);
  /* For compatibility with netCDF readers, only use dataset as dimension */
  if (( H5LTset_attribute_string(ncid_mpi, IT_NAME, "NAME",
                                 NETCDF_COMPAT) ) < 0) HERR(routineName);
  /* variables*/
  dims[0] = mpi.nx;
  if (( H5LTmake_dataset(ncid_mpi, XNUM_NAME, 1, dims,
                H5T_NATIVE_INT, mpi.xnum) ) < 0) HERR(routineName);
  if (( id_tmp = H5Dopen2(ncid_mpi, XNUM_NAME,
                          H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(id_tmp, id_x, 0)) < 0) HERR(routineName);
  if (( H5Dclose(id_tmp) ) < 0) HERR(routineName);
  dims[0] = mpi.ny;
  if (( H5LTmake_dataset(ncid_mpi, YNUM_NAME, 1, dims,
                H5T_NATIVE_INT, mpi.ynum) ) < 0) HERR(routineName);
  if (( id_tmp = H5Dopen2(ncid_mpi, YNUM_NAME,
                          H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(id_tmp, id_y, 0)) < 0) HERR(routineName);
  if (( H5Dclose(id_tmp) ) < 0) HERR(routineName);
  dims[0] = mpi.nx;
  dims[1] = mpi.ny;
  if (( file_dspace = H5Screate_simple(2, dims, NULL) ) < 0) HERR(routineName);
  if (( plist = H5Pcreate(H5P_DATASET_CREATE) ) < 0) HERR(routineName);
  if (( H5Pset_fill_value(plist, H5T_NATIVE_FLOAT, &FILLVALUE) ) < 0)
    HERR(routineName);
  if (( H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY) ) < 0) HERR(routineName);
  if (( H5Pset_fill_time(plist, H5D_FILL_TIME_ALLOC) ) < 0) HERR(routineName);
  if (( io.in_mpi_tm = H5Dcreate(ncid_mpi, TASK_MAP, H5T_NATIVE_LONG,
         file_dspace, H5P_DEFAULT, plist, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_tm, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_tm, id_y, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid_mpi, TASK_MAP, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);
  if (( io.in_mpi_tn = H5Dcreate(ncid_mpi, TASK_NUMBER, H5T_NATIVE_LONG,
         file_dspace, H5P_DEFAULT, plist, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_tn, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_tn, id_y, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid_mpi, TASK_NUMBER, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);
  if (( io.in_mpi_it = H5Dcreate(ncid_mpi, ITER_NAME, H5T_NATIVE_LONG,
         file_dspace, H5P_DEFAULT, plist, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_it, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_it, id_y, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid_mpi, ITER_NAME, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);
  if (( io.in_mpi_conv = H5Dcreate(ncid_mpi, CONV_NAME, H5T_NATIVE_LONG,
         file_dspace, H5P_DEFAULT, plist, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_conv, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_conv, id_y, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid_mpi, CONV_NAME, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);
  if (( io.in_mpi_dm = H5Dcreate(ncid_mpi, DM_NAME, H5T_NATIVE_FLOAT,
         file_dspace, H5P_DEFAULT, plist, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_dm, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_dm, id_y, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid_mpi, DM_NAME, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);
  if (( io.in_mpi_zc = H5Dcreate(ncid_mpi, ZC_NAME, H5T_NATIVE_INT,
         file_dspace, H5P_DEFAULT, plist, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_zc, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_zc, id_y, 1)) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid_mpi, ZC_NAME, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);
  if (( H5Pclose(plist) ) < 0) HERR(routineName);
  if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);

  dims[0] = mpi.nx;
  dims[1] = mpi.ny;
  dims[2] = NMaxIter;
  if (( file_dspace = H5Screate_simple(3, dims, NULL) ) < 0) HERR(routineName);
  if (( plist = H5Pcreate(H5P_DATASET_CREATE) ) < 0) HERR(routineName);
  if (( H5Pset_fill_value(plist, H5T_NATIVE_FLOAT, &FILLVALUE) ) < 0)
    HERR(routineName);
  if (( H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY) ) < 0) HERR(routineName);
  if (( H5Pset_fill_time(plist, H5D_FILL_TIME_ALLOC) ) < 0) HERR(routineName);
  if (( io.in_mpi_dmh = H5Dcreate(ncid_mpi, DMH_NAME, H5T_NATIVE_FLOAT,
         file_dspace, H5P_DEFAULT, plist, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( id_tmp = H5Dopen2(ncid_mpi, IT_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_dmh, id_x, 0)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_dmh, id_y, 1)) < 0) HERR(routineName);
  if (( H5DSattach_scale(io.in_mpi_dmh, id_tmp, 2)) < 0) HERR(routineName);
  if (( H5LTset_attribute_float(ncid_mpi, DMH_NAME, "_FillValue",
                                &FILLVALUE, 1) ) < 0) HERR(routineName);
  if (( H5Dclose(id_tmp) ) < 0) HERR(routineName);
  if (( H5Pclose(plist) ) < 0) HERR(routineName);
  if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);

  /* attributes */
  if (( H5LTset_attribute_int(ncid_mpi, ".", "x_start",
                              &input.p15d_x0, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_mpi, ".", "x_end",
                              &input.p15d_x1, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_mpi, ".", "x_step",
                              &input.p15d_xst, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_mpi, ".", "y_start",
                              &input.p15d_y0, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_mpi, ".", "y_end",
                              &input.p15d_y1, 1) ) < 0) HERR(routineName);
  if (( H5LTset_attribute_int(ncid_mpi, ".", "y_step",
                              &input.p15d_yst, 1) ) < 0) HERR(routineName);

  /* Tiago: most of the arrays involving Ntasks or rank as index are not
            currently being written. They should eventually be migrated into
            arrays of [ix, iy] and be written for each task. This is to
            avoid causing problems with pool mode, where these quantities are
            not known from the start.
  */
  if (( H5Dclose(id_x) ) < 0) HERR(routineName);
  if (( H5Dclose(id_y) ) < 0) HERR(routineName);
  if (( H5Dclose(id_z) ) < 0) HERR(routineName);
  /* Flush ensures file is created in case of crash */
  if (( H5Fflush(ncid, H5F_SCOPE_LOCAL) ) < 0) HERR(routineName);
  /* --- Copy stuff to the IO data struct --- */
  io.in_ncid       = ncid;
  io.in_input_ncid = ncid_input;
  io.in_atmos_ncid = ncid_atmos;
  io.in_mpi_ncid   = ncid_mpi;
  return;
}
/* ------- end   --------------------------   init_hdf5_indata.c  --- */

/* ------- begin --------------------------   init_hdf5_indata_old.c  --- */
void init_hdf5_indata_existing(void)
/* Opens an existing input data file, loads structures and ids */
{
  const char routineName[] = "init_hdf5_indata_existing";
  size_t  attr_size;
  hid_t   ncid, plist;
  char   *atmosID;
  int     NMaxIter;
  H5T_class_t type_class;

  /* Open the file with parallel MPI-IO access */
  if (( plist = H5Pcreate(H5P_FILE_ACCESS )) < 0) HERR(routineName);
  if (( H5Pset_fapl_mpio(plist, mpi.comm, mpi.info) ) < 0) HERR(routineName);
  if (( ncid = H5Fopen(INPUTDATA_FILE, H5F_ACC_RDWR, plist) ) < 0)
    HERR(routineName);
  if (( H5Pclose(plist) ) < 0) HERR(routineName);
  io.in_ncid = ncid;
  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if (( H5LTget_attribute_info(ncid, "/", "atmosID", NULL, &type_class,
                               &attr_size) ) < 0) HERR(routineName);
  atmosID = (char *) malloc(attr_size + 1);
  if (( H5LTget_attribute_string(ncid, "/", "atmosID", atmosID) ) < 0)
    HERR(routineName);
  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,
       "Indata file was calculated for different atmosphere (%s) than current",
	     atmosID);
    Error(WARNING, routineName, messageStr);
    }
  free(atmosID);
  /* Check that NMaxIter is enough to cover N_MAX_ITER from keyword.input */
  if (( H5LTget_attribute_int(ncid, "/mpi", "niterations", &NMaxIter) ) < 0)
    HERR(routineName);
  if (input.NmaxIter > NMaxIter) {  /* Raise error to avoid hdf5 writing error */
      sprintf(messageStr,  
              "N_MAX_ITER is %d, while maximum size in file is %d. Reduce N_MAX_ITER or start a new file with larger NMaxIter.\n",
              input.NmaxIter, NMaxIter);
      Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* Get group IDs */
  if (( io.in_input_ncid = H5Gopen(ncid, "input", H5P_DEFAULT) ) < 0)
      HERR(routineName);
  if (( io.in_atmos_ncid = H5Gopen(ncid, "atmos", H5P_DEFAULT) ) < 0)
      HERR(routineName);
  if (( io.in_mpi_ncid = H5Gopen(ncid, "mpi", H5P_DEFAULT) ) < 0)
      HERR(routineName);
  /* --- Open datasets collectively ---*/
  if (( io.in_atmos_T = H5Dopen(io.in_atmos_ncid, TEMP_NAME,
                                H5P_DEFAULT) ) < 0) HERR(routineName);
  if (( io.in_atmos_vz = H5Dopen(io.in_atmos_ncid, VZ_NAME,
                                H5P_DEFAULT) ) < 0) HERR(routineName);
  if (( io.in_atmos_z = H5Dopen(io.in_atmos_ncid, ZH_NAME,
                                H5P_DEFAULT) ) < 0) HERR(routineName);
  if (( io.in_atmos_ne = H5Dopen(io.in_atmos_ncid, NE_NAME,
                                H5P_DEFAULT) ) < 0) HERR(routineName);

  if (( io.in_mpi_tm = H5Dopen(io.in_mpi_ncid, TASK_MAP,
                               H5P_DEFAULT) ) < 0) HERR(routineName);
  if (( io.in_mpi_tn = H5Dopen(io.in_mpi_ncid, TASK_NUMBER,
                               H5P_DEFAULT) ) < 0) HERR(routineName);
  if (( io.in_mpi_it = H5Dopen(io.in_mpi_ncid, ITER_NAME,
                               H5P_DEFAULT) ) < 0) HERR(routineName);
  if (( io.in_mpi_conv = H5Dopen(io.in_mpi_ncid, CONV_NAME,
                               H5P_DEFAULT) ) < 0) HERR(routineName);
  if (( io.in_mpi_dm = H5Dopen(io.in_mpi_ncid, DM_NAME,
                               H5P_DEFAULT) ) < 0) HERR(routineName);
  if (( io.in_mpi_dmh = H5Dopen(io.in_mpi_ncid, DMH_NAME,
                               H5P_DEFAULT) ) < 0) HERR(routineName);
  if (( io.in_mpi_zc = H5Dopen(io.in_mpi_ncid, ZC_NAME,
                               H5P_DEFAULT) ) < 0) HERR(routineName);
  if (input.wavetable != NULL) free(input.wavetable);
  return;
}
/* ------- end   --------------------------   init_hdf5_indata_old.c  --- */


/* ------- begin --------------------------   close_hdf5_indata.c --- */
void close_hdf5_indata(void)
/* Closes the indata file */
{
  const char routineName[] = "close_hdf5_indata";
  int n;
  Atom *atom;

  /* Close all datasets */
  if (( H5Dclose(io.in_atmos_T) ) < 0) HERR(routineName);
  if (( H5Dclose(io.in_atmos_vz) ) < 0) HERR(routineName);
  if (( H5Dclose(io.in_atmos_z) ) < 0) HERR(routineName);
  if (( H5Dclose(io.in_atmos_ne) ) < 0) HERR(routineName);
  if (( H5Dclose(io.in_mpi_tm) ) < 0) HERR(routineName);
  if (( H5Dclose(io.in_mpi_tn) ) < 0) HERR(routineName);
  if (( H5Dclose(io.in_mpi_it) ) < 0) HERR(routineName);
  if (( H5Dclose(io.in_mpi_conv) ) < 0) HERR(routineName);
  if (( H5Dclose(io.in_mpi_dm) ) < 0) HERR(routineName);
  if (( H5Dclose(io.in_mpi_dmh) ) < 0) HERR(routineName);
  if (( H5Dclose(io.in_mpi_zc) ) < 0) HERR(routineName);

  /* Close all groups */
  if (( H5Gclose(io.in_input_ncid) ) < 0) HERR(routineName);
  if (( H5Gclose(io.in_atmos_ncid) ) < 0) HERR(routineName);
  if (( H5Gclose(io.in_mpi_ncid) ) < 0) HERR(routineName);

  /* Close file */
  if (( H5Fclose(io.in_ncid) ) < 0) HERR(routineName);

  /* Free other resources */
  if (input.keyword_file_contents != NULL) free(input.keyword_file_contents);
  if (input.atomic_file_contents != NULL) {  /* For reruns */
      for (n=0; n < atmos.Natom; n++)
          free(input.atomic_file_contents[n]);
      free(input.atomic_file_contents);
  } else {
      for (n=0; n < atmos.Natom; n++) {   /* Outside of reruns */
          atom = &atmos.atoms[n];
          free(atom->fp_input);
      }
  }
  if (input.kurucz_file_contents != NULL)
      free(input.kurucz_file_contents);
  if (input.kurucz_line_file_contents != NULL) {
      for (n = 0;  n < input.Nkurucz_files;  n++) {
          free(input.kurucz_line_file_contents[n]);
      }
      free(input.kurucz_line_file_contents);
  }

  return;
}
/* ------- end   --------------------------   close_hdf5_indata.c --- */

/* ------- begin --------------------------   writeAtmos_p.c --- */
void writeAtmos_p(void)
{
  /* Write atmos arrays. This has now been modified and writes the interpolated
     arrays, from depth_refine. With that, now this is the only viable option
     to write the atmos data, as there is no option to save in memory and
     writeAtmos_all used to write from the input file, not the interpolated
     quantities

     IMPORTANT: at the moment this is a trimmed version, only writing z to save
                space and computational time.

     */
  const char routineName[] = "writeAtmos_p";
  hsize_t  offset[] = {0, 0, 0, 0};
  hsize_t  count[] = {1, 1, 1, 1};
  hsize_t  dims[4];
  hid_t    file_dspace, mem_dspace;

  /* Memory dataspace */
  dims[0] = atmos.Nspace;
  if (( mem_dspace = H5Screate_simple(1, dims, NULL) ) < 0)
    HERR(routineName);
  /* File dataspace */
  offset[0] = mpi.ix;
  offset[1] = mpi.iy;
  offset[2] = mpi.zcut;
  count[2] = atmos.Nspace;
  if (( file_dspace = H5Dget_space(io.in_atmos_T) ) < 0) HERR(routineName);
  if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                            NULL, count, NULL) ) < 0) HERR(routineName);
  if (( H5Dwrite(io.in_atmos_T, H5T_NATIVE_DOUBLE, mem_dspace,
                file_dspace, H5P_DEFAULT, atmos.T) ) < 0) HERR(routineName);
  if (( H5Dwrite(io.in_atmos_vz, H5T_NATIVE_DOUBLE, mem_dspace,
           file_dspace, H5P_DEFAULT, geometry.vel) ) < 0) HERR(routineName);
  if (( H5Dwrite(io.in_atmos_z, H5T_NATIVE_DOUBLE, mem_dspace,
        file_dspace, H5P_DEFAULT, geometry.height) ) < 0) HERR(routineName);
  if (( H5Dwrite(io.in_atmos_ne, H5T_NATIVE_DOUBLE, mem_dspace,
               file_dspace, H5P_DEFAULT, atmos.ne) ) < 0) HERR(routineName);
  /* release dataspace resources */
  if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
  if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
  return;
}
/* ------- end   --------------------------   writeAtmos_p.c --- */

/* ------- begin --------------------------   writeMPI_all.c --- */
void writeMPI_all(void) {
/* Writes output on indata file, MPI group, all tasks at once */
  const char routineName[] = "writeMPI_all";
  int      task;
  hsize_t  offset[] = {0, 0, 0, 0};
  hsize_t  count[] = {1, 1, 1, 1};
  hsize_t  dims[4];
  hid_t    file_dspace, mem_dspace;

  /* Write single values of Ntasks, one value at a time */
  dims[0] = 1;
  if (( mem_dspace = H5Screate_simple(1, dims, NULL) ) < 0)
    HERR(routineName);

  for (task = 0; task < mpi.Ntasks; task++) {
    offset[0] = mpi.taskmap[task + mpi.my_start][0];
    offset[1] = mpi.taskmap[task + mpi.my_start][1];
    if (( file_dspace = H5Dget_space(io.in_mpi_it) ) < 0) HERR(routineName);
    if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                              NULL, count, NULL) ) < 0) HERR(routineName);
    if (( H5Dwrite(io.in_mpi_it, H5T_NATIVE_INT, mem_dspace, file_dspace,
                   H5P_DEFAULT, &mpi.niter[task]) ) < 0) HERR(routineName);
    if (( H5Dwrite(io.in_mpi_conv, H5T_NATIVE_INT, mem_dspace, file_dspace,
                H5P_DEFAULT, &mpi.convergence[task]) ) < 0) HERR(routineName);
    if (( H5Dwrite(io.in_mpi_zc, H5T_NATIVE_INT, mem_dspace, file_dspace,
                  H5P_DEFAULT, &mpi.zcut_hist[task]) ) < 0) HERR(routineName);
    if (( H5Dwrite(io.in_mpi_dm, H5T_NATIVE_DOUBLE, mem_dspace, file_dspace,
                   H5P_DEFAULT, &mpi.dpopsmax[task]) ) < 0) HERR(routineName);
    if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
  }
  if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);

  /* Write array with multiple values */
  for (task = 0; task < mpi.Ntasks; task++) {
    dims[0] = mpi.niter[task];
    if (( mem_dspace = H5Screate_simple(1, dims, NULL) ) < 0)
      HERR(routineName);
    offset[0] = mpi.taskmap[task + mpi.my_start][0];
    offset[1] = mpi.taskmap[task + mpi.my_start][1];
    count[2] = mpi.niter[task];
    if (( file_dspace = H5Dget_space(io.in_mpi_dmh) ) < 0) HERR(routineName);
    if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                              NULL, count, NULL) ) < 0) HERR(routineName);
    if (( H5Dwrite(io.in_mpi_dmh, H5T_NATIVE_DOUBLE, mem_dspace, file_dspace,
               H5P_DEFAULT, mpi.dpopsmax_hist[task]) ) < 0) HERR(routineName);
    if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
    if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
  }
  return;
}
/* ------- end   --------------------------   writeMPI_all.c --- */


/* ------- begin --------------------------   writeMPI_p.c ----- */
void writeMPI_p(int task) {
/* Writes output on indata file, MPI group, one task at once */
  const char routineName[] = "writeMPI_p";
  hsize_t  offset[] = {0, 0, 0, 0};
  hsize_t  count[] = {1, 1, 1, 1};
  hsize_t  dims[4];
  hid_t    file_dspace, mem_dspace;

  dims[0] = 1;
  if (( mem_dspace = H5Screate_simple(1, dims, NULL) ) < 0) HERR(routineName);
  offset[0] = mpi.ix;
  offset[1] = mpi.iy;
  if (( file_dspace = H5Dget_space(io.in_mpi_tm) ) < 0) HERR(routineName);
  if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                            NULL, count, NULL) ) < 0) HERR(routineName);
  if (( H5Dwrite(io.in_mpi_tm, H5T_NATIVE_INT, mem_dspace, file_dspace,
                   H5P_DEFAULT, &mpi.rank) ) < 0) HERR(routineName);
  if (( H5Dwrite(io.in_mpi_tn, H5T_NATIVE_INT, mem_dspace, file_dspace,
                   H5P_DEFAULT, &task) ) < 0) HERR(routineName);
  if (( H5Dwrite(io.in_mpi_it, H5T_NATIVE_INT, mem_dspace, file_dspace,
                   H5P_DEFAULT, &mpi.niter[0]) ) < 0) HERR(routineName);
  if (( H5Dwrite(io.in_mpi_conv, H5T_NATIVE_INT, mem_dspace, file_dspace,
                   H5P_DEFAULT, &mpi.convergence[0]) ) < 0) HERR(routineName);
  if (( H5Dwrite(io.in_mpi_zc, H5T_NATIVE_INT, mem_dspace, file_dspace,
                   H5P_DEFAULT, &mpi.zcut_hist[0]) ) < 0) HERR(routineName);
  if (( H5Dwrite(io.in_mpi_dm, H5T_NATIVE_DOUBLE, mem_dspace, file_dspace,
                   H5P_DEFAULT, &mpi.dpopsmax[0]) ) < 0) HERR(routineName);
  if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
  if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);

  dims[0] = mpi.niter[0];
  if (( mem_dspace = H5Screate_simple(1, dims, NULL) ) < 0) HERR(routineName);
  offset[0] = mpi.ix;
  offset[1] = mpi.iy;
  count[2] = mpi.niter[0];
  if (( file_dspace = H5Dget_space(io.in_mpi_dmh) ) < 0) HERR(routineName);
  if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                            NULL, count, NULL) ) < 0) HERR(routineName);
  if (( H5Dwrite(io.in_mpi_dmh, H5T_NATIVE_DOUBLE, mem_dspace, file_dspace,
                 H5P_DEFAULT, mpi.dpopsmax_hist[0]) ) < 0) HERR(routineName);
  if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
  if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
  return;
}
/* ------- end   --------------------------   writeMPI_p.c ------- */


void readConvergence(void) {
  /* This is a self-contained function to read the convergence matrix,
     written by RH. */
  const char routineName[] = "readConvergence";
  char *atmosID;
  int nx, ny;
  size_t attr_size;
  hid_t ncid, ncid_mpi, plist;
  H5T_class_t type_class;

  mpi.rh_converged = matrix_int(mpi.nx, mpi.ny);

  /* --- Open the inputdata file --- */
  if (( plist = H5Pcreate(H5P_FILE_ACCESS )) < 0) HERR(routineName);
  if (( H5Pset_fapl_mpio(plist, mpi.comm, mpi.info) ) < 0) HERR(routineName);
  if (( ncid = H5Fopen(INPUTDATA_FILE, H5F_ACC_RDONLY, plist) ) < 0)
    HERR(routineName);
  if (( H5Pclose(plist) ) < 0) HERR(routineName);
  /* Get ncid of the MPI group */
  if (( ncid_mpi = H5Gopen(ncid, "mpi", H5P_DEFAULT) ) < 0) HERR(routineName);

  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if (( H5LTget_attribute_info(ncid, "/", "atmosID", NULL, &type_class,
                               &attr_size) ) < 0) HERR(routineName);
  atmosID = (char *) malloc(attr_size + 1);
  if (( H5LTget_attribute_string(ncid, "/", "atmosID", atmosID) ) < 0)
    HERR(routineName);
  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr, "Indata file was calculated for different "
            "atmosphere (%s) than current (%s)", atmosID, atmos.ID);
    Error(WARNING, routineName, messageStr);
    }
  free(atmosID);
  /* Check that dimension sizes match */
  if (( H5LTget_attribute_int(ncid, "/", "nx", &nx) ) < 0) HERR(routineName);
  if (nx != mpi.nx) {
    sprintf(messageStr,
	    "Number of x points mismatch: expected %d, found %d.",
	    mpi.nx, (int)nx);
    Error(WARNING, routineName, messageStr);
  }
  if (( H5LTget_attribute_int(ncid, "/", "ny", &ny) ) < 0) HERR(routineName);
  if (ny != mpi.ny) {
    sprintf(messageStr,
	    "Number of y points mismatch: expected %d, found %d.",
	    mpi.ny, (int)ny);
    Error(WARNING, routineName, messageStr);
  }
  /* --- Read variable --- */
  if (( H5LTread_dataset_int(ncid_mpi, CONV_NAME,
                             mpi.rh_converged[0]) ) < 0) HERR(routineName);
  /* --- Close inputdata file --- */
  if (( H5Gclose(ncid_mpi) ) < 0) HERR(routineName);
  if (( H5Fclose(ncid) ) < 0) HERR(routineName);
  return;
}


void readSavedKeywords(void) {
  /* Reads keyword.input saved into output_indata file */
  const char routineName[] = "readSavedKeywords";
  char *atmosID;
  size_t attr_size, str_size = 0;
  hid_t ncid, ncid_input, plist;
  H5T_class_t type_class;
  bool_t saved_p15d_refine, saved_p15d_zcut, saved_accelerate_mols;
  int saved_NpescIter, saved_Ngdelay, saved_Ngorder, saved_Ngperiod;
  int saved_NmaxScatter, saved_NmaxIter, saved_PRD_NmaxIter;
  int saved_PRD_Ngdelay, saved_PRD_Ngorder, saved_PRD_Ngperiod;
  enum S_interpol saved_S_interpolation;
  enum S_interpol_stokes saved_S_interpolation_stokes;
  double saved_crsw, saved_crsw_ini, saved_prdswitch, saved_prdsw;
  double saved_p15d_tmax, saved_iterLimit, saved_PRDiterLimit;

  /* --- Open the inputdata file --- */
  if (( plist = H5Pcreate(H5P_FILE_ACCESS )) < 0) HERR(routineName);
  if (( H5Pset_fapl_mpio(plist, mpi.comm, mpi.info) ) < 0) HERR(routineName);
  if (( ncid = H5Fopen(INPUTDATA_FILE, H5F_ACC_RDONLY, plist) ) < 0)
    HERR(routineName);
  if (( H5Pclose(plist) ) < 0) HERR(routineName);

  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if (( H5LTget_attribute_info(ncid, "/", "atmosID", NULL, &type_class,
                               &attr_size) ) < 0) HERR(routineName);
  atmosID = (char *) malloc(attr_size + 1);
  if (( H5LTget_attribute_string(ncid, "/", "atmosID", atmosID) ) < 0)
    HERR(routineName);
  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,  "Indata file was calculated for different "
           "atmosphere (%s) than current (%s)", atmosID, atmos.ID);
    Error(WARNING, routineName, messageStr);
    }
  free(atmosID);

  /* --- Read input files --- */
  if (( ncid_input = H5Gopen(ncid, "input", H5P_DEFAULT) ) < 0) HERR(routineName);
  if (H5LTfind_dataset(ncid_input, "keyword_file_contents")) {
    /* For H5T_STRING datasets, size of string is saved under last argument */
    if ((H5LTget_dataset_info(ncid_input, "keyword_file_contents", NULL,
			                  NULL, &str_size)) < 0) HERR(routineName);
  } else {
    sprintf(messageStr, "Could not read keyword.input file in indata file, "
	                    "no rerun is possible. Aborting.\n");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  input.keyword_file_contents = (char *) malloc(str_size + 1);
  if (( H5LTread_dataset_string(ncid_input, "keyword_file_contents",
                         input.keyword_file_contents) ) < 0) HERR(routineName);
  /* --- Close inputdata file --- */
  if (( H5Gclose(ncid_input) ) < 0) HERR(routineName);
  if (( H5Fclose(ncid) ) < 0) HERR(routineName);

  /* Save keywords that can change in a rerun */
  saved_p15d_zcut = input.p15d_zcut;
  saved_p15d_tmax = input.p15d_tmax;
  saved_p15d_refine = input.p15d_refine;
  saved_crsw = input.crsw;
  saved_crsw_ini = input.crsw_ini;
  saved_prdsw = input.prdsw;
  saved_prdswitch = input.prdswitch;
  saved_NpescIter = input.NpescIter;
  saved_p15d_tmax = input.p15d_tmax;
  saved_NmaxScatter = input.NmaxScatter;
  saved_NmaxIter = input.NmaxIter;
  saved_iterLimit = input.iterLimit;
  saved_Ngdelay = input.Ngdelay;
  saved_Ngorder = input.Ngorder;
  saved_Ngperiod = input.Ngperiod;
  saved_accelerate_mols = input.accelerate_mols;
  saved_PRD_NmaxIter = input.PRD_NmaxIter;
  saved_PRDiterLimit = input.PRDiterLimit;
  saved_PRD_Ngdelay = input.PRD_Ngdelay;
  saved_PRD_Ngorder = input.PRD_Ngorder;
  saved_PRD_Ngperiod = input.PRD_Ngperiod;
  saved_S_interpolation = input.S_interpolation;
  saved_S_interpolation_stokes = input.S_interpolation_stokes;
  /* Overwrite non-changeable keyword options with saved ones */
  readInput(input.keyword_file_contents);
  /* Put back changeable keyword options */
  input.p15d_rerun = TRUE;  /* This is only called for reruns */
  input.p15d_zcut = saved_p15d_zcut;
  input.p15d_tmax = saved_p15d_tmax;
  input.p15d_refine = saved_p15d_refine;
  input.crsw = saved_crsw;
  input.crsw_ini = saved_crsw_ini;
  input.prdsw = saved_prdsw;
  input.prdswitch = saved_prdswitch;
  input.NpescIter = saved_NpescIter;
  input.p15d_tmax = saved_p15d_tmax;
  input.NmaxScatter = saved_NmaxScatter;
  input.NmaxIter = saved_NmaxIter;
  input.iterLimit = saved_iterLimit;
  input.Ngdelay = saved_Ngdelay;
  input.Ngorder = saved_Ngorder;
  input.Ngperiod = saved_Ngperiod;
  input.accelerate_mols = saved_accelerate_mols;
  input.PRD_NmaxIter = saved_PRD_NmaxIter;
  input.PRDiterLimit = saved_PRDiterLimit;
  input.PRD_Ngdelay = saved_PRD_Ngdelay;
  input.PRD_Ngorder = saved_PRD_Ngorder;
  input.PRD_Ngperiod = saved_PRD_Ngperiod;
  input.S_interpolation = saved_S_interpolation;
  input.S_interpolation_stokes = saved_S_interpolation_stokes;
  return;
}


void readSavedInput(void) {
  /* Reads saved input and convergence info for rerun */
  const char routineName[] = "readSavedInput";
  char *atmosID, **atom_names, group_name[ARR_STRLEN];
  int n, atom_name_size;
  size_t attr_size, str_size = 0;
  hid_t ncid, ncid_input, ncid_tmp, plist;
  hsize_t dims[5];
  H5T_class_t type_class;


  /* --- Open the inputdata file --- */
  if (( plist = H5Pcreate(H5P_FILE_ACCESS )) < 0) HERR(routineName);
  if (( H5Pset_fapl_mpio(plist, mpi.comm, mpi.info) ) < 0) HERR(routineName);
  if (( ncid = H5Fopen(INPUTDATA_FILE, H5F_ACC_RDONLY, plist) ) < 0)
    HERR(routineName);
  if (( H5Pclose(plist) ) < 0) HERR(routineName);

  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if (( H5LTget_attribute_info(ncid, "/", "atmosID", NULL, &type_class,
                               &attr_size) ) < 0) HERR(routineName);
  atmosID = (char *) malloc(attr_size + 1);
  if (( H5LTget_attribute_string(ncid, "/", "atmosID", atmosID) ) < 0)
    HERR(routineName);
  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,  "Indata file was calculated for different "
           "atmosphere (%s) than current (%s)", atmosID, atmos.ID);
    Error(WARNING, routineName, messageStr);
    }
  free(atmosID);

  /* --- Read input files --- */
  if (( ncid_input = H5Gopen(ncid, "input", H5P_DEFAULT) ) < 0) HERR(routineName);
  if (H5LTfind_dataset(ncid_input, "atoms_file_contents")) {
    /* For H5T_STRING datasets, size of string is saved under last argument */
    if ((H5LTget_dataset_info(ncid_input, "atoms_file_contents", NULL,
			                  NULL, &str_size)) < 0) HERR(routineName);
  } else {
    sprintf(messageStr, "Could not read atoms file in indata file, "
	                    "no rerun is possible. Aborting.\n");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  input.atoms_file_contents = (char *) malloc(str_size + 1);
  if (( H5LTread_dataset_string(ncid_input, "atoms_file_contents",
                          input.atoms_file_contents) ) < 0) HERR(routineName);
  if (( H5LTget_attribute_int(ncid_input, ".", "natoms",
                              &input.Natoms) ) < 0) HERR(routineName);
  /* Load atom files */
  if ((H5LTget_dataset_info(ncid_input, "atom_groups", dims, NULL, NULL)) < 0)
    HERR(routineName);
  input.Natoms = dims[0];
  atom_name_size = dims[1];
  atom_names = matrix_char(input.Natoms, atom_name_size);
  if (( H5LTread_dataset(ncid_input, "atom_groups", H5T_C_S1,
                         atom_names[0]) ) < 0) HERR(routineName);
  input.atomic_file_contents = (char **) malloc(input.Natoms * sizeof(char *));
  for (n = 0;  n < input.Natoms;  n++) {
      if (( ncid_tmp = H5Gopen(ncid_input, atom_names[n], H5P_DEFAULT) ) < 0)
        HERR(routineName);
      if (H5LTfind_dataset(ncid_tmp, "file_contents")) {
          /* H5T_STRING datasets: size of string is saved on last argument */
          if ((H5LTget_dataset_info(ncid_tmp, "file_contents", NULL,
                                    NULL, &str_size)) < 0) HERR(routineName);
      } else {
          sprintf(messageStr, "Could not read atomic file in indata file, "
                              "no rerun is possible. Aborting.\n");
          Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      input.atomic_file_contents[n] = (char *) malloc(str_size + 1);
      if (( H5LTread_dataset_string(ncid_tmp, "file_contents",
                       input.atomic_file_contents[n]) ) < 0) HERR(routineName);
      if (( H5Gclose(ncid_tmp) ) < 0) HERR(routineName);
  }
  freeMatrix((void **) atom_names);

  /* Load data from ray.input */
  if (( H5LTget_attribute_double(ncid_input, ".", INPUT_MU,
                                 &io.ray_muz) ) < 0) HERR(routineName);
  if (( H5LTget_attribute_uint(ncid_input, ".", WAVE_SEL,
                  (unsigned int *) &io.ray_nwave_sel) ) < 0) HERR(routineName);
  io.ray_wave_idx = NULL;
  if (io.ray_nwave_sel > 0) {
      io.ray_wave_idx = (int *) malloc(io.ray_nwave_sel * sizeof(int));
      if ((H5LTread_dataset_int(ncid_input, WAVE_SEL_IDX, io.ray_wave_idx)) < 0)
          HERR(routineName);
  }
  /* --- Save geometry values to change back after --    ------------ */
  geometry.save_Nrays = atmos.Nrays;
  geometry.save_wmu = geometry.wmu[0];
  geometry.save_muz = geometry.muz[0];
  geometry.save_mux = geometry.mux[0];
  geometry.save_muy = geometry.muy[0];
  /* Data from wavetable */
  if (H5LTfind_dataset(ncid_input, WAVETABLE)) {
      if (( H5LTget_attribute_uint(ncid_input, ".", NXWAVE,
                      (unsigned int *) &input.Nxwave) ) < 0) HERR(routineName);
      input.wavetable = (double *) malloc(input.Nxwave * sizeof(double));
      if ((H5LTread_dataset_double(ncid_input, WAVETABLE,
                                   input.wavetable)) < 0) HERR(routineName);
  }
  /* Data from Kurucz line file, if existing */
  if (( H5LTget_attribute_uint(ncid_input, ".", NKURUCZ,
     (unsigned int *) &input.Nkurucz_files) ) < 0) HERR(routineName);
  if (input.Nkurucz_files > 0) {
      if (H5LTfind_dataset(ncid_input, "kurucz_file_contents")) {
          if ((H5LTget_dataset_info(ncid_input, "kurucz_file_contents", NULL,
                                    NULL, &str_size)) < 0) HERR(routineName);
      } else {
          sprintf(messageStr, "Could not read Kurucz file contents in indata "
                              "file, no rerun is possible. Aborting.\n");
          Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      input.kurucz_file_contents = (char *) malloc(str_size + 1);
      if (( H5LTread_dataset_string(ncid_input, "kurucz_file_contents",
                       input.kurucz_file_contents) ) < 0) HERR(routineName);
      input.kurucz_line_file_contents = (char **) malloc(input.Nkurucz_files *
                                                         sizeof(char *));
      for (n = 0;  n < input.Nkurucz_files;  n++) {
          sprintf(group_name, KURUCZ_LINE_FILE, n);
          if (( ncid_tmp = H5Gopen(ncid_input, group_name, H5P_DEFAULT) ) < 0)
              HERR(routineName);
          if (H5LTfind_dataset(ncid_tmp, "file_contents")) {
              if ((H5LTget_dataset_info(ncid_tmp, "file_contents", NULL,
                                     NULL, &str_size)) < 0) HERR(routineName);
          } else {
              sprintf(messageStr, "Could not read Kurucz line file in indata "
                                  "file, no rerun is possible. Aborting.\n");
              Error(ERROR_LEVEL_2, routineName, messageStr);
          }
          input.kurucz_line_file_contents[n] = (char *) malloc(str_size + 1);
          if (( H5LTread_dataset_string(ncid_tmp, "file_contents",
                  input.kurucz_line_file_contents[n]) ) < 0) HERR(routineName);
          if (( H5Gclose(ncid_tmp) ) < 0) HERR(routineName);
      }
  }

  /* --- Close group and file --- */
  if (( H5Gclose(ncid_input) ) < 0) HERR(routineName);
  if (( H5Fclose(ncid) ) < 0) HERR(routineName);
  return;
}
