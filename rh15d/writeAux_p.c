/* ------- file: -------------------------- writeAux_p.c ---------

       Version:       rh2.0, 1.5-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Mon Jan 03 14:28:25 2011 --

       --------------------------                      -----------RH-- */

/* --- Writes auxiliary data to output file            --------------- */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"
#include "io.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];
extern Input_Atmos_file infile;
extern MPI_data  mpi;
extern IO_data   io;
extern IO_buffer iobuf;

/* ------- begin --------------------------   init_hdf5_aux.c     --- */
void init_hdf5_aux(void) {
  /* Wrapper to find out if we should use old file or create new one */
  int     nact;
  bool_t  use_old = FALSE;
  Atom   *atom;


  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];

    if (atom->initial_solution == OLD_POPULATIONS) {
      use_old = TRUE;
      break;
    }
  }

  if ((use_old) || (input.p15d_rerun)) init_aux_existing(); else init_aux_new();

  return;
}
/* ------- end --------------------------   init_hdf5_aux.c     --- */

/* ------- begin --------------------------   init_aux_new.c --   --- */
void init_aux_new(void) {
  /* Creates the HDF5 file for the auxiliary data */
  const char routineName[] = "init_aux_new";
  unsigned int *tmp;
  double   *tmp_double;
  int       i;
  hid_t     plist, ncid, file_dspace, ncid_atom, ncid_mol;
  hid_t     id_x, id_y, id_z, id_n, id_tmp;
  hsize_t   dims[5];
  char      group_name[ARR_STRLEN];
  Atom     *atom;
  Molecule *molecule;

  /* Create the file  */
  if (( plist = H5Pcreate(H5P_FILE_ACCESS) ) < 0) HERR(routineName);
  if (( H5Pset_fapl_mpio(plist, mpi.comm, mpi.info) ) < 0) HERR(routineName);
  if (( ncid = H5Fcreate(AUX_FILE, H5F_ACC_TRUNC, H5P_DEFAULT,
                         plist) ) < 0) HERR(routineName);
  if (( H5Pclose(plist) ) < 0) HERR(routineName);


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

  /* Create arrays for multiple-atom/molecule output */
  io.aux_atom_ncid   = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
  io.aux_mol_ncid    = (hid_t *) malloc(atmos.Nactivemol  * sizeof(hid_t));
  if (input.p15d_wpop) {
      io.aux_atom_pop    = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
      io.aux_atom_poplte = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
      io.aux_mol_pop     = (hid_t *) malloc(atmos.Nactivemol  * sizeof(hid_t));
      io.aux_mol_poplte  = (hid_t *) malloc(atmos.Nactivemol  * sizeof(hid_t));
  }
  if (input.p15d_wrates) {
      io.aux_atom_RijL   = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
      io.aux_atom_RjiL   = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
      io.aux_atom_RijC   = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
      io.aux_atom_RjiC   = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
  }
  if (input.p15d_wcrates) {
      io.aux_atom_Cij    = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
  }


  /* Fill value */
  if (( plist = H5Pcreate(H5P_DATASET_CREATE) ) < 0) HERR(routineName);
  if (( H5Pset_fill_value(plist, H5T_NATIVE_FLOAT, &FILLVALUE) ) < 0)
      HERR(routineName);
  if (( H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY) ) < 0) HERR(routineName);
  if (( H5Pset_fill_time(plist, H5D_FILL_TIME_ALLOC) ) < 0) HERR(routineName);

  /* --- Group loop over active ATOMS --- */
  for (i=0; i < atmos.Nactiveatom; i++) {
    atom = atmos.activeatoms[i];
    /* Get group name */
    sprintf(group_name, (atom->ID[1] == ' ') ? "atom_%.1s" : "atom_%.2s",
            atom->ID);
    if (( ncid_atom = H5Gcreate(ncid, group_name, H5P_DEFAULT, H5P_DEFAULT,
                                H5P_DEFAULT) ) < 0) HERR(routineName);
    io.aux_atom_ncid[i] = ncid_atom;
    /* --- dimensions as attributes --- */
    if (( H5LTset_attribute_int(ncid_atom, ".", "nlevel",
                                &atom->Nlevel, 1)) < 0) HERR(routineName);
    if (( H5LTset_attribute_int(ncid_atom, ".", "nline",
                                &atom->Nline, 1)) < 0) HERR(routineName);
    if (( H5LTset_attribute_int(ncid_atom, ".", "ncontinuum",
                                &atom->Ncont, 1)) < 0) HERR(routineName);
    /* --- dimension datasets --- */
    dims[0] = mpi.nx;
    if (( H5LTmake_dataset(ncid_atom, X_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                           geometry.xscale) ) < 0)  HERR(routineName);
    if (( id_x = H5Dopen2(ncid_atom, X_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
    dims[0] = mpi.ny;
    if (( H5LTmake_dataset(ncid_atom, Y_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                           geometry.yscale) ) < 0)  HERR(routineName);
    if (( id_y = H5Dopen2(ncid_atom, Y_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
    dims[0] = infile.nz;
    tmp_double = (double *) calloc(infile.nz , sizeof(double));
    if (( H5LTmake_dataset(ncid_atom, ZOUT_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                           tmp_double) ) < 0)  HERR(routineName);
    free(tmp_double);
    if (( id_z = H5Dopen2(ncid_atom, ZOUT_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
    dims[0] = atom->Nlevel;
    tmp = (unsigned int *) calloc(atom->Nlevel , sizeof(unsigned int));
    if (( H5LTmake_dataset(ncid_atom, LEVEL_NAME, 1, dims, H5T_NATIVE_UINT,
                           tmp) ) < 0)  HERR(routineName);
    free(tmp);
    dims[0] = atom->Nline;
    tmp = (unsigned int *) calloc(atom->Nline , sizeof(unsigned int));
    if (( H5LTmake_dataset(ncid_atom, LINE_NAME, 1, dims, H5T_NATIVE_UINT,
                           tmp) ) < 0)  HERR(routineName);
    free(tmp);
    if (atom->Ncont > 0) {
        dims[0] = atom->Ncont;
        tmp = (unsigned int *) calloc(atom->Ncont , sizeof(unsigned int));
        if (( H5LTmake_dataset(ncid_atom, CONT_NAME, 1, dims, H5T_NATIVE_UINT,
                               tmp) ) < 0)  HERR(routineName);
        free(tmp);
    }
    /* For compatibility with netCDF readers, only use dataset as dimension */
    if (( H5LTset_attribute_string(ncid_atom, ZOUT_NAME, "NAME",
                                   NETCDF_COMPAT) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid_atom, LEVEL_NAME, "NAME",
                                   NETCDF_COMPAT) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid_atom, LINE_NAME, "NAME",
                                   NETCDF_COMPAT) ) < 0) HERR(routineName);
    if (atom->Ncont > 0) {
        if (( H5LTset_attribute_string(ncid_atom, CONT_NAME, "NAME",
                                       NETCDF_COMPAT) ) < 0) HERR(routineName);
    }
    /* --- variables --- */
    dims[0] = atom->Nlevel;
    dims[1] = mpi.nx;
    dims[2] = mpi.ny;
    dims[3] = infile.nz;
    /* Populations */
    if (input.p15d_wpop) {
      if (( file_dspace = H5Screate_simple(4, dims, NULL) ) < 0)
        HERR(routineName);
      if (( id_n = H5Dopen2(ncid_atom, LEVEL_NAME,
                            H5P_DEFAULT)) < 0) HERR(routineName);
      if (atom->n != NULL) {
        if (( id_tmp = H5Dcreate(ncid_atom, POP_NAME, H5T_NATIVE_FLOAT,
                                 file_dspace, H5P_DEFAULT,
                                 plist, H5P_DEFAULT)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_x, 1)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_y, 2)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_z, 3)) < 0) HERR(routineName);
        if (( H5LTset_attribute_float(ncid_atom, POP_NAME, "_FillValue",
                                      &FILLVALUE, 1) ) < 0) HERR(routineName);
        io.aux_atom_pop[i] = id_tmp;
      }
      if (atom->nstar != NULL) {
        if (( id_tmp = H5Dcreate(ncid_atom, POPLTE_NAME, H5T_NATIVE_FLOAT,
                                 file_dspace, H5P_DEFAULT,
                                 plist, H5P_DEFAULT)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_x, 1)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_y, 2)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_z, 3)) < 0) HERR(routineName);
        if (( H5LTset_attribute_float(ncid_atom, POPLTE_NAME, "_FillValue",
                                      &FILLVALUE, 1) ) < 0) HERR(routineName);
        io.aux_atom_poplte[i] = id_tmp;
      }
      if (( H5Dclose(id_n) ) < 0) HERR(routineName);
      if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
    }
    if (input.p15d_wrates) {
      /* Radiative rates */
      dims[0] = atom->Nline;
      dims[1] = mpi.nx;
      dims[2] = mpi.ny;
      dims[3] = infile.nz;
      if (( file_dspace = H5Screate_simple(4, dims, NULL) ) < 0) HERR(routineName);
      if (( id_n = H5Dopen2(ncid_atom, LINE_NAME,
                            H5P_DEFAULT)) < 0) HERR(routineName);
      if (( id_tmp = H5Dcreate(ncid_atom, RIJ_L_NAME, H5T_NATIVE_FLOAT,
                               file_dspace, H5P_DEFAULT, plist,
                               H5P_DEFAULT)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_x, 1)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_y, 2)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_z, 3)) < 0) HERR(routineName);
      if (( H5LTset_attribute_float(ncid_atom, RIJ_L_NAME, "_FillValue",
                                    &FILLVALUE, 1) ) < 0) HERR(routineName);
      io.aux_atom_RijL[i] = id_tmp;
      if (( id_tmp = H5Dcreate(ncid_atom, RJI_L_NAME, H5T_NATIVE_FLOAT,
                               file_dspace, H5P_DEFAULT, plist,
                               H5P_DEFAULT)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_x, 1)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_y, 2)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_z, 3)) < 0) HERR(routineName);
      if (( H5LTset_attribute_float(ncid_atom, RJI_L_NAME, "_FillValue",
                                    &FILLVALUE, 1) ) < 0) HERR(routineName);
      io.aux_atom_RjiL[i] = id_tmp;
      if (( H5Dclose(id_n) ) < 0) HERR(routineName);
      if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
      if (atom->Ncont > 0) {
          dims[0] = atom->Ncont;
          if (( file_dspace = H5Screate_simple(4, dims, NULL) ) < 0)
            HERR(routineName);
          if (( id_n = H5Dopen2(ncid_atom, CONT_NAME,
                                H5P_DEFAULT)) < 0) HERR(routineName);
          if (( id_tmp = H5Dcreate(ncid_atom, RIJ_C_NAME, H5T_NATIVE_FLOAT,
                                   file_dspace, H5P_DEFAULT,  plist,
                                   H5P_DEFAULT)) < 0) HERR(routineName);
          if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
          if (( H5DSattach_scale(id_tmp, id_x, 1)) < 0) HERR(routineName);
          if (( H5DSattach_scale(id_tmp, id_y, 2)) < 0) HERR(routineName);
          if (( H5DSattach_scale(id_tmp, id_z, 3)) < 0) HERR(routineName);
          if (( H5LTset_attribute_float(ncid_atom, RIJ_C_NAME, "_FillValue",
                                        &FILLVALUE, 1) ) < 0) HERR(routineName);
          io.aux_atom_RijC[i] = id_tmp;
          if (( id_tmp = H5Dcreate(ncid_atom, RJI_C_NAME, H5T_NATIVE_FLOAT,
                                   file_dspace, H5P_DEFAULT, plist,
                                   H5P_DEFAULT)) < 0) HERR(routineName);
          if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
          if (( H5DSattach_scale(id_tmp, id_x, 1)) < 0) HERR(routineName);
          if (( H5DSattach_scale(id_tmp, id_y, 2)) < 0) HERR(routineName);
          if (( H5DSattach_scale(id_tmp, id_z, 3)) < 0) HERR(routineName);
          if (( H5LTset_attribute_float(ncid_atom, RJI_C_NAME, "_FillValue",
                                        &FILLVALUE, 1) ) < 0) HERR(routineName);
          io.aux_atom_RjiC[i] = id_tmp;
          if (( H5Dclose(id_n) ) < 0) HERR(routineName);
          if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
      }
    }
    if (input.p15d_wcrates) {  
      /* Collisional rates */
      dims[0] = atom->Nlevel;
      dims[1] = atom->Nlevel;
      dims[2] = mpi.nx;
      dims[3] = mpi.ny;
      dims[4] = infile.nz;
      if (( file_dspace = H5Screate_simple(5, dims, NULL) ) < 0) HERR(routineName);
      if (( id_n = H5Dopen2(ncid_atom, LEVEL_NAME,
                            H5P_DEFAULT)) < 0) HERR(routineName);
      if (( id_tmp = H5Dcreate(ncid_atom, COLL_NAME, H5T_NATIVE_FLOAT,
                               file_dspace, H5P_DEFAULT, plist,
                               H5P_DEFAULT)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_n, 1)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_x, 2)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_y, 3)) < 0) HERR(routineName);
      if (( H5DSattach_scale(id_tmp, id_z, 4)) < 0) HERR(routineName);
      if (( H5LTset_attribute_float(ncid_atom, COLL_NAME, "_FillValue",
                                    &FILLVALUE, 1) ) < 0) HERR(routineName);
      io.aux_atom_Cij[i] = id_tmp;
      if (( H5Dclose(id_n) ) < 0) HERR(routineName);
      if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
    }
  if (( H5Dclose(id_x) ) < 0) HERR(routineName);
  if (( H5Dclose(id_y) ) < 0) HERR(routineName);
  if (( H5Dclose(id_z) ) < 0) HERR(routineName);
  }   /* end active ATOMS loop */

  /* --- Group loop over active MOLECULES --- */
  for (i=0; i < atmos.Nactivemol; i++) {
    molecule = atmos.activemols[i];
    /* Get group name */
    sprintf( group_name, "molecule_%s", molecule->ID);
    if (( ncid_mol = H5Gcreate(ncid, group_name, H5P_DEFAULT, H5P_DEFAULT,
                                H5P_DEFAULT) ) < 0) HERR(routineName);
    io.aux_mol_ncid[i] = ncid_mol;
    /* --- dimensions as attributes --- */
    if (( H5LTset_attribute_int(ncid_mol, ".", "nlevel_vibr",
                                &molecule->Nv, 1)) < 0) HERR(routineName);
    if (( H5LTset_attribute_int(ncid_mol, ".", "nline_molecule",
                                &molecule->Nrt, 1)) < 0) HERR(routineName);
    if (( H5LTset_attribute_int(ncid_mol, ".", "nJ",
                                &molecule->NJ, 1)) < 0) HERR(routineName);
    /* --- dimension datasets --- */
    dims[0] = mpi.nx;
    if (( H5LTmake_dataset(ncid_mol, X_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                           geometry.xscale) ) < 0)  HERR(routineName);
    if (( id_x = H5Dopen2(ncid_mol, X_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
    dims[0] = mpi.ny;
    if (( H5LTmake_dataset(ncid_mol, Y_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                           geometry.yscale) ) < 0)  HERR(routineName);
    if (( id_y = H5Dopen2(ncid_mol, Y_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
    dims[0] = infile.nz;
    tmp_double = (double *) calloc(infile.nz , sizeof(double));
    if (( H5LTmake_dataset(ncid_mol, ZOUT_NAME, 1, dims, H5T_NATIVE_DOUBLE,
                           tmp_double) ) < 0)  HERR(routineName);
    free(tmp_double);
    if (( id_z = H5Dopen2(ncid_mol, ZOUT_NAME, H5P_DEFAULT)) < 0) HERR(routineName);
    dims[0] = molecule->Nv;
    tmp = (unsigned int *) calloc(molecule->Nv, sizeof(unsigned int));
    if (( H5LTmake_dataset(ncid_mol, VLEVEL_NAME, 1, dims, H5T_NATIVE_UINT,
                           tmp) ) < 0)  HERR(routineName);
    free(tmp);
    dims[0] = molecule->Nrt;
    tmp = (unsigned int *) calloc(molecule->Nrt, sizeof(unsigned int));
    if (( H5LTmake_dataset(ncid_mol, VLINE_NAME, 1, dims, H5T_NATIVE_UINT,
                           tmp) ) < 0)  HERR(routineName);
    free(tmp);
    dims[0] = molecule->NJ;
    tmp = (unsigned int *) calloc(molecule->NJ, sizeof(unsigned int));
    if (( H5LTmake_dataset(ncid_mol, NJ_NAME, 1, dims, H5T_NATIVE_UINT,
                           tmp) ) < 0)  HERR(routineName);
    free(tmp);
    /* For compatibility with netCDF readers, only use dataset as dimension */
    if (( H5LTset_attribute_string(ncid_mol, ZOUT_NAME, "NAME",
                                   NETCDF_COMPAT) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid_mol, VLEVEL_NAME, "NAME",
                                   NETCDF_COMPAT) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid_mol, VLINE_NAME, "NAME",
                                   NETCDF_COMPAT) ) < 0) HERR(routineName);
    if (( H5LTset_attribute_string(ncid_mol, NJ_NAME, "NAME",
                                   NETCDF_COMPAT) ) < 0) HERR(routineName);
    /* --- variables --- */
    dims[0] = molecule->Nv;
    dims[1] = mpi.nx;
    dims[2] = mpi.ny;
    dims[3] = infile.nz;
    /* Populations */
    if (input.p15d_wpop) {
      if (( file_dspace = H5Screate_simple(4, dims, NULL) ) < 0)
        HERR(routineName);
      if (( id_n = H5Dopen2(ncid_mol, VLEVEL_NAME,
                            H5P_DEFAULT)) < 0) HERR(routineName);
      if (molecule->nv != NULL) {
        if (( id_tmp = H5Dcreate(ncid_mol, POP_NAME, H5T_NATIVE_FLOAT,
                                 file_dspace, H5P_DEFAULT, plist,
                                 H5P_DEFAULT)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_x, 1)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_y, 2)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_z, 3)) < 0) HERR(routineName);
        io.aux_mol_pop[i] = id_tmp;
      }
      if (molecule->nvstar != NULL) {
        if (( id_tmp = H5Dcreate(ncid_mol, POPLTE_NAME, H5T_NATIVE_FLOAT,
                                 file_dspace, H5P_DEFAULT, plist,
                                 H5P_DEFAULT)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_n, 0)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_x, 1)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_y, 2)) < 0) HERR(routineName);
        if (( H5DSattach_scale(id_tmp, id_z, 3)) < 0) HERR(routineName);
        io.aux_mol_poplte[i] = id_tmp;
      }
      if (( H5Dclose(id_n) ) < 0) HERR(routineName);
      if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
      // TODO:  molecule->Ediss, molecule->Tmin, molecule->Tmax
    }
  if (( H5Dclose(id_x) ) < 0) HERR(routineName);
  if (( H5Dclose(id_y) ) < 0) HERR(routineName);
  if (( H5Dclose(id_z) ) < 0) HERR(routineName);
  } /* end active MOLECULES loop */
  io.aux_ncid = ncid;   /* Copy stuff to the IO data struct */
  if (( H5Pclose(plist) ) < 0) HERR(routineName);  /* Free hdf5 resources */
  /* Flush ensures file is created in case of crash */
  if (( H5Fflush(ncid, H5F_SCOPE_LOCAL) ) < 0) HERR(routineName);
  return;
}
/* ------- end   --------------------------   init_aux_new.c  --- */

/* ------- begin --------------------   init_aux_existing.c   --- */
void init_aux_existing(void) {
  const char routineName[] = "init_aux_existing";
  int     i, nlevel, nline, ncont;
  size_t  attr_size;
  hid_t   ncid, ncid_atom, ncid_mol, plist;
  H5T_class_t type_class;
  char    group_name[ARR_STRLEN], *atmosID;
  Atom   *atom;
  Molecule *molecule;

  /* Open the file with parallel MPI-IO access */
  if (( plist = H5Pcreate(H5P_FILE_ACCESS )) < 0) HERR(routineName);
  if (( H5Pset_fapl_mpio(plist, mpi.comm, mpi.info) ) < 0) HERR(routineName);
  if (( ncid = H5Fopen(AUX_FILE, H5F_ACC_RDWR, plist) ) < 0)
    HERR(routineName);
  if (( H5Pclose(plist) ) < 0) HERR(routineName);
  io.aux_ncid = ncid;
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

  /* Create arrays for multiple-atom/molecule output */
  io.aux_atom_ncid   = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
  io.aux_mol_ncid    = (hid_t *) malloc(atmos.Nactivemol  * sizeof(hid_t));
  if (input.p15d_wpop) {
      io.aux_atom_pop    = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
      io.aux_atom_poplte = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
      io.aux_mol_pop     = (hid_t *) malloc(atmos.Nactivemol  * sizeof(hid_t));
      io.aux_mol_poplte  = (hid_t *) malloc(atmos.Nactivemol  * sizeof(hid_t));
  }
  if (input.p15d_wrates) {
      io.aux_atom_RijL   = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
      io.aux_atom_RjiL   = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
      io.aux_atom_RijC   = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
      io.aux_atom_RjiC   = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
  }
  if (input.p15d_wcrates) {
      io.aux_atom_Cij    = (hid_t *) malloc(atmos.Nactiveatom * sizeof(hid_t));
  }
  /* --- Group loop over active ATOMS --- */
  for (i=0; i < atmos.Nactiveatom; i++) {
    atom = atmos.activeatoms[i];

    /* --- Get ncid of the atom group --- */
    sprintf(group_name,(atom->ID[1] == ' ') ? "atom_%.1s" : "atom_%.2s", atom->ID);
    if (( io.aux_atom_ncid[i] = H5Gopen(io.aux_ncid, group_name,
                                        H5P_DEFAULT) ) < 0)  HERR(routineName);
    ncid_atom = io.aux_atom_ncid[i];

    /* Check that dimension sizes match */
    if (( H5LTget_attribute_int(ncid_atom, ".", "nlevel", &nlevel) ) < 0)
      HERR(routineName);
    if (nlevel != atom->Nlevel) {
        sprintf(messageStr,
  	      "Number of levels mismatch: expected %d, found %d.",
  	      atom->Nlevel, (int)nlevel);
        Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    if (( H5LTget_attribute_int(ncid_atom, ".", "nline", &nline) ) < 0)
      HERR(routineName);
    if (nline != atom->Nline) {
      sprintf(messageStr,
          "Number of lines mismatch: expected %d, found %d.",
          atom->Nline, (int)nline);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    if (( H5LTget_attribute_int(ncid_atom, ".", "ncontinuum", &ncont) ) < 0)
      HERR(routineName);
    if (ncont != atom->Ncont) {
      sprintf(messageStr,
          "Number of continua mismatch: expected %d, found %d.",
          atom->Ncont, (int)ncont);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }

    /* --- Open datasets collectively ---*/
    if (input.p15d_wpop) {
      if (( io.aux_atom_pop[i] = H5Dopen(ncid_atom, POP_NAME,
                                         H5P_DEFAULT) )  < 0) HERR(routineName);
      if (( io.aux_atom_poplte[i] = H5Dopen(ncid_atom, POPLTE_NAME,
                                          H5P_DEFAULT) ) < 0) HERR(routineName);
    }
    if (input.p15d_wrates) {
      if (( io.aux_atom_RijL[i] = H5Dopen(ncid_atom, RIJ_L_NAME,
                                         H5P_DEFAULT) )  < 0) HERR(routineName);
      if (( io.aux_atom_RjiL[i] = H5Dopen(ncid_atom, RJI_L_NAME,
                                         H5P_DEFAULT) )  < 0) HERR(routineName);
      if (atom->Ncont > 0) {
          if (( io.aux_atom_RijC[i] = H5Dopen(ncid_atom, RIJ_C_NAME,
                                         H5P_DEFAULT) )  < 0) HERR(routineName);
          if (( io.aux_atom_RjiC[i] = H5Dopen(ncid_atom, RJI_C_NAME,
                                         H5P_DEFAULT) )  < 0) HERR(routineName);
      }
    }
    if (input.p15d_wcrates) {
      if (( io.aux_atom_Cij[i] = H5Dopen(ncid_atom, COLL_NAME,
                                         H5P_DEFAULT) )  < 0) HERR(routineName);
    }
  } /* end active ATOMS loop */


  /* --- Group loop over active MOLECULES --- */
  for (i=0; i < atmos.Nactivemol; i++) {
    molecule = atmos.activemols[i];

    /* --- Get ncid of the molecule group --- */
    sprintf( group_name, "molecule_%s", molecule->ID);
    if (( io.aux_mol_ncid[i] = H5Gopen(io.aux_ncid, group_name,
                                        H5P_DEFAULT) ) < 0)  HERR(routineName);
    ncid_mol = io.aux_mol_ncid[i];

    /* Check that dimension sizes match */
    if (( H5LTget_attribute_int(ncid_mol, ".", "nlevel_vibr", &nlevel) ) < 0)
      HERR(routineName);
    if (nlevel != molecule->Nv) {
      sprintf(messageStr,
          "Number of levels mismatch: expected %d, found %d.",
          molecule->Nv, (int)nlevel);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    if (( H5LTget_attribute_int(ncid_mol, ".", "nline_molecule", &nline) ) < 0)
      HERR(routineName);
    if (nline != molecule->Nrt) {
      sprintf(messageStr,
          "Number of lines mismatch: expected %d, found %d.",
          molecule->Nrt, (int)nline);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    if (( H5LTget_attribute_int(ncid_mol, ".", "nJ", &ncont) ) < 0)
      HERR(routineName);
    if (ncont != molecule->NJ) {
      sprintf(messageStr,
          "Number of J mismatch: expected %d, found %d.",
          molecule->NJ, (int)ncont);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }

    /* --- Open datasets collectively ---*/
    if (input.p15d_wpop) {
      if (( io.aux_mol_pop[i] = H5Dopen(ncid_mol, POP_NAME,
                                         H5P_DEFAULT) )  < 0) HERR(routineName);
      if (( io.aux_mol_poplte[i] = H5Dopen(ncid_mol, POPLTE_NAME,
                                         H5P_DEFAULT) )  < 0) HERR(routineName);
    }
  } /* end active MOLECULES loop */

  return;
}
/* ------- end   ---------------------   init_aux_existing.c   --- */

/* ------- begin --------------------------   close_hdf5_aux.c --- */
void close_hdf5_aux(void)
/* Closes the aux file */
{
  const char routineName[] = "close_hdf5_aux";
  Atom   *atom;
  int i;

  /* Close all datasets and groups */
  for (i=0; i < atmos.Nactiveatom; i++) {
    atom = atmos.activeatoms[i];
    if (input.p15d_wpop) {
      if (( H5Dclose(io.aux_atom_pop[i]) ) < 0) HERR(routineName);
      if (( H5Dclose(io.aux_atom_poplte[i]) ) < 0) HERR(routineName);
    }
    if (input.p15d_wrates) {
      if (( H5Dclose(io.aux_atom_RijL[i]) ) < 0) HERR(routineName);
      if (( H5Dclose(io.aux_atom_RjiL[i]) ) < 0) HERR(routineName);
      if (atom->Ncont > 0) {
          if (( H5Dclose(io.aux_atom_RijC[i]) ) < 0) HERR(routineName);
          if (( H5Dclose(io.aux_atom_RjiC[i]) ) < 0) HERR(routineName);
      }
    }
    if (input.p15d_wcrates) {
      if (( H5Dclose(io.aux_atom_Cij[i]) ) < 0) HERR(routineName);
    }
    if (( H5Gclose(io.aux_atom_ncid[i]) ) < 0) HERR(routineName);
  }
  for (i=0; i < atmos.Nactivemol; i++) {
    if (input.p15d_wpop) {
      if (( H5Dclose(io.aux_mol_pop[i]) ) < 0) HERR(routineName);
      if (( H5Dclose(io.aux_mol_poplte[i]) ) < 0) HERR(routineName);
    }
    if (( H5Gclose(io.aux_mol_ncid[i]) ) < 0) HERR(routineName);
  }

  /* Close file */
  if (( H5Fclose(io.aux_ncid) ) < 0) HERR(routineName);

  /* Free memory */
  free(io.aux_atom_ncid);
  free(io.aux_mol_ncid);
  if (input.p15d_wpop) {
    free(io.aux_atom_pop);
    free(io.aux_atom_poplte);
    free(io.aux_mol_pop);
    free(io.aux_mol_poplte);
  }
  if (input.p15d_wrates) {
    free(io.aux_atom_RijL);
    free(io.aux_atom_RjiL);
    free(io.aux_atom_RijC);
    free(io.aux_atom_RjiC);
  }
  if (input.p15d_wcrates) {
    free(io.aux_atom_Cij);
  }

  return;
}
/* ------- end   --------------------------   close_hdf5_aux.c --- */

/* ------- begin --------------------------   writeAux_all.c   --- */
void writeAux_all(void) {
  const char routineName[] = "writeAux_all";
  hsize_t  offset[] = {0, 0, 0, 0};
  hsize_t  count[] = {1, 1, 1, 1};
  hsize_t  dims[4];
  hid_t    file_dspace, mem_dspace;
  Atom      *atom;
  Molecule  *molecule;
  int        nact, task;
  long       ind = 0;

  /* ATOM loop */
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    ind = 0;   /* Index of each task inside variables */
    for (task = 0; task < mpi.Ntasks; task++) {
      /* If there was a crash, no data were written into buffer variables */
      if (mpi.convergence[task] < 0) continue;
      /* Memory dataspace */
      dims[0] = atom->Nlevel;
      dims[1] = infile.nz - mpi.zcut_hist[task];
      if (( mem_dspace = H5Screate_simple(2, dims, NULL) ) < 0)
        HERR(routineName);
      /* File dataspace */
      offset[0] = 0;
      offset[1] = mpi.taskmap[task + mpi.my_start][0];
      offset[2] = mpi.taskmap[task + mpi.my_start][1];
      offset[3] = mpi.zcut_hist[task];
      count[0] = atom->Nlevel;
      count[1] = 1;
      count[2] = 1;
      count[3] = infile.nz - mpi.zcut_hist[task];
      if (input.p15d_wpop) {
        if (( file_dspace = H5Dget_space(io.aux_atom_pop[nact]) ) < 0)
          HERR(routineName);
        if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                                    NULL, count, NULL) ) < 0) HERR(routineName);
        if (( H5Dwrite(io.aux_atom_pop[nact], H5T_NATIVE_DOUBLE,
               mem_dspace, file_dspace, H5P_DEFAULT,
               &iobuf.n[nact][ind * atom->Nlevel]) ) < 0) HERR(routineName);
        if (( H5Dwrite(io.aux_atom_poplte[nact], H5T_NATIVE_DOUBLE,
               mem_dspace, file_dspace, H5P_DEFAULT,
               &iobuf.nstar[nact][ind * atom->Nlevel]) ) < 0) HERR(routineName);
        /* release dataspace resources */
        if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
        if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
      }
      if (input.p15d_wrates) {
        /* Memory dataspace */
        dims[0] = atom->Nline;
        if (( mem_dspace = H5Screate_simple(2, dims, NULL) ) < 0)
          HERR(routineName);
        count[0] = atom->Nline;
        if (( file_dspace = H5Dget_space(io.aux_atom_RijL[nact]) ) < 0)
          HERR(routineName);
        if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                                  NULL, count, NULL) ) < 0) HERR(routineName);
        if (( H5Dwrite(io.aux_atom_RijL[nact], H5T_NATIVE_DOUBLE,
                 mem_dspace, file_dspace, H5P_DEFAULT,
                 &iobuf.RijL[nact][ind * atom->Nline]) ) < 0) HERR(routineName);
        if (( H5Dwrite(io.aux_atom_RjiL[nact], H5T_NATIVE_DOUBLE,
                 mem_dspace, file_dspace, H5P_DEFAULT,
                 &iobuf.RjiL[nact][ind * atom->Nline]) ) < 0) HERR(routineName);
        /* release dataspace resources */
        if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
        if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);

        if (atom->Ncont > 0) {
            /* Memory dataspace */
            dims[0] = atom->Ncont;
            if (( mem_dspace = H5Screate_simple(2, dims, NULL) ) < 0)
              HERR(routineName);
            count[0] = atom->Ncont;
            if (( file_dspace = H5Dget_space(io.aux_atom_RijC[nact]) ) < 0)
              HERR(routineName);
            if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                                    NULL, count, NULL) ) < 0) HERR(routineName);
            if (( H5Dwrite(io.aux_atom_RijC[nact], H5T_NATIVE_DOUBLE,
                 mem_dspace, file_dspace, H5P_DEFAULT,
                 &iobuf.RijC[nact][ind * atom->Ncont]) ) < 0) HERR(routineName);
            if (( H5Dwrite(io.aux_atom_RjiC[nact], H5T_NATIVE_DOUBLE,
                 mem_dspace, file_dspace, H5P_DEFAULT,
                 &iobuf.RjiC[nact][ind * atom->Ncont]) ) < 0) HERR(routineName);
            /* release dataspace resources */
            if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
            if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
        }
      }
      ind += count[3];
    }
  }

  /* MOLECULE loop */
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    ind = 0;   /* Index of each task inside variables */
    for (task = 0; task < mpi.Ntasks; task++) {
      /* If there was a crash, no data were written into buffer variables */
      if (mpi.convergence[task] < 0) continue;
      if (input.p15d_wpop) {
        /* Memory dataspace */
        dims[0] = molecule->Nv;
        dims[1] = infile.nz - mpi.zcut_hist[task];
        if (( mem_dspace = H5Screate_simple(2, dims, NULL) ) < 0)
          HERR(routineName);
        /* File dataspace */
        offset[0] = 0;
        offset[1] = mpi.taskmap[task + mpi.my_start][0];
        offset[2] = mpi.taskmap[task + mpi.my_start][1];
        offset[3] = mpi.zcut_hist[task];
        count[0] = molecule->Nv;
        count[3] = infile.nz - mpi.zcut_hist[task];
        if (( file_dspace = H5Dget_space(io.aux_mol_pop[nact]) ) < 0)
          HERR(routineName);
        if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                                  NULL, count, NULL) ) < 0) HERR(routineName);
        if (( H5Dwrite(io.aux_mol_pop[nact], H5T_NATIVE_DOUBLE,
                mem_dspace, file_dspace, H5P_DEFAULT,
                &iobuf.nv[nact][ind * molecule->Nv]) ) < 0) HERR(routineName);
        if (( H5Dwrite(io.aux_mol_poplte[nact], H5T_NATIVE_DOUBLE,
               mem_dspace, file_dspace, H5P_DEFAULT,
              &iobuf.nvstar[nact][ind * molecule->Nv]) ) < 0) HERR(routineName);
        /* release dataspace resources */
        if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
        if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
      }
      ind += count[3];
    }
  }
  return;
}
/* ------- end --------------------------   writeAux_all.c   --- */


/* ------- begin --------------------------   writeAux_p.c     --- */
void writeAux_p(void) {
  /* this will write: populations, radrates, coll, damping */
  const char routineName[] = "writeAux_p";
  int      nact, kr, i, j, k, idx;
  hsize_t  offset[] = {0, 0, 0, 0, 0};
  hsize_t  count[] = {1, 1, 1, 1, 1};
  hsize_t  dims[5];
  double  *coll_rates;
  hid_t    file_dspace, mem_dspace;
  Atom      *atom;
  Molecule  *molecule;
  AtomicLine      *line;
  AtomicContinuum *continuum;

  /* --- ATOMS ---  */
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    /* Write populations */
    /* Memory dataspace */
    dims[0] = atom->Nlevel;
    dims[1] = atmos.Nspace;
    /* File dataspace */
    offset[0] = 0;
    offset[1] = mpi.ix;
    offset[2] = mpi.iy;
    offset[3] = mpi.zcut;
    count[0] = atom->Nlevel;
    count[1] = 1;
    count[2] = 1;
    count[3] = atmos.Nspace;
    if (input.p15d_wpop) {
      if (( mem_dspace = H5Screate_simple(2, dims, NULL) ) < 0)
          HERR(routineName);
      if (( file_dspace = H5Dget_space(io.aux_atom_pop[nact]) ) < 0)
          HERR(routineName);
      if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                                NULL, count, NULL) ) < 0) HERR(routineName);
      if (( H5Dwrite(io.aux_atom_pop[nact], H5T_NATIVE_DOUBLE, mem_dspace,
                  file_dspace, H5P_DEFAULT, atom->n[0]) ) < 0) HERR(routineName);
      if (( H5Dwrite(io.aux_atom_poplte[nact], H5T_NATIVE_DOUBLE, mem_dspace,
             file_dspace, H5P_DEFAULT, atom->nstar[0]) ) < 0) HERR(routineName);
      if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
      if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
    }
    if (input.p15d_wrates) {
      /* Write radiative rates */
      dims[0] = atmos.Nspace;
      count[0] = 1;
      if (( mem_dspace = H5Screate_simple(1, dims, NULL) ) < 0)
          HERR(routineName);
      if (( file_dspace = H5Dget_space(io.aux_atom_RijL[nact]) ) < 0)
          HERR(routineName);
      for (kr=0; kr < atom->Nline; kr++) {
          offset[0] = kr;
          line = &atom->line[kr];
          /* Convention for output, both R and C, ij means transition i -> j */
          if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                                    NULL, count, NULL) ) < 0) HERR(routineName);
          if (( H5Dwrite(io.aux_atom_RijL[nact], H5T_NATIVE_DOUBLE, mem_dspace,
                  file_dspace, H5P_DEFAULT, line->Rij) ) < 0) HERR(routineName);
          if (( H5Dwrite(io.aux_atom_RjiL[nact], H5T_NATIVE_DOUBLE, mem_dspace,
                  file_dspace, H5P_DEFAULT, line->Rji) ) < 0) HERR(routineName);
      }
      for (kr=0; kr < atom->Ncont; kr++) {
          offset[0] = kr;
          continuum = &atom->continuum[kr];
          if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                                    NULL, count, NULL) ) < 0) HERR(routineName);
          if (( H5Dwrite(io.aux_atom_RijC[nact], H5T_NATIVE_DOUBLE, mem_dspace,
             file_dspace, H5P_DEFAULT, continuum->Rij) ) < 0) HERR(routineName);
          if (( H5Dwrite(io.aux_atom_RjiC[nact], H5T_NATIVE_DOUBLE, mem_dspace,
             file_dspace, H5P_DEFAULT, continuum->Rji) ) < 0) HERR(routineName);
      }
      if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
      if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
    }
    if (input.p15d_wcrates) {
      /* Write collisional rates */
      /* Copy data to better-aligned memory structure */
      coll_rates = (double *) calloc(SQ(atom->Nlevel) * atmos.Nspace, sizeof(double));
      for (i=0; i < atom->Nlevel; i++) {
        for (j=0; j < atom->Nlevel; j++) {
          if (i != j) {
            for (k=0; k < atmos.Nspace; k++) {
              idx = (i * atom->Nlevel + j) * atmos.Nspace + k;
              coll_rates[idx] = atom->C[j * atom->Nlevel + i][k];
            }
          }
        }
      }
      /* Memory dataspace */
      dims[0] = atom->Nlevel;
      dims[1] = atom->Nlevel;
      dims[2] = atmos.Nspace;
      /* File dataspace */
      offset[0] = 0;
      offset[1] = 0;
      offset[2] = mpi.ix;
      offset[3] = mpi.iy;
      offset[4] = mpi.zcut;
      count[0] = atom->Nlevel;
      count[1] = atom->Nlevel;
      count[2] = 1;
      count[3] = 1;
      count[4] = atmos.Nspace;
      if (( mem_dspace = H5Screate_simple(3, dims, NULL) ) < 0)
          HERR(routineName);
      if (( file_dspace = H5Dget_space(io.aux_atom_Cij[nact]) ) < 0)
          HERR(routineName);
      if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                                NULL, count, NULL) ) < 0) HERR(routineName);
      if (( H5Dwrite(io.aux_atom_Cij[nact], H5T_NATIVE_DOUBLE, mem_dspace,
                   file_dspace, H5P_DEFAULT, coll_rates) ) < 0) HERR(routineName);
      if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
      if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
      free(coll_rates);
    }
  }

  /* --- MOLECULES ---  */
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    if (input.p15d_wpop) {
      /* Write populations */
      dims[0] = molecule->Nv;
      dims[1] = atmos.Nspace;
      if (( mem_dspace = H5Screate_simple(2, dims, NULL) ) < 0)
          HERR(routineName);
      offset[0] = 0;
      offset[1] = mpi.ix;
      offset[2] = mpi.iy;
      offset[3] = mpi.zcut;
      count[0] = molecule->Nv;
      count[3] = atmos.Nspace;
      if (( file_dspace = H5Dget_space(io.aux_mol_pop[nact]) ) < 0)
          HERR(routineName);
      if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                                NULL, count, NULL) ) < 0) HERR(routineName);
      if (( H5Dwrite(io.aux_mol_pop[nact], H5T_NATIVE_DOUBLE, mem_dspace,
            file_dspace, H5P_DEFAULT, molecule->nv[0]) ) < 0) HERR(routineName);
      if (( H5Dwrite(io.aux_mol_poplte[nact], H5T_NATIVE_DOUBLE, mem_dspace,
        file_dspace, H5P_DEFAULT, molecule->nvstar[0]) ) < 0) HERR(routineName);
      if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
      if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
    }
  }
  return;
}
/* ------- end   --------------------------   writeAux_p.c     --- */

/* ------- begin -------------------------- readPopulations_p.c -- */
void readPopulations(Atom *atom) {

  /* --- Read populations from file.

   Note: readPopulations only reads the true populations and not
         the LTE populations.
   Note2: This routine will not work yet in HFD5. Access to the aux file
          needs to be better regulated in OLD_POPULATIONS cases (ie, append
          and read/write acess). The code below assumes the file is already
          open, which may not be the case.
         --                                            -------------- */
  const char routineName[] = "readPopulations_p";
  char    group_name[ARR_STRLEN], *atmosID;
  int nz, nlevel;
  H5T_class_t type_class;
  size_t attr_size;
  hsize_t offset[] = {0, 0, 0, 0};
  hsize_t count[] = {1, 1, 1, 1};
  hsize_t dims[4];
  hid_t ncid, file_dspace, mem_dspace, pop_var;

  /* --- Open atom group --- */
  sprintf(group_name,(atom->ID[1] == ' ') ? "atom_%.1s" : "atom_%.2s", atom->ID);
  if (( ncid = H5Gopen(io.aux_ncid, group_name, H5P_DEFAULT) ) < 0)
    HERR(routineName);
  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if (( H5LTget_attribute_info(io.aux_ncid, "/", "atmosID", NULL, &type_class,
                               &attr_size) ) < 0) HERR(routineName);
  atmosID = (char *) malloc(attr_size + 1);
  if (( H5LTget_attribute_string(io.aux_ncid, "/", "atmosID", atmosID) ) < 0)
    HERR(routineName);
  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,
       "Indata file was calculated for different atmosphere (%s) than current",
       atmosID);
    Error(WARNING, routineName, messageStr);
  }
  free(atmosID);
  /* Check that dimension sizes match */
  if (( H5LTget_attribute_int(ncid, ".", "nlevel", &nlevel) ) < 0)
    HERR(routineName);
  if (nlevel != atom->Nlevel) {
      sprintf(messageStr,
          "Number of levels mismatch: expected %d, found %d.",
          atom->Nlevel, (int)nlevel);
      Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  if (( H5LTget_attribute_int(io.aux_ncid, ".", "nz", &nz) ) < 0)
    HERR(routineName);
  if (nz < atmos.Nspace) {
    sprintf(messageStr,
      "Number of depth points mismatch: expected %ld, found %d.",
      atmos.Nspace, (int)nz);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* --- Read data --- */
  if (( pop_var = H5Dopen2(ncid, POP_NAME, H5P_DEFAULT)) < 0)
    HERR(routineName);
  dims[0] = atom->Nlevel;
  dims[1] = atmos.Nspace;
  if (( mem_dspace = H5Screate_simple(2, dims, NULL) ) < 0) HERR(routineName);
  /* File dataspace */
  offset[0] = 0;
  offset[1] = mpi.ix;
  offset[2] = mpi.iy;
  offset[3] = mpi.zcut;
  count[0] = atom->Nlevel;
  count[3] = atmos.Nspace;
  if (( file_dspace = H5Dget_space(pop_var) ) < 0) HERR(routineName);
  if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                            NULL, count, NULL) ) < 0) HERR(routineName);
  if (( H5Dread(pop_var, H5T_NATIVE_DOUBLE, mem_dspace,
                file_dspace, H5P_DEFAULT, atom->n[0]) ) < 0) HERR(routineName);
  if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
  if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
  if (( H5Dclose(pop_var) ) < 0) HERR(routineName);

  return;
}
/* ------- end   -------------------------- readPopulations_p.c -- */


/* ------- begin -------------------------- readMolPops_p.c -- */
void readMolPops(Molecule *molecule) {

  /* --- Read populations from file.

   Note: readPopulations only reads the true populations and not
         the LTE populations.
         --                                            -------------- */

  const char routineName[] = "readMolPops";
  char    group_name[ARR_STRLEN], *atmosID;

  int nz, nlevel;
  H5T_class_t type_class;
  size_t attr_size;
  hsize_t offset[] = {0, 0, 0, 0};
  hsize_t count[] = {1, 1, 1, 1};
  hsize_t dims[4];
  hid_t ncid, file_dspace, mem_dspace, pop_var;

  /* --- Open molecule group --- */
  sprintf(group_name, "molecule_%s", molecule->ID);
  if (( ncid = H5Gopen(io.aux_ncid, group_name, H5P_DEFAULT) ) < 0)
    HERR(routineName);
  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if (( H5LTget_attribute_info(io.aux_ncid, "/", "atmosID", NULL, &type_class,
                               &attr_size) ) < 0) HERR(routineName);
  atmosID = (char *) malloc(attr_size + 1);
  if (( H5LTget_attribute_string(io.aux_ncid, "/", "atmosID", atmosID) ) < 0)
    HERR(routineName);
  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,
       "Indata file was calculated for different atmosphere (%s) than current",
       atmosID);
    Error(WARNING, routineName, messageStr);
  }
  free(atmosID);
  /* Check that dimension sizes match */
  if (( H5LTget_attribute_int(ncid, ".", "nlevel_vibr", &nlevel) ) < 0)
    HERR(routineName);
  if (nlevel != molecule->Nv) {
    sprintf(messageStr,
      "Number of levels mismatch: expected %d, found %d.",
      molecule->Nv, (int)nlevel);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  if (( H5LTget_attribute_int(io.aux_ncid, ".", "nz", &nz) ) < 0)
    HERR(routineName);
  if (nz < atmos.Nspace) {
    sprintf(messageStr,
      "Number of depth points mismatch: expected %ld, found %d.",
      atmos.Nspace, (int)nz);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* --- Read data --- */
  if (( pop_var = H5Dopen2(ncid, POP_NAME, H5P_DEFAULT)) < 0)
    HERR(routineName);
  dims[0] = molecule->Nv;
  dims[1] = atmos.Nspace;
  if (( mem_dspace = H5Screate_simple(2, dims, NULL) ) < 0) HERR(routineName);
  /* File dataspace */
  offset[0] = 0;
  offset[1] = mpi.ix;
  offset[2] = mpi.iy;
  offset[3] = mpi.zcut;
  count[0] = molecule->Nv;
  count[3] = atmos.Nspace;
  if (( file_dspace = H5Dget_space(pop_var) ) < 0) HERR(routineName);
  if (( H5Sselect_hyperslab(file_dspace, H5S_SELECT_SET, offset,
                            NULL, count, NULL) ) < 0) HERR(routineName);
  if (( H5Dread(pop_var, H5T_NATIVE_DOUBLE, mem_dspace,
            file_dspace, H5P_DEFAULT, molecule->nv[0]) ) < 0) HERR(routineName);
  if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
  if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
  if (( H5Dclose(pop_var) ) < 0) HERR(routineName);

  return;
}
/* ------- end   -------------------------- readMolPops_p.c -- */

/* ------- begin -------------------------- writeMolPops.c ---------- */

void writeMolPops(struct Molecule *molecule)
{
  /* Blank function to avoid unresolved symbols.
     (not to worry, molecular pops are written by writeAux_p)
   */
  return;
}
/* ------- end ---------------------------- writeMolPops.c ---------- */




// LOW PRIORITY:
/* ------- begin -------------------------- writeAtom_metadata_p.c */
// this will write the same as writeAtom() for active atoms
/* ------- begin -------------------------- writeAtom_metadata_p.c */
