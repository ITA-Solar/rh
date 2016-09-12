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
/* Creates the netCDF file for the auxiliary data */
  const char routineName[] = "init_aux_new";
  int     i, nx_id, ny_id, nspace_id,
          nrays_id, nwad_id, nwai_id, nlevel_id, nline_id, ncont_id, dimids[5];
  int     nwave, *ai_idx, *ad_idx, nai, nad;
  hid_t   plist, ncid, file_dspace, ncid_atom, ncid_mol, ncid_op;
  hsize_t dims[4];
  float   fillval = 9.96921e+36;  /* from netcdf */ 
  double *wave, *lambda_air;
  char    group_name[ARR_STRLEN];
  Atom   *atom;
  Molecule *molecule;
  ActiveSet *as;
  
  /* Create the file  */
  if (( plist = H5Pcreate(H5P_FILE_ACCESS )) < 0) HERR(routineName);
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
  io.aux_atom_ncid   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_pop    = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_poplte = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RijL   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RjiL   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RijC   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RjiC   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  /* And for multiple molecule output */
  io.aux_mol_ncid    = (int *) malloc(atmos.Nactivemol  * sizeof(int));
  io.aux_mol_pop     = (int *) malloc(atmos.Nactivemol  * sizeof(int));
  io.aux_mol_poplte  = (int *) malloc(atmos.Nactivemol  * sizeof(int));
  
  /* Fill value */
  if (( plist = H5Pcreate(H5P_DATASET_CREATE) ) < 0) HERR(routineName);
  if (( H5Pset_fill_value(plist, H5T_NATIVE_FLOAT, &fillval) ) < 0)
      HERR(routineName);  
  
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
                                atom->Nlevel, 1)) < 0) HERR(routineName);
    if (( H5LTset_attribute_int(ncid_atom, ".", "nline",
                                atom->Nline, 1)) < 0) HERR(routineName);
    if (( H5LTset_attribute_int(ncid_atom, ".", "ncontinuum",
                                atom->Ncont, 1)) < 0) HERR(routineName);
    /* --- variables --- */
    dims[0] = atom->Nlevel;
    dims[1] = mpi.nx;
    dims[2] = mpi.ny;
    dims[3] = infile.nz;
    if (( file_dspace = H5Screate_simple(4, dims, NULL) ) < 0)
      HERR(routineName);      
    /* Populations */
    if (atom->n != NULL) 
      if (( io.aux_atom_pop[i] = H5Dcreate(ncid_atom, POP_NAME,
                                   H5T_NATIVE_FLOAT, file_dspace, H5P_DEFAULT,
                                   plist, H5P_DEFAULT)) < 0) HERR(routineName);
    if (atom->nstar != NULL) 
      if (( io.aux_atom_poplte[i] = H5Dcreate(ncid_atom, POPLTE_NAME,
                                   H5T_NATIVE_FLOAT, file_dspace, H5P_DEFAULT,
                                   plist, H5P_DEFAULT)) < 0) HERR(routineName);
    if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
    /* Radiative rates */
    dims[0] = atom->Nline;
    dims[1] = mpi.nx;
    dims[2] = mpi.ny;
    dims[3] = infile.nz;
    if (( file_dspace = H5Screate_simple(4, dims, NULL) ) < 0)
      HERR(routineName);
    if (( io.aux_atom_RijL[i] = H5Dcreate(ncid_atom, RIJ_L_NAME,
                                   H5T_NATIVE_FLOAT, file_dspace, H5P_DEFAULT,
                                   plist, H5P_DEFAULT)) < 0) HERR(routineName);
    if (( io.aux_atom_RjiL[i] = H5Dcreate(ncid_atom, RJI_L_NAME,
                                   H5T_NATIVE_FLOAT, file_dspace, H5P_DEFAULT,
                                   plist, H5P_DEFAULT)) < 0) HERR(routineName);
    if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
    dims[0] = atom->Ncont;
    if (( file_dspace = H5Screate_simple(4, dims, NULL) ) < 0)
      HERR(routineName);
    if (( io.aux_atom_RijC[i] = H5Dcreate(ncid_atom, RIJ_C_NAME,
                                   H5T_NATIVE_FLOAT, file_dspace, H5P_DEFAULT,
                                   plist, H5P_DEFAULT)) < 0) HERR(routineName);
    if (( io.aux_atom_RjiC[i] = H5Dcreate(ncid_atom, RJI_C_NAME,
                                   H5T_NATIVE_FLOAT, file_dspace, H5P_DEFAULT,
                                   plist, H5P_DEFAULT)) < 0) HERR(routineName);
    if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
  } /* end active ATOMS loop */

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
                                molecule->Nv, 1)) < 0) HERR(routineName);
    if (( H5LTset_attribute_int(ncid_mol, ".", "nline_molecule",
                                molecule->Nrt, 1)) < 0) HERR(routineName);
    if (( H5LTset_attribute_int(ncid_mol, ".", "nJ",
                                molecule->NJ, 1)) < 0) HERR(routineName);
    /* --- variables --- */
    dims[0] = molecule->Nv;
    dims[1] = mpi.nx;
    dims[2] = mpi.ny;
    dims[3] = infile.nz;
    if (( file_dspace = H5Screate_simple(4, dims, NULL) ) < 0)
      HERR(routineName);      
    /* Populations */
    if (molecule->nv != NULL) 
      if (( io.aux_mol_pop[i] = H5Dcreate(ncid_mol, POP_NAME,
                                   H5T_NATIVE_FLOAT, file_dspace, H5P_DEFAULT,
                                   plist, H5P_DEFAULT)) < 0) HERR(routineName);
    if (molecule->nvstar != NULL) 
      if (( io.aux_mol_poplte[i] = H5Dcreate(ncid_mol, POPLTE_NAME,
                                   H5T_NATIVE_FLOAT, file_dspace, H5P_DEFAULT,
                                   plist, H5P_DEFAULT)) < 0) HERR(routineName);
    if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
    // TODO:  molecule->Ediss, molecule->Tmin, molecule->Tmax
  } /* end active MOLECULES loop */
  io.aux_ncid = ncid;   /* Copy stuff to the IO data struct */
  if (( H5Pclose(plist) ) < 0) HERR(routineName);  /* Free hdf5 resources */
  return;
}
/* ------- end   --------------------------   init_aux_new.c  --- */

/* ------- begin --------------------   init_aux_existing.c   --- */
void init_aux_existing(void) {
  const char routineName[] = "init_aux_existing";
  int     ierror, ncid, i, dimid;
  size_t  len_id, nlevel, nline, ncont;
  char    group_name[ARR_STRLEN], *atmosID;
  Atom   *atom;
  Molecule *molecule;

  /* --- Open the file --- */
  if ((ierror = nc_open_par(AUX_FILE, NC_WRITE | NC_MPIPOSIX, 
                  mpi.comm, mpi.info, &io.aux_ncid))) ERR(ierror,routineName);

  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if ((ierror = nc_inq_attlen(io.aux_ncid, NC_GLOBAL, "atmosID", &len_id ))) 
    ERR(ierror,routineName);

  atmosID = (char *) malloc(len_id+1);

  if ((ierror = nc_get_att_text(io.aux_ncid, NC_GLOBAL, "atmosID", atmosID ))) 
    ERR(ierror,routineName);

  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,
         "Populations were derived from different atmosphere (%s) than current",
	    atmosID);
    Error(WARNING, routineName, messageStr);
    }
  free(atmosID);

  /* Check that dimension sizes match 
  if ((ierror = nc_inq_dimid(io.aux_ncid, "nz", &dimid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_inq_dimlen(io.aux_ncid, dimid, &nz ))) 
    ERR(ierror,routineName);    

  if (nz < atmos.Nspace) {
    sprintf(messageStr,
	    "Number of depth points mismatch: expected %ld, found %d.",
	    atmos.Nspace, (int)nz);
    Error(WARNING, routineName, messageStr);
  }
  */

  /* Create arrays for multiple-atom/molecule output */
  io.aux_atom_ncid   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_pop    = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_poplte = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RijL   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RjiL   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RijC   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RjiC   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_mol_ncid    = (int *) malloc(atmos.Nactivemol  * sizeof(int));
  io.aux_mol_pop     = (int *) malloc(atmos.Nactivemol  * sizeof(int));
  io.aux_mol_poplte  = (int *) malloc(atmos.Nactivemol  * sizeof(int));

  /* --- Group loop over active ATOMS --- */
  for (i=0; i < atmos.Nactiveatom; i++) {
    atom = atmos.activeatoms[i];
    
    /* --- Get ncid of the atom group --- */
    sprintf(group_name,(atom->ID[1] == ' ') ? "atom_%.1s" : "atom_%.2s", atom->ID);

    if ((ierror = nc_inq_ncid(io.aux_ncid, group_name, &ncid)))
      ERR(ierror,routineName);

    io.aux_atom_ncid[i] = ncid;
    
    /* Check that dimension sizes match */
    if ((ierror = nc_inq_dimid(ncid, "nlevel", &dimid ))) 
      ERR(ierror,routineName);  
    if ((ierror = nc_inq_dimlen(ncid, dimid, &nlevel ))) 
      ERR(ierror,routineName);   
    
    if (nlevel != atom->Nlevel) {
      sprintf(messageStr,
	      "Number of levels mismatch: expected %d, found %d.",
	      atom->Nlevel, (int)nlevel);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }

    if ((ierror = nc_inq_dimid(ncid, "nline", &dimid ))) 
      ERR(ierror,routineName);  
    if ((ierror = nc_inq_dimlen(ncid, dimid, &nline ))) 
      ERR(ierror,routineName);   
    
    if (nline != atom->Nline) {
      sprintf(messageStr,
	      "Number of lines mismatch: expected %d, found %d.",
	      atom->Nline, (int)nline);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    
    if ((ierror = nc_inq_dimid(ncid, "ncontinuum", &dimid ))) 
      ERR(ierror,routineName);  
    if ((ierror = nc_inq_dimlen(ncid, dimid, &ncont ))) 
      ERR(ierror,routineName);   
    
    if (ncont != atom->Ncont) {
      sprintf(messageStr,
	      "Number of continua mismatch: expected %d, found %d.",
	      atom->Ncont, (int)ncont);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }

    /* Get var ids */
  
    if ((ierror = nc_inq_varid(ncid, POP_NAME,    &io.aux_atom_pop[i]    ))) 
      ERR(ierror,routineName);  
    if ((ierror = nc_inq_varid(ncid, POPLTE_NAME, &io.aux_atom_poplte[i] ))) 
      ERR(ierror,routineName);  
    if ((ierror = nc_inq_varid(ncid, RIJ_L_NAME,  &io.aux_atom_RijL[i]   ))) 
      ERR(ierror,routineName);  
    if ((ierror = nc_inq_varid(ncid, RJI_L_NAME,  &io.aux_atom_RjiL[i]   ))) 
      ERR(ierror,routineName);
    if ((ierror = nc_inq_varid(ncid, RIJ_C_NAME,  &io.aux_atom_RijC[i]   ))) 
      ERR(ierror,routineName);
    if ((ierror = nc_inq_varid(ncid, RJI_C_NAME,  &io.aux_atom_RjiC[i]   ))) 
      ERR(ierror,routineName);  
  } /* end active ATOMS loop */


  /* --- Group loop over active MOLECULES --- */
  for (i=0; i < atmos.Nactivemol; i++) {
    molecule = atmos.activemols[i];
    
    /* --- Get ncid of the atom group --- */
    /* Get group name */
    sprintf( group_name, "molecule_%s", molecule->ID);
    if ((ierror = nc_inq_ncid(io.aux_ncid, group_name, &ncid)))
      ERR(ierror,routineName);

    io.aux_mol_ncid[i] = ncid;
    
    /* Check that dimension sizes match */
    if ((ierror = nc_inq_dimid(ncid, "nlevel_vibr", &dimid ))) 
      ERR(ierror,routineName);  
    if ((ierror = nc_inq_dimlen(ncid, dimid, &nlevel ))) 
      ERR(ierror,routineName);   
    
    if (nlevel != molecule->Nv) {
      sprintf(messageStr,
	      "Number of levels mismatch: expected %d, found %d.",
	      molecule->Nv, (int)nlevel);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }

    if ((ierror = nc_inq_dimid(ncid, "nline_molecule", &dimid ))) 
      ERR(ierror,routineName);  
    if ((ierror = nc_inq_dimlen(ncid, dimid, &nline ))) 
      ERR(ierror,routineName);   
    
    if (nline != molecule->Nrt) {
      sprintf(messageStr,
	      "Number of lines mismatch: expected %d, found %d.",
	      molecule->Nrt, (int)nline);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    
    if ((ierror = nc_inq_dimid(ncid, "nJ", &dimid ))) 
      ERR(ierror,routineName);  
    if ((ierror = nc_inq_dimlen(ncid, dimid, &ncont ))) 
      ERR(ierror,routineName);   
    
    if (ncont != molecule->NJ) {
      sprintf(messageStr,
	      "Number of J mismatch: expected %d, found %d.",
	      molecule->NJ, (int)ncont);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }

    /* Get var ids */
  
    if ((ierror = nc_inq_varid(ncid, POP_NAME,    &io.aux_mol_pop[i]    ))) 
      ERR(ierror,routineName);  
    if ((ierror = nc_inq_varid(ncid, POPLTE_NAME, &io.aux_mol_poplte[i] ))) 
      ERR(ierror,routineName); 
  } /* end active MOLECULES loop */

  if (input.p15d_wxtra) {
    /* Now for the opacity group */
    if ((ierror = nc_inq_ncid(io.aux_ncid, "opacity", &ncid)))
      ERR(ierror,routineName);
    
    io.aux_op_ncid = ncid;

    /* Get var ids */
    if ((ierror = nc_inq_varid(ncid, CHI_AI_NAME, &io.aux_op_chi_ai ))) 
      ERR(ierror,routineName); 
    if ((ierror = nc_inq_varid(ncid, ETA_AI_NAME, &io.aux_op_eta_ai ))) 
      ERR(ierror,routineName); 
    if ((ierror = nc_inq_varid(ncid, CHI_AD_NAME, &io.aux_op_chi_ad ))) 
      ERR(ierror,routineName); 
    if ((ierror = nc_inq_varid(ncid, ETA_AD_NAME, &io.aux_op_eta_ad ))) 
      ERR(ierror,routineName);  
  }
  
  
  return;
}
/* ------- end   ---------------------   init_aux_existing.c   --- */



/* ------- begin --------------------------   close_hdf5_aux.c --- */
void close_hdf5_aux(void)
/* Closes the spec netCDF file */ 
{
  const char routineName[] = "close_hdf5_aux";
  int        ierror;

  if ((ierror = nc_close(io.aux_ncid))) ERR(ierror,routineName);

  free(io.aux_atom_ncid);
  free(io.aux_atom_pop);
  free(io.aux_atom_poplte);
  free(io.aux_atom_RijL);
  free(io.aux_atom_RjiL);
  free(io.aux_atom_RijC);
  free(io.aux_atom_RjiC);
  free(io.aux_mol_ncid);
  free(io.aux_mol_pop);
  free(io.aux_mol_poplte);

  return; 
}
/* ------- end   --------------------------   close_hdf5_aux.c --- */

/* ------- begin --------------------------   writeAux_all.c   --- */
void writeAux_all(void) {
  const char routineName[] = "writeAux_all";
  int        ncid, ierror, nact, task;
  long       ind = 0;
  size_t     start[] = {0, 0, 0, 0};
  size_t     count[] = {1, 1, 1, 1};
  Atom      *atom;
  Molecule  *molecule;



  /* --- Task loop --- */
  for (task = 0; task < mpi.Ntasks; task++) {
    
    /* If there was a crash, no data were written into buffer variables */
    if (mpi.convergence[task] < 0) continue;

    start[1] = mpi.taskmap[task + mpi.my_start][0];  count[2] = 1;
    start[2] = mpi.taskmap[task + mpi.my_start][1];  count[2] = 1;
    start[3] = mpi.zcut_hist[task];   count[3] = infile.nz - start[3];
    

    /* ATOM loop */
    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      atom = atmos.activeatoms[nact];
      ncid = io.aux_atom_ncid[nact];

      /* n, nstar */
      start[0] = 0; count[0] = atom->Nlevel;
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_pop[nact], start, count,
	    	    &iobuf.n[nact][ind*atom->Nlevel] ))) ERR(ierror,routineName);

      if ((ierror=nc_put_vara_double(ncid,io.aux_atom_poplte[nact], start, count,
                &iobuf.nstar[nact][ind*atom->Nlevel] ))) ERR(ierror,routineName); 

      /* Rij, Rji for lines */
      start[0] = 0; count[0] = atom->Nline;
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RijL[nact], start, count,
                 &iobuf.RijL[nact][ind*atom->Nline] ))) ERR(ierror,routineName);  
  
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RjiL[nact], start, count,
                 &iobuf.RjiL[nact][ind*atom->Nline] ))) ERR(ierror,routineName);  

      /* Rij, Rji for continua */
      start[0] = 0; count[0] = atom->Ncont;
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RijC[nact], start, count,
                 &iobuf.RijC[nact][ind*atom->Ncont] ))) ERR(ierror,routineName);
      
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RjiC[nact], start, count,
                 &iobuf.RjiC[nact][ind*atom->Ncont] ))) ERR(ierror,routineName);
    }

    /* MOLECULE loop */
    for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
      molecule = atmos.activemols[nact];
      
      /* nv, nvstar */
      start[0] = 0; count[0] = molecule->Nv;
      if ((ierror=nc_put_vara_double(ncid, io.aux_mol_pop[nact],   start, count,
       	          &iobuf.nv[nact][ind*molecule->Nv] ))) ERR(ierror,routineName);

      if ((ierror=nc_put_vara_double(ncid,io.aux_mol_poplte[nact], start, count,
              &iobuf.nvstar[nact][ind*molecule->Nv] ))) ERR(ierror,routineName);  
    }

    ind += count[3];
  } 
  /* --- End task loop --- */


  return;
}
/* ------- end --------------------------   writeAux_all.c   --- */


/* ------- begin --------------------------   writeAux_p.c     --- */
void writeAux_p(void) {
/* this will write: populations, radrates, coll, damping */
  const char routineName[] = "writeAux_p";
  int        ierror, ncid, nact, kr;
  size_t     start[] = {0, 0, 0, 0, 0};
  size_t     count[] = {1, 1, 1, 1, 1};
  double   **adamp, **E;
  Atom      *atom;
  Molecule  *molecule;
  AtomicLine      *line;
  AtomicContinuum *continuum;
  MolecularLine   *mrt;

  /* --- Main ATOM loop --- */
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    ncid = io.aux_atom_ncid[nact];

    /* --- write populations --- */
    start[1] = mpi.ix;
    start[2] = mpi.iy;
    start[3] = mpi.zcut;
    count[0] = atom->Nlevel; 
    count[3] = atmos.Nspace;
    if (atom->n != NULL) 
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_pop[nact], start, count,
     			     atom->n[0] )))     ERR(ierror,"writeAux_p n");
    if (atom->nstar != NULL)
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_poplte[nact], start, count,
				     atom->nstar[0] ))) ERR(ierror,"writeAux_p nstar");
      
    if (input.p15d_wxtra) {
      /* --- write damping --- */
      adamp = matrix_double(atom->Nline, atmos.Nspace);
      
      for (kr=0; kr < atom->Nline; kr++) {
	if (atom->line[kr].Voigt)
	  Damping(&atom->line[kr], adamp[kr]);
      }
      count[0] = atom->Nline;
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_damp[nact], start, count,
				     adamp[0] ))) ERR(ierror,"writeAux_p adamp");
      freeMatrix((void **) adamp);
    }


    /* --- write radiative rates --- */
    count[0] = 1;
    /* for bound-bound transitions */
    for (kr=0; kr < atom->Nline; kr++) {
      start[0] = kr;
      line = &atom->line[kr];
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RijL[nact], start, count,
				     line->Rij ))) ERR(ierror,"writeAux_p RijL");
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RjiL[nact], start, count,
				     line->Rji ))) ERR(ierror,"writeAux_p RjiL");
    }
    /* for bound-free transitions */
    for (kr=0; kr < atom->Ncont; kr++) {
      start[0] = kr;
      continuum = &atom->continuum[kr];
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RijC[nact], start, count,
				     continuum->Rij ))) ERR(ierror,"writeAux_p RijC");
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RjiC[nact], start, count,
				     continuum->Rji ))) ERR(ierror,"writeAux_p RjiC");
    }
  } /* End ATOM loop */


  /* --- Main MOLECULE loop --- */
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    ncid     = io.aux_mol_ncid[nact];

    /* --- write populations --- */
    start[1] = mpi.ix;
    start[2] = mpi.iy;
    start[3] = mpi.zcut;
    count[0] = molecule->Nv; 
    count[3] = atmos.Nspace;
    if (molecule->nv != NULL) 
      if ((ierror=nc_put_vara_double(ncid, io.aux_mol_pop[nact], start, count,
     			     molecule->nv[0] )))     ERR(ierror,routineName);
    if (molecule->nvstar != NULL)
      if ((ierror=nc_put_vara_double(ncid, io.aux_mol_poplte[nact], start, count,
			     molecule->nvstar[0] ))) ERR(ierror,routineName);
  } /* End MOLECULE loop */


  return;
}
/* ------- end   --------------------------   writeAux_p.c     --- */

/* ------- begin -------------------------- readPopulations_p.c -- */
void readPopulations(Atom *atom) {

  /* --- Read populations from file.

   Note: readPopulations only reads the true populations and not
         the LTE populations. 
         --                                            -------------- */

  const char routineName[] = "readPopulations_p";
  char    group_name[ARR_STRLEN], *atmosID;

  int ierror, ncid, varid, dimid;
  size_t nlevel, nz, len_id;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};



  /* --- Get ncid of the atom group --- */
  sprintf(group_name,(atom->ID[1] == ' ') ? "atom_%.1s" : "atom_%.2s", atom->ID);


  if ((ierror = nc_inq_ncid(io.aux_ncid, group_name, &ncid)))
    ERR(ierror,routineName);


  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if ((ierror = nc_inq_attlen(io.aux_ncid, NC_GLOBAL, "atmosID", &len_id ))) 
    ERR(ierror,routineName);

  atmosID = (char *) malloc(len_id+1);

  if ((ierror = nc_get_att_text(io.aux_ncid, NC_GLOBAL, "atmosID", atmosID ))) 
    ERR(ierror,routineName);

  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,
         "Populations were derived from different atmosphere (%s) than current",
	    atmosID);
    Error(WARNING, routineName, messageStr);
    }
  free(atmosID);

  /* Check that dimension sizes match */
  if ((ierror = nc_inq_dimid(ncid, "nlevel", &dimid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_inq_dimlen(ncid, dimid, &nlevel ))) 
    ERR(ierror,routineName);   

  if (nlevel != atom->Nlevel) {
    sprintf(messageStr,
	    "Number of levels mismatch: expected %d, found %d.",
	    atom->Nlevel, (int)nlevel);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  if ((ierror = nc_inq_dimid(io.aux_ncid, "nz", &dimid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_inq_dimlen(io.aux_ncid, dimid, &nz ))) 
    ERR(ierror,routineName);    

  if (nz < atmos.Nspace) {
    sprintf(messageStr,
	    "Number of depth points mismatch: expected %ld, found %d.",
	    atmos.Nspace, (int)nz);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
 

  /* --- Get variable id for populations --- */
  if ((ierror = nc_inq_varid(ncid, POP_NAME, &varid ))) 
    ERR(ierror,routineName); 

  /* --- Read data --- */
  start[1] = mpi.ix;
  start[2] = mpi.iy;
  start[3] = mpi.zcut;
  
  count[0] = atom->Nlevel;
  count[3] = atmos.Nspace;

  /* read as double, although it is written as float */
  if ((ierror = nc_get_vara_double(ncid, varid, start, count, atom->n[0] )))
    ERR(ierror,routineName);


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

  int ierror, ncid, varid, dimid;
  size_t nlevel, nz, len_id, nline, ncont;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};


  /* --- Get ncid of the atom group --- */
  sprintf(group_name, "molecule_%s", molecule->ID);


  if ((ierror = nc_inq_ncid(io.aux_ncid, group_name, &ncid)))
    ERR(ierror,routineName);


  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if ((ierror = nc_inq_attlen(io.aux_ncid, NC_GLOBAL, "atmosID", &len_id ))) 
    ERR(ierror,routineName);

  atmosID = (char *) malloc(len_id+1);

  if ((ierror = nc_get_att_text(io.aux_ncid, NC_GLOBAL, "atmosID", atmosID ))) 
    ERR(ierror,routineName);

  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,
         "Populations were derived from different atmosphere (%s) than current",
	    atmosID);
    Error(WARNING, routineName, messageStr);
    }
  free(atmosID);


  /* Check that dimension sizes match */
  if ((ierror = nc_inq_dimid(io.aux_ncid, "nz", &dimid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_inq_dimlen(io.aux_ncid, dimid, &nz ))) 
    ERR(ierror,routineName);    

  if (nz < atmos.Nspace) {
    sprintf(messageStr,
	    "Number of depth points mismatch: expected %ld, found %d.",
	    atmos.Nspace, (int)nz);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
 
  if ((ierror = nc_inq_dimid(ncid, "nlevel_vibr", &dimid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_inq_dimlen(ncid, dimid, &nlevel ))) 
    ERR(ierror,routineName);   
  
  if (nlevel != molecule->Nv) {
    sprintf(messageStr,
	    "Number of levels mismatch: expected %d, found %d.",
	    molecule->Nv, (int)nlevel);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  if ((ierror = nc_inq_dimid(ncid, "nline_molecule", &dimid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_inq_dimlen(ncid, dimid, &nline ))) 
    ERR(ierror,routineName);   
  
  if (nline != molecule->Nrt) {
    sprintf(messageStr,
	    "Number of lines mismatch: expected %d, found %d.",
	    molecule->Nrt, (int)nline);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  
  if ((ierror = nc_inq_dimid(ncid, "nJ", &dimid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_inq_dimlen(ncid, dimid, &ncont ))) 
    ERR(ierror,routineName);   
  
  if (ncont != molecule->NJ) {
    sprintf(messageStr,
	    "Number of J mismatch: expected %d, found %d.",
	    molecule->NJ, (int)ncont);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  

  /* --- Get variable id for populations --- */
  if ((ierror = nc_inq_varid(ncid, POP_NAME, &varid ))) 
    ERR(ierror,routineName); 

  /* --- Read data --- */
  start[1] = mpi.ix;
  start[2] = mpi.iy;
  start[3] = mpi.zcut;
  count[0] = molecule->Nv;
  count[3] = atmos.Nspace;

  /* read as double, although it is written as float */
  if ((ierror = nc_get_vara_double(ncid, varid, start, count, molecule->nv[0] )))
    ERR(ierror,routineName);


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
