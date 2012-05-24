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
extern NCDF_Atmos_file infile;
extern MPI_data  mpi;
extern IO_data   io; 
extern IO_buffer iobuf; 

/* ------- begin --------------------------   init_ncdf_aux.c     --- */
void init_ncdf_aux(void) {
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

  if ((use_old) || (input.p15d_rerun)) init_aux_old(); else init_aux_new();

  return;
}
/* ------- end --------------------------   init_ncdf_aux.c     --- */

/* ------- begin --------------------------   init_aux_new.c --   --- */
void init_aux_new(void) {
/* Creates the netCDF file for the auxiliary data */
  const char routineName[] = "init_ncdf_aux";
  int     ierror, ncid, ncid_atom, ncid_mol, ncid_op, i, nx_id, ny_id, nspace_id,
          nrays_id, nwad_id, nwai_id, nlevel_id, nline_id, ncont_id, dimids[5];
  int     nwave, *ai_idx, *ad_idx, nai, nad;
  double *wave, *lambda_air;
  char    group_name[ARR_STRLEN];
  FILE   *test;
  Atom   *atom;
  Molecule *molecule;
  ActiveSet *as;
  

  /* Create the file  */
  if ((ierror = nc_create_par(AUX_FILE, NC_NETCDF4 | NC_CLOBBER | NC_MPIPOSIX, 
  //if ((ierror = nc_create_par(AUX_FILE, NC_NETCDF4 | NC_CLOBBER | NC_MPIIO, 
			      mpi.comm, mpi.info, &ncid))) ERR(ierror,routineName);
  
  
  /* --- Definitions for the root group --- */
  /* dimensions */
  if ((ierror = nc_def_dim(ncid, "nx",       mpi.nx,           &nx_id    ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid, "ny",       mpi.ny,           &ny_id    ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid, "nz",       infile.nz,        &nspace_id))) 
    ERR(ierror,routineName);

  /* attributes */
  if ((ierror = nc_put_att_text(ncid, NC_GLOBAL, "atmosID", strlen(atmos.ID),
				atmos.ID ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text(ncid, NC_GLOBAL, "svn_id", strlen(mpi.svn_id),
				mpi.svn_id ))) ERR(ierror,routineName);


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
  /* Allocate only if we're writting the extra stuff */
  if (input.p15d_wxtra) {
    io.aux_atom_coll   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
    io.aux_atom_damp   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
    io.aux_atom_vbroad = (int *) malloc(atmos.Nactiveatom * sizeof(int));
    io.aux_mol_E       = (int *) malloc(atmos.Nactivemol  * sizeof(int));
    io.aux_mol_vbroad  = (int *) malloc(atmos.Nactivemol  * sizeof(int));
  }
  
  /* --- Group loop over active ATOMS --- */
  for (i=0; i < atmos.Nactiveatom; i++) {
    atom = atmos.activeatoms[i];

    /* Get group name */
    sprintf( group_name, (atom->ID[1] == ' ') ? "atom_%.1s" : "atom_%.2s", atom->ID);
    if ((ierror = nc_def_grp(ncid, group_name, &io.aux_atom_ncid[i]))) 
      ERR(ierror,routineName);

    ncid_atom = io.aux_atom_ncid[i];

    /* --- dimensions --- */
    if ((ierror = nc_def_dim(ncid_atom, "nlevel",     atom->Nlevel, &nlevel_id)))
      ERR(ierror,routineName);  
    if ((ierror = nc_def_dim(ncid_atom, "nline",      atom->Nline,  &nline_id )))
      ERR(ierror,routineName);  
    if ((ierror = nc_def_dim(ncid_atom, "ncontinuum", atom->Ncont,  &ncont_id )))
      ERR(ierror,routineName);  

    /* --- variables --- */
    dimids[0] = nlevel_id;
    dimids[1] = nx_id;
    dimids[2] = ny_id;
    dimids[3] = nspace_id;
    /* Populations */
    if (atom->n != NULL) 
      if ((ierror = nc_def_var(ncid_atom, POP_NAME,    NC_FLOAT, 4, dimids,
			       &io.aux_atom_pop[i]))) ERR(ierror,routineName);
    if (atom->nstar != NULL)
      if ((ierror = nc_def_var(ncid_atom, POPLTE_NAME, NC_FLOAT, 4, dimids,
		            &io.aux_atom_poplte[i]))) ERR(ierror,routineName);

    /* Radiative rates */
    dimids[0] = nline_id;
    dimids[1] = nx_id;
    dimids[2] = ny_id;
    dimids[3] = nspace_id;
    if ((ierror = nc_def_var(ncid_atom, RIJ_L_NAME, NC_FLOAT, 4, dimids,
			     &io.aux_atom_RijL[i]))) ERR(ierror,routineName);
    if ((ierror = nc_def_var(ncid_atom, RJI_L_NAME, NC_FLOAT, 4, dimids,
			     &io.aux_atom_RjiL[i]))) ERR(ierror,routineName); 
    dimids[0] = ncont_id;
    if ((ierror = nc_def_var(ncid_atom, RIJ_C_NAME, NC_FLOAT, 4, dimids,
			     &io.aux_atom_RijC[i]))) ERR(ierror,routineName);    
    if ((ierror = nc_def_var(ncid_atom, RJI_C_NAME, NC_FLOAT, 4, dimids,
			     &io.aux_atom_RjiC[i]))) ERR(ierror,routineName); 

    /* Write extra output? */
    if (input.p15d_wxtra) {
      /* Collision rates */
      dimids[0] = nlevel_id;
      dimids[1] = nlevel_id;
      dimids[2] = nx_id;
      dimids[3] = ny_id;
      dimids[4] = nspace_id;
      if ((ierror = nc_def_var(ncid_atom, COLL_NAME, NC_FLOAT, 5, dimids,
			       &io.aux_atom_coll[i]))) ERR(ierror,routineName);

      /* Damping */
      dimids[0] = nline_id;
      dimids[1] = nx_id;
      dimids[2] = ny_id;
      dimids[3] = nspace_id;
      if ((ierror = nc_def_var(ncid_atom, DAMP_NAME, NC_FLOAT, 4, dimids,
			       &io.aux_atom_damp[i]))) ERR(ierror,routineName);

      
      /* Broadening velocity */
      dimids[0] = nx_id;
      dimids[1] = ny_id;
      dimids[2] = nspace_id;
      if ((ierror = nc_def_var(ncid_atom, VBROAD_NAME, NC_FLOAT, 3, dimids,
			       &io.aux_atom_vbroad[i]))) ERR(ierror,routineName);
    } 

  } /* end active ATOMS loop */


  /* --- Group loop over active MOLECULES --- */
  for (i=0; i < atmos.Nactivemol; i++) {
    molecule = atmos.activemols[i];

    /* Get group name */
    sprintf( group_name, "molecule_%s", molecule->ID);
    if ((ierror = nc_def_grp(ncid, group_name, &io.aux_mol_ncid[i]))) 
      ERR(ierror,routineName);

    ncid_mol = io.aux_mol_ncid[i];

    /* --- dimensions --- */
    if ((ierror = nc_def_dim(ncid_mol, "nlevel_vibr",    molecule->Nv, 
      &nlevel_id)))  ERR(ierror,routineName);  
    if ((ierror = nc_def_dim(ncid_mol, "nline_molecule", molecule->Nrt,  
      &nline_id )))  ERR(ierror,routineName);  
    if ((ierror = nc_def_dim(ncid_mol, "nJ",             molecule->NJ,
      &ncont_id )))  ERR(ierror,routineName);  

    /* --- variables --- */
    dimids[0] = nlevel_id;
    dimids[1] = nx_id;
    dimids[2] = ny_id;
    dimids[3] = nspace_id;
    /* Populations */
    if (molecule->nv != NULL) 
      if ((ierror = nc_def_var(ncid_mol, POP_NAME,    NC_FLOAT, 4, dimids,
			       &io.aux_mol_pop[i]))) ERR(ierror,routineName);
    if (molecule->nvstar != NULL)
      if ((ierror = nc_def_var(ncid_mol, POPLTE_NAME, NC_FLOAT, 4, dimids,
		            &io.aux_mol_poplte[i]))) ERR(ierror,routineName);

    if (input.p15d_wxtra) {
      /* Energy matrix */
      dimids[0] = nx_id;
      dimids[1] = ny_id;
      dimids[2] = nlevel_id;
      dimids[3] = ncont_id;
      if ((ierror = nc_def_var(ncid_mol, EM_NAME,     NC_FLOAT, 4, dimids,
			       &io.aux_mol_E[i]))) ERR(ierror,routineName);
      
      /* Broadening velocity */
      dimids[0] = nx_id;
      dimids[1] = ny_id;
      dimids[2] = nspace_id;
      if ((ierror = nc_def_var(ncid_mol, VBROAD_NAME, NC_FLOAT, 3, dimids,
			       &io.aux_mol_vbroad[i]))) ERR(ierror,routineName); 
    }
    /* --- attributes --- */
    // TODO:  molecule->Ediss, molecule->Tmin, molecule->Tmax


  } /* end active MOLECULES loop */


  /* Opacity for the quadrature intensities is only written if write extra */ 
  if (input.p15d_wxtra) {
    /* --- Create opacity group --- */
    if ((ierror = nc_def_grp(ncid, "opacity", &ncid_op))) ERR(ierror,routineName);

    /* Find out which wavelengths contain angle dependent opacities, which do not, 
       and get corresponding indices
     
       spectrum.Nspect = total number of wavelengths
       nwave           = wavelengths where there are active transitions (includes H)
       nai             = wavelengths where there are NO bound-bound, polarised or 
                         angle-dep transitions
       nad             = wavelengths where there ARE bound-bound, polarised or 
                         angle-dep transitions

    */
    nwave = 0; nai = 0; nad = 0;
    /* Indices relative to nwave */
    ai_idx  = (int *)    malloc(spectrum.Nspect * sizeof(int));
    ad_idx  = (int *)    malloc(spectrum.Nspect * sizeof(int));
    wave    = (double *) malloc(spectrum.Nspect * sizeof(double));

    for (i=0; i < spectrum.Nspect; i++) {
      as = &spectrum.as[i];
      if (containsActive(as)) {
	if ( containsPolarized(as) || (containsPRDline(as) && input.PRD_angle_dep != PRD_ANGLE_INDEP) 
	     || (atmos.moving && containsBoundBound(as))) {
	  ad_idx[nad++] = nwave;
	} else {
	  ai_idx[nai++] = nwave;
	}
	wave[nwave] = spectrum.lambda[i];
	++nwave;
      }
    }

    /* --- Definitions for the OPACITY group --- */
    /* dimensions */
    if ((ierror = nc_def_dim(ncid_op, "nrays",    geometry.Nrays,  &nrays_id ))) 
      ERR(ierror,routineName);
    if ((ierror = nc_def_dim(ncid_op, NW_AD_NAME, nad,             &nwad_id  ))) 
      ERR(ierror,routineName); 
    if ((ierror = nc_def_dim(ncid_op, NW_AI_NAME, nai,             &nwai_id  ))) 
      ERR(ierror,routineName); 
    

    /* variables*/
    dimids[0] = nwai_id;
    dimids[1] = nx_id;
    dimids[2] = ny_id;
    dimids[3] = nspace_id;
    if ((ierror = nc_def_var(ncid_op, CHI_AI_NAME, NC_FLOAT, 4, dimids,
			     &io.aux_op_chi_ai))) ERR(ierror,routineName);
    if ((ierror = nc_def_var(ncid_op, ETA_AI_NAME, NC_FLOAT, 4, dimids,
			     &io.aux_op_eta_ai))) ERR(ierror,routineName);
    
    dimids[0] = nwad_id;
    dimids[1] = nrays_id;
    dimids[2] = nx_id;
    dimids[3] = ny_id;
    dimids[4] = nspace_id;
    if ((ierror = nc_def_var(ncid_op, CHI_AD_NAME, NC_FLOAT, 5, dimids,
			     &io.aux_op_chi_ad))) ERR(ierror,routineName);
    if ((ierror = nc_def_var(ncid_op, ETA_AD_NAME, NC_FLOAT, 5, dimids,
			     &io.aux_op_eta_ad))) ERR(ierror,routineName);


    /* attributes */
    /* write active set wavelengths as attribute, not variable! */
    if (spectrum.vacuum_to_air) {
      lambda_air = (double *) malloc(nwave * sizeof(double));
      vacuum_to_air(nwave, wave, lambda_air);
      if ((ierror=nc_put_att_double(ncid_op, NC_GLOBAL, WAVET_NAME, NC_DOUBLE, 
				    nwave, lambda_air ))) ERR(ierror,routineName);
      free(lambda_air);
      free(wave);
    } else {
      if ((ierror=nc_put_att_double(ncid_op, NC_GLOBAL, WAVET_NAME, NC_DOUBLE, 
				    nwave, wave )))       ERR(ierror,routineName);
      free(wave);
    }
  
    /* wavelength angle indices */
    if ((ierror=nc_put_att_int(ncid_op, NC_GLOBAL, WAVE_AD_IDX_NAME,
			       NC_INT, nad, ad_idx ))) ERR(ierror,routineName);
    if ((ierror=nc_put_att_int(ncid_op, NC_GLOBAL, WAVE_AI_IDX_NAME,
			       NC_INT, nai, ai_idx ))) ERR(ierror,routineName);
    free(ad_idx);
    free(ai_idx);

  } /* end write extra condition */


  /* End define mode */
  if ((ierror = nc_enddef(ncid))) ERR(ierror,routineName);

  /* Copy stuff to the IO data struct */
  io.aux_ncid    = ncid;
  io.aux_op_ncid = ncid_op;


  return;
}
/* ------- end   --------------------------   init_aux_new.c  --- */

/* ------- begin --------------------------   init_aux_old.c   --- */
void init_aux_old(void) {
  const char routineName[] = "init_aux_old";
  int     ierror, ncid, i, dimid;
  size_t  len_id, nlevel, nline, ncont, nz;
  char    group_name[ARR_STRLEN], *atmosID;
  FILE   *test;
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
  if (input.p15d_wxtra) {
    io.aux_atom_coll   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
    io.aux_atom_damp   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
    io.aux_atom_vbroad = (int *) malloc(atmos.Nactiveatom * sizeof(int));
    io.aux_mol_E       = (int *) malloc(atmos.Nactivemol  * sizeof(int));
    io.aux_mol_vbroad  = (int *) malloc(atmos.Nactivemol  * sizeof(int));
  }

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

    if (input.p15d_wxtra) {
      if ((ierror = nc_inq_varid(ncid, COLL_NAME,   &io.aux_atom_coll[i]   ))) 
	ERR(ierror,routineName); 
      if ((ierror = nc_inq_varid(ncid, DAMP_NAME,   &io.aux_atom_damp[i]   ))) 
	ERR(ierror,routineName);  
      if ((ierror = nc_inq_varid(ncid, VBROAD_NAME, &io.aux_atom_vbroad[i] )))
	ERR(ierror,routineName);  
    }

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
 
    if (input.p15d_wxtra) {
      if ((ierror = nc_inq_varid(ncid, EM_NAME,     &io.aux_mol_E[i]      ))) 
	ERR(ierror,routineName);  
      if ((ierror = nc_inq_varid(ncid, VBROAD_NAME, &io.aux_mol_vbroad[i] )))
	ERR(ierror,routineName);  
    }

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
/* ------- end   --------------------------   init_aux_old.c   --- */



/* ------- begin --------------------------   close_ncdf_aux.c --- */
void close_ncdf_aux(void)
/* Closes the spec netCDF file */ 
{
  const char routineName[] = "close_ncdf_aux";
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

  if (input.p15d_wxtra) {
    free(io.aux_atom_coll);
    free(io.aux_atom_damp);
    free(io.aux_atom_vbroad);
    free(io.aux_mol_E);
    free(io.aux_mol_vbroad);
  }

  return; 
}
/* ------- end   --------------------------   close_ncdf_aux.c --- */

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

      /* Set collective access for variables 
      if ((ierror = nc_var_par_access(ncid, io.aux_atom_pop[nact], NC_COLLECTIVE)))
	ERR(ierror,routineName);
      if ((ierror = nc_var_par_access(ncid, io.aux_atom_poplte[nact], NC_COLLECTIVE)))
	ERR(ierror,routineName);
      if ((ierror = nc_var_par_access(ncid, io.aux_atom_RijL[nact], NC_COLLECTIVE)))
	ERR(ierror,routineName);
      if ((ierror = nc_var_par_access(ncid, io.aux_atom_RjiL[nact], NC_COLLECTIVE)))
	ERR(ierror,routineName);
      if ((ierror = nc_var_par_access(ncid, io.aux_atom_RijC[nact], NC_COLLECTIVE)))
	ERR(ierror,routineName);
      if ((ierror = nc_var_par_access(ncid, io.aux_atom_RjiC[nact], NC_COLLECTIVE)))
	ERR(ierror,routineName);
      // end collective set 
      */


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

    if (input.p15d_wxtra) {
      /* --- write broadening velocity --- */
      start[0] = mpi.ix;
      start[1] = mpi.iy;
      start[2] = mpi.zcut;
      count[0] = 1;
      count[1] = 1;
      count[2] = atmos.Nspace;
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_vbroad[nact], start, count,
				     atom->vbroad ))) ERR(ierror,"writeAux_p vbroad");

      /* --- write collisional rates --- */
      start[0] = 0;
      start[1] = 0;
      start[2] = mpi.ix;
      start[3] = mpi.iy;
      start[4] = mpi.zcut;
      count[0] = atom->Nlevel;
      count[1] = atom->Nlevel;
      count[2] = 1;
      count[3] = 1;
      count[4] = atmos.Nspace;
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_coll[nact], start, count,
				     atom->C[0] ))) ERR(ierror,"writeAux_p C");
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

    if (input.p15d_wxtra) {
      /* --- write energy matrix --- */
      if (molecule->Nv > 0  &&  molecule->NJ > 0) {
	E = matrix_double(molecule->Nv, molecule->NJ);
	
	for (kr = 0;  kr < molecule->Nrt;  kr++) {
	  mrt = molecule->mrt + kr;
	  
	  E[mrt->vi][(int) (mrt->gi - 1)/2] = mrt->Ei;
	  E[mrt->vj][(int) (mrt->gj - 1)/2] = mrt->Ej;
	}
	
	start[0] = mpi.ix;
	start[1] = mpi.iy;
	start[2] = 0;
	start[3] = 0;
	count[0] = 1;
	count[1] = 1;
	count[2] = molecule->Nv;
	count[3] = molecule->NJ;
	if ((ierror=nc_put_vara_double(ncid, io.aux_mol_E[nact], start, count,
				       E[0] ))) ERR(ierror,routineName);
	freeMatrix((void **) E);
      }

      /* --- write broadening velocity --- */
      start[0] = mpi.ix;
      start[1] = mpi.iy;
      start[2] = mpi.zcut;
      count[0] = 1;
      count[1] = 1;
      count[2] = atmos.Nspace;
      if ((ierror=nc_put_vara_double(ncid, io.aux_mol_vbroad[nact], start, count,
				     molecule->vbroad ))) ERR(ierror,routineName);
    } 

  } /* End MOLECULE loop */


  /* Now write opacity (also goes on aux file) */
  if (input.p15d_wxtra) writeOpacity_p();

  return;
}
/* ------- end   --------------------------   writeAux_p.c     --- */

/* ------- begin --------------------------   writeOpacity_p.c --- */
void writeOpacity_p(void) {
  const char routineName[] = "writeOpacity_p";
  int      ierror, ncid, mu, nspect;
  size_t   start_ai[] = {0, 0, 0, 0};
  size_t   count_ai[] = {1, 1, 1, 1};
  size_t   start_ad[] = {0, 0, 0, 0, 0};
  size_t   count_ad[] = {1, 1, 1, 1, 1};
  bool_t   boundbound, PRD_angle_dep, polarized, crosscoupling, initialize, 
           to_obs;
  ActiveSet *as;

  ncid = io.aux_op_ncid;

  start_ai[1] = mpi.ix;
  start_ai[2] = mpi.iy;
  start_ai[3] = mpi.zcut;
  start_ad[2] = mpi.ix;
  start_ad[3] = mpi.iy;
  start_ad[4] = mpi.zcut;
  count_ai[3] = atmos.Nspace;
  count_ad[4] = atmos.Nspace;

  for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
      as = &spectrum.as[nspect];

      if (containsActive(as)) {
	alloc_as(nspect, crosscoupling=FALSE);

      /* --- Check whether current active set includes a bound-bound
             and/or polarized transition and/or angledependent PRD
             transition. Otherwise, only angle-independent opacity and
             source functions are needed --            -------------- */ 

	boundbound    = containsBoundBound(as);
	PRD_angle_dep = (containsPRDline(as) && input.PRD_angle_dep != PRD_ANGLE_INDEP);
	polarized     = containsPolarized(as);

	/* --- Case of angle-dependent opacity and source function -- - */
	if (polarized || PRD_angle_dep || (atmos.moving && boundbound)) {
	  /* using wave angle dependent indices */
	  for (mu = 0;  mu < atmos.Nrays;  mu++) {
	    start_ad[1] = mu;

	    initialize = (mu == 0);
	    Opacity(nspect, mu, to_obs=TRUE, initialize);

	    if ((ierror=nc_put_vara_double(ncid, io.aux_op_chi_ad, start_ad,
			      count_ad, as->chi ))) ERR(ierror,routineName);
	    if ((ierror=nc_put_vara_double(ncid, io.aux_op_eta_ad, start_ad, 
                              count_ad, as->eta ))) ERR(ierror,routineName);
	  }
	  ++start_ad[0];
	} else {
	  /* using wave angle independent indices */
	  Opacity(nspect, 0, to_obs=TRUE, initialize=TRUE);
	  
	  if ((ierror=nc_put_vara_double(ncid, io.aux_op_chi_ai, start_ai,
                            count_ai, as->chi ))) ERR(ierror,routineName);
	  if ((ierror=nc_put_vara_double(ncid, io.aux_op_eta_ai, start_ai,
			    count_ai, as->eta ))) ERR(ierror,routineName);
	  ++start_ai[0];
	}
	
	free_as(nspect, crosscoupling=FALSE);
      }
  }

  return;
}
/* ------- end   --------------------------   writeOpacity_p.c --- */

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

/* ------- begin -------------------------- readRadRates_p.c   --- */
bool_t readRadRate(Atom *atom) {
  const char routineName[] = "readRadRates_p";
  char    group_name[ARR_STRLEN], *atmosID;

  int ierror, ncid, kr, dimid, RijL_id, RjiL_id, RijC_id, RjiC_id;
  size_t nline, ncont, nz, len_id;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};
  AtomicLine *line;
  AtomicContinuum *continuum;
  
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
  if ((ierror = nc_inq_dimid(ncid, "nline", &dimid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_inq_dimlen(ncid, dimid, &nline ))) 
    ERR(ierror,routineName);   

  if (nline != atom->Nline) {
    sprintf(messageStr,
	    "Number of lines mismatch: expected %d, found %d.",
	    atom->Nline, (int)nline);
    Error(WARNING, routineName, messageStr);
  }

  if ((ierror = nc_inq_dimid(ncid, "ncontinuum", &dimid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_inq_dimlen(ncid, dimid, &ncont ))) 
    ERR(ierror,routineName);   

  if (ncont != atom->Ncont) {
    sprintf(messageStr,
	    "Number of continua mismatch: expected %d, found %d.",
	    atom->Ncont, (int)ncont);
    Error(WARNING, routineName, messageStr);
  }

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
 

  /* --- Get variable ids  --- */
  if ((ierror = nc_inq_varid(ncid, RIJ_L_NAME, &RijL_id ))) 
    ERR(ierror,routineName); 
  if ((ierror = nc_inq_varid(ncid, RJI_L_NAME, &RjiL_id ))) 
    ERR(ierror,routineName); 
  if ((ierror = nc_inq_varid(ncid, RIJ_C_NAME, &RijC_id ))) 
    ERR(ierror,routineName); 
  if ((ierror = nc_inq_varid(ncid, RJI_C_NAME, &RjiC_id ))) 
    ERR(ierror,routineName); 


  /* --- Read data --- */
  start[1] = mpi.ix;
  start[2] = mpi.iy;
  start[3] = mpi.zcut;
  count[3] = atmos.Nspace;

  for (kr = 0; kr < atom->Nline; kr++) {
    line = &atom->line[kr];
    
    start[0] = kr;

    if ((ierror = nc_get_vara_double(ncid, RijL_id, start, count, line->Rij )))
      ERR(ierror,routineName);
    if ((ierror = nc_get_vara_double(ncid, RjiL_id, start, count, line->Rji )))
      ERR(ierror,routineName);
  }

  for (kr = 0; kr < atom->Ncont; kr++) {
    continuum = &atom->continuum[kr];
    
    start[0] = kr;

    if ((ierror = nc_get_vara_double(ncid,RijC_id,start,count,continuum->Rij)))
      ERR(ierror,routineName);
    if ((ierror = nc_get_vara_double(ncid,RjiC_id,start,count,continuum->Rji)))
      ERR(ierror,routineName);
  }


  return TRUE;
}
/* ------- end   -------------------------- readRadRates_p.c   --- */

/* ------- begin -------------------------- readCollRates_p.c  --- */
void readCollRates_p(Atom *atom) {
  const char routineName[] = "readCollRates_p";
  char    group_name[ARR_STRLEN], *atmosID;

  int ierror, ncid, varid, dimid;
  size_t nlevel, nz, len_id;
  size_t start[] = {0, 0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1, 1};



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
    Error(WARNING, routineName, messageStr);
  }

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
 

  /* --- Get variable id for collision rates --- */
  if ((ierror = nc_inq_varid(ncid, COLL_NAME, &varid ))) 
    ERR(ierror,routineName); 

  /* --- Read data --- */
  start[2] = mpi.ix;
  start[3] = mpi.iy;
  start[4] = mpi.zcut;
  count[0] = atom->Nlevel;
  count[1] = atom->Nlevel;
  count[4] = atmos.Nspace;

  /* read as double, although it is written as float */
  if ((ierror = nc_get_vara_double(ncid, varid, start, count, atom->C[0] )))
    ERR(ierror,routineName);


  return;
}
/* ------- end   -------------------------- readCollRates_p.c  --- */

/* ------- begin -------------------------- readDamping_p.c    --- */
void readDamping_p(Atom *atom) {

  const char routineName[] = "readDamping_p";
  char    group_name[ARR_STRLEN], *atmosID;

  int ierror, ncid, varid, dimid;
  size_t nz, len_id;
  size_t start[] = {0, 0, 0};
  size_t count[] = {1, 1, 1};

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
 

  /* --- Get variable id for vbroad --- */
  /* Reads only vbroad, as damp is not stored in any structure */
  if ((ierror = nc_inq_varid(ncid, COLL_NAME, &varid ))) 
    ERR(ierror,routineName); 

  /* --- Read data --- */
  start[0] = mpi.ix;
  start[1] = mpi.iy;
  start[2] = mpi.zcut;
  count[2] = atmos.Nspace;

  /* read as double, although it is written as float */
  if ((ierror = nc_get_vara_double(ncid, varid, start, count, atom->vbroad )))
    ERR(ierror,routineName);


  return;
}
/* ------- end   -------------------------- readDamping_p.c    --- */


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
