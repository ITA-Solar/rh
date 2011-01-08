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
extern MPI_data mpi;
extern IO_data io; 

/* ------- begin --------------------------   init_ncdf_aux.c     --- */
void init_ncdf_aux(void)
/* Creates the netCDF file for the auxiliary data */
{
  const char routineName[] = "init_ncdf_aux";
  int     ierror, ncid, ncid_atom, ncid_op, i, nx_id, ny_id, nspace_id, 
          nrays_id, nwad_id, nwai_id, nlevel_id, nline_id, ncont_id, dimids[5];
  int     nwave, *ai_idx, *ad_idx, nai, nad;
  double *wave, *lambda_air;
  char    group_name[ARR_STRLEN];
  FILE   *test;
  Atom   *atom;
  ActiveSet *as;

  /* Check if we can open the file */
  if ((test = fopen(AUX_FILE, "w")) == NULL) {
    sprintf(messageStr, "Unable to open spectrum output file %s", AUX_FILE);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(test);
  }

  /* Create the file  */
  if ((ierror = nc_create_par(AUX_FILE, NC_NETCDF4 | NC_CLOBBER | NC_MPIIO, 
			      mpi.comm, mpi.info, &ncid))) ERR(ierror,routineName);
  
  
  /* --- Definitions for the root group --- */
  /* dimensions */
  if ((ierror = nc_def_dim(ncid, "nx",       mpi.nx,           &nx_id    ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid, "ny",       mpi.ny,           &ny_id    ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid, "nz",       atmos.Nspace,     &nspace_id))) 
    ERR(ierror,routineName);

  /* variables*/
  // None?

  /* attributes */
  if ((ierror = nc_put_att_text(ncid, NC_GLOBAL, "atmosID", strlen(atmos.ID),
				atmos.ID ))) ERR(ierror,routineName);


  /* Create arrays for multiple-atom output */
  io.aux_atom_ncid   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_pop    = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_poplte = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RijL   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RjiL   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RijC   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_RjiC   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_coll   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  io.aux_atom_damp   = (int *) malloc(atmos.Nactiveatom * sizeof(int));
  
  /* --- Group loop over active atoms --- */
  // Should create groups for active molecules as well!
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
    
    /* --- attributes --- */
    // ?


  } /* end active atoms loop */


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
      if ( containsPolarized(as) || (containsPRDline(as) && input.PRD_angle_dep) 
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

  // MUST CREATE INPUT OPTION TO WRITE (OR NOT) THE OPACITY

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


  /* End define mode */
  if ((ierror = nc_enddef(ncid))) ERR(ierror,routineName);

  /* Copy stuff to the IO data struct */
  io.aux_ncid    = ncid;
  io.aux_op_ncid = ncid_op;


  return;
}
/* ------- end   --------------------------   init_ncdf_aux.c     --- */

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
  free(io.aux_atom_coll);
  free(io.aux_atom_damp);

  return; 
}
/* ------- end   --------------------------   close_ncdf_aux.c --- */

/* ------- begin --------------------------   writeAux_p.c     --- */
void writeAux_p(void) {
// this will write: populations, radrates, coll, damping
  const char routineName[] = "writeAux_p";
  int      ierror, ncid, nact, kr;
  size_t   start[] = {0, 0, 0, 0, 0};
  size_t   count[] = {1, 1, 1, 1, 1};
  double **adamp;
  Atom    *atom;
  AtomicLine      *line;
  AtomicContinuum *continuum;

  /* --- Main atom loop --- */
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    ncid = io.aux_atom_ncid[nact];

    /* --- write populations --- */
    start[1] = mpi.ix;
    start[2] = mpi.iy;
    count[0] = atom->Nlevel; 
    count[3] = atmos.Nspace;
    if (atom->n != NULL) 
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_pop[nact], start, count,
     			     atom->n[0] )))     ERR(ierror,routineName);
    if (atom->nstar != NULL)
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_poplte[nact], start, count,
			     atom->nstar[0] ))) ERR(ierror,routineName);
      
    /* --- write damping --- */
    adamp = matrix_double(atom->Nline, atmos.Nspace);
    
    for (kr=0; kr < atom->Nline; kr++) {
      if (atom->line[kr].Voigt)
	Damping(&atom->line[kr], adamp[kr]);
    }
    count[0] = atom->Nline;
    if ((ierror=nc_put_vara_double(ncid, io.aux_atom_damp[nact], start, count,
				   adamp[0] ))) ERR(ierror,routineName);
    freeMatrix((void **) adamp);

    /* --- write radiative rates --- */
    count[0] = 1;
    /* for bound-bound transitions */
    for (kr=0; kr < atom->Nline; kr++) {
      start[0] = kr;
      line = &atom->line[kr];
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RijL[nact], start, count,
				     line->Rij ))) ERR(ierror,routineName);
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RjiL[nact], start, count,
				     line->Rji ))) ERR(ierror,routineName);
    }
    /* for bound-free transitions */
    for (kr=0; kr < atom->Ncont; kr++) {
      start[0] = kr;
      continuum = &atom->continuum[kr];
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RijC[nact], start, count,
				     continuum->Rij ))) ERR(ierror,routineName);
      if ((ierror=nc_put_vara_double(ncid, io.aux_atom_RjiC[nact], start, count,
				     continuum->Rji ))) ERR(ierror,routineName);
    }

    /* --- write collisional rates --- */
    start[0] = 0;
    start[1] = 0;
    start[2] = mpi.ix;
    start[3] = mpi.iy;
    count[0] = atom->Nlevel;
    count[1] = atom->Nlevel;
    count[2] = 1;
    count[3] = 1;
    count[4] = atmos.Nspace;
    if ((ierror=nc_put_vara_double(ncid, io.aux_atom_coll[nact], start, count,
				   atom->C[0] ))) ERR(ierror,routineName);

    /* free atom (it will be reallocated for each task) */
    //freeAtom(atom);

  } /* End atom loop */

  /* Now write opacity (also goes on aux file) */
  writeOpacity_p();

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
  start_ad[2] = mpi.ix;
  start_ad[3] = mpi.iy;

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
	PRD_angle_dep = (containsPRDline(as) && input.PRD_angle_dep);
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




// LOW PRIORITY:
/* ------- begin -------------------------- writeAtom_metadata_p.c */
// this will write the same as writeAtom() for active atoms
/* ------- begin -------------------------- writeAtom_metadata_p.c */
