/* ------- file: -------------------------- writeray.c ----- --------

       Version:       rh2.0
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Thu Jan 13 16:35:11 2011 --

       --------------------------                      ----------RH-- */

/* --- Writes data from solveray into output file  --  -------------- */

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

void init_ncdf_ray_new(void);
void init_ncdf_ray_old(void);

/* --- Global variables --                             -------------- */

extern enum Topology topology;

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern Geometry geometry;
extern char messageStr[];
extern NCDF_Atmos_file infile;
extern MPI_data mpi;
extern IO_data io; 


/* ------- begin --------------------------   init_ncdf_ray.c  ------ */
void init_ncdf_ray(void) {
  /* Wrapper to find out if we should use old file or create new one. */
  
  if (input.p15d_rerun) init_ncdf_ray_old(); else init_ncdf_ray_new();
  
  return;
}
/* ------- end   --------------------------   init_ncdf_ray.c  ------ */

/* ------- begin -------------------------- init_ncdf_ray_new.c ----- */
void init_ncdf_ray_new(void)
/* Creates the netCDF file for the ray */
{
  const char routineName[] = "init_ncdf_ray_new";
  int     ierror, ncid, nx_id, ny_id, nspect_id, wave_var, wave_sel_id, 
          nspace_id, intensity_var, stokes_u_var, stokes_q_var, stokes_v_var,
          //chi_l_var, eta_l_var, chi_c_var, eta_c_var, sca_c_var, 
          chi_var, S_var, wave_idx_var, dimids[4];
  bool_t  write_xtra;
  FILE   *test;
  double *lambda_air;
  char    timestr[ARR_STRLEN];
  time_t  curtime;
  struct tm *loctime;

  
  write_xtra = (io.ray_nwave_sel > 0);

  /* Check if we can open the file */
  if ((test = fopen(RAY_FILE, "w")) == NULL) {
    sprintf(messageStr, "Unable to open spectrum output file %s", SPEC_FILE);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(test);
  }

  /* Create the file  */
  if ((ierror = nc_create_par(RAY_FILE, NC_NETCDF4 | NC_CLOBBER | NC_MPIPOSIX, 
			      mpi.comm, mpi.info, &ncid))) ERR(ierror,routineName);
  
  /* Global attributes */
  if ((ierror = nc_put_att_text(ncid, NC_GLOBAL, "atmosID", strlen(atmos.ID),
				atmos.ID ))) ERR(ierror,routineName);

  if ((ierror = nc_put_att_int(ncid, NC_GLOBAL, "snapshot_number", NC_USHORT, 1,
			       &mpi.snap_number))) ERR(ierror,routineName);
  
  if ((ierror = nc_put_att_text(ncid, NC_GLOBAL, "svn_id", strlen(mpi.svn_id),
				mpi.svn_id ))) ERR(ierror,routineName);

  

  /* Create dimensions */ 
  if ((ierror = nc_def_dim( ncid, "nx",     mpi.nx,          &nx_id      ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim( ncid, "ny",     mpi.ny,          &ny_id      ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim( ncid, "nz",     infile.nz,       &nspace_id  ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim( ncid, "nwave",  spectrum.Nspect, &nspect_id  ))) 
    ERR(ierror,routineName);
  if (write_xtra) {
    if ((ierror = nc_def_dim( ncid, WAVE_SEL, io.ray_nwave_sel, &wave_sel_id)))
      ERR(ierror,routineName);
  }
  
  /* Create variables */
  dimids[0] = nx_id;
  dimids[1] = ny_id;
  dimids[2] = nspect_id;

  /* Intensity */
  if ((ierror = nc_def_var( ncid,  INT_NAME, NC_FLOAT, 3, dimids, &intensity_var)))
    ERR(ierror,routineName);

  /* Other Stokes parameters, if available */
  if (atmos.Stokes || input.backgr_pol) {
    if ((ierror = nc_def_var( ncid, STOKES_Q, NC_FLOAT, 3, dimids, &stokes_q_var)))
      ERR(ierror,routineName); 
    if ((ierror = nc_def_var( ncid, STOKES_U, NC_FLOAT, 3, dimids, &stokes_u_var)))
      ERR(ierror,routineName); 
    if ((ierror = nc_def_var( ncid, STOKES_V, NC_FLOAT, 3, dimids, &stokes_v_var)))
      ERR(ierror,routineName); 
  }


  if (write_xtra) {

    dimids[0] = nx_id;
    dimids[1] = ny_id;
    dimids[2] = nspace_id;
    dimids[3] = wave_sel_id;
    /* Source function, opacity and emissivity, line and continuum */
    if ((ierror = nc_def_var( ncid, CHI_NAME,  NC_FLOAT, 4, dimids, &chi_var)))
      ERR(ierror,routineName); 
    if ((ierror = nc_def_var( ncid, S_NAME,    NC_FLOAT, 4, dimids, &S_var  )))
      ERR(ierror,routineName); 


    /* Tiago: these take too much space, so not writing them at the moment
    if ((ierror = nc_def_var( ncid, CHI_L_NAME, NC_FLOAT, 4, dimids, &chi_l_var)))
      ERR(ierror,routineName); 
    if ((ierror = nc_def_var( ncid, ETA_L_NAME, NC_FLOAT, 4, dimids, &eta_l_var)))
      ERR(ierror,routineName); 
    if ((ierror = nc_def_var( ncid, CHI_C_NAME, NC_FLOAT, 4, dimids, &chi_c_var)))
      ERR(ierror,routineName); 
    if ((ierror = nc_def_var( ncid, ETA_C_NAME, NC_FLOAT, 4, dimids, &eta_c_var)))
      ERR(ierror,routineName);
    if ((ierror = nc_def_var( ncid, SCA_C_NAME, NC_FLOAT, 4, dimids, &sca_c_var)))
      ERR(ierror,routineName); 
    */
  }


  /* Array with wavelengths */
  dimids[0] = nspect_id;
  if ((ierror = nc_def_var(ncid, WAVE_NAME, NC_DOUBLE, 1, dimids,&wave_var)))
    ERR(ierror,routineName);

  if (write_xtra) {
    /* Array with selected wavelength indices */
    dimids[0] = wave_sel_id;
    if ((ierror = nc_def_var(ncid, WAVE_SEL_IDX, NC_INT, 1, dimids,&wave_idx_var)))
      ERR(ierror,routineName);
  }

  /* --- Write attributes --- */
  /* Time of creation in ISO 8601 */
  curtime = time(NULL);
  loctime = localtime(&curtime);
  strftime(timestr, ARR_STRLEN, "%Y-%m-%dT%H:%M:%S%z", loctime);
  
  if ((ierror = nc_put_att_text(ncid, NC_GLOBAL, "creation_time", strlen(timestr),
      (const char *) &timestr )))     ERR(ierror,routineName); 

  /*  units  */
  if ((ierror = nc_put_att_text( ncid, wave_var,        "units",  2,
                         "nm" )))                      ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid, intensity_var,   "units",  23,
                         "J s^-1 m^-2 Hz^-1 sr^-1" ))) ERR(ierror,routineName);

  if (atmos.Stokes || input.backgr_pol) {
    if ((ierror = nc_put_att_text( ncid, stokes_q_var,  "units",  23,
                         "J s^-1 m^-2 Hz^-1 sr^-1" ))) ERR(ierror,routineName);
    if ((ierror = nc_put_att_text( ncid, stokes_u_var,  "units",  23,
                         "J s^-1 m^-2 Hz^-1 sr^-1" ))) ERR(ierror,routineName);
    if ((ierror = nc_put_att_text( ncid, stokes_v_var,  "units",  23,
                         "J s^-1 m^-2 Hz^-1 sr^-1" ))) ERR(ierror,routineName);
  }
  if (write_xtra) {
    if ((ierror = nc_put_att_text( ncid, S_var,         "units",  23,
			 "J s^-1 m^-2 Hz^-1 sr^-1" ))) ERR(ierror,routineName);
    if ((ierror = nc_put_att_text( ncid, S_var,       "description",    40,
        "Total source function (line + continuum)" ))) ERR(ierror,routineName);
    if ((ierror = nc_put_att_text( ncid, chi_var,       "units",  4,
			 "m^-1" ))) ERR(ierror,routineName);
    if ((ierror = nc_put_att_text( ncid, chi_var,       "description",  35,
             "Total absorption (line + continuum)" ))) ERR(ierror,routineName);
  }

  /* End define mode */
  if ((ierror = nc_enddef(ncid))) ERR(ierror,routineName);


  /* Copy stuff to the IO data struct */
  io.ray_ncid         = ncid;
  io.ray_wave_var     = wave_var;
  io.ray_wave_idx_var = wave_idx_var;
  io.ray_int_var      = intensity_var;
  io.ray_stokes_q_var = stokes_q_var;
  io.ray_stokes_u_var = stokes_u_var;
  io.ray_stokes_v_var = stokes_v_var;
  io.ray_chi_var      = chi_var;
  io.ray_S_var        = S_var;
  /*
  io.ray_chi_l_var    = chi_l_var;
  io.ray_eta_l_var    = eta_l_var;
  io.ray_chi_c_var    = chi_c_var;
  io.ray_sca_c_var    = sca_c_var;
  */


  /* --- Write wavelength and wavelength indices --- */

  /* Write wavelength */
  if (spectrum.vacuum_to_air) {
    lambda_air = (double *) malloc(spectrum.Nspect * sizeof(double));
    vacuum_to_air(spectrum.Nspect, spectrum.lambda, lambda_air);
    if ((ierror = nc_put_var_double(io.ray_ncid, io.ray_wave_var, lambda_air )))
      ERR(ierror,routineName);
    free(lambda_air);
  } else
    if ((ierror = nc_put_var_double(io.ray_ncid, io.ray_wave_var, 
				    spectrum.lambda ))) ERR(ierror,routineName);


  /* Write wavelength indices */
  if (write_xtra) 
    if ((ierror = nc_put_var_int(io.ray_ncid, io.ray_wave_idx_var, 
				 io.ray_wave_idx ))) ERR(ierror,routineName);
    
  return; 
}

/* ------- end   -------------------------- init_ncdf_ray_new.c ------- */


/* ------- begin -------------------------- init_ncdf_ray_old.c ------ */
void init_ncdf_ray_old(void)
/* Opens an existing ray netCDF file */
{
  const char routineName[] = "init_ncdf_ray_old";
  int     ierror, ncid;
  bool_t  write_xtra;
  size_t  len_id;
  FILE   *test;
  char   *atmosID;

  
  write_xtra = (io.ray_nwave_sel > 0);

  /* Check if we can open the file */
  if ((test = fopen(RAY_FILE, "a")) == NULL) {
    sprintf(messageStr, "Unable to open spectrum output file %s", SPEC_FILE);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(test);
  }

  /* Open the file  */
  if ((ierror = nc_open_par(RAY_FILE, NC_WRITE | NC_MPIPOSIX, 
			      mpi.comm, mpi.info, &ncid))) ERR(ierror,routineName);
  io.ray_ncid = ncid;
  
  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if ((ierror = nc_inq_attlen(ncid, NC_GLOBAL, "atmosID", &len_id ))) 
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

  /* --- Get variable IDs ---*/
  if ((ierror = nc_inq_varid(ncid, INT_NAME,  &io.ray_int_var ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(ncid, WAVE_NAME, &io.ray_wave_var))) 
      ERR(ierror,routineName);
      
  if (atmos.Stokes || input.backgr_pol) {
    if ((ierror = nc_inq_varid(ncid, STOKES_Q, &io.ray_stokes_q_var))) 
        ERR(ierror,routineName);
    if ((ierror = nc_inq_varid(ncid, STOKES_U, &io.ray_stokes_u_var))) 
        ERR(ierror,routineName);
    if ((ierror = nc_inq_varid(ncid, STOKES_V, &io.ray_stokes_v_var))) 
        ERR(ierror,routineName);
  }
  
  if (write_xtra) {
    if ((ierror = nc_inq_varid(ncid, CHI_NAME,     &io.ray_chi_var     ))) 
        ERR(ierror,routineName);
    if ((ierror = nc_inq_varid(ncid, S_NAME,       &io.ray_S_var       ))) 
        ERR(ierror,routineName);
    if ((ierror = nc_inq_varid(ncid, WAVE_SEL_IDX, &io.ray_wave_idx_var))) 
        ERR(ierror,routineName);
  }
  
  return; 
}

/* ------- end   -------------------------- init_ncdf_ray_old.c ------ */

/* ------- begin -------------------------- close_ncdf_ray.c --------- */
void close_ncdf_ray(void)
/* Closes the spec netCDF file */ 
{
  const char routineName[] = "close_ncdf_ray";
  int        ierror;

  if ((ierror = nc_close(io.ray_ncid))) ERR(ierror,routineName);
  return; 
}
/* ------- end   -------------------------- close_ncdf_ray.c --------- */


/* ------- begin -------------------------- writeRay.c --------------- */
void writeRay(void)
/* Writes ray data to netCDF file. */
{
  const char routineName[] = "writeRay";
  int        ierror, idx, ncid, k, l, nspect;
  double    *sca, *J;
  float    **chi, **S;
  size_t     start[] = {0, 0, 0, 0};
  size_t     count[] = {1, 1, 1, 1};
  bool_t     write_xtra, crosscoupling, to_obs, initialize;
  ActiveSet *as;

  write_xtra = (io.ray_nwave_sel > 0);
  ncid = io.ray_ncid;

  start[0] = mpi.ix; 
  start[1] = mpi.iy;

  /* Write intensity */
  count[2] = spectrum.Nspect;
  if ((ierror = nc_put_vara_double(ncid, io.ray_int_var, start, count,
				   spectrum.I[0] ))) ERR(ierror,routineName);
  
  if (atmos.Stokes || input.backgr_pol) { /* Write rest of Stokes vector */
    if ((ierror = nc_put_vara_double(ncid, io.ray_stokes_q_var, start,
		     count, spectrum.Stokes_Q[0] ))) ERR(ierror,routineName);
    if ((ierror = nc_put_vara_double(ncid, io.ray_stokes_u_var, start,
		     count, spectrum.Stokes_U[0] ))) ERR(ierror,routineName);
    if ((ierror = nc_put_vara_double(ncid, io.ray_stokes_v_var, start,
		     count, spectrum.Stokes_V[0] ))) ERR(ierror,routineName);

  }


  if (write_xtra) {
    /* Write opacity and emissivity for line and continuum */

    /* compute scattering emissivity (from continuum) */
    // make chi and S matrix float arrays, write into netcdf in one step
    // and not for each wavelength!!!
    sca = (double *) malloc(atmos.Nspace * sizeof(double));
    chi = matrix_float(infile.nz, io.ray_nwave_sel);
    S   = matrix_float(infile.nz, io.ray_nwave_sel);

    if (input.limit_memory) J = (double *) malloc(atmos.Nspace * sizeof(double));

    for (nspect = 0; nspect < io.ray_nwave_sel; nspect++) {
      idx = io.ray_wave_idx[nspect];
      as  = &spectrum.as[idx];
      
      alloc_as(idx, crosscoupling=FALSE);
      Opacity(idx, 0, to_obs=TRUE, initialize=TRUE);
      readBackground(idx, 0, to_obs=TRUE);
      
      if (input.limit_memory) {
	readJlambda_single(idx, J);
      } else
	J = spectrum.J[idx];

      /* Zero S and chi */
      for (k = 0; k < infile.nz; k++) {
	S[k][nspect]   = 0.0;
	chi[k][nspect] = 0.0;
      }

      for (k = 0;  k < atmos.Nspace;  k++) {
	l = k + mpi.zcut;
	sca[k] = as->sca_c[k]*J[k];
	chi[l][nspect] = (float) (as->chi[k] + as->chi_c[k]);
	S[l][nspect]   = (float) ((as->eta[k] + as->eta_c[k] + sca[k]) /
				  (as->chi[k] + as->chi_c[k]));
      }

      /* Write variables 
      if ((ierror=nc_put_vara_double(ncid, io.ray_chi_var, start, count,
				     chi ))) ERR(ierror,routineName);
      if ((ierror=nc_put_vara_double(ncid, io.ray_S_var,   start, count,
				     S )))   ERR(ierror,routineName);
      */
      /* Tiago: these take too much space, so not writing them at the moment
      if ((ierror=nc_put_vara_double(ncid, io.ray_chi_l_var, start, count,
				     as->chi )))   ERR(ierror,routineName);
      if ((ierror=nc_put_vara_double(ncid, io.ray_eta_l_var, start, count,
				     as->eta )))   ERR(ierror,routineName);
      if ((ierror=nc_put_vara_double(ncid, io.ray_chi_c_var, start, count,
				     as->chi_c ))) ERR(ierror,routineName);
      if ((ierror=nc_put_vara_double(ncid, io.ray_eta_c_var, start, count,
				     as->eta_c ))) ERR(ierror,routineName);
      if ((ierror=nc_put_vara_double(ncid, io.ray_sca_c_var, start, count,
				     sca )))       ERR(ierror,routineName);     
      */
    }
    
    /* Write variables */
    start[0] = mpi.ix;  count[0] = 1;
    start[1] = mpi.iy;  count[1] = 1;
    start[2] = 0;       count[2] = infile.nz;
    start[3] = 0;       count[3] = io.ray_nwave_sel;


    if ((ierror=nc_put_vara_float(ncid, io.ray_chi_var, start, count,
				  chi[0] ))) ERR(ierror,routineName);
    if ((ierror=nc_put_vara_float(ncid, io.ray_S_var,   start, count,
				  S[0] )))   ERR(ierror,routineName);
    
    
    free(sca); freeMatrix((void **) chi); freeMatrix((void **) S);
    if (input.limit_memory) free(J);
  }


  return;
}
/* ------- end   -------------------------- writeRay.c -------------- */

/* ------- begin -------------------------- calculate_ray.c --------- */
void calculate_ray(void) {
  /* performs necessary reallocations and inits, and solves for ray */
  
      int i, nact, ierror;
      bool_t analyze_output, equilibria_only;
      Atom *atom;
      AtomicLine *line;
      
      
          /*
    atom = atmos.activeatoms[0];
    for (i=0; i < atmos.Nspace; i++) printf("%e  %e  %e  %e\n", atom->nstar[0][i], atom->nstar[1][i],
					    atom->nstar[2][i], atom->nstar[3][i]);  
    */
      
  
      /* Must calculate background opacity for new ray, need some prep first */
      for (nact = 0; nact < atmos.Nactiveatom; nact++) {
	atom = atmos.activeatoms[nact];
	
        /* Rewind atom files to point just before collisional data */
        if ((ierror = fseek(atom->fp_input, io.atom_file_pos[nact], SEEK_SET))) {
          sprintf(messageStr, "Unable to rewind atom file for %s", atom->ID);
          Error(ERROR_LEVEL_2, "rh15d_ray", messageStr);
	}
	
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
      
      
      Background_p(analyze_output=FALSE, equilibria_only=FALSE);
        
      /*
      as = &spectrum.as[176];
      alloc_as(176, FALSE);
      readBackground(176, 0, 0);
      Opacity(176, 0, 0, TRUE);
      
      for (i=0; i < atmos.Nspace; i++) {
        printf("%e   %e   %e   %e   %e   %e\n", as->chi_c[i], as->eta_c[i],
       as->chi[i], as->eta[i], as->sca_c[i], spectrum.J[176][i]);	     
      } */

      /* --- Solve radiative transfer for ray --           -------------- */
      solveSpectrum(FALSE, FALSE);

      /*
      for (i=0; i < spectrum.Nspect; i++) {
	printf("%i   %e\n", i, spectrum.I[i][0]);
      } */
    
  
}

/* ------- end   -------------------------- calculate_ray.c --------- */
