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
extern NCDF_Atmos_file infile;
extern MPI_data mpi;
extern IO_data io; 

/* ------- begin --------------------------   init_ncdf_indata.c  --- */
void init_ncdf_indata(void) {
  /* Wrapper to find out if we should use old file or create new one */
  
  if (input.p15d_rerun) init_ncdf_indata_old(); else init_ncdf_indata_new();
  
  return;
}
/* ------- end   --------------------------   init_ncdf_indata.c  --- */


/* ------- begin --------------------------   init_ncdf_indata.c  --- */
void init_ncdf_indata_new(void)
/* Creates the netCDF file for the input data */
{
  const char routineName[] = "init_ncdf_indata_new";
  int     i, ierror, ncid, ncid_input, ncid_atmos, ncid_mpi, nx_id, ny_id, 
          nspace_id, nhydr_id, nelem_id, nrays_id, nproc_id, atstr_id, niter_id,
          temp_var, ne_var, vz_var, vturb_var, Bx_var, By_var, Bz_var, nh_var, 
          /* B_var, gB_var, chiB_var, */
          ew_var, ab_var, eid_var, mu_var, wmu_var, height_var,x_var, y_var,
          xnum_var, ynum_var, tm_var, tn_var, it_var, conv_var, dm_var, dmh_var, 
          zch_var, ntsk_var, host_var, st_var, ft_var, z_varid, dimids[4];
  /* This value is harcoded for efficiency. Maximum number of iterations ever needed */
  int     NMaxIter = 500; 
  long    task;
  size_t  start[] = {0, 0, 0};
  size_t  count[] = {1, 1, 1};
  bool_t  PRD_angle_dep, XRD;
  double *height;
  char    startJ[MAX_LINE_SIZE], StokesMode[MAX_LINE_SIZE], angleSet[MAX_LINE_SIZE],
          hostname[ARR_STRLEN], timestr[ARR_STRLEN];
  FILE   *test;
  time_t  curtime;
  struct tm *loctime;

  /* Check if we can open the file */
  if ((test = fopen(INPUTDATA_FILE, "w")) == NULL) {
    sprintf(messageStr, "Unable to open spectrum output file %s", INPUTDATA_FILE);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(test);
  }

  /* Create the file  */
  if ((ierror = nc_create_par(INPUTDATA_FILE, NC_NETCDF4 | NC_CLOBBER | NC_MPIPOSIX, 
  //if ((ierror = nc_create_par(INPUTDATA_FILE, NC_NETCDF4 | NC_CLOBBER | NC_MPIIO, 
			      mpi.comm, mpi.info, &ncid))) ERR(ierror,routineName);

  /* Create groups */
  if ((ierror = nc_def_grp(ncid, "input", &ncid_input))) ERR(ierror,routineName);
  if ((ierror = nc_def_grp(ncid, "atmos", &ncid_atmos))) ERR(ierror,routineName);
  if ((ierror = nc_def_grp(ncid, "mpi",   &ncid_mpi)))   ERR(ierror,routineName);

  /* --- Definitions for the root group --- */ 
  /* dimensions */
  if ((ierror = nc_def_dim(ncid, "nx",     mpi.nx,       &nx_id    ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid, "ny",     mpi.ny,       &ny_id    ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid, "nz",     infile.nz,   &nspace_id))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid, "strlen", ARR_STRLEN,   &atstr_id ))) 
    ERR(ierror,routineName);

  /* attributes */
  if ((ierror = nc_put_att_text(ncid, NC_GLOBAL, "atmosID", strlen(atmos.ID),
				atmos.ID ))) ERR(ierror,routineName);

  /* --- Definitions for the INPUT group --- */
  /* attributes */
  PRD_angle_dep = (input.PRD_angle_dep  &&  atmos.NPRDactive > 0);
  XRD           = (input.XRD  &&  atmos.NPRDactive > 0);

  if ((ierror=nc_put_att_ubyte(ncid_input,NC_GLOBAL,"Magneto_optical",NC_UBYTE, 1,
               (unsigned char *) &input.magneto_optical ))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_ubyte(ncid_input,NC_GLOBAL, "PRD_angle_dep", NC_UBYTE, 1,
               (unsigned char *) &PRD_angle_dep )))    ERR(ierror,routineName);
  if ((ierror=nc_put_att_ubyte(ncid_input,NC_GLOBAL, "XRD",           NC_UBYTE, 1,
               (unsigned char *) &XRD )))              ERR(ierror,routineName);
  if ((ierror=nc_put_att_ubyte(ncid_input,NC_GLOBAL, "Background_polarization",
    NC_UBYTE,1,(unsigned char *) &input.backgr_pol ))) ERR(ierror,routineName);
  
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
  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL, "Start_J", strlen(startJ),
			      startJ )))  ERR(ierror,routineName);

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
  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL,"Stokes_mode",strlen(StokesMode),
			      StokesMode )))  ERR(ierror,routineName);

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
  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL,"Angle_set",strlen(angleSet),
			      angleSet )))  ERR(ierror,routineName);

  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL,"Atmos_file",
       strlen(input.atmos_input), input.atmos_input))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL,"Abundances_file",
       strlen(input.abund_input), input.abund_input))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_text(ncid_input,NC_GLOBAL,"Kurucz_PF_data",
       strlen(input.pfData),      input.pfData)))  ERR(ierror,routineName);

  if ((ierror=nc_put_att_double( ncid_input, NC_GLOBAL, "Iteration_limit", 
                NC_DOUBLE, 1, &input.iterLimit )))    ERR(ierror,routineName);
  if ((ierror=nc_put_att_double( ncid_input, NC_GLOBAL, "PRD_Iteration_limit", 
                NC_DOUBLE, 1, &input.PRDiterLimit ))) ERR(ierror,routineName);

  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "N_max_iter",    NC_INT, 1,
			      &input.NmaxIter )))     ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "Ng_delay",      NC_INT, 1,
			      &input.Ngdelay )))      ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "Ng_order",      NC_INT, 1,
			      &input.Ngorder )))      ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "Ng_period",     NC_INT, 1,
			      &input.Ngperiod )))     ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "PRD_N_max_iter",NC_INT, 1,
			      &input.PRD_NmaxIter ))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "PRD_Ng_delay",  NC_INT, 1,
			      &input.PRD_Ngdelay )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "PRD_Ng_order",  NC_INT, 1,
			      &input.PRD_Ngorder )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_input, NC_GLOBAL, "PRD_Ng_period", NC_INT, 1,
			      &input.PRD_Ngperiod ))) ERR(ierror,routineName);

  if ((ierror=nc_put_att_double(ncid_input, NC_GLOBAL, "Metallicity", NC_DOUBLE, 1,
				&input.metallicity ))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_double(ncid_input, NC_GLOBAL, "Lambda_reference",
                   NC_DOUBLE, 1, &atmos.lambda_ref ))) ERR(ierror,routineName);


  /* --- Definitions for the ATMOS group --- */
  /* dimensions */
  if ((ierror = nc_def_dim(ncid_atmos, "nhydr",     atmos.H->Nlevel, &nhydr_id )))
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid_atmos, "nelements", atmos.Nelem,     &nelem_id ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid_atmos, "nrays",     geometry.Nrays,  &nrays_id ))) 
    ERR(ierror,routineName);

  /* variables*/
  dimids[0] = nx_id;
  dimids[1] = ny_id;
  dimids[2] = nspace_id;
  if ((ierror = nc_def_var(ncid_atmos, "temperature",        NC_FLOAT, 3, dimids,
			   &temp_var)))  ERR(ierror,routineName);
  /*
  if ((ierror = nc_def_var(ncid_atmos, "electron_density",   NC_DOUBLE, 3, dimids,
			   &ne_var)))    ERR(ierror,routineName);
  */
  if ((ierror = nc_def_var(ncid_atmos, "velocity_z",         NC_FLOAT, 3, dimids,
			   &vz_var)))    ERR(ierror,routineName);
  /* No need to write vturb, at least for now... 
  if ((ierror = nc_def_var(ncid_atmos, "velocity_turbulent", NC_FLOAT, 3, dimids,
			   &vturb_var))) ERR(ierror,routineName);
  */
  if ((ierror = nc_def_var(ncid_atmos, "height",              NC_FLOAT, 3, dimids,
			   &height_var))) ERR(ierror,routineName);
  if (atmos.Stokes) {
    /* New definitions (Bx, By, Bz) 
    if ((ierror=nc_def_var(ncid_atmos, "Bx", NC_FLOAT, 3, dimids,
			   &Bx_var)))  ERR(ierror,routineName); 
    if ((ierror=nc_put_att_text( ncid_atmos, Bx_var, "units", 1, "T" ))) 
                                         ERR(ierror,routineName);
    if ((ierror=nc_def_var(ncid_atmos, "By", NC_FLOAT, 3, dimids,
			   &By_var)))  ERR(ierror,routineName); 
    if ((ierror=nc_put_att_text( ncid_atmos, By_var, "units", 1, "T" ))) 
                                         ERR(ierror,routineName);
    if ((ierror=nc_def_var(ncid_atmos, "Bz", NC_FLOAT, 3, dimids,
			   &Bz_var)))  ERR(ierror,routineName); 
    if ((ierror=nc_put_att_text( ncid_atmos, Bz_var, "units", 1, "T" ))) 
                                         ERR(ierror,routineName);
    */

    /* Old definitions (B, gamma_B, chi_B) 
    if ((ierror=nc_def_var(ncid_atmos, "B",                  NC_FLOAT, 3, dimids,
			   &B_var)))     ERR(ierror,routineName); 
    if ((ierror=nc_put_att_text( ncid_atmos, B_var,    "units", 1, "T"   ))) 
                                         ERR(ierror,routineName);
    if ((ierror=nc_def_var(ncid_atmos, "gamma_B",            NC_FLOAT, 3, dimids,
			   &gB_var)))    ERR(ierror,routineName); 
    if ((ierror=nc_put_att_text( ncid_atmos, gB_var,   "units", 3, "rad" ))) 
                                         ERR(ierror,routineName);
    if ((ierror=nc_def_var(ncid_atmos, "chi_B",              NC_FLOAT, 3, dimids,
			   &chiB_var)))  ERR(ierror,routineName); 
    if ((ierror=nc_put_att_text( ncid_atmos, chiB_var, "units", 3, "rad" ))) 
                                         ERR(ierror,routineName);
    */
  }
  dimids[0] = nhydr_id;
  dimids[1] = nx_id;
  dimids[2] = ny_id;
  dimids[3] = nspace_id;
  /*
  if ((ierror = nc_def_var(ncid_atmos, "hydrogen_populations",NC_FLOAT, 4, dimids,
			   &nh_var)))    ERR(ierror,routineName);
  */
  dimids[0] = nelem_id;
  if ((ierror = nc_def_var(ncid_atmos, "element_weight",      NC_DOUBLE,1, dimids,
			   &ew_var)))    ERR(ierror,routineName);  
  if ((ierror = nc_def_var(ncid_atmos, "element_abundance",   NC_DOUBLE,1, dimids,
			   &ab_var)))    ERR(ierror,routineName); 
  dimids[1] = atstr_id;
  if ((ierror = nc_def_var(ncid_atmos, "element_id",          NC_CHAR,  2, dimids,
                           &eid_var)))   ERR(ierror,routineName); 
  /* next 3 came from geometry */
  dimids[0] = nrays_id;
  if ((ierror = nc_def_var(ncid_atmos, "muz",                 NC_DOUBLE,1, dimids,
			   &mu_var)))    ERR(ierror,routineName);  
  if ((ierror = nc_def_var(ncid_atmos, "wmu",                 NC_DOUBLE,1, dimids,
			   &wmu_var)))   ERR(ierror,routineName);
  /*
  dimids[0] = nspace_id;
  if ((ierror = nc_def_var(ncid_atmos, "height",              NC_FLOAT, 1, dimids,
			   &height_var))) ERR(ierror,routineName);
  */
  dimids[0] = nx_id;
  if ((ierror = nc_def_var(ncid_atmos, "x",                   NC_FLOAT, 1, dimids,
			   &x_var)))      ERR(ierror,routineName);
  dimids[0] = ny_id;
  if ((ierror = nc_def_var(ncid_atmos, "y",                   NC_FLOAT, 1, dimids,
			   &y_var)))      ERR(ierror,routineName);

  /* attributes */
  if ((ierror = nc_put_att_ubyte(ncid_atmos, NC_GLOBAL, "moving", NC_UBYTE, 1,
		     (unsigned char *) &atmos.moving ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_ubyte(ncid_atmos, NC_GLOBAL, "stokes", NC_UBYTE, 1,
	             (unsigned char *) &atmos.Stokes ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, temp_var,  "units",  1,
				 "K" ))) ERR(ierror,routineName);
  /*
  if ((ierror = nc_put_att_text( ncid_atmos, ne_var,    "units",  4,
				 "m^-3" ))) ERR(ierror,routineName);
  */
  if ((ierror = nc_put_att_text( ncid_atmos, vz_var,    "units",  6,
				 "m s^-1" ))) ERR(ierror,routineName);
  /*
  if ((ierror = nc_put_att_text( ncid_atmos, vturb_var, "units",  6,
				 "m s^-1" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, nh_var,    "units",  4,
				 "m^-3" ))) ERR(ierror,routineName);
  */
  if ((ierror = nc_put_att_text( ncid_atmos, ew_var,    "units", 17,
				 "atomic_mass_units" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, height_var,"units",  1,
				 "m" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, x_var,     "units",  1,
				 "m" ))) ERR(ierror,routineName);
  if ((ierror = nc_put_att_text( ncid_atmos, y_var,     "units",  1,
				 "m" ))) ERR(ierror,routineName);
  
  /* --- Definitions for the MPI group --- */
  /* dimensions */
  if ((ierror = nc_def_dim(ncid_mpi, "nprocesses",  mpi.size, &nproc_id))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_def_dim(ncid_mpi, "niterations", NMaxIter, &niter_id))) 
    ERR(ierror,routineName);  

  /* variables*/
  dimids[0] = nx_id;
  if ((ierror = nc_def_var(ncid_mpi, XNUM_NAME,   NC_USHORT, 1, dimids,
			   &xnum_var))) ERR(ierror,routineName);
  dimids[0] = ny_id;
  if ((ierror = nc_def_var(ncid_mpi, YNUM_NAME,   NC_USHORT, 1, dimids,
			   &ynum_var))) ERR(ierror,routineName);

  dimids[0] = nx_id;
  dimids[1] = ny_id;
  if ((ierror = nc_def_var(ncid_mpi, TASK_MAP,    NC_LONG, 2, dimids,
			   &tm_var)))   ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, TASK_NUMBER, NC_LONG, 2, dimids,
			   &tn_var)))   ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, ITER_NAME,   NC_LONG, 2, dimids,
			   &it_var)))   ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, CONV_NAME,   NC_LONG,  2, dimids,
			   &conv_var))) ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, DM_NAME,     NC_FLOAT,  2, dimids,
			   &dm_var)))   ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, ZC_NAME,     NC_INT,    2, dimids,
			   &zch_var))) ERR(ierror,routineName);

  dimids[0] = nproc_id;
  dimids[1] = atstr_id;
  if ((ierror = nc_def_var(ncid_mpi, NTASKS,      NC_LONG,   1, dimids,
			   &ntsk_var))) ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, HOSTNAME,    NC_CHAR,   2, dimids,
			   &host_var))) ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, START_TIME,  NC_CHAR,   2, dimids,
			   &st_var))) ERR(ierror,routineName);
  if ((ierror = nc_def_var(ncid_mpi, FINISH_TIME, NC_CHAR,   2, dimids,
			   &ft_var))) ERR(ierror,routineName);

  dimids[0] = nx_id;
  dimids[1] = ny_id;
  dimids[2] = niter_id;
  if ((ierror = nc_def_var(ncid_mpi, DMH_NAME,    NC_FLOAT,  3, dimids,
			   &dmh_var))) ERR(ierror,routineName);


  /* attributes */
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "x_start", NC_INT, 1,
			      &input.p15d_x0 )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "x_end",   NC_INT, 1,
			      &input.p15d_x1 )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "x_step",  NC_INT, 1,
			      &input.p15d_xst ))) ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "y_start", NC_INT, 1,
			      &input.p15d_y0 )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "y_end",   NC_INT, 1,
			      &input.p15d_y1 )))  ERR(ierror,routineName);
  if ((ierror=nc_put_att_int( ncid_mpi, NC_GLOBAL, "y_step",  NC_INT, 1,
			      &input.p15d_yst ))) ERR(ierror,routineName);

  /* --- End define mode --- */
  if ((ierror = nc_enddef(ncid))) ERR(ierror,routineName);


  if (mpi.rank == 0) printf("writeindata end of define\n");
  /* --- Write some data that does not depend on xi, yi, ATMOS group --- */
  /* only the first process needs to write the constant things */
  if (mpi.rank == 0) {
    /* arrays of number of elements */
    count[0] = 1;
    start[1] = 0;   count[1] = ATOM_ID_WIDTH+1;
    for (i=0; i < atmos.Nelem; i++) {
      start[0] = (size_t) i;
      if ((ierror = nc_put_var1_double(ncid_atmos, ew_var,  &start[0], 
		      &atmos.elements[i].weight ))) ERR(ierror,routineName);  
      if ((ierror = nc_put_var1_double(ncid_atmos, ab_var,  &start[0], 
		      &atmos.elements[i].abund )))  ERR(ierror,routineName); 
      if ((ierror = nc_put_vara_text(ncid_atmos,  eid_var, start, count,
       (const char *)&atmos.elements[i].ID )))     ERR(ierror,routineName); 
    }

    /* arrays from geometry */
    if ((ierror = nc_put_var_double(ncid_atmos, mu_var,      geometry.muz )))
      ERR(ierror,routineName);   
    if ((ierror = nc_put_var_double(ncid_atmos, wmu_var,     geometry.wmu )))
      ERR(ierror,routineName);   

    /* Must read full z from NetCDF file */
    /* Not writing z anymore, with dept_refine this is now a 3D variable and is
       written later.
       
    start[0] = input.p15d_nt; count[0] = 1;
    start[1] = 0;             count[1] = infile.nz;
    height = (double *) malloc(infile.nz * sizeof(double));

    if ((ierror=nc_inq_varid(infile.ncid, "z",  &z_varid)))   
      ERR(ierror,routineName);
    if ((ierror = nc_get_vara_double(infile.ncid, z_varid, start, count, height))) 
      ERR(ierror,routineName);
    
    if ((ierror = nc_put_var_double(ncid_atmos, height_var,  height)))
      ERR(ierror,routineName);   

    free(height);
    */

    /* write x and y, convert from MPI indices to actual values */
    for (i=0; i < mpi.nx; i++) {
      start[0] = (size_t) i;
      if ((ierror = nc_put_var1_double(ncid_atmos, x_var,  &start[0], 
	       &infile.x[mpi.xnum[i]]))) ERR(ierror,routineName);  
    }
    
    for (i=0; i < mpi.ny; i++) {
      start[0] = (size_t) i;
      if ((ierror = nc_put_var1_double(ncid_atmos, y_var,  &start[0], 
	       &infile.y[mpi.ynum[i]]))) ERR(ierror,routineName);  
    }
  }

  /* --- Write some data that does not depend on xi, yi, MPI group --- */

  if (mpi.rank == 0) {
  /* xnum, ynum */
    if ((ierror = nc_put_var_int(ncid_mpi, xnum_var, mpi.xnum))) 
      ERR(ierror,routineName);
    if ((ierror = nc_put_var_int(ncid_mpi, ynum_var, mpi.ynum))) 
      ERR(ierror,routineName);
  }

  if (mpi.rank == 0) printf("writeindata end of mpi.rank=0 writes\n");

  /* Number of tasks */
  start[0] = mpi.rank;
  if ((ierror = nc_put_var1_long(ncid_mpi, ntsk_var, start, &mpi.Ntasks )))
    ERR(ierror,routineName);

  /* Hostname of each process */
  if ((ierror = gethostname((char *) &hostname, ARR_STRLEN)) != 0 )
    printf("(EEE) %s: error getting hostname.\n",routineName);

  start[0] = mpi.rank; count[0] = 1;
  start[1] = 0;        count[1] = strlen(hostname);
  if ((ierror = nc_put_vara_text(ncid_mpi, host_var, start, count,
      (const char *) &hostname )))     ERR(ierror,routineName); 

  if (mpi.rank == 0) printf("writeindata end of hostname\n");

  /* Get time in ISO 8601 */
  curtime = time(NULL);
  loctime = localtime(&curtime);
  strftime(timestr, ARR_STRLEN, "%Y-%m-%dT%H:%M:%S%z", loctime);
  
  start[1] = 0; count[1] = strlen(timestr);
  if ((ierror = nc_put_vara_text(ncid_mpi, st_var, start, count,
      (const char *) &timestr )))     ERR(ierror,routineName); 

  if (mpi.rank == 0) printf("writeindata end of time\n");
  /* Write arrays of Ntasks, one value at a time */
  for (task = 0; task < mpi.Ntasks; task++) {

    start[0] = mpi.taskmap[task + mpi.my_start][0];  count[0] = 1;
    start[1] = mpi.taskmap[task + mpi.my_start][1];  count[1] = 1;
    
    /* Task map */
    if ((ierror = nc_put_var1_int(ncid_mpi,  tm_var, start, &mpi.rank )))
      ERR(ierror,routineName);
    /* Task number */
    if ((ierror = nc_put_var1_long(ncid_mpi, tn_var, start, &task )))
      ERR(ierror,routineName);
  }
  
  if (mpi.rank == 0) printf("writeindata end of task maps\n");

  /* --- Copy stuff to the IO data struct --- */
  io.in_ncid       = ncid;
  io.in_input_ncid = ncid_input;
  io.in_atmos_ncid = ncid_atmos;
  io.in_mpi_ncid   = ncid_mpi;
  
  io.in_atmos_T    = temp_var;
  io.in_atmos_ne   = ne_var;
  io.in_atmos_vz   = vz_var;
  io.in_atmos_vt   = vturb_var;
  io.in_atmos_Bx   = Bx_var;
  io.in_atmos_By   = By_var;
  io.in_atmos_Bz   = Bz_var;
  /*
  io.in_atmos_B    = B_var;
  io.in_atmos_gB   = gB_var;
  io.in_atmos_chiB = chiB_var;
  */
  io.in_atmos_nh   = nh_var;
  io.in_atmos_ew   = ew_var;
  io.in_atmos_ab   = ab_var;
  io.in_atmos_eid  = eid_var;
  io.in_atmos_mu   = mu_var;
  io.in_atmos_wmu  = wmu_var;
  io.in_atmos_z    = height_var;
  io.in_atmos_x    = x_var;
  io.in_atmos_y    = y_var;
  
  io.in_mpi_xnum   = xnum_var;
  io.in_mpi_ynum   = ynum_var;
  io.in_mpi_tm     = tm_var;
  io.in_mpi_tn     = tn_var;
  io.in_mpi_it     = it_var;
  io.in_mpi_st     = st_var;
  io.in_mpi_conv   = conv_var;
  io.in_mpi_dm     = dm_var;
  io.in_mpi_dmh    = dmh_var;
  io.in_mpi_zc     = zch_var;
  io.in_mpi_ntsk   = ntsk_var;
  io.in_mpi_host   = host_var;
  io.in_mpi_ft     = ft_var;


  return;
}
/* ------- end   --------------------------   init_ncdf_indata.c  --- */

/* ------- begin --------------------------   init_ncdf_indata_old.c  --- */
void init_ncdf_indata_old(void)
/* Opens an existing NetCDF input data file, loads structures and ids */
{
  const char routineName[] = "init_ncdf_indata_old";
  int     i, ierror, ncid, ncid_input, ncid_atmos, ncid_mpi, nx_id, ny_id, 
          nspace_id, nhydr_id, nelem_id, nrays_id, nproc_id, atstr_id, niter_id,
          temp_var, ne_var, vz_var, vturb_var, Bx_var, By_var, Bz_var, nh_var, 
          /* B_var, gB_var, chiB_var, */
          ew_var, ab_var, eid_var, mu_var, wmu_var, height_var,x_var, y_var,
          xnum_var, ynum_var, tm_var, tn_var, it_var, conv_var, dm_var, dmh_var, 
          zch_var, ntsk_var, host_var, st_var, ft_var, z_varid, dimids[4];
  long    task;
  size_t  start[] = {0, 0, 0};
  size_t  count[] = {1, 1, 1};
  size_t  len_id;
  bool_t  PRD_angle_dep, XRD;
  double *height;
  char    startJ[MAX_LINE_SIZE], StokesMode[MAX_LINE_SIZE], angleSet[MAX_LINE_SIZE],
          hostname[ARR_STRLEN], timestr[ARR_STRLEN], *atmosID;
  FILE   *test;
  time_t  curtime;
  struct tm *loctime;

  /* Check if we can open the file */
  if ((test = fopen(INPUTDATA_FILE, "a")) == NULL) {
    sprintf(messageStr, "Unable to open spectrum output file %s", INPUTDATA_FILE);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(test);
  }

  /* Opwn the file  */
  if ((ierror = nc_open_par(INPUTDATA_FILE, NC_WRITE | NC_MPIPOSIX, 
			      mpi.comm, mpi.info, &ncid))) ERR(ierror,routineName);
  io.in_ncid = ncid;

  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if ((ierror = nc_inq_attlen(ncid, NC_GLOBAL, "atmosID", &len_id ))) 
    ERR(ierror,routineName);

  atmosID = (char *) malloc(len_id+1);

  if ((ierror = nc_get_att_text(ncid, NC_GLOBAL, "atmosID", atmosID ))) 
    ERR(ierror,routineName);

  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,
         "Populations were derived from different atmosphere (%s) than current",
	    atmosID);
    Error(WARNING, routineName, messageStr);
    }
  free(atmosID);
    
  /* Get group IDs */
  if ((ierror = nc_inq_ncid(ncid, "input", &io.in_input_ncid))) ERR(ierror,routineName);
  if ((ierror = nc_inq_ncid(ncid, "atmos", &io.in_atmos_ncid))) ERR(ierror,routineName);
  if ((ierror = nc_inq_ncid(ncid, "mpi",   &io.in_mpi_ncid  ))) ERR(ierror,routineName);

  /* --- Definitions for the ATMOS group --- */
  /* get variable IDs */
  if ((ierror = nc_inq_varid(io.in_atmos_ncid, "temperature",       &io.in_atmos_T  ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_atmos_ncid, "velocity_z",        &io.in_atmos_vz ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_atmos_ncid, "height",            &io.in_atmos_z  ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_atmos_ncid, "element_weight",    &io.in_atmos_ew ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_atmos_ncid, "element_abundance", &io.in_atmos_ab ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_atmos_ncid, "element_id",        &io.in_atmos_eid))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_atmos_ncid, "muz",               &io.in_atmos_mu ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_atmos_ncid, "wmu",               &io.in_atmos_wmu))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_atmos_ncid, "x",                 &io.in_atmos_x  ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_atmos_ncid, "y",                 &io.in_atmos_y  ))) 
      ERR(ierror,routineName);

  
  /* --- Definitions for the MPI group --- */
  /* get variable IDs */
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, XNUM_NAME,   &io.in_mpi_xnum))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, YNUM_NAME,   &io.in_mpi_ynum))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, TASK_MAP,    &io.in_mpi_tm  ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, TASK_NUMBER, &io.in_mpi_tn  ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, ITER_NAME,   &io.in_mpi_it  ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, CONV_NAME,   &io.in_mpi_conv))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, DM_NAME,     &io.in_mpi_dm  ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, ZC_NAME,     &io.in_mpi_zc  ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, NTASKS,      &io.in_mpi_ntsk))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, HOSTNAME,    &io.in_mpi_host))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, START_TIME,  &io.in_mpi_st  ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, FINISH_TIME, &io.in_mpi_ft  ))) 
      ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(io.in_mpi_ncid, DMH_NAME,    &io.in_mpi_dmh ))) 
      ERR(ierror,routineName);

  return;
}
/* ------- end   --------------------------   init_ncdf_indata_old.c  --- */


/* ------- begin --------------------------   close_ncdf_indata.c --- */
void close_ncdf_indata(void)
/* Closes the indata netCDF file */ 
{
  const char routineName[] = "close_ncdf_indata";
  int        ierror;

  if ((ierror = nc_close(io.in_ncid))) ERR(ierror,routineName);
  return; 
}
/* ------- end   --------------------------   close_ncdf_indata.c --- */

/* ------- begin --------------------------   writeAtmos_p.c --- */
void writeAtmos_all(void) {
/* Reads the NCDF atmos file, and writes data into indata file, all tasks at once */

  const char routineName[] = "writeAtmos_all";
  int        ierror, ncid_in, ncid_out, task;
  size_t     start[]  = {0, 0, 0, 0};
  size_t     count[]  = {1, 1, 1, 1};
  size_t     starti[] = {0, 0, 0, 0, 0};
  size_t     counti[] = {1, 1, 1, 1, 1};
  double    *tmp, **mtmp, *zeros;
  size_t    *st, *ct;


  ncid_in  = infile.ncid;
  ncid_out = io.in_atmos_ncid;

  /* set collective access for variables 
  if ((ierror = nc_var_par_access(ncid_out, io.in_atmos_T, NC_COLLECTIVE)))
    ERR(ierror,routineName);
  if ((ierror = nc_var_par_access(ncid_out, io.in_atmos_ne, NC_COLLECTIVE)))
    ERR(ierror,routineName);
  if ((ierror = nc_var_par_access(ncid_out, io.in_atmos_vz, NC_COLLECTIVE)))
    ERR(ierror,routineName);
  if ((ierror = nc_var_par_access(ncid_out, io.in_atmos_vt, NC_COLLECTIVE)))
    ERR(ierror,routineName);
  if ((ierror = nc_var_par_access(ncid_out, io.in_atmos_nh, NC_COLLECTIVE)))
    ERR(ierror,routineName);
  */
  // end collective set


  tmp   = (double *) calloc(infile.nz,  sizeof(double));
  zeros = (double *) calloc(infile.nz,  sizeof(double)); 
  mtmp  = matrix_double(atmos.NHydr, infile.nz);
  
  starti[0] = input.p15d_nt; counti[0] = 1;

  for (task = 0; task < mpi.Ntasks; task++) {
    start[0] = 0;   count[0] = atmos.NHydr;
    start[1] = mpi.taskmap[task + mpi.my_start][0];  count[2] = 1;
    start[2] = mpi.taskmap[task + mpi.my_start][1];  count[2] = 1;
    start[3] = mpi.zcut_hist[task];   count[3] = infile.nz - start[3];

    /* convert ix, iy into xnum, ynum of the original file */
    //starti[0] = start[0];           starti[3] = start[3];
    starti[0] = input.p15d_nt;      counti[0] = 1;
    starti[1] = mpi.xnum[start[1]]; counti[1] = 1;
    starti[2] = mpi.ynum[start[2]]; counti[2] = 1;
    starti[3] = start[3];           counti[3] = count[3];

    st = &start[1];  ct = &count[1];  

    /* Temperature */
    if ((ierror = nc_get_vara_double(ncid_in,  infile.T_varid, starti, counti, tmp)))
      ERR(ierror,routineName);
    if ((ierror = nc_put_vara_double(ncid_out, io.in_atmos_T,  st,  ct, tmp)))
      ERR(ierror,routineName);

    /* Electron density 
    if ((ierror = nc_get_vara_double(ncid_in,  infile.ne_varid, starti, counti, tmp)))
      ERR(ierror,routineName);
    if ((ierror = nc_put_vara_double(ncid_out, io.in_atmos_ne,  st,  ct, tmp)))
      ERR(ierror,routineName);
    */

    /* Vz */
    if ((ierror = nc_get_vara_double(ncid_in,  infile.vz_varid, starti, counti, tmp)))
      ERR(ierror,routineName);
    if ((ierror = nc_put_vara_double(ncid_out, io.in_atmos_vz,  st,  ct, tmp)))
      ERR(ierror,routineName);

    /* Vturb is always zero for now, not writing anything 
    if ((ierror = nc_put_vara_double(ncid_out, io.in_atmos_vt,  st, ct, zeros)))
      ERR(ierror,routineName);
    */

    if (atmos.Stokes) {
      /* Bx, By, Bz 
      if ((ierror = nc_get_vara_double(ncid_in,  infile.Bx_varid, starti, counti, tmp)))
	ERR(ierror,routineName);
      if ((ierror = nc_put_vara_double(ncid_out, io.in_atmos_Bx,  st,  ct, tmp)))
	ERR(ierror,routineName);

      if ((ierror = nc_get_vara_double(ncid_in,  infile.By_varid, starti, counti, tmp)))
	ERR(ierror,routineName);
      if ((ierror = nc_put_vara_double(ncid_out, io.in_atmos_By,  st,  ct, tmp)))
	ERR(ierror,routineName);

      if ((ierror = nc_get_vara_double(ncid_in,  infile.Bz_varid, starti, counti, tmp)))
	ERR(ierror,routineName);
      if ((ierror = nc_put_vara_double(ncid_out, io.in_atmos_Bz,  st,  ct, tmp)))
	ERR(ierror,routineName);
      */
    }

    /* Hydrogen populations 
    starti[0] = input.p15d_nt;      counti[0] = 1;
    starti[1] = 0;                  counti[1] = atmos.NHydr;
    starti[2] = mpi.xnum[start[1]]; counti[2] = 1;
    starti[3] = mpi.ynum[start[2]]; counti[3] = 1;
    starti[4] = start[3];           counti[4] = count[3];
    if ((ierror=nc_get_vara_double(ncid_in,  infile.nh_varid, starti, counti, 
				   mtmp[0]))) ERR(ierror,routineName);
    if ((ierror=nc_put_vara_double(ncid_out, io.in_atmos_nh,  start,  count,
				   mtmp[0]))) ERR(ierror,routineName);
    */
  }

  free(tmp);
  free(zeros);
  freeMatrix((void **) mtmp);

  return;
}

/* ------- begin --------------------------   writeAtmos_p.c --- */
void writeAtmos_p(void)
{
  /* Write atmos arrays. This has now been modified and writes the interpolated
     arrays, from depth_refine. With that, now this is the only viable option
     to write the atmos data, as there is no option to save in memory and
     writeAtmos_all just writes from the input file, not the interpolated
     quantities
     
     IMPORTANT: at the moment this is a trimmed version, only writing z to save
                space and computational time.
     
     */
  const char routineName[] = "writeAtmos_p";
  int     ierror, ncid;
  size_t  start[] = {0, 0, 0, 0};
  size_t  count[] = {1, 1, 1, 1};

  ncid = io.in_atmos_ncid;

  /* put atmosphere variables */
  start[0] = mpi.ix;   count[0] = 1;
  start[1] = mpi.iy;   count[1] = 1;
  start[2] = mpi.zcut; count[2] = atmos.Nspace;


  /* Tiago: modified, now only writing T, vz, and z */
  
  if ((ierror = nc_put_vara_double(ncid, io.in_atmos_T, start, count,
				   atmos.T ))) ERR(ierror,routineName);
  //if ((ierror = nc_put_vara_double(ncid, io.in_atmos_ne, start, count,
  //				   atmos.ne ))) ERR(ierror,routineName);
  if ((ierror = nc_put_vara_double(ncid, io.in_atmos_vz, start, count,
				   geometry.vel ))) ERR(ierror,routineName);
  
  if ((ierror = nc_put_vara_double(ncid, io.in_atmos_z, start, count,
				   geometry.height ))) ERR(ierror,routineName);
  
  /*
  if ((ierror = nc_put_vara_double(ncid, io.in_atmos_vt, start, count,
				   atmos.vturb ))) ERR(ierror,routineName);
  */
  if (atmos.Stokes) {
    /* These are commented out, because now we're writing Bx, By, Bz,
       in writeAtmos_all. 
    if ((ierror = nc_put_vara_double(ncid, io.in_atmos_B, start, count,
				     atmos.B ))) ERR(ierror,routineName); 
    if ((ierror = nc_put_vara_double(ncid, io.in_atmos_gB, start, count,
				     atmos.gamma_B ))) ERR(ierror,routineName);
    if ((ierror = nc_put_vara_double(ncid, io.in_atmos_chiB, start, count,
				     atmos.chi_B ))) ERR(ierror,routineName);
    */
  }

  /* put hydrogen populations */
  /*
  start[0] = 0;        count[0] = atmos.H->Nlevel;
  start[1] = mpi.ix;   count[1] = 1;
  start[2] = mpi.iy;   count[2] = 1;
  start[3] = 0;        count[3] = atmos.Nspace;

  if ((ierror=nc_put_vara_double(ncid, io.in_atmos_nh, start, count,
				 atmos.H->n[0] )))   ERR(ierror,routineName);  
  */

  return;
}
/* ------- end   --------------------------   writeAtmos_p.c --- */

/* ------- begin --------------------------   writeMPI_p.c --- */
void writeMPI_all(void) {
/* Writes output on indata file, MPI group, all tasks at once */ 
  const char routineName[] = "writeMPI_p";
  int     ierror, task;
  size_t  start[] = {0, 0, 0, 0};
  size_t  count[] = {1, 1, 1, 1};
  char       timestr[ARR_STRLEN];
  time_t     curtime;
  struct tm *loctime;


  /* set collective access for variables 
  if ((ierror = nc_var_par_access(io.in_mpi_ncid, io.in_mpi_ft, NC_COLLECTIVE)))
    ERR(ierror,routineName);
  if ((ierror = nc_var_par_access(io.in_mpi_ncid, io.in_mpi_it, NC_COLLECTIVE)))
    ERR(ierror,routineName);
  if ((ierror = nc_var_par_access(io.in_mpi_ncid, io.in_mpi_conv, NC_COLLECTIVE)))
    ERR(ierror,routineName);
  if ((ierror = nc_var_par_access(io.in_mpi_ncid, io.in_mpi_dm, NC_COLLECTIVE)))
    ERR(ierror,routineName);
  if ((ierror = nc_var_par_access(io.in_mpi_ncid, io.in_mpi_dmh, NC_COLLECTIVE)))
    ERR(ierror,routineName);
    end set collective */



  /* Get finish time in ISO 8601 */
  curtime = time(NULL);
  loctime = localtime(&curtime);
  strftime(timestr, ARR_STRLEN, "%Y-%m-%dT%H:%M:%S%z", loctime);
  
  start[0] = mpi.rank; count[0] = 1;
  start[1] = 0;        count[1] = strlen(timestr);
  if ((ierror = nc_put_vara_text(io.in_mpi_ncid, io.in_mpi_ft, start, count,
      (const char *) &timestr )))     ERR(ierror,routineName); 

  /* Write arrays of Ntasks, one value at a time */
  for (task = 0; task < mpi.Ntasks; task++) {

    start[0] = mpi.taskmap[task + mpi.my_start][0];  count[0] = 1;
    start[1] = mpi.taskmap[task + mpi.my_start][1];  count[1] = 1;
    start[2] = 0;                                    count[2] = mpi.niter[task];

    /* number of iterations */
    if ((ierror = nc_put_var1_int(io.in_mpi_ncid,    io.in_mpi_it, start, 
			  &mpi.niter[task])))        ERR(ierror,routineName);
    /* convergence */
    if ((ierror = nc_put_var1_int(io.in_mpi_ncid,    io.in_mpi_conv, start, 
			  &mpi.convergence[task])))  ERR(ierror,routineName);  
    /* zcut hist */
    if ((ierror = nc_put_var1_int(io.in_mpi_ncid,    io.in_mpi_zc, start,
		          &mpi.zcut_hist[task])))  ERR(ierror,routineName);
    /* dpopsmax */
    if ((ierror = nc_put_var1_double(io.in_mpi_ncid, io.in_mpi_dm, start, 
			  &mpi.dpopsmax[task])))     ERR(ierror,routineName); 
    /* dpopsmax hist */
    if ((ierror = nc_put_vara_double(io.in_mpi_ncid, io.in_mpi_dmh, start,
		     count, mpi.dpopsmax_hist[task])))  ERR(ierror,routineName);
  }
  


  /*
  ncid = io.in_mpi_ncid;

  start[0] = mpi.ix;
  start[1] = mpi.iy;
  if ((ierror = nc_put_var1_int(ncid,    io.in_mpi_it, start, &mpi.niter)))
    ERR(ierror,routineName);
  if ((ierror = nc_put_var1_double(ncid, io.in_mpi_dm, start, &mpi.dpopsmax)))
    ERR(ierror,routineName);
  if ((ierror = nc_put_var1_int(ncid,  io.in_mpi_conv, start, &mpi.convergence)))
    ERR(ierror,routineName);

  start[2] = 0; count[2] = mpi.niter;
  if ((ierror = nc_put_vara_double(ncid, io.in_mpi_dmh, start, count, 
				   mpi.dpopsmax_hist))) ERR(ierror,routineName);


  */
  return;
}

/* ------- end   --------------------------   writeMPI_p.c --- */


/* ------- begin -------------------------- readConvergence.c  --- */
void readConvergence(void) {
  /* This is a self-contained function to read the convergence matrix,
     written by RH. */
  const char routineName[] = "readConvergence";
  size_t len_id, nx, ny;
  int    ncid, ncid_mpi, ierror, dimid, varid;
  char  *atmosID;
  FILE  *test;

  mpi.rh_converged = matrix_int(mpi.nx, mpi.ny);
  
  /* --- Check if we can open the file --- */
  if ((test = fopen(INPUTDATA_FILE, "r")) == NULL) {
    sprintf(messageStr, "Unable to read input data file %s", AUX_FILE);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(test);
  }

  /* --- Open the inputdata file --- */
  if ((ierror=nc_open(INPUTDATA_FILE,NC_NOWRITE,&ncid))) ERR(ierror,routineName);
  
  /* Get ncid of the MPI group */
  if ((ierror = nc_inq_ncid(ncid, "mpi", &ncid_mpi))) ERR(ierror,routineName);


  /* --- Consistency checks --- */
  /* Check that atmosID is the same */
  if ((ierror = nc_inq_attlen(ncid, NC_GLOBAL, "atmosID", &len_id ))) 
    ERR(ierror,routineName);

  atmosID = (char *) malloc(len_id+1);

  if ((ierror = nc_get_att_text(ncid, NC_GLOBAL, "atmosID", atmosID ))) 
    ERR(ierror,routineName);

  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,
         "Input data were derived from different atmosphere (%s) than current",
	    atmosID);
    Error(WARNING, routineName, messageStr);
    }
  free(atmosID);

  /* Check that dimension sizes match */
  if ((ierror = nc_inq_dimid(ncid, "nx", &dimid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_inq_dimlen(ncid, dimid, &nx ))) 
    ERR(ierror,routineName);    

  if (nx != mpi.nx) {
    sprintf(messageStr,
	    "Number of x points mismatch: expected %d, found %d.",
	    mpi.nx, (int)nx);
    Error(WARNING, routineName, messageStr);
  }

  if ((ierror = nc_inq_dimid(ncid, "ny", &dimid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_inq_dimlen(ncid, dimid, &ny ))) 
    ERR(ierror,routineName);    

  if (ny != mpi.ny) {
    sprintf(messageStr,
	    "Number of y points mismatch: expected %d, found %d.",
	    mpi.ny, (int)ny);
    Error(WARNING, routineName, messageStr);
  }
  

  /* --- Read variable --- */
  if ((ierror = nc_inq_varid(ncid_mpi, CONV_NAME, &varid ))) 
    ERR(ierror,routineName);  
  if ((ierror = nc_get_var_int(ncid_mpi, varid, mpi.rh_converged[0]))) 
    ERR(ierror,routineName);  


  /* --- Close inputdata file --- */
  if ((ierror = nc_close(ncid))) ERR(ierror,routineName);

  return;
}
/* ------- end   -------------------------- readConvergence.c  --- */

