/* ------- file: -------------------------- readatmos_ncdf.c ------------

       Version:       rh2.0, 1.5-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       --------------------------                      ----------RH-- */

/* --- Reads atmospheric model in NetCDF format. --     -------------- */



// Must check if all of these are really needed!
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "constant.h"
#include "background.h"
#include "error.h"
#include "inputs.h"
#include "statistics.h"
#include "parallel.h"
#include "io.h"



/* --- Function prototypes --                          -------------- */
void setTcut(Atmosphere *atmos, Geometry *geometry, double Tmax);
void realloc_ndep(Atmosphere *atmos, Geometry *geometry);
void depth_refine(Atmosphere *atmos, Geometry *geometry, double tmax);

/* --- Global variables --                             -------------- */

extern MPI_data mpi;
extern InputData input;
extern char messageStr[];

/* ------- begin --------------------------   init_ncdf_atmos   ----- */

void init_ncdf_atmos(Atmosphere *atmos, Geometry *geometry, NCDF_Atmos_file *infile)
/* Initialises the input atmosphere file, gets dimensions and variable ids. 
   Also performs other basic RH initialisations like readAbundance       */
{ 
  const char routineName[] = "init_ncdf";
  struct  stat statBuffer;
  int ierror, ncid, x_varid, y_varid, z_varid, has_B, sn_varid;
  size_t nn;
  size_t start[] = {0, 0};
  size_t count[] = {1, 1};
  char *filename;


  /* --- Get abundances of background elements --      -------------- */

  readAbundance(atmos);

  /* --- Open input file for model atmosphere --       -------------- */
  if ((ierror = nc_open_par(input.atmos_input ,NC_NOWRITE | NC_MPIIO, mpi.comm, 
			    mpi.info, &infile->ncid))) ERR(ierror,routineName);

  ncid = infile->ncid;

  /* Is magnetic field included? */
  if ((ierror = nc_get_att_int( ncid, NC_GLOBAL, "has_B", &has_B))) 
    ERR(ierror,routineName);

  atmos->Stokes = FALSE;
  if ((has_B) && (strcmp(input.Stokes_input, "none"))) atmos->Stokes = TRUE;


  /* Get the dimids and values */
  if ((ierror = nc_inq_dimid( ncid, "nx", &infile->nx_id)))       ERR(ierror,routineName);
  if ((ierror = nc_inq_dimlen(ncid, infile->nx_id, &infile->nx))) ERR(ierror,routineName);

  if ((ierror = nc_inq_dimid( ncid, "ny", &infile->ny_id)))       ERR(ierror,routineName);
  if ((ierror = nc_inq_dimlen(ncid, infile->ny_id, &infile->ny))) ERR(ierror,routineName);

  if ((ierror = nc_inq_dimid( ncid, "nz", &infile->nz_id)))       ERR(ierror,routineName);
  if ((ierror = nc_inq_dimlen(ncid, infile->nz_id, &infile->nz))) ERR(ierror,routineName);

  if ((ierror = nc_inq_dimid( ncid, "nhydr", &infile->nhyd_id)))  ERR(ierror,routineName);
  if ((ierror = nc_inq_dimlen(ncid, infile->nhyd_id, &nn)))       ERR(ierror,routineName);

  /* get some values in atmos/geometry structures */
  geometry->Ndep = (int)  infile->nz;
  atmos->Nspace  = (long) infile->nz;
  atmos->NHydr   = (int)  nn;


  if (atmos->NHydr < 2  &&  !atmos->H_LTE) {
    sprintf(messageStr, "NHydr has to be at least 2, not %d to run with"
	    " NLTE hydrogen populations in background\n", atmos->NHydr);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* Get the varids */
  if ((ierror=nc_inq_varid(ncid, TEMP_NAME, &infile->T_varid)))  ERR(ierror,routineName);
  if ((ierror=nc_inq_varid(ncid, NE_NAME,   &infile->ne_varid))) ERR(ierror,routineName);
  if ((ierror=nc_inq_varid(ncid, VZ_NAME,   &infile->vz_varid))) ERR(ierror,routineName);
  if ((ierror=nc_inq_varid(ncid, NH_NAME,   &infile->nh_varid))) ERR(ierror,routineName);
  if ((ierror=nc_inq_varid(ncid, "z",       &z_varid)))          ERR(ierror,routineName);
  if ((ierror=nc_inq_varid(ncid, "y",       &y_varid)))          ERR(ierror,routineName);
  if ((ierror=nc_inq_varid(ncid, "x",       &x_varid)))          ERR(ierror,routineName);
  if ((ierror=nc_inq_varid(ncid, SNAPNAME,  &sn_varid)))         ERR(ierror,routineName);
  
  /* Microturbulence, get ID if variable found */
  if ((ierror=nc_inq_varid(ncid, VTURB_NAME, &infile->vturb_varid))) {
    /* exception for variable not found (errcode -49)*/
    (ierror == -49) ? infile->vturb_varid = -1 : ERR(ierror,routineName);
  }

  if (atmos->Stokes) {
      if ((ierror = nc_inq_varid(ncid, BX_NAME, &infile->Bx_varid)))
	ERR(ierror,routineName);
      if ((ierror = nc_inq_varid(ncid, BY_NAME, &infile->By_varid)))
	ERR(ierror,routineName);
      if ((ierror = nc_inq_varid(ncid, BZ_NAME, &infile->Bz_varid)))
	ERR(ierror,routineName);
  }

  /* read things that don't depend on x, y */
  start[0] = input.p15d_nt; count[0] = 1;
  start[1] = 0;             count[1] = infile->nz;

  geometry->height = (double *) malloc(infile->nz * sizeof(double));
  if ((ierror = nc_get_vara_double(ncid, z_varid, start, count, geometry->height))) 
    ERR(ierror,routineName);

  infile->y   = (double *) malloc(infile->ny * sizeof(double));
  if ((ierror = nc_get_var_double(ncid, y_varid, infile->y))) ERR(ierror,routineName);
  infile->x   = (double *) malloc(infile->nx * sizeof(double));
  if ((ierror = nc_get_var_double(ncid, x_varid, infile->x))) ERR(ierror,routineName);

  if ((ierror = nc_get_var1_int(ncid, sn_varid, start, &mpi.snap_number)))
    ERR(ierror,routineName);

  /* allocate arrays */
  geometry->vel = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->T      = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->ne     = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->vturb  = (double *) calloc(atmos->Nspace , sizeof(double)); /* default zero */
  atmos->nHtot  = (double *) malloc(atmos->Nspace * sizeof(double));

  if (atmos->Stokes) {
    atmos->B       = (double *) malloc(atmos->Nspace * sizeof(double));
    atmos->gamma_B = (double *) malloc(atmos->Nspace * sizeof(double));
    atmos->chi_B   = (double *) malloc(atmos->Nspace * sizeof(double));
  }

  /* some other housekeeping */ 
  geometry->vboundary[TOP]    = ZERO;
  geometry->vboundary[BOTTOM] = THERMALIZED;
  geometry->scale             = GEOMETRIC;

  /* --- Construct atmosID from filename and last modification date - */
  // NOTE: perhaps later this should be built into the NetCDF file (description attr?)
  stat(input.atmos_input, &statBuffer);
  if ((filename = strrchr(input.atmos_input, '/')) != NULL)
    filename++;
  else
    filename = input.atmos_input;
  sprintf(atmos->ID, "%s (%.24s)", filename,
	  asctime(localtime(&statBuffer.st_mtime)));


  /* --- Get angle-quadrature and copy geometry independent quantity
         wmu to atmos structure. --                    -------------- */

  getAngleQuad(geometry);
  atmos->wmu = geometry->wmu;


  /* --- set up pointers for background opacities --- */

  atmos->chi_b = NULL;
  atmos->eta_b = NULL;
  atmos->sca_b = NULL;

  return;

}
/* ------- end ---------------------------- init_ncdf_atmos  -------- */

/* ------- begin -------------------------- readAtmos_ncdf  --------- */

void readAtmos_ncdf(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
		    NCDF_Atmos_file *infile)
/* Reads the variables T, ne, vel, nh for a given (xi,yi) pair */ 
{
  const char routineName[] = "readAtmos_ncdf";
  size_t     start[]    = {0, 0, 0, 0}; /* starting values */
  size_t     count[]    = {1, 1, 1, 1};
  size_t     start_nh[] = {0, 0, 0, 0, 0};
  size_t     count_nh[] = {1, 1, 1, 1, 1};
  int        ncid, ierror, i, j, z_varid, imin, i50k;
  bool_t     old_moving;
  double     Tmin, diff;
  double    *Bx, *By, *Bz;

  ncid = infile->ncid;

  atmos->Nspace = geometry->Ndep = infile->nz;

  /* read full T column, to see where to zcut */
  start[0] = input.p15d_nt; count[0] = 1;
  start[1] = (size_t) xi;   count[1] = 1;
  start[2] = (size_t) yi;   count[2] = 1;
  start[3] = 0;             count[3] = infile->nz;
  
  atmos->T = (double *) realloc(atmos->T, infile->nz * sizeof(double)); 

  if ((ierror = nc_get_vara_double(ncid, infile->T_varid,  start, count, atmos->T)))
    ERR(ierror,routineName);

  /* Finds z value for Tmax cut, redefines Nspace, reallocates arrays */
  /* Tiago: not using this at the moment, only z cut in depth_refine */
  if (input.p15d_zcut) {
    setTcut(atmos, geometry, input.p15d_tmax);
  } else {
    mpi.zcut = 0;
  }

  /* Get z again */
  start[0] = input.p15d_nt; count[0] = 1;
  start[1] = mpi.zcut;      count[1] = atmos->Nspace;

  if ((ierror=nc_inq_varid(ncid, "z",  &z_varid)))          
    ERR(ierror,routineName);
  if ((ierror = nc_get_vara_double(ncid, z_varid, start, count, geometry->height))) 
    ERR(ierror,routineName);
 

  start[0] = input.p15d_nt; count[0] = 1;
  start[1] = (size_t) xi;   count[1] = 1;
  start[2] = (size_t) yi;   count[2] = 1;
  start[3] = mpi.zcut;      count[3] = atmos->Nspace;

   /* read variables */
  if ((ierror = nc_get_vara_double(ncid, infile->T_varid,  start, count, atmos->T)))
    ERR(ierror,routineName);
  if ((ierror = nc_get_vara_double(ncid, infile->ne_varid, start, count, atmos->ne)))
    ERR(ierror,routineName);
  if ((ierror = nc_get_vara_double(ncid, infile->vz_varid, start, count, geometry->vel)))
    ERR(ierror,routineName);
  /* vturb, if available */   
  if (infile->vturb_varid != -1) {
    if ((ierror = nc_get_vara_double(ncid, infile->vturb_varid, &start[3], &count[3],
				     atmos->vturb))) ERR(ierror,routineName);
  }

  /* Read magnetic field */
  if (atmos->Stokes) {
    Bx = (double *) malloc(atmos->Nspace * sizeof(double));
    By = (double *) malloc(atmos->Nspace * sizeof(double));
    Bz = (double *) malloc(atmos->Nspace * sizeof(double));
    /* Read in cartesian coordinates */
    if ((ierror = nc_get_vara_double(ncid, infile->Bx_varid,  start, count,
				     Bx))) ERR(ierror,routineName);
    if ((ierror = nc_get_vara_double(ncid, infile->By_varid,  start, count,
				     By))) ERR(ierror,routineName);
    if ((ierror = nc_get_vara_double(ncid, infile->Bz_varid,  start, count,
				     Bz))) ERR(ierror,routineName);
    
    /* Convert to spherical coordinates */
    for (j = 0; j < atmos->Nspace; j++) {
      atmos->B[j]       = sqrt(SQ(Bx[j]) + SQ(By[j]) + SQ(Bz[j]));
      atmos->gamma_B[j] = acos(Bz[j]/atmos->B[j]);
      atmos->chi_B[j]   = atan(By[j]/Bx[j]);
      
      /* Protect from undefined cases */
      if ((Bx[j] == 0) && (By[j] == 0) && (Bz[j] == 0))
	atmos->gamma_B[j] = 0.0;
      
      if ((Bx[j] == 0) && (By[j] == 0))
	atmos->chi_B[j]   = 1.0;
    }

    free(Bx); free(By); free(Bz);
  }
  


  /* allocate and zero nHtot */
  atmos->nH = matrix_double(atmos->NHydr, atmos->Nspace);  
  for (j = 0; j < atmos->Nspace; j++) atmos->nHtot[j] = 0.0; 

  /* read nH, all at once */
  start_nh[0] = input.p15d_nt; count_nh[0] = 1;
  start_nh[1] = 0;             count_nh[1] = atmos->NHydr;
  start_nh[2] = (size_t) xi;   count_nh[2] = 1;
  start_nh[3] = (size_t) yi;   count_nh[3] = 1;
  start_nh[4] = mpi.zcut;      count_nh[4] = atmos->Nspace;
  if ((ierror = nc_get_vara_double(ncid, infile->nh_varid, start_nh, count_nh, 
				   atmos->nH[0]))) ERR(ierror,routineName);

  /* Depth grid refinement */
  if (input.p15d_refine)
    depth_refine(atmos, geometry, input.p15d_tmax);
  
  /* Fix vturb: remove zeros, use multiplier and add */
  for (i = 0; i < atmos->Nspace; i++) {
    if (atmos->vturb[i] < 0.0) atmos->vturb[i] = 0.0;
    atmos->vturb[i] = atmos->vturb[i] * input.vturb_mult + input.vturb_add;
  }

  /* Sum to get nHtot */
  for (i = 0; i < atmos->NHydr; i++){
    for (j = 0; j < atmos->Nspace; j++) atmos->nHtot[j] += atmos->nH[i][j];
  }
  
  /* Some other housekeeping */
  old_moving = atmos->moving;
  atmos->moving = FALSE;
  for (i = 0;  i < atmos->Nspace;  i++) {
    if (fabs(geometry->vel[i]) >= atmos->vmacro_tresh) {
      atmos->moving = TRUE;
      /* old_moving should only be false*/
      if ((old_moving == FALSE) & (atmos->moving == TRUE)) {
	sprintf(messageStr,
		"Moving atmosphere detected when the previous column\n"
		" (or column [0,0] in file) was not. This will cause problems\n"
		" and the code will abort.\n"
		" To prevent this situation one can force all columns\n"
		" to be moving by setting VMACRO_TRESH = 0 in keyword.input\n");
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      break;
    }
  }

  return;
}


/* ------- end ---------------------------- readAtmos_ncdf  --------- */

/* ------- begin -------------------------- close_ncdf  ------------- */
void close_ncdf_atmos(Atmosphere *atmos, Geometry *geometry, NCDF_Atmos_file *infile)
/* Closes the NCDF file and frees memory */
{
  const char routineName[] = "close_atmos_ncdf";
  int ierror;
  
  /* Close the file. */
  if ((ierror = nc_close(infile->ncid))) ERR(ierror,routineName);

  /* Free stuff */
  free(atmos->T);
  free(atmos->ne);
  free(atmos->vturb);
  free(atmos->nHtot);
  free(geometry->vel);
  free(geometry->height);
  free(infile->y);
  free(infile->x);

  return; 

}

/* ------- end ---------------------------- close_ncdf   ------------- */

/* ------- begin--------------------------- setTcut      ------------- */
void setTcut(Atmosphere *atmos, Geometry *geometry, double Tmax) {
  /* Find the point where temperature (in TR region) gets below Tmax,
     set atmos.Nspace to remaining points, repoint all the atmospheric
     quantities to start at that point.  */

  const char routineName[] = "setTcut";
  int   i, cutpoint;

  mpi.zcut =  0;

  if (Tmax < 0) return; /* Do nothing when Tmax < 0 */

  cutpoint = -1;

  for (i=0; i < atmos->Nspace; i++) {
    if (atmos->T[i] <= Tmax) {
      cutpoint = i;
      break;
    }
  }

  if (cutpoint < 0) {
    sprintf(messageStr,
	    "\n-Could not find temperature cut point! Aborting.\n");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  mpi.zcut = cutpoint;

  /* subtract from total number of points */
  atmos->Nspace  -= mpi.zcut;
  geometry->Ndep -= mpi.zcut;

  /* Reallocate arrays of Nspace */ 
  realloc_ndep(atmos, geometry);
   
  return;
}
/* ------- end  --------------------------- setTcut      ----------- */

/* ------- begin--------------------------- realloc_ndep ----------- */
void realloc_ndep(Atmosphere *atmos, Geometry *geometry) {
  /* Reallocates the arrays of Nspace */

  

  atmos->T     = (double *) realloc(atmos->T,     atmos->Nspace * sizeof(double));
  atmos->ne    = (double *) realloc(atmos->ne,    atmos->Nspace * sizeof(double));
  atmos->nHtot = (double *) realloc(atmos->nHtot, atmos->Nspace * sizeof(double));
  geometry->vel    = (double *) realloc(geometry->vel, 
					atmos->Nspace * sizeof(double));  
  geometry->height = (double *) realloc(geometry->height, 
					atmos->Nspace * sizeof(double));
  if (atmos->nHmin != NULL) 
      atmos->nHmin = (double *) realloc(atmos->nHmin,
					atmos->Nspace * sizeof(double));
  /* default zero */
  free(atmos->vturb);
  atmos->vturb  = (double *) calloc(atmos->Nspace , sizeof(double)); 

  if (atmos->Stokes) {
    atmos->B       = (double *) realloc(atmos->B, 
					atmos->Nspace * sizeof(double));
    atmos->gamma_B = (double *) realloc(atmos->gamma_B,
					atmos->Nspace * sizeof(double));
    atmos->chi_B   = (double *) realloc(atmos->chi_B,
					atmos->Nspace * sizeof(double));
  }

  return;
}
/* ------- end  --------------------------- realloc_ndep ----------- */

/* ------- begin--------------------------- depth_refine ----------- */
void depth_refine(Atmosphere *atmos, Geometry *geometry, double Tmax) {
  /* Performs depth refinement to optimise for gradients in temperature,
     density and optical depth. In the same fashion as multi23's ipol_dscal
     (in fact, shamelessly copied). 
  */
  
  bool_t nhm_flag, hunt;
  long    i, k, k0=0, k1=0;
  size_t  bufsize;
  double  CI, PhiHmin, *chi, *eta, *tau, tdiv, rdiv, taudiv, *aind, *xpt;
  double *new_height, *buf;
  const double taumax = 100.0, lg1 = log10(1.1);
  
  if (Tmax < 0) Tmax = 1.e20; /* No zcut if Tmax is negative */ 
  
  chi = (double *) malloc(atmos->Nspace * sizeof(double));
  eta = (double *) malloc(atmos->Nspace * sizeof(double));
  tau = (double *) calloc(atmos->Nspace , sizeof(double));
  xpt = (double *) calloc(atmos->Nspace , sizeof(double));
  buf = (double *) malloc(atmos->Nspace * sizeof(double));
  aind = (double *) calloc(atmos->Nspace , sizeof(double));
  new_height = (double *) calloc(atmos->Nspace , sizeof(double));
  
  bufsize = atmos->Nspace * sizeof(double);
  
  nhm_flag = FALSE;
  CI = (HPLANCK/(2.0*PI*M_ELECTRON)) * (HPLANCK/KBOLTZMANN);
  k1 = atmos->Nspace;
  
  /* --- Calculate tau500 scale from H-bf opacity alone. --- */
  /* First, build nHmin because ChemicalEquilibrium has not been called yet.
     Unless H level pops are given in ncdf file, all H is assumed neutral.  */
  if (atmos->nHmin == NULL) {
    atmos->nHmin = (double *) malloc(atmos->Nspace * sizeof(double));
    nhm_flag = TRUE;
  }
  
  for (k = 0; k < atmos->Nspace; k++) {
    PhiHmin = 0.25*pow(CI/atmos->T[k], 1.5) *
        exp(0.754 * EV / (KBOLTZMANN * atmos->T[k]));
    atmos->nHmin[k] = atmos->ne[k] * atmos->nH[0][k] * PhiHmin;
  }
  
  Hminus_bf(500.0, chi, eta);
  
  /* integrate for optical depth */
  for (k = 1; k < atmos->Nspace; k++) {
    if ((atmos->T[k] > Tmax) && (k0 == k - 1)) k0 = k;
    tau[k] = tau[k-1] + 0.5 * (chi[k] + chi[k-1]) *
      (fabs(geometry->height[k-1] - geometry->height[k]));
    if (tau[k] < taumax) k1 = k;
    xpt[k] = (double) k;
  }
  
  tau[0] = tau[1];
  
  /* --- Compute log variations --- */
  for (k = k0+1; k <= k1; k++) {  
    tdiv = fabs(log10(atmos->T[k]) - log10(atmos->T[k-1]))/lg1;
    /* rho is not available, so nH[0] used instead */
    rdiv = fabs(log10(atmos->nH[0][k]) - log10(atmos->nH[0][k-1]))/lg1;
    taudiv = fabs(log10(tau[k]) - log10(tau[k-1]))/0.1;
    aind[k] = aind[k-1] + MAX(MAX(tdiv,rdiv),taudiv);
  }
  
  for (k = 1; k <= k1; k++)  
    aind[k] *= (atmos->Nspace-1)/aind[k1];
    
  /* --- Create new height scale --- */
  splineCoef(k1-k0+1, &aind[k0], &geometry->height[k0]);
  splineEval(atmos->Nspace, xpt, new_height, hunt=FALSE);
    
  /* --- Interpolate quantities to new scale --- */
  /* Take logs of densities to avoid interpolation to negatives */
  for (k = 0; k < atmos->Nspace; k++) {
    atmos->ne[k] = log(atmos->ne[k]);
    for (i = 0; i < atmos->NHydr; i++)
      atmos->nH[i][k] = log(atmos->nH[i][k]);    
  }

  splineCoef(atmos->Nspace, geometry->height, atmos->T);
  splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
  memcpy((void *) atmos->T, (void *) buf, bufsize);
  /* condition for T < 1000 K ? */
 
  splineCoef(atmos->Nspace, geometry->height, atmos->ne);
  splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
  memcpy((void *) atmos->ne, (void *) buf, bufsize);
  
  for (i = 0; i < atmos->NHydr; i++) {
    splineCoef(atmos->Nspace, geometry->height, atmos->nH[i]);
    splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
    memcpy((void *) atmos->nH[i], (void *) buf, bufsize);
  }
  
  splineCoef(atmos->Nspace, geometry->height, geometry->vel);
  splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
  memcpy((void *) geometry->vel, (void *) buf, bufsize);
  
  splineCoef(atmos->Nspace, geometry->height, atmos->vturb);
  splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
  memcpy((void *) atmos->vturb, (void *) buf, bufsize);
  
  if (atmos->Stokes) {
    splineCoef(atmos->Nspace, geometry->height, atmos->B);
    splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
    memcpy((void *) atmos->B, (void *) buf, bufsize);

    splineCoef(atmos->Nspace, geometry->height, atmos->gamma_B);
    splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
    memcpy((void *) atmos->gamma_B, (void *) buf, bufsize);

    splineCoef(atmos->Nspace, geometry->height, atmos->chi_B);
    splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
    memcpy((void *) atmos->chi_B, (void *) buf, bufsize);
  }
  
  /* Take back logs */
  for (k = 0; k < atmos->Nspace; k++) {
    atmos->ne[k] = exp(atmos->ne[k]);
    for (i = 0; i < atmos->NHydr; i++)
      atmos->nH[i][k] = exp(atmos->nH[i][k]);    
  }
  
  memcpy((void *) geometry->height, (void *) new_height, bufsize);
    
  if (nhm_flag) {
    free(atmos->nHmin);
    atmos->nHmin = NULL;
  }
  free(buf); free(new_height);
  free(chi); free(eta);
  free(tau); free(aind); free(xpt);
  
  return; 
}
/* ------- end  --------------------------- depth_refine ----------- */
