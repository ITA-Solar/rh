/* ------- file: -------------------------- multiatmos.c ------------

       Version:       rh2.0, 1.5-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Tue Nov 23 16:03:23 2010 --

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
#include <netcdf.h>
#include <mpi.h>


#include "rh.h"


#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "error.h"
#include "inputs.h"
#include "statistics.h"
#include "parallel.h"


//#define FILE_NAME "/Users/tiago/data/bifrost/bifrost_s20_slice.ncdf"
#define TEMP_NAME "temperature"
#define VZ_NAME   "velocity_z"
#define NE_NAME   "electron_density"
#define NH_NAME   "hydrogen_populations"
#define Z_NAME    "z"



/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern MPI_data mpi;
extern InputData input;
extern char messageStr[];

/* ------- begin --------------------------   init_ncdf   ------------ */

void init_ncdf(Atmosphere *atmos, Geometry *geometry, NCDF_Atmos_file *infile)
/* Initialises the input atmosphere file, gets dimensions and variable ids. 
   Also performs other basic RH initialisations like readAbundance       */
{ 
  const char routineName[] = "init_ncdf";
  struct  stat statBuffer;
  int ierror, ncid, z_varid, has_B;
  size_t nn;
  char *filename;


  /* --- Get abundances of background elements --      -------------- */

  readAbundance(atmos);

  /* --- Open input file for model atmosphere --       -------------- */
  
  if ((atmos->fp_atmos = fopen(input.atmos_input, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", input.atmos_input);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(atmos->fp_atmos);
  }
  

  /* Open the file. */
  if ((ierror = nc_open_par(input.atmos_input ,NC_NOWRITE | NC_MPIIO, mpi.comm, 
			    mpi.info, &infile->ncid))) ERR(ierror,routineName);

  ncid = infile->ncid;

  /* Is magnetic field included? */
  if ((ierror = nc_get_att_int( ncid, NC_GLOBAL, "has_B", &has_B))) ERR(ierror,routineName);

  if (has_B) printf("NCDF has B, but we are not doing anything about it for now.\n");
  atmos->Stokes = FALSE;

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


  /* Get the varids */
  if ((ierror = nc_inq_varid(ncid, TEMP_NAME, &infile->T_varid)))  ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(ncid, NE_NAME,   &infile->ne_varid))) ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(ncid, VZ_NAME,   &infile->vz_varid))) ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(ncid, NH_NAME,   &infile->nh_varid))) ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(ncid, "z",       &z_varid)))          ERR(ierror,routineName);

  /* read things that don't depend on x, y */
  geometry->height = (double *) malloc(atmos->Nspace * sizeof(double));
  if ((ierror = nc_get_var_double(ncid, z_varid, geometry->height))) ERR(ierror,routineName);

  /* allocate arrays */
  geometry->vel = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->T      = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->ne     = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->vturb  = (double *) calloc(atmos->Nspace , sizeof(double)); /* default zero */
  atmos->nHtot  = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->nH     = matrix_double(atmos->NHydr, atmos->Nspace);


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


  /* For now, not reading B (could be extended in future to place it ncdf file */



  /* 
     TODO: 
        - Put some attribute on NCDF file to say if atmos->moving = True --NOPE
	- atmos->Stokes = readB(atmos); (placeholder, not for now...) --DONE
	- geometry->Nrays stuff --DONE
	- more flags needed on netcdf file:
	  * H_LTE? --NOPE
 	  * B available?
     	  * ?
	- nHtot --DONE

   */


  return;

}
/* ------- end ---------------------------- init_ncdf   ------------- */

/* ------- begin -------------------------- readAtmos_ncdf  --------- */

void readAtmos_ncdf(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
		    NCDF_Atmos_file *infile)
/* Reads the variables T, ne, vel, nh for a given (xi,yi) pair */ 
{
  const char routineName[] = "readAtmos_ncdf";
  size_t start[]    = {0, 0, 0}; /* starting values */
  size_t count[]    = {1, 1, 1};
  size_t start_nh[] = {0, 0, 0, 0};
  size_t count_nh[] = {1, 1, 1, 1};
  int    ncid, ierror, i, j;

  ncid = infile->ncid;

  count[2] = infile->nz;
  start[0] = (size_t) xi;
  start[1] = (size_t) yi;
  
   /* read variables */
  if ((ierror = nc_get_vara_double(ncid, infile->T_varid,  start, count, atmos->T)))
    ERR(ierror,routineName);
  if ((ierror = nc_get_vara_double(ncid, infile->ne_varid, start, count, atmos->ne)))
    ERR(ierror,routineName);
  if ((ierror = nc_get_vara_double(ncid, infile->vz_varid, start, count, geometry->vel)))
    ERR(ierror,routineName);

  /* zero nHtot */
  for (j = 0; j < atmos->Nspace; j++) atmos->nHtot[j] = 0.0; 

  /* read nH */ 
  count_nh[3] = infile->nz;
  start_nh[1] = (size_t) xi;
  start_nh[2] = (size_t) yi;

  for (i = 0; i < atmos->NHydr; i++){
    start_nh[0] = i;
    if ((ierror = nc_get_vara_double(ncid, infile->nh_varid, start_nh, count_nh, 
				     atmos->nH[i]))) ERR(ierror,routineName);

    /* Sum to get nHtot */
    for (j = 0; j < atmos->Nspace; j++) atmos->nHtot[j] += atmos->nH[i][j];
  }
  

  /* Some other housekeeping */
  // Should this be necessary? Perhaps atmos should always be moving...
  atmos->moving = FALSE;
  for (i = 0;  i < atmos->Nspace;  i++) {
    if (fabs(geometry->vel[i]) >= atmos->vmacro_tresh) {
      atmos->moving = TRUE;
      break;
    }
  }



  return;

}


/* ------- end ---------------------------- readAtmos_ncdf  --------- */

/* ------- begin -------------------------- close_ncdf  ------------- */
void close_atmos_ncdf(Atmosphere *atmos, Geometry *geometry, NCDF_Atmos_file *infile)
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
  //freeMatrix((void **) atmos->nH); // this is free by distribute_nH

  return; 

}

/* ------- end ---------------------------- close_ncdf   ------------- */
