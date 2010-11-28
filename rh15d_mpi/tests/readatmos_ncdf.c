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

//#include "parallel.h"
#include "../../rh.h"


#include "../../atom.h"
#include "../../atmos.h"
#include "geometry.h"
//#include "spectrum.h"
//#include "background.h"
//#include "constant.h"
#include "../../error.h"
#include "../../inputs.h"
//#include "statistics.h"
//#include "xdr.h"


#define FILE_NAME "/Users/tiago/data/bifrost/bifrost_s20_slice.ncdf"
#define TEMP_NAME "temperature"
#define VZ_NAME   "velocity_z"
#define NE_NAME   "electron_density"
#define NH_NAME   "hydrogen_populations"
#define Z_NAME    "z"

#define ERR(e) {printf("(EEE) NetCDF: %s\n", nc_strerror(e)); exit(EXIT_FAILURE);}


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

Atmosphere atmos;
Geometry geometry;

MPI_Comm comm = MPI_COMM_WORLD;
MPI_Info info = MPI_INFO_NULL;

CommandLine commandline;
char messageStr[256];



/*
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[]
*/


/* ------- begin --------------------------   init_ncdf   ------------ */

void init_ncdf(Atmosphere *atmos, Geometry *geometry, NCDF_Atmos_file *infile)
/* Initialises the input atmosphere file, gets dimensions and variable ids. 
   Also performs other basic RH initialisations like readAbundance       */
{ 
  int ierror, ncid;
  int z_varid;
  size_t nn;


  /* --- Get abundances of background elements --      -------------- */

  //readAbundance(atmos);

  /* --- Open input file for model atmosphere --       -------------- */
  /*
  if ((atmos->fp_atmos = fopen(input.atmos_input, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", input.atmos_input);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(fp_atmos);
  }
  */

  /* Open the file. */
  if ((ierror = nc_open_par(infile->file_name,NC_NOWRITE, comm, info, &infile->ncid)))
    ERR(ierror);

  ncid = infile->ncid;

  /* Get the dimids and values */
  if ((ierror = nc_inq_dimid( ncid, "nx", &infile->nx_id)))       ERR(ierror);
  if ((ierror = nc_inq_dimlen(ncid, infile->nx_id, &infile->nx))) ERR(ierror);

  if ((ierror = nc_inq_dimid( ncid, "ny", &infile->ny_id)))       ERR(ierror);
  if ((ierror = nc_inq_dimlen(ncid, infile->ny_id, &infile->ny))) ERR(ierror);

  if ((ierror = nc_inq_dimid( ncid, "nz", &infile->nz_id)))       ERR(ierror);
  if ((ierror = nc_inq_dimlen(ncid, infile->nz_id, &infile->nz))) ERR(ierror);

  if ((ierror = nc_inq_dimid( ncid, "nhydr", &infile->nhyd_id)))  ERR(ierror);
  if ((ierror = nc_inq_dimlen(ncid, infile->nhyd_id, &nn)))       ERR(ierror);

  /* get some values in atmos/geometry structures */
  atmos->Nspace = geometry->Ndep = (int) infile->nz;
  atmos->NHydr  = (int) nn;


  /* Get the varids */
  if ((ierror = nc_inq_varid(ncid, TEMP_NAME, &infile->T_varid)))  ERR(ierror);
  if ((ierror = nc_inq_varid(ncid, NE_NAME,   &infile->ne_varid))) ERR(ierror);
  if ((ierror = nc_inq_varid(ncid, VZ_NAME,   &infile->vz_varid))) ERR(ierror);
  if ((ierror = nc_inq_varid(ncid, NH_NAME,   &infile->nh_varid))) ERR(ierror);
  if ((ierror = nc_inq_varid(ncid, "z",       &z_varid)))          ERR(ierror);

  /* read things that don't depend on x, y */
  geometry->height = (double *) malloc(atmos->Nspace * sizeof(double));
  if ((ierror = nc_get_var_double(ncid, z_varid, geometry->height))) ERR(ierror);

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

  /* vmacro_tres comes from input, wait until later */
  /*
  atmos->moving = FALSE;
  for (k = 0;  k < atmos->Nspace;  k++) {
    if (fabs(geometry->vel[k]) >= atmos->vmacro_tresh) {
      atmos->moving = TRUE;
      break;
    }
  */
  atmos->moving = TRUE;

  /* --- Construct atmosID from filename and last modification date - */
  /*
  stat(input.atmos_input, &statBuffer);
  if ((filename = strrchr(input.atmos_input, '/')) != NULL)
    filename++;
  else
    filename = input.atmos_input;
  sprintf(atmos->ID, "%s (%.24s)", filename,
	  asctime(localtime(&statBuffer.st_mtime)));
  */

  /* For now, not reading B (could be extended in future to place it ncdf file */
  atmos->Stokes = FALSE;




  /* 
     TODO: 
        - Put some attribute on NCDF file to say if atmos->moving = True
	- atmos->Stokes = readB(atmos); (placeholder, not for now...)
	- geometry->Nrays stuff
	- more flags needed on netcdf file:
	  * H_LTE?
 	  * B available?
     	  * ?
	- nHtot

   */


  return;

}
/* ------- end ---------------------------- init_ncdf   ------------- */

/* ------- begin -------------------------- readAtmos_ncdf  --------- */

void readAtmos_ncdf(int xi, int yi, Atmosphere *atmos, Geometry *geometry, NCDF_Atmos_file *infile)
/* Reads the variables T, ne, vel, nh for a given (xi,yi) pair */ 
{
  
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
    ERR(ierror);
  if ((ierror = nc_get_vara_double(ncid, infile->ne_varid, start, count, atmos->ne)))
    ERR(ierror);
  if ((ierror = nc_get_vara_double(ncid, infile->vz_varid, start, count, geometry->vel)))
    ERR(ierror);

  /* zero nHtot */
  for (j = 0; j < atmos->Nspace; j++) atmos->nHtot[j] = 0.0; 

  /* read nH */ 
  count_nh[3] = infile->nz;
  start_nh[1] = (size_t) xi;
  start_nh[2] = (size_t) yi;

  for (i = 0; i < atmos->NHydr; i++){
    start_nh[0] = i;
    if ((ierror = nc_get_vara_double(ncid, infile->nh_varid, start_nh, count_nh, atmos->nH[i])))
      ERR(ierror);

    /* Sum to get nHtot */
    for (j = 0; j < atmos->Nspace; j++) atmos->nHtot[j] += atmos->nH[i][j];
  }
  


  return;

}


/* ------- end ---------------------------- readAtmos_ncdf  --------- */

/* ------- begin -------------------------- close_ncdf  ------------- */
void close_ncdf(Atmosphere *atmos, Geometry *geometry, NCDF_Atmos_file *infile)
/* Closes the NCDF file and frees memory */
{
  int ierror;
  
  /* Close the file. */
  if ((ierror = nc_close(infile->ncid))) ERR(ierror);

  /* Free stuff */
  free(atmos->T);
  free(atmos->ne);
  free(atmos->vturb);
  free(atmos->nHtot);
  free(geometry->vel);
  free(geometry->height);
  freeMatrix((void **) atmos->nH);

  return; 

}

/* ------- end ---------------------------- close_ncdf   ------------- */

int main(int argc, char **argv)
{


  /* MPI stuff. */
  int mpi_namelen;		
  char mpi_name[MPI_MAX_PROCESSOR_NAME];
  int mpi_size, mpi_rank;

  NCDF_Atmos_file infile;
  int i;


  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Get_processor_name(mpi_name, &mpi_namelen);


  infile.file_name = FILE_NAME;
  
  printf("-- Start\n");

  init_ncdf(&atmos, &geometry, &infile);
  readAtmos_ncdf(mpi_rank, mpi_rank, &atmos, &geometry, &infile);

  for(i=0; i<10; i++){
    //printf("Process %d: temp[%d] = %e\n",mpi_rank,i,atmos.T[i]);
    //printf("Process %d: vz[%d] = %e\n",mpi_rank,i,geometry.vel[i]);
    printf("P %d: nh[0,%d]=%10.3e nh[1,%d]=%10.3e nh[2,%d]=%10.3e nh[3,%d]=%10.3e nh[4,%d]=%10.3e nh[5,%d]=%10.3e nhTot[%d]=%10.3e\n",
    	   mpi_rank,i,atmos.nH[0][i], i,atmos.nH[1][i], i,atmos.nH[2][i],
                    i,atmos.nH[3][i], i,atmos.nH[4][i], i,atmos.nH[5][i],
	            i,atmos.nHtot[i]);
  }
  

  close_ncdf(&atmos, &geometry, &infile);
  
  printf("-- End\n");
  MPI_Finalize();

  return 0;
}
