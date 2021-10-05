/*
    Functions for HDF5 reading of input atmosphere
*/

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
#include "parallel.h"
#include "io.h"

#define FAIL -1
#define MULTI_COMMENT_CHAR  "*"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */
extern MPI_data mpi;
extern InputData input;
extern char messageStr[];


/* ------- begin --------------------------   init_hdf5_atmos   ----- */
void init_hdf5_atmos(Atmosphere *atmos, Geometry *geometry,
    Input_Atmos_file *infile) {
  /* Initialises the input atmosphere file, gets dimensions and variable ids.
   Also performs other basic RH initialisations like readAbundance       */
  const char routineName[] = "init_hdf5_atmos";
  struct  stat statBuffer;
  int has_B;
  hid_t plist_id, ncid;
  hsize_t dims[5];
  size_t nn = 0;
  size_t start[] = {0, 0};
  size_t count[] = {1, 1};
  char *filename;

  /* --- Open input file for model atmosphere --       -------------- */
  if ((plist_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) HERR(routineName);
  if ((H5Pset_fapl_mpio(plist_id, mpi.comm, mpi.info)) < 0) HERR(routineName);
  if ((ncid = H5Fopen(input.atmos_input, H5F_ACC_RDONLY, plist_id)) < 0)
    HERR(routineName);
  infile->ncid = ncid;
  if ((H5Pclose(plist_id)) < 0) HERR(routineName); /* plist no longer needed */

  /* Is magnetic field included? */
  if ((H5LTget_attribute_int(ncid, "/", "has_B", &has_B)) < 0)
    HERR(routineName);

  atmos->Stokes = FALSE;
  if ((has_B) && (strcmp(input.Stokes_input, "none"))) atmos->Stokes = TRUE;

  /* Get dimensions from hydrogen_populations (if present) */
  if (H5LTfind_dataset(ncid, "hydrogen_populations")) {
    if ((H5LTget_dataset_info(ncid, "hydrogen_populations", dims,
			      NULL, NULL)) < 0)
      HERR(routineName);
    nn = dims[1];
    infile->nx = dims[2];
    infile->ny = dims[3];
    infile->nz = dims[4];
  } else {
    sprintf(messageStr, "Atmosphere file missing hydrogen_populations"
	    " array, aborting.\n");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* get some values in atmos/geometry structures */
  geometry->Ndep = (int)  infile->nz;
  atmos->Nspace  = (long) infile->nz;
  atmos->NHydr   = (int)  nn;

  if (atmos->NHydr < 2  &&  !atmos->H_LTE) {
    sprintf(messageStr, "NHydr has to be at least 2, not %d to run with"
	    " NLTE hydrogen populations in background\n", atmos->NHydr);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* Get the varids (collective operation!) */
  if (input.solve_ne == NONE) {
      if (H5LTfind_dataset(ncid, NE_NAME)) {
        if ((infile->ne_varid = H5Dopen2(ncid, NE_NAME, H5P_DEFAULT)) < 0)
          HERR(routineName);
      } else {
          sprintf(messageStr, "Electron density not found in atmosphere file,"
                  " calculating it assuming Saha-Boltzmann ionisation.\n");
          Error(WARNING, routineName, messageStr);
          input.solve_ne = ONCE;
      }
  }
  if ((infile->vz_varid = H5Dopen2(ncid, VZ_NAME, H5P_DEFAULT)) < 0)
    HERR(routineName);
  if ((infile->nh_varid = H5Dopen2(ncid, NH_NAME, H5P_DEFAULT)) < 0)
    HERR(routineName);
  if ((infile->T_varid = H5Dopen2(ncid, TEMP_NAME, H5P_DEFAULT)) < 0)
    HERR(routineName);
  if ((infile->z_varid = H5Dopen2(ncid, "z", H5P_DEFAULT)) < 0)
    HERR(routineName);

  /* Microturbulence, get ID if variable found */
  if (H5LTfind_dataset(ncid, VTURB_NAME)) {
    if ((infile->vturb_varid = H5Dopen2(ncid, VTURB_NAME, H5P_DEFAULT)) < 0)
      HERR(routineName);
  } else {
    infile->vturb_varid = -1;
  }

  /* Find out if z is written for each (x, y) column */
  if ((H5LTget_dataset_ndims(ncid, "z", &mpi.ndims_z)) < 0) HERR(routineName);

  if (atmos->Stokes) {
    if ((infile->Bx_varid = H5Dopen2(ncid, BX_NAME, H5P_DEFAULT)) < 0)
      HERR(routineName);
    if ((infile->By_varid = H5Dopen2(ncid, BY_NAME, H5P_DEFAULT)) < 0)
      HERR(routineName);
    if ((infile->Bz_varid = H5Dopen2(ncid, BZ_NAME, H5P_DEFAULT)) < 0)
      HERR(routineName);
  }

  /* read things that don't depend on x, y */
  start[0] = input.p15d_nt;
  count[0] = 1;
  start[1] = 0;
  count[1] = infile->nz;

  geometry->height = (double *) malloc(infile->nz * sizeof(double));

  infile->y   = (double *) malloc(infile->ny * sizeof(double));
  if ((H5LTread_dataset_double(ncid, "y", infile->y)) < 0) HERR(routineName);
  infile->x   = (double *) malloc(infile->nx * sizeof(double));
  if ((H5LTread_dataset_double(ncid, "x", infile->x)) < 0) HERR(routineName);
  if ((H5LTread_dataset_int(ncid, SNAPNAME, &mpi.snap_number)) < 0)
    HERR(routineName);

  /* allocate arrays */
  geometry->vel = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->T      = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->ne     = (double *) malloc(atmos->Nspace * sizeof(double));
  atmos->vturb  = (double *) calloc(atmos->Nspace , sizeof(double));
  atmos->nHtot  = (double *) malloc(atmos->Nspace * sizeof(double));

  if (atmos->Stokes) {
    atmos->B       = (double *) malloc(atmos->Nspace * sizeof(double));
    atmos->gamma_B = (double *) malloc(atmos->Nspace * sizeof(double));
    atmos->chi_B   = (double *) malloc(atmos->Nspace * sizeof(double));
  }

  /* get boundary conditions from file, assume zero at top and
   * thermalised at bottom if attribute not found */
  if (H5LTfind_attribute(ncid, BTOP_NAME)) {
    if ((H5LTget_attribute_uint(ncid, "/", BTOP_NAME,
				&geometry->vboundary[TOP])) < 0)
      HERR(routineName);
  } else {
    geometry->vboundary[TOP] = ZERO;
  }
  if (H5LTfind_attribute(ncid, BBOT_NAME)) {
    if ((H5LTget_attribute_uint(ncid, "/", BBOT_NAME,
				&geometry->vboundary[BOTTOM])) < 0)
      HERR(routineName);
  } else {
    geometry->vboundary[BOTTOM] = THERMALIZED;
  }
  /* some other housekeeping */
  geometry->scale = GEOMETRIC;
  /* --- Construct atmosID from filename and last modification date - */
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
}
/* ------- end ---------------------------- init_hdf5_atmos  -------- */

/* ------- begin -------------------------- readAtmos_hdf5  --------- */
void readAtmos_hdf5(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
		    Input_Atmos_file *infile) {
  /* Reads the variables T, ne, vel, nh for a given (xi,yi) pair */
  const char routineName[] = "readAtmos_hdf5";
  hsize_t     start[]    = {0, 0, 0, 0}; /* starting values */
  hsize_t     count[]    = {1, 1, 1, 1};
  hsize_t     start_nh[] = {0, 0, 0, 0, 0};
  hsize_t     count_nh[] = {1, 1, 1, 1, 1};
  hsize_t    dims_memory[2];
  hid_t      ncid, dataspace_id, memspace_id;
  int        ierror, i, j;
  bool_t     old_moving;
  double    *Bx, *By, *Bz;

  ncid = infile->ncid;
  atmos->Nspace = geometry->Ndep = infile->nz;

  /* read full T column, to see where to zcut */
  /* Memory dataspace */
  dims_memory[0] = infile->nz;
  if ((memspace_id = H5Screate_simple(1, dims_memory, NULL)) < 0)
    HERR(routineName);
  atmos->T = (double *) realloc(atmos->T, infile->nz * sizeof(double));
  /* File dataspace */
  start[0] = input.p15d_nt; count[0] = 1;
  start[1] = (size_t) xi;   count[1] = 1;
  start[2] = (size_t) yi;   count[2] = 1;
  start[3] = 0;             count[3] = infile->nz;
  if ((dataspace_id = H5Dget_space(infile->T_varid)) < 0) HERR(routineName);
  if ((H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start,
                               NULL, count, NULL)) < 0) HERR(routineName);
  if ((H5Dread(infile->T_varid, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id,
	      H5P_DEFAULT, atmos->T)) < 0) HERR(routineName);
  if (( H5Sclose(dataspace_id) ) < 0) HERR(routineName);
  if (( H5Sclose(memspace_id) ) < 0) HERR(routineName);
  /* Finds z value for Tmax cut, redefines Nspace, reallocates arrays */
  /* Tiago: not using this at the moment, only z cut in depth_refine */
  if (input.p15d_zcut) {
    setTcut(atmos, geometry, input.p15d_tmax);
  } else {
    mpi.zcut = 0;
    /*Free vturb in order to avoid successive addition of vturb_add. 
      This is done in realloc_ndep (called by setTcut) if 15D_DEFINE_ZCUT = TRUE*/
    free(atmos->vturb);
    atmos->vturb  = (double *) calloc(atmos->Nspace , sizeof(double));
  }

  /* Memory dataspace, redefine for Nspace */
  dims_memory[0] = atmos->Nspace;
  memspace_id = H5Screate_simple(1, dims_memory, NULL);
  /* Get z again */
  if ((mpi.ndims_z) == 2) {   /* z scale is specified only once per snapshot */
    start[0] = input.p15d_nt; count[0] = 1;
    start[1] = mpi.zcut;      count[1] = atmos->Nspace;
  } else if ((mpi.ndims_z) == 4) {  /* z scale is specified for every column */
    start[0] = input.p15d_nt; count[0] = 1;
    start[1] = (size_t) xi;   count[1] = 1;
    start[2] = (size_t) yi;   count[2] = 1;
    start[3] = mpi.zcut;      count[3] = atmos->Nspace;
  }
  dataspace_id = H5Dget_space(infile->z_varid);
  ierror = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start,
                               NULL, count, NULL);
  if ((ierror = H5Dread(infile->z_varid, H5T_NATIVE_DOUBLE, memspace_id,
		        dataspace_id, H5P_DEFAULT, geometry->height)) < 0)
    HERR(routineName);
  if (( H5Sclose(dataspace_id) ) < 0) HERR(routineName);
  /* redefine dataspace for 4-D variables */
  start[0] = input.p15d_nt; count[0] = 1;
  start[1] = (size_t) xi;   count[1] = 1;
  start[2] = (size_t) yi;   count[2] = 1;
  start[3] = mpi.zcut;      count[3] = atmos->Nspace;
  dataspace_id = H5Dget_space(infile->T_varid);
  if ((ierror = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start,
                               NULL, count, NULL)) < 0) HERR(routineName);
   /* read variables */
  if ((ierror = H5Dread(infile->T_varid, H5T_NATIVE_DOUBLE, memspace_id,
	          dataspace_id, H5P_DEFAULT, atmos->T)) < 0) HERR(routineName);
  if (input.solve_ne == NONE) {
      if ((ierror = H5Dread(infile->ne_varid, H5T_NATIVE_DOUBLE, memspace_id,
      	          dataspace_id, H5P_DEFAULT, atmos->ne)) < 0) HERR(routineName);
  }
  if ((ierror = H5Dread(infile->vz_varid, H5T_NATIVE_DOUBLE, memspace_id,
	     dataspace_id, H5P_DEFAULT, geometry->vel)) < 0) HERR(routineName);
  /* vturb, if available */
  if (infile->vturb_varid != -1) {
    ierror = H5Dread(infile->vturb_varid, H5T_NATIVE_DOUBLE, memspace_id,
		   dataspace_id, H5P_DEFAULT, atmos->vturb);
  }
  /* Read magnetic field */
  if (atmos->Stokes) {
    Bx = (double *) malloc(atmos->Nspace * sizeof(double));
    By = (double *) malloc(atmos->Nspace * sizeof(double));
    Bz = (double *) malloc(atmos->Nspace * sizeof(double));
    /* Read in cartesian coordinates */
    ierror = H5Dread(infile->Bx_varid, H5T_NATIVE_DOUBLE, memspace_id,
		   dataspace_id, H5P_DEFAULT, Bx);
    ierror = H5Dread(infile->By_varid, H5T_NATIVE_DOUBLE, memspace_id,
		   dataspace_id, H5P_DEFAULT, By);
    ierror = H5Dread(infile->Bz_varid, H5T_NATIVE_DOUBLE, memspace_id,
		   dataspace_id, H5P_DEFAULT, Bz);
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
  if (( H5Sclose(dataspace_id) ) < 0) HERR(routineName);
  if (( H5Sclose(memspace_id) ) < 0) HERR(routineName);
  /* allocate and zero nHtot */
  atmos->nH = matrix_double(atmos->NHydr, atmos->Nspace);
  for (j = 0; j < atmos->Nspace; j++) atmos->nHtot[j] = 0.0;

  /* read nH, all at once */
  start_nh[0] = input.p15d_nt; count_nh[0] = 1;
  start_nh[1] = 0;             count_nh[1] = atmos->NHydr;
  start_nh[2] = (size_t) xi;   count_nh[2] = 1;
  start_nh[3] = (size_t) yi;   count_nh[3] = 1;
  start_nh[4] = mpi.zcut;      count_nh[4] = atmos->Nspace;
  dataspace_id = H5Dget_space(infile->nh_varid);
  ierror = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_nh,
                               NULL, count_nh, NULL);
  dims_memory[0] = atmos->NHydr;
  dims_memory[1] = atmos->Nspace;
  memspace_id = H5Screate_simple(2, dims_memory, NULL);
  ierror = H5Dread(infile->nh_varid, H5T_NATIVE_DOUBLE,
		   memspace_id, dataspace_id, H5P_DEFAULT, atmos->nH[0]);
  if (( H5Sclose(dataspace_id) ) < 0) HERR(routineName);
  if (( H5Sclose(memspace_id) ) < 0) HERR(routineName);
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
}
/* ------- end ---------------------------- readAtmos_hdf5  --------- */

/* ------- begin -------------------------- close_hdf5  ------------- */
void close_hdf5_atmos(Atmosphere *atmos, Geometry *geometry,
    Input_Atmos_file *infile) {
  /* Closes the HDF5 file and frees memory */
  int ierror;
  /* Close the file. */
  ierror = H5Dclose(infile->z_varid);
  ierror = H5Dclose(infile->T_varid);
  if (input.solve_ne == NONE) ierror = H5Dclose(infile->ne_varid);
  ierror = H5Dclose(infile->vz_varid);
  ierror = H5Dclose(infile->nh_varid);
  if (atmos->Stokes) {
    ierror = H5Dclose(infile->Bx_varid);
    ierror = H5Dclose(infile->By_varid);
    ierror = H5Dclose(infile->Bz_varid);
  }
  if (infile->vturb_varid != -1) ierror = H5Dclose(infile->vturb_varid);
  ierror = H5Fclose(infile->ncid);
  /* Free stuff */
  free(atmos->T);
  free(atmos->ne);
  free(atmos->vturb);
  free(atmos->nHtot);
  free(geometry->vel);
  free(geometry->height);
  free(infile->y);
  free(infile->x);
  free(geometry->xscale);
  free(geometry->yscale);
}
