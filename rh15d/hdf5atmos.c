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

/* --- Cached file dataspace handles for readAtmos_hdf5 --------------
   Opening the file dataspace via H5Dget_space and closing it again on
   every column was a major source of HDF5 metadata overhead in pool
   mode (~5+ get/select/close per column × thousands of columns).  We
   cache one handle per dataset on the first call and just re-call
   H5Sselect_hyperslab on each subsequent column.  Cleaned up in
   close_hdf5_atmos. --------------------------------------------- */
static int   rfa_initialized   = 0;
static hid_t rfa_T_fspace      = -1;
static hid_t rfa_z_fspace      = -1;
static hid_t rfa_ne_fspace     = -1;
static hid_t rfa_vz_fspace     = -1;
static hid_t rfa_vturb_fspace  = -1;
static hid_t rfa_nh_fspace     = -1;
static hid_t rfa_Bx_fspace     = -1;
static hid_t rfa_By_fspace     = -1;
static hid_t rfa_Bz_fspace     = -1;

static void rfa_init_cached_dataspaces(Atmosphere *atmos,
                                       Input_Atmos_file *infile) {
  const char routineName[] = "rfa_init_cached_dataspaces";
  if ((rfa_T_fspace  = H5Dget_space(infile->T_varid))  < 0) HERR(routineName);
  if ((rfa_z_fspace  = H5Dget_space(infile->z_varid))  < 0) HERR(routineName);
  if ((rfa_vz_fspace = H5Dget_space(infile->vz_varid)) < 0) HERR(routineName);
  if ((rfa_nh_fspace = H5Dget_space(infile->nh_varid)) < 0) HERR(routineName);
  if (input.solve_ne == NONE) {
    if ((rfa_ne_fspace = H5Dget_space(infile->ne_varid)) < 0) HERR(routineName);
  }
  if (infile->vturb_varid != -1) {
    if ((rfa_vturb_fspace = H5Dget_space(infile->vturb_varid)) < 0)
      HERR(routineName);
  }
  if (atmos->Stokes) {
    if ((rfa_Bx_fspace = H5Dget_space(infile->Bx_varid)) < 0) HERR(routineName);
    if ((rfa_By_fspace = H5Dget_space(infile->By_varid)) < 0) HERR(routineName);
    if ((rfa_Bz_fspace = H5Dget_space(infile->Bz_varid)) < 0) HERR(routineName);
  }
  rfa_initialized = 1;
}

static void rfa_close_cached_dataspaces(void) {
  if (!rfa_initialized) return;
  if (rfa_T_fspace     >= 0) H5Sclose(rfa_T_fspace);
  if (rfa_z_fspace     >= 0) H5Sclose(rfa_z_fspace);
  if (rfa_ne_fspace    >= 0) H5Sclose(rfa_ne_fspace);
  if (rfa_vz_fspace    >= 0) H5Sclose(rfa_vz_fspace);
  if (rfa_vturb_fspace >= 0) H5Sclose(rfa_vturb_fspace);
  if (rfa_nh_fspace    >= 0) H5Sclose(rfa_nh_fspace);
  if (rfa_Bx_fspace    >= 0) H5Sclose(rfa_Bx_fspace);
  if (rfa_By_fspace    >= 0) H5Sclose(rfa_By_fspace);
  if (rfa_Bz_fspace    >= 0) H5Sclose(rfa_Bz_fspace);
  rfa_T_fspace = rfa_z_fspace = rfa_ne_fspace = rfa_vz_fspace = -1;
  rfa_vturb_fspace = rfa_nh_fspace = -1;
  rfa_Bx_fspace = rfa_By_fspace = rfa_Bz_fspace = -1;
  rfa_initialized = 0;
}


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

  /* --- Open input file for model atmosphere (Lustre-optimised) ---
     Use independent metadata: in pool mode drones read columns
     asynchronously, so collective metadata ops would deadlock. */
  plist_id = create_hdf5_fapl_indep();

  /* --- Sets IO cache enviction for beter perfomance on parallel io */
  if (IO_CACHE_EVICTION) { // disable metadata cache eviction
    H5AC_cache_config_t mdc_config;
    mdc_config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    if ((H5Pget_mdc_config(plist_id, &mdc_config)) < 0) HERR(routineName);
    mdc_config.evictions_enabled = false;
    mdc_config.incr_mode = H5C_incr__off;
    mdc_config.decr_mode = H5C_decr__off;
    mdc_config.flash_incr_mode = H5C_flash_incr__off;
    H5Pset_mdc_config(plist_id, &mdc_config);
  }

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

/* ------- begin -------------------------- readAtmos_hdf5  ---------
   Reads T, ne, vz, vturb, nH (and B if Stokes) for a given (xi,yi).

   Tier 2 I/O optimisations relative to the original version:
   (a) Eliminated the redundant second T read.  T is read once at the
       full nz extent; if zcut > 0 the data is shifted in place with
       memmove and the per-column arrays are realloced.
   (b) File dataspaces are cached on first call (rfa_*_fspace) so
       H5Dget_space / H5Sclose are not called per-column.
   (c) The remaining 5–8 reads are issued in a single H5Dread_multi
       call, collapsing 5–8 collective sync points into one.       */
void readAtmos_hdf5(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
                    Input_Atmos_file *infile) {
  const char routineName[] = "readAtmos_hdf5";
  hsize_t  start4[]   = {0, 0, 0, 0};
  hsize_t  count4[]   = {1, 1, 1, 0};
  hsize_t  start_z[]  = {0, 0, 0, 0};
  hsize_t  count_z[]  = {1, 1, 1, 0};
  hsize_t  start_nh[] = {0, 0, 0, 0, 0};
  hsize_t  count_nh[] = {1, 1, 1, 1, 1};
  hsize_t  dims_mem[2];
  hid_t    plist_id;
  hid_t    mem_T, mem_1d, mem_nh;
  bool_t   old_moving;
  double  *Bx = NULL, *By = NULL, *Bz = NULL;
  int      i, j;

  /* --- Lazy init of cached file dataspaces --- */
  if (!rfa_initialized) rfa_init_cached_dataspaces(atmos, infile);

  /* --- Transfer property list (collective vs independent) --- */
  if (COLLECTIVE_IO_R && mpi.isbalanced) {
    if ((plist_id = H5Pcreate(H5P_DATASET_XFER)) < 0) HERR(routineName);
    if (H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE) < 0)
      HERR(routineName);
  } else {
    plist_id = H5P_DEFAULT;
  }

  atmos->Nspace = geometry->Ndep = infile->nz;

  /* === Step 1: read full T column (single read, used to find zcut) === */
  start4[0] = input.p15d_nt; count4[0] = 1;
  start4[1] = (size_t) xi;   count4[1] = 1;
  start4[2] = (size_t) yi;   count4[2] = 1;
  start4[3] = 0;             count4[3] = infile->nz;

  if (H5Sselect_hyperslab(rfa_T_fspace, H5S_SELECT_SET, start4,
                          NULL, count4, NULL) < 0) HERR(routineName);
  dims_mem[0] = infile->nz;
  if ((mem_T = H5Screate_simple(1, dims_mem, NULL)) < 0) HERR(routineName);
  atmos->T = (double *) realloc(atmos->T, infile->nz * sizeof(double));
  if (H5Dread(infile->T_varid, H5T_NATIVE_DOUBLE, mem_T, rfa_T_fspace,
              plist_id, atmos->T) < 0) HERR(routineName);
  if (H5Sclose(mem_T) < 0) HERR(routineName);

  /* === Step 2: determine zcut, shift T in place, resize columns === */
  mpi.zcut = 0;
  if (input.p15d_zcut && input.p15d_tmax >= 0) {
    int cut = -1;
    for (i = 0; i < infile->nz; i++) {
      if (atmos->T[i] <= input.p15d_tmax) { cut = i; break; }
    }
    if (cut < 0) {
      sprintf(messageStr,
              "\n-Could not find temperature cut point! Aborting.\n");
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    mpi.zcut = cut;
  }
  if (mpi.zcut > 0) {
    long keep = (long)infile->nz - mpi.zcut;
    memmove(atmos->T, atmos->T + mpi.zcut, keep * sizeof(double));
  }
  atmos->Nspace  = (long)infile->nz - mpi.zcut;
  geometry->Ndep = (int) atmos->Nspace;
  realloc_ndep(atmos, geometry);

  /* === Step 3: build the multi-read for everything except T === */
  hid_t   dsets  [10];
  hid_t   types  [10];
  hid_t   fspaces[10];
  hid_t   mspaces[10];
  void   *bufs   [10];
  size_t  ndset = 0;

  /* 1D memory dataspace, shared by all 4D variables */
  dims_mem[0] = atmos->Nspace;
  if ((mem_1d = H5Screate_simple(1, dims_mem, NULL)) < 0) HERR(routineName);

  /* Update 4D hyperslab to apply zcut */
  start4[3] = mpi.zcut;
  count4[3] = atmos->Nspace;

  /* z (different shape: 2D or 4D depending on file layout) */
  if (mpi.ndims_z == 2) {
    start_z[0] = input.p15d_nt; count_z[0] = 1;
    start_z[1] = mpi.zcut;      count_z[1] = atmos->Nspace;
  } else if (mpi.ndims_z == 4) {
    start_z[0] = input.p15d_nt; count_z[0] = 1;
    start_z[1] = (size_t) xi;   count_z[1] = 1;
    start_z[2] = (size_t) yi;   count_z[2] = 1;
    start_z[3] = mpi.zcut;      count_z[3] = atmos->Nspace;
  }
  if (H5Sselect_hyperslab(rfa_z_fspace, H5S_SELECT_SET, start_z,
                          NULL, count_z, NULL) < 0) HERR(routineName);
  dsets  [ndset] = infile->z_varid;
  types  [ndset] = H5T_NATIVE_DOUBLE;
  fspaces[ndset] = rfa_z_fspace;
  mspaces[ndset] = mem_1d;
  bufs   [ndset] = geometry->height;
  ndset++;

  /* ne */
  if (input.solve_ne == NONE) {
    if (H5Sselect_hyperslab(rfa_ne_fspace, H5S_SELECT_SET, start4,
                            NULL, count4, NULL) < 0) HERR(routineName);
    dsets  [ndset] = infile->ne_varid;
    types  [ndset] = H5T_NATIVE_DOUBLE;
    fspaces[ndset] = rfa_ne_fspace;
    mspaces[ndset] = mem_1d;
    bufs   [ndset] = atmos->ne;
    ndset++;
  }

  /* vz */
  if (H5Sselect_hyperslab(rfa_vz_fspace, H5S_SELECT_SET, start4,
                          NULL, count4, NULL) < 0) HERR(routineName);
  dsets  [ndset] = infile->vz_varid;
  types  [ndset] = H5T_NATIVE_DOUBLE;
  fspaces[ndset] = rfa_vz_fspace;
  mspaces[ndset] = mem_1d;
  bufs   [ndset] = geometry->vel;
  ndset++;

  /* vturb (if present in file) */
  if (infile->vturb_varid != -1) {
    if (H5Sselect_hyperslab(rfa_vturb_fspace, H5S_SELECT_SET, start4,
                            NULL, count4, NULL) < 0) HERR(routineName);
    dsets  [ndset] = infile->vturb_varid;
    types  [ndset] = H5T_NATIVE_DOUBLE;
    fspaces[ndset] = rfa_vturb_fspace;
    mspaces[ndset] = mem_1d;
    bufs   [ndset] = atmos->vturb;
    ndset++;
  }

  /* Magnetic field (Cartesian, converted to spherical after read) */
  if (atmos->Stokes) {
    Bx = (double *) malloc(atmos->Nspace * sizeof(double));
    By = (double *) malloc(atmos->Nspace * sizeof(double));
    Bz = (double *) malloc(atmos->Nspace * sizeof(double));
    if (H5Sselect_hyperslab(rfa_Bx_fspace, H5S_SELECT_SET, start4,
                            NULL, count4, NULL) < 0) HERR(routineName);
    if (H5Sselect_hyperslab(rfa_By_fspace, H5S_SELECT_SET, start4,
                            NULL, count4, NULL) < 0) HERR(routineName);
    if (H5Sselect_hyperslab(rfa_Bz_fspace, H5S_SELECT_SET, start4,
                            NULL, count4, NULL) < 0) HERR(routineName);
    dsets  [ndset] = infile->Bx_varid; types[ndset] = H5T_NATIVE_DOUBLE;
    fspaces[ndset] = rfa_Bx_fspace;    mspaces[ndset] = mem_1d;
    bufs   [ndset] = Bx; ndset++;
    dsets  [ndset] = infile->By_varid; types[ndset] = H5T_NATIVE_DOUBLE;
    fspaces[ndset] = rfa_By_fspace;    mspaces[ndset] = mem_1d;
    bufs   [ndset] = By; ndset++;
    dsets  [ndset] = infile->Bz_varid; types[ndset] = H5T_NATIVE_DOUBLE;
    fspaces[ndset] = rfa_Bz_fspace;    mspaces[ndset] = mem_1d;
    bufs   [ndset] = Bz; ndset++;
  }

  /* nH (5D, separate memspace) */
  start_nh[0] = input.p15d_nt; count_nh[0] = 1;
  start_nh[1] = 0;             count_nh[1] = atmos->NHydr;
  start_nh[2] = (size_t) xi;   count_nh[2] = 1;
  start_nh[3] = (size_t) yi;   count_nh[3] = 1;
  start_nh[4] = mpi.zcut;      count_nh[4] = atmos->Nspace;
  if (H5Sselect_hyperslab(rfa_nh_fspace, H5S_SELECT_SET, start_nh,
                          NULL, count_nh, NULL) < 0) HERR(routineName);
  dims_mem[0] = atmos->NHydr;
  dims_mem[1] = atmos->Nspace;
  if ((mem_nh = H5Screate_simple(2, dims_mem, NULL)) < 0) HERR(routineName);
  /* (existing behaviour: nH matrix is reallocated each call) */
  atmos->nH = matrix_double(atmos->NHydr, atmos->Nspace);
  dsets  [ndset] = infile->nh_varid;
  types  [ndset] = H5T_NATIVE_DOUBLE;
  fspaces[ndset] = rfa_nh_fspace;
  mspaces[ndset] = mem_nh;
  bufs   [ndset] = atmos->nH[0];
  ndset++;

  /* === Step 4: single multi-read for all remaining datasets === */
  if (H5Dread_multi(ndset, dsets, types, mspaces, fspaces,
                    plist_id, bufs) < 0) HERR(routineName);

  if (H5Sclose(mem_1d) < 0) HERR(routineName);
  if (H5Sclose(mem_nh) < 0) HERR(routineName);
  if (plist_id != H5P_DEFAULT) {
    if (H5Pclose(plist_id) < 0) HERR(routineName);
  }

  /* === Step 5: post-processing === */
  if (atmos->Stokes) {
    for (j = 0; j < atmos->Nspace; j++) {
      atmos->B[j]       = sqrt(SQ(Bx[j]) + SQ(By[j]) + SQ(Bz[j]));
      atmos->gamma_B[j] = acos(Bz[j]/atmos->B[j]);
      atmos->chi_B[j]   = atan(By[j]/Bx[j]);
      if ((Bx[j] == 0) && (By[j] == 0) && (Bz[j] == 0))
        atmos->gamma_B[j] = 0.0;
      if ((Bx[j] == 0) && (By[j] == 0))
        atmos->chi_B[j]   = 1.0;
    }
    free(Bx); free(By); free(Bz);
  }

  for (j = 0; j < atmos->Nspace; j++) atmos->nHtot[j] = 0.0;

  /* Depth grid refinement */
  if (input.p15d_refine)
    depth_refine(atmos, geometry, input.p15d_tmax);

  /* Fix vturb: remove zeros, apply multiplier and offset */
  for (i = 0; i < atmos->Nspace; i++) {
    if (atmos->vturb[i] < 0.0) atmos->vturb[i] = 0.0;
    atmos->vturb[i] = atmos->vturb[i] * input.vturb_mult + input.vturb_add;
  }

  /* Sum to get nHtot */
  for (i = 0; i < atmos->NHydr; i++)
    for (j = 0; j < atmos->Nspace; j++)
      atmos->nHtot[j] += atmos->nH[i][j];

  /* Moving-atmosphere check */
  old_moving = atmos->moving;
  atmos->moving = FALSE;
  for (i = 0; i < atmos->Nspace; i++) {
    if (fabs(geometry->vel[i]) >= atmos->vmacro_tresh) {
      atmos->moving = TRUE;
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
  /* Release cached file dataspaces from readAtmos_hdf5 */
  rfa_close_cached_dataspaces();
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
