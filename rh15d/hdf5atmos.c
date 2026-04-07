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
static void node_cache_init(Atmosphere *atmos, Input_Atmos_file *infile);
static void node_cache_free(void);

/* --- Global variables --                             -------------- */
extern MPI_data mpi;
extern InputData input;
extern char messageStr[];

/* ==================================================================
   Node-shared atmosphere cache
   ------------------------------------------------------------------
   When input.use_node_atmos_cache is TRUE, each node reads its
   assigned slab of the atmosphere ONCE at startup into a set of
   MPI-3 shared-memory windows (one per variable).  readAtmos then
   serves columns from this cache instead of touching the file — no
   per-column HDF5 I/O in the hot path.

   Layout per variable (reduced coordinates):
       T, ne, vz, vturb, [Bx, By, Bz]:
           [n_local_rows × mpi.ny × nz] of doubles, row-major
       nH:
           [NHydr × n_local_rows × mpi.ny × nz] of doubles
       z:  [n_local_rows × mpi.ny × nz] (4D file) or [nz] (2D file)

   A node owns reduced rows ix ∈ [mpi.node_ix0, mpi.node_ix1).
   Local index for reduced (ix, iy): (ix - node_ix0) * mpi.ny + iy
   File → reduced: xi = p15d_x0 + ix*p15d_xst (and analogously for y).
   ================================================================== */

typedef struct {
  bool_t    initialized;
  int       nrows;                /* node's owned reduced row count    */
  int       ix0;                  /* first reduced row (=node_ix0)     */
  int       ny;                   /* reduced ny                        */
  long      nz;                   /* full z extent (same as infile.nz) */
  int       NHydr;
  bool_t    stokes;
  bool_t    has_vturb;
  int       ndims_z;              /* 2 or 4                            */

  /* Pointers into shared-memory windows */
  double   *T;
  double   *ne;
  double   *vz;
  double   *vturb;
  double   *z;       /* size depends on ndims_z */
  double   *nH;      /* NHydr-major */
  double   *Bx;
  double   *By;
  double   *Bz;

  /* Shared-memory window handles (only node-rank 0 does the alloc,
     but every rank keeps the handle so MPI_Win_free can be called). */
  MPI_Win   win_T;
  MPI_Win   win_ne;
  MPI_Win   win_vz;
  MPI_Win   win_vturb;
  MPI_Win   win_z;
  MPI_Win   win_nH;
  MPI_Win   win_Bx;
  MPI_Win   win_By;
  MPI_Win   win_Bz;
} AtmosNodeCache;

static AtmosNodeCache node_cache = { .initialized = FALSE };

/* --- Allocate one shared-memory window on node_comm and return
       the pointer that every rank on the node can use.  Only
       node-rank 0 contributes non-zero bytes; other ranks allocate
       zero and immediately MPI_Win_shared_query to get the mapping. */
static double *alloc_shared_double(MPI_Comm node_comm, int node_rank,
                                   size_t n_doubles, MPI_Win *out_win) {
  const char routineName[] = "alloc_shared_double";
  MPI_Aint  alloc_bytes = (node_rank == 0)
                          ? (MPI_Aint)(n_doubles * sizeof(double)) : 0;
  double   *base_ptr = NULL;

  if (MPI_Win_allocate_shared(alloc_bytes, (int)sizeof(double),
                              MPI_INFO_NULL, node_comm,
                              &base_ptr, out_win) != MPI_SUCCESS) {
    sprintf(messageStr,
            "MPI_Win_allocate_shared failed for %.3f MiB\n",
            (double)alloc_bytes / (1024.0 * 1024.0));
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* Non-zero ranks map the memory allocated by rank 0 */
  if (node_rank != 0) {
    MPI_Aint sz;
    int      disp_unit;
    if (MPI_Win_shared_query(*out_win, 0, &sz, &disp_unit,
                             &base_ptr) != MPI_SUCCESS) {
      Error(ERROR_LEVEL_2, routineName,
            "MPI_Win_shared_query failed\n");
    }
  }
  return base_ptr;
}


/* --- Strided-hyperslab read of one 4D variable [nt, nx, ny, nz] ---
   Reads this node's row slab from the file into `dst` (already sized
   nrows * ny_red * nz doubles).  Only called by node-rank 0.       */
static void read_slab_4d(hid_t dset, hid_t file_type,
                         int node_ix0, int node_ix1,
                         long nz, double *dst) {
  const char routineName[] = "read_slab_4d";
  hid_t    fspace, mspace;
  hsize_t  start[4], stride[4], count[4];
  hsize_t  mem_dims[3];
  int      nrows = node_ix1 - node_ix0;

  if ((fspace = H5Dget_space(dset)) < 0) HERR(routineName);

  start[0]  = input.p15d_nt;
  start[1]  = (hsize_t)(input.p15d_x0 + node_ix0 * input.p15d_xst);
  start[2]  = (hsize_t) input.p15d_y0;
  start[3]  = 0;
  stride[0] = 1;
  stride[1] = (hsize_t) input.p15d_xst;
  stride[2] = (hsize_t) input.p15d_yst;
  stride[3] = 1;
  count[0]  = 1;
  count[1]  = (hsize_t) nrows;
  count[2]  = (hsize_t) mpi.ny;
  count[3]  = (hsize_t) nz;

  if (H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, stride,
                          count, NULL) < 0) HERR(routineName);

  mem_dims[0] = (hsize_t) nrows;
  mem_dims[1] = (hsize_t) mpi.ny;
  mem_dims[2] = (hsize_t) nz;
  if ((mspace = H5Screate_simple(3, mem_dims, NULL)) < 0) HERR(routineName);

  if (H5Dread(dset, file_type, mspace, fspace, H5P_DEFAULT, dst) < 0)
    HERR(routineName);

  if (H5Sclose(mspace) < 0) HERR(routineName);
  if (H5Sclose(fspace) < 0) HERR(routineName);
}


/* --- Strided read of nH [nt, NHydr, nx, ny, nz] --- */
static void read_slab_nh(hid_t dset, int node_ix0, int node_ix1,
                         int NHydr, long nz, double *dst) {
  const char routineName[] = "read_slab_nh";
  hid_t    fspace, mspace;
  hsize_t  start[5], stride[5], count[5];
  hsize_t  mem_dims[4];
  int      nrows = node_ix1 - node_ix0;

  if ((fspace = H5Dget_space(dset)) < 0) HERR(routineName);

  start[0] = input.p15d_nt; start[1] = 0;
  start[2] = (hsize_t)(input.p15d_x0 + node_ix0 * input.p15d_xst);
  start[3] = (hsize_t) input.p15d_y0;
  start[4] = 0;
  stride[0] = 1;                          stride[1] = 1;
  stride[2] = (hsize_t) input.p15d_xst;   stride[3] = (hsize_t) input.p15d_yst;
  stride[4] = 1;
  count[0] = 1;                           count[1] = (hsize_t) NHydr;
  count[2] = (hsize_t) nrows;             count[3] = (hsize_t) mpi.ny;
  count[4] = (hsize_t) nz;

  if (H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, stride,
                          count, NULL) < 0) HERR(routineName);

  mem_dims[0] = (hsize_t) NHydr;
  mem_dims[1] = (hsize_t) nrows;
  mem_dims[2] = (hsize_t) mpi.ny;
  mem_dims[3] = (hsize_t) nz;
  if ((mspace = H5Screate_simple(4, mem_dims, NULL)) < 0) HERR(routineName);

  if (H5Dread(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, dst) < 0)
    HERR(routineName);

  if (H5Sclose(mspace) < 0) HERR(routineName);
  if (H5Sclose(fspace) < 0) HERR(routineName);
}


/* --- Read z depending on its dimensionality (2D or 4D in file) --- */
static void read_slab_z(hid_t dset, int node_ix0, int node_ix1,
                        long nz, int ndims_z, double *dst) {
  const char routineName[] = "read_slab_z";
  hid_t    fspace, mspace;

  if ((fspace = H5Dget_space(dset)) < 0) HERR(routineName);

  if (ndims_z == 2) {
    /* z[nt, nz] — shared for the whole snapshot */
    hsize_t start[2] = {(hsize_t)input.p15d_nt, 0};
    hsize_t count[2] = {1, (hsize_t)nz};
    hsize_t mem_dims[1] = {(hsize_t)nz};

    if (H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL,
                            count, NULL) < 0) HERR(routineName);
    if ((mspace = H5Screate_simple(1, mem_dims, NULL)) < 0) HERR(routineName);
    if (H5Dread(dset, H5T_NATIVE_DOUBLE, mspace, fspace,
                H5P_DEFAULT, dst) < 0) HERR(routineName);
    if (H5Sclose(mspace) < 0) HERR(routineName);
  } else {
    /* z[nt, nx, ny, nz] — strided like the 4D variables */
    read_slab_4d(dset, H5T_NATIVE_DOUBLE, node_ix0, node_ix1, nz, dst);
    if (H5Sclose(fspace) < 0) HERR(routineName);
    return;
  }
  if (H5Sclose(fspace) < 0) HERR(routineName);
}


/* --- Initialize the node-shared cache.  Called from init_hdf5_atmos
       when input.use_node_atmos_cache is TRUE. --------------------- */
static void node_cache_init(Atmosphere *atmos, Input_Atmos_file *infile) {
  const char routineName[] = "node_cache_init";
  AtmosNodeCache *c = &node_cache;
  int    nrows;
  size_t n_1d, n_nh, n_z;

  if (c->initialized) return;

  nrows = mpi.node_ix1 - mpi.node_ix0;
  if (nrows <= 0) {
    sprintf(messageStr,
            "node_cache_init: node %d has zero rows (ix0=%d ix1=%d)\n",
            mpi.node_id, mpi.node_ix0, mpi.node_ix1);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  c->nrows     = nrows;
  c->ix0       = mpi.node_ix0;
  c->ny        = mpi.ny;
  c->nz        = (long) infile->nz;
  c->NHydr     = atmos->NHydr;
  c->stokes    = atmos->Stokes;
  c->has_vturb = (infile->vturb_varid != -1);
  c->ndims_z   = mpi.ndims_z;

  n_1d = (size_t)nrows * (size_t)mpi.ny * (size_t)c->nz;
  n_nh = (size_t)c->NHydr * n_1d;
  n_z  = (c->ndims_z == 2) ? (size_t)c->nz : n_1d;

  /* --- Allocate shared windows --- */
  c->T  = alloc_shared_double(mpi.node_comm, mpi.node_rank, n_1d, &c->win_T);
  c->ne = alloc_shared_double(mpi.node_comm, mpi.node_rank, n_1d, &c->win_ne);
  c->vz = alloc_shared_double(mpi.node_comm, mpi.node_rank, n_1d, &c->win_vz);
  if (c->has_vturb) {
    c->vturb = alloc_shared_double(mpi.node_comm, mpi.node_rank, n_1d,
                                   &c->win_vturb);
  }
  c->z  = alloc_shared_double(mpi.node_comm, mpi.node_rank, n_z,  &c->win_z);
  c->nH = alloc_shared_double(mpi.node_comm, mpi.node_rank, n_nh, &c->win_nH);
  if (c->stokes) {
    c->Bx = alloc_shared_double(mpi.node_comm, mpi.node_rank, n_1d, &c->win_Bx);
    c->By = alloc_shared_double(mpi.node_comm, mpi.node_rank, n_1d, &c->win_By);
    c->Bz = alloc_shared_double(mpi.node_comm, mpi.node_rank, n_1d, &c->win_Bz);
  }

  /* --- Only node-rank 0 reads from the HDF5 file; others wait --- */
  if (mpi.node_rank == 0) {
    read_slab_4d(infile->T_varid,  H5T_NATIVE_DOUBLE,
                 mpi.node_ix0, mpi.node_ix1, c->nz, c->T);
    if (input.solve_ne == NONE) {
      read_slab_4d(infile->ne_varid, H5T_NATIVE_DOUBLE,
                   mpi.node_ix0, mpi.node_ix1, c->nz, c->ne);
    }
    read_slab_4d(infile->vz_varid, H5T_NATIVE_DOUBLE,
                 mpi.node_ix0, mpi.node_ix1, c->nz, c->vz);
    if (c->has_vturb) {
      read_slab_4d(infile->vturb_varid, H5T_NATIVE_DOUBLE,
                   mpi.node_ix0, mpi.node_ix1, c->nz, c->vturb);
    }
    read_slab_z(infile->z_varid, mpi.node_ix0, mpi.node_ix1, c->nz,
                c->ndims_z, c->z);
    read_slab_nh(infile->nh_varid, mpi.node_ix0, mpi.node_ix1,
                 c->NHydr, c->nz, c->nH);
    if (c->stokes) {
      read_slab_4d(infile->Bx_varid, H5T_NATIVE_DOUBLE,
                   mpi.node_ix0, mpi.node_ix1, c->nz, c->Bx);
      read_slab_4d(infile->By_varid, H5T_NATIVE_DOUBLE,
                   mpi.node_ix0, mpi.node_ix1, c->nz, c->By);
      read_slab_4d(infile->Bz_varid, H5T_NATIVE_DOUBLE,
                   mpi.node_ix0, mpi.node_ix1, c->nz, c->Bz);
    }
  }
  /* Barrier so other ranks only read memory after rank 0 has filled it. */
  MPI_Barrier(mpi.node_comm);

  c->initialized = TRUE;

  if (mpi.rank == 0) {
    double gib_per_node = ((double)n_1d * (c->stokes ? 7 : 4)  /* T,ne,vz,vturb?,B* */
                          + (double)n_nh
                          + (double)n_z)
                         * sizeof(double) / (1024.0 * 1024.0 * 1024.0);
    sprintf(messageStr,
            "node cache: %d nodes, %d rows/node (first), %d ny, %ld nz, "
            "%.3f GiB per node\n",
            mpi.n_nodes, nrows, mpi.ny, c->nz, gib_per_node);
    Error(MESSAGE, routineName, messageStr);
  }
}


/* --- Free shared windows --- */
static void node_cache_free(void) {
  AtmosNodeCache *c = &node_cache;
  if (!c->initialized) return;
  MPI_Win_free(&c->win_T);
  MPI_Win_free(&c->win_ne);
  MPI_Win_free(&c->win_vz);
  if (c->has_vturb) MPI_Win_free(&c->win_vturb);
  MPI_Win_free(&c->win_z);
  MPI_Win_free(&c->win_nH);
  if (c->stokes) {
    MPI_Win_free(&c->win_Bx);
    MPI_Win_free(&c->win_By);
    MPI_Win_free(&c->win_Bz);
  }
  memset(c, 0, sizeof(*c));
  c->initialized = FALSE;
}


/* --- Does a reduced (ix, iy) fall in this node's cache? --- */
static inline bool_t node_cache_owns(int ix, int iy) {
  (void) iy;
  return (node_cache.initialized &&
          ix >= node_cache.ix0 &&
          ix <  node_cache.ix0 + node_cache.nrows);
}

/* --- Local byte offset of a column in the 1D shared buffers --- */
static inline long node_cache_col_offset(int ix, int iy) {
  long local = (long)(ix - node_cache.ix0) * (long)node_cache.ny + (long)iy;
  return local * node_cache.nz;
}

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


/* ------- begin -------------------------- init_atmos_node_cache ---
   Public wrapper that populates the node-shared atmosphere cache.
   Must be called AFTER distribute_jobs has run (needs
   mpi.node_ix0/ix1, mpi.ny populated) but BEFORE the first
   readAtmos call in the hot loop.

   No-op when input.use_node_atmos_cache is FALSE or if this node
   happens to own zero reduced rows.  Safe to call multiple times;
   subsequent calls are no-ops.                                     */
void init_atmos_node_cache(Atmosphere *atmos, Input_Atmos_file *infile) {
  if (!input.use_node_atmos_cache) return;
  if (mpi.n_nodes <= 0 || mpi.ny <= 0) return;
  if (mpi.node_ix1 <= mpi.node_ix0)    return;
  if (node_cache.initialized)          return;
  node_cache_init(atmos, infile);
}
/* ------- end ---------------------------- init_atmos_node_cache --- */

/* ------- begin -------------------------- readAtmos_hdf5  ---------
   Reads T, ne, vz, vturb, nH (and B if Stokes) for a given (xi,yi)
   in a SINGLE H5Dread_multi call.

   Strategy: read every variable at the FULL nz extent in one
   collective operation, then if temperature-cutoff is requested,
   determine zcut from the in-memory T array and compact each buffer
   in place via memmove + realloc.  Trades a small amount of extra
   bytes (the cut-off region) for a single sync point per column.   */
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
  hid_t    mem_1d, mem_nh;
  bool_t   old_moving;
  double  *Bx = NULL, *By = NULL, *Bz = NULL;
  double  *nh_scratch;
  long     full_nz, Nspace;
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

  full_nz       = (long) infile->nz;
  atmos->Nspace = full_nz;
  geometry->Ndep = (int) full_nz;

  /* --- Realloc all per-column 1D buffers to the full nz extent so the
     multi-read can fill them in one shot. --- */
  atmos->T         = (double *) realloc(atmos->T,         full_nz * sizeof(double));
  atmos->ne        = (double *) realloc(atmos->ne,        full_nz * sizeof(double));
  geometry->vel    = (double *) realloc(geometry->vel,    full_nz * sizeof(double));
  atmos->vturb     = (double *) realloc(atmos->vturb,     full_nz * sizeof(double));
  geometry->height = (double *) realloc(geometry->height, full_nz * sizeof(double));
  atmos->nHtot     = (double *) realloc(atmos->nHtot,     full_nz * sizeof(double));

  /* B-field scratch buffers (only if Stokes) */
  if (atmos->Stokes) {
    Bx = (double *) malloc(full_nz * sizeof(double));
    By = (double *) malloc(full_nz * sizeof(double));
    Bz = (double *) malloc(full_nz * sizeof(double));
  }

  /* nH scratch — flat NHydr*nz buffer; final atmos->nH built post-zcut. */
  nh_scratch = (double *) malloc((long)atmos->NHydr * full_nz * sizeof(double));

  /* === Build the multi-read at full nz extent === */
  hid_t   dsets  [10];
  hid_t   types  [10];
  hid_t   fspaces[10];
  hid_t   mspaces[10];
  void   *bufs   [10];
  size_t  ndset = 0;

  /* 1D memory dataspace at full nz, shared by all 4D variables */
  dims_mem[0] = full_nz;
  if ((mem_1d = H5Screate_simple(1, dims_mem, NULL)) < 0) HERR(routineName);

  /* 4D hyperslab template for column (xi, yi) at full nz */
  start4[0] = input.p15d_nt; count4[0] = 1;
  start4[1] = (size_t) xi;   count4[1] = 1;
  start4[2] = (size_t) yi;   count4[2] = 1;
  start4[3] = 0;             count4[3] = full_nz;

  /* T */
  if (H5Sselect_hyperslab(rfa_T_fspace, H5S_SELECT_SET, start4,
                          NULL, count4, NULL) < 0) HERR(routineName);
  dsets  [ndset] = infile->T_varid;
  types  [ndset] = H5T_NATIVE_DOUBLE;
  fspaces[ndset] = rfa_T_fspace;
  mspaces[ndset] = mem_1d;
  bufs   [ndset] = atmos->T;
  ndset++;

  /* z (2D or 4D layout) */
  if (mpi.ndims_z == 2) {
    start_z[0] = input.p15d_nt; count_z[0] = 1;
    start_z[1] = 0;             count_z[1] = full_nz;
  } else if (mpi.ndims_z == 4) {
    start_z[0] = input.p15d_nt; count_z[0] = 1;
    start_z[1] = (size_t) xi;   count_z[1] = 1;
    start_z[2] = (size_t) yi;   count_z[2] = 1;
    start_z[3] = 0;             count_z[3] = full_nz;
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

  /* vturb (optional) */
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

  /* Magnetic field */
  if (atmos->Stokes) {
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

  /* nH (5D, separate memspace [NHydr, nz]) */
  start_nh[0] = input.p15d_nt; count_nh[0] = 1;
  start_nh[1] = 0;             count_nh[1] = atmos->NHydr;
  start_nh[2] = (size_t) xi;   count_nh[2] = 1;
  start_nh[3] = (size_t) yi;   count_nh[3] = 1;
  start_nh[4] = 0;             count_nh[4] = full_nz;
  if (H5Sselect_hyperslab(rfa_nh_fspace, H5S_SELECT_SET, start_nh,
                          NULL, count_nh, NULL) < 0) HERR(routineName);
  dims_mem[0] = atmos->NHydr;
  dims_mem[1] = full_nz;
  if ((mem_nh = H5Screate_simple(2, dims_mem, NULL)) < 0) HERR(routineName);
  dsets  [ndset] = infile->nh_varid;
  types  [ndset] = H5T_NATIVE_DOUBLE;
  fspaces[ndset] = rfa_nh_fspace;
  mspaces[ndset] = mem_nh;
  bufs   [ndset] = nh_scratch;
  ndset++;

  /* === The single H5Dread_multi call: ONE I/O op per column === */
  if (H5Dread_multi(ndset, dsets, types, mspaces, fspaces,
                    plist_id, bufs) < 0) HERR(routineName);

  if (H5Sclose(mem_1d) < 0) HERR(routineName);
  if (H5Sclose(mem_nh) < 0) HERR(routineName);
  if (plist_id != H5P_DEFAULT) {
    if (H5Pclose(plist_id) < 0) HERR(routineName);
  }

  /* === Determine zcut from T (in memory, no I/O) === */
  mpi.zcut = 0;
  if (input.p15d_zcut && input.p15d_tmax >= 0) {
    int cut = -1;
    for (i = 0; i < full_nz; i++) {
      if (atmos->T[i] <= input.p15d_tmax) { cut = i; break; }
    }
    if (cut < 0) {
      sprintf(messageStr,
              "\n-Could not find temperature cut point! Aborting.\n");
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    mpi.zcut = cut;
  }
  Nspace = full_nz - mpi.zcut;

  /* === Compact 1D buffers in place if zcut > 0 === */
  if (mpi.zcut > 0) {
    memmove(atmos->T,         atmos->T         + mpi.zcut,
            Nspace * sizeof(double));
    memmove(atmos->ne,        atmos->ne        + mpi.zcut,
            Nspace * sizeof(double));
    memmove(geometry->vel,    geometry->vel    + mpi.zcut,
            Nspace * sizeof(double));
    memmove(atmos->vturb,     atmos->vturb     + mpi.zcut,
            Nspace * sizeof(double));
    memmove(geometry->height, geometry->height + mpi.zcut,
            Nspace * sizeof(double));
    if (atmos->Stokes) {
      memmove(Bx, Bx + mpi.zcut, Nspace * sizeof(double));
      memmove(By, By + mpi.zcut, Nspace * sizeof(double));
      memmove(Bz, Bz + mpi.zcut, Nspace * sizeof(double));
    }
  }

  /* Shrink 1D buffers down to Nspace */
  atmos->T         = (double *) realloc(atmos->T,         Nspace * sizeof(double));
  atmos->ne        = (double *) realloc(atmos->ne,        Nspace * sizeof(double));
  geometry->vel    = (double *) realloc(geometry->vel,    Nspace * sizeof(double));
  atmos->vturb     = (double *) realloc(atmos->vturb,     Nspace * sizeof(double));
  geometry->height = (double *) realloc(geometry->height, Nspace * sizeof(double));
  atmos->nHtot     = (double *) realloc(atmos->nHtot,     Nspace * sizeof(double));
  if (atmos->Stokes) {
    atmos->B       = (double *) realloc(atmos->B,       Nspace * sizeof(double));
    atmos->gamma_B = (double *) realloc(atmos->gamma_B, Nspace * sizeof(double));
    atmos->chi_B   = (double *) realloc(atmos->chi_B,   Nspace * sizeof(double));
  }

  atmos->Nspace  = Nspace;
  geometry->Ndep = (int) Nspace;

  /* === Build atmos->nH at the post-zcut size from scratch buffer === */
  /* (existing behaviour: previous nH is leaked, not fixing here) */
  atmos->nH = matrix_double(atmos->NHydr, Nspace);
  for (i = 0; i < atmos->NHydr; i++) {
    memcpy(atmos->nH[i],
           nh_scratch + (long)i * full_nz + mpi.zcut,
           Nspace * sizeof(double));
  }
  free(nh_scratch);

  /* === Post-processing === */
  if (atmos->Stokes) {
    for (j = 0; j < Nspace; j++) {
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

  for (j = 0; j < Nspace; j++) atmos->nHtot[j] = 0.0;

  /* Depth grid refinement */
  if (input.p15d_refine)
    depth_refine(atmos, geometry, input.p15d_tmax);

  /* Fix vturb: remove zeros, apply multiplier and offset */
  for (i = 0; i < Nspace; i++) {
    if (atmos->vturb[i] < 0.0) atmos->vturb[i] = 0.0;
    atmos->vturb[i] = atmos->vturb[i] * input.vturb_mult + input.vturb_add;
  }

  /* Sum to get nHtot */
  for (i = 0; i < atmos->NHydr; i++)
    for (j = 0; j < Nspace; j++)
      atmos->nHtot[j] += atmos->nH[i][j];

  /* Moving-atmosphere check */
  old_moving = atmos->moving;
  atmos->moving = FALSE;
  for (i = 0; i < Nspace; i++) {
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
  /* Release the node-shared atmosphere cache (no-op if not populated) */
  node_cache_free();
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
