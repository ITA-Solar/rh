/* ------- file: -------------------------- pool_buffer.c ----------------

       Two-phase collective I/O buffer for pool mode.

       During the pool compute phase, each drone stores its completed
       columns in a PoolOutputBuf.  After the pool loop, all MPI ranks
       (including the idle overlord) call writeCollective_pool() which
       performs a single collective HDF5 write per dataset.

       --------------------------                      ----------RH-- */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"
#include "io.h"

/* --- Function prototypes -- */
void loadBackground(int la, int mu, bool_t to_obs);

/* --- Global variables -- */
extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;
extern Input_Atmos_file infile;
extern MPI_data  mpi;
extern IO_data   io;
extern char messageStr[];

/* ======================================================================
   Buffer management
   ====================================================================== */

/* ------- begin --------------------------   poolbuf_init  ------------- */
void poolbuf_init(PoolOutputBuf *buf, int initial_capacity) {
  buf->ncols    = 0;
  buf->capacity = initial_capacity;
  buf->cols     = (PoolColumnBuf *) calloc(initial_capacity,
                                           sizeof(PoolColumnBuf));
}
/* ------- end   --------------------------   poolbuf_init  ------------- */

/* ------- begin --------------------------   poolbuf_grow  ------------- */
static void poolbuf_grow(PoolOutputBuf *buf) {
  buf->capacity *= 2;
  buf->cols = (PoolColumnBuf *) realloc(buf->cols,
                  buf->capacity * sizeof(PoolColumnBuf));
  /* Zero-initialise the new slots */
  memset(&buf->cols[buf->ncols], 0,
         (buf->capacity - buf->ncols) * sizeof(PoolColumnBuf));
}
/* ------- end   --------------------------   poolbuf_grow  ------------- */

/* ------- begin --------------------------   poolbuf_store_mpi  -------- */
void poolbuf_store_mpi(PoolOutputBuf *buf, long task_id) {
/* First call for every column.  Allocates a new slot and records
   the grid coordinates, convergence metadata and MPI bookkeeping.
   Must be called while mpi.ix/iy/task/niter/convergence are current. */

  if (buf->ncols >= buf->capacity) poolbuf_grow(buf);

  PoolColumnBuf *c = &buf->cols[buf->ncols];
  memset(c, 0, sizeof(PoolColumnBuf));

  c->ix          = mpi.ix;
  c->iy          = mpi.iy;
  c->task_id     = task_id;
  c->zcut        = mpi.zcut;
  c->Nspace      = atmos.Nspace;
  c->niter       = mpi.niter[0];
  c->convergence = mpi.convergence[0];
  c->dpopsmax    = mpi.dpopsmax[0];
  c->rank        = mpi.rank;

  /* Copy dpopsmax history, padded to NmaxIter */
  c->dpopsmax_hist = (double *) calloc(input.NmaxIter, sizeof(double));
  memcpy(c->dpopsmax_hist, mpi.dpopsmax_hist[0],
         c->niter * sizeof(double));

  /* ncols is incremented here so that store_aux/store_ray can find
     this slot via buf->cols[buf->ncols - 1] */
  buf->ncols++;
}
/* ------- end   --------------------------   poolbuf_store_mpi  -------- */

/* ------- begin --------------------------   poolbuf_store_aux_atmos --- */
void poolbuf_store_aux_atmos(PoolOutputBuf *buf) {
/* Second call for converged columns.  Captures atmosphere arrays and
   atomic/molecular populations and rates.  Must be called BEFORE the
   geometry is redefined for the output ray. */

  PoolColumnBuf *c = &buf->cols[buf->ncols - 1];
  int nact, kr, i, j, k, idx;
  Atom *atom;
  Molecule *molecule;
  AtomicLine      *line;
  AtomicContinuum *continuum;

  /* --- Atmosphere data --- */
  c->atmos_T  = (double *) malloc(c->Nspace * sizeof(double));
  c->atmos_vz = (double *) malloc(c->Nspace * sizeof(double));
  c->atmos_z  = (double *) malloc(c->Nspace * sizeof(double));
  c->atmos_ne = (double *) malloc(c->Nspace * sizeof(double));
  memcpy(c->atmos_T,  atmos.T,         c->Nspace * sizeof(double));
  memcpy(c->atmos_vz, geometry.vel,    c->Nspace * sizeof(double));
  memcpy(c->atmos_z,  geometry.height, c->Nspace * sizeof(double));
  memcpy(c->atmos_ne, atmos.ne,        c->Nspace * sizeof(double));

  /* --- Aux: atoms --- */
  if (atmos.Nactiveatom > 0) {
    c->atom_n     = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
    c->atom_nstar = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
    c->atom_RijL  = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
    c->atom_RjiL  = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
    c->atom_RijC  = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
    c->atom_RjiC  = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
    c->atom_Cij   = (double **) malloc(atmos.Nactiveatom * sizeof(double *));

    for (nact = 0; nact < atmos.Nactiveatom; nact++) {
      atom = atmos.activeatoms[nact];

      if (input.p15d_wpop) {
        long psize = atom->Nlevel * (long)c->Nspace * sizeof(double);
        c->atom_n[nact]     = (double *) malloc(psize);
        c->atom_nstar[nact] = (double *) malloc(psize);
        memcpy(c->atom_n[nact],     atom->n[0],     psize);
        memcpy(c->atom_nstar[nact], atom->nstar[0], psize);
      } else {
        c->atom_n[nact] = c->atom_nstar[nact] = NULL;
      }

      if (input.p15d_wrates) {
        long Lsize = atom->Nline * (long)c->Nspace * sizeof(double);
        long Csize = atom->Ncont * (long)c->Nspace * sizeof(double);

        c->atom_RijL[nact] = (double *) malloc(Lsize);
        c->atom_RjiL[nact] = (double *) malloc(Lsize);

        for (kr = 0; kr < atom->Nline; kr++) {
          line = &atom->line[kr];
          memcpy(c->atom_RijL[nact] + kr * c->Nspace, line->Rij,
                 c->Nspace * sizeof(double));
          memcpy(c->atom_RjiL[nact] + kr * c->Nspace, line->Rji,
                 c->Nspace * sizeof(double));
        }

        c->atom_RijC[nact] = (double *) malloc(Csize);
        c->atom_RjiC[nact] = (double *) malloc(Csize);

        for (kr = 0; kr < atom->Ncont; kr++) {
          continuum = &atom->continuum[kr];
          memcpy(c->atom_RijC[nact] + kr * c->Nspace, continuum->Rij,
                 c->Nspace * sizeof(double));
          memcpy(c->atom_RjiC[nact] + kr * c->Nspace, continuum->Rji,
                 c->Nspace * sizeof(double));
        }
      } else {
        c->atom_RijL[nact] = c->atom_RjiL[nact] = NULL;
        c->atom_RijC[nact] = c->atom_RjiC[nact] = NULL;
      }

      if (input.p15d_wcrates) {
        /* Capture full collision-rate matrix, layout matches writeAux_p
           for direct hyperslab write: [i*Nlevel+j][k] = C[j*Nlevel+i][k]. */
        long Ksize = (long)atom->Nlevel * atom->Nlevel * c->Nspace;
        c->atom_Cij[nact] = (double *) calloc(Ksize, sizeof(double));
        for (i = 0; i < atom->Nlevel; i++) {
          for (j = 0; j < atom->Nlevel; j++) {
            if (i == j) continue;
            for (k = 0; k < c->Nspace; k++) {
              idx = (i * atom->Nlevel + j) * c->Nspace + k;
              c->atom_Cij[nact][idx] = atom->C[j * atom->Nlevel + i][k];
            }
          }
        }
      } else {
        c->atom_Cij[nact] = NULL;
      }
    }
  }

  /* --- Aux: molecules --- */
  if (atmos.Nactivemol > 0) {
    c->mol_nv     = (double **) malloc(atmos.Nactivemol * sizeof(double *));
    c->mol_nvstar = (double **) malloc(atmos.Nactivemol * sizeof(double *));

    for (nact = 0; nact < atmos.Nactivemol; nact++) {
      molecule = atmos.activemols[nact];
      if (input.p15d_wpop) {
        long msize = molecule->Nv * (long)c->Nspace * sizeof(double);
        c->mol_nv[nact]     = (double *) malloc(msize);
        c->mol_nvstar[nact] = (double *) malloc(msize);
        memcpy(c->mol_nv[nact],     molecule->nv[0],     msize);
        memcpy(c->mol_nvstar[nact], molecule->nvstar[0], msize);
      } else {
        c->mol_nv[nact] = c->mol_nvstar[nact] = NULL;
      }
    }
  }
}
/* ------- end   --------------------------   poolbuf_store_aux_atmos --- */

/* ------- begin --------------------------   poolbuf_store_ray  -------- */
void poolbuf_store_ray(PoolOutputBuf *buf) {
/* Third call for converged columns.  Captures ray output (intensity,
   Stokes, and optionally chi/S/Jnu/sca/tau_one).  Must be called
   AFTER calculate_ray() so that spectrum.I is populated. */

  PoolColumnBuf *c = &buf->cols[buf->ncols - 1];
  int            nspect, k, l, idx;
  double        *J;
  bool_t         write_xtra, crosscoupling, to_obs, initialize;
  bool_t         prdh_limit_mem_save;
  ActiveSet     *as;

  write_xtra = (io.ray_nwave_sel > 0);

  /* --- Intensity --- */
  c->intensity = (double *) malloc(spectrum.Nspect * sizeof(double));
  memcpy(c->intensity, spectrum.I[0], spectrum.Nspect * sizeof(double));

  /* --- Stokes --- */
  if (atmos.Stokes || input.backgr_pol) {
    c->stokes_Q = (double *) malloc(spectrum.Nspect * sizeof(double));
    c->stokes_U = (double *) malloc(spectrum.Nspect * sizeof(double));
    c->stokes_V = (double *) malloc(spectrum.Nspect * sizeof(double));
    memcpy(c->stokes_Q, spectrum.Stokes_Q[0],
           spectrum.Nspect * sizeof(double));
    memcpy(c->stokes_U, spectrum.Stokes_U[0],
           spectrum.Nspect * sizeof(double));
    memcpy(c->stokes_V, spectrum.Stokes_V[0],
           spectrum.Nspect * sizeof(double));
  }

  /* --- tau=1 height --- */
  if (input.p15d_wtau) {
    float tau_cur, tau_prev, tmp;
    float *chi_tmp;

    c->tau_one = (float *) calloc(spectrum.Nspect, sizeof(float));

    prdh_limit_mem_save = FALSE;
    if (input.prdh_limit_mem) prdh_limit_mem_save = TRUE;
    input.prdh_limit_mem = TRUE;

    for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
      as = &spectrum.as[nspect];
      alloc_as(nspect, crosscoupling=FALSE);
      Opacity(nspect, 0, to_obs=TRUE, initialize=TRUE);
      chi_tmp = (float *) calloc(infile.nz, sizeof(float));

      if (input.backgr_in_mem)
        loadBackground(nspect, 0, to_obs=TRUE);
      else
        readBackground(nspect, 0, to_obs=TRUE);

      tau_prev = 0.0;
      tau_cur  = 0.0;

      for (k = 0; k < atmos.Nspace; k++) {
        l = k + mpi.zcut;
        chi_tmp[l] = (float)(as->chi[k] + as->chi_c[k]);

        if (k > 0) {
          tau_prev = tau_cur;
          tmp = 0.5 * (geometry.height[k-1] - geometry.height[k]);
          tau_cur = tau_prev + tmp * (chi_tmp[l] + chi_tmp[l-1]);

          if (((tau_cur > 1.) && (k == 1)) ||
              ((tau_cur < 1.) && (k == atmos.Nspace - 1))) {
            c->tau_one[nspect] = geometry.height[k];
          } else if ((tau_cur > 1.) && (tau_prev < 1.)) {
            tmp = (1. - tau_prev) / (tau_cur - tau_prev);
            c->tau_one[nspect] = geometry.height[k-1] +
              tmp * (geometry.height[k] - geometry.height[k-1]);
          }
        }
      }
      free_as(nspect, crosscoupling=FALSE);
      free(chi_tmp);
    }

    if (input.PRD_angle_dep == PRD_ANGLE_APPROX && atmos.NPRDactive > 0)
      input.prdh_limit_mem = prdh_limit_mem_save;
  }

  /* --- Extra quantities (chi, S, Jnu, sca) --- */
  if (write_xtra) {
    long xsize = (long)infile.nz * io.ray_nwave_sel * sizeof(float);
    c->chi   = (float *) malloc(xsize);
    c->S_ray = (float *) malloc(xsize);
    c->Jnu   = (float *) malloc(xsize);
    c->sca   = (float *) malloc(xsize);

    /* Initialise with fill value */
    for (k = 0; k < infile.nz * io.ray_nwave_sel; k++) {
      c->chi[k]   = FILLVALUE;
      c->S_ray[k] = FILLVALUE;
      c->Jnu[k]   = FILLVALUE;
      c->sca[k]   = FILLVALUE;
    }

    if (input.limit_memory)
      J = (double *) malloc(atmos.Nspace * sizeof(double));

    prdh_limit_mem_save = FALSE;
    if (input.prdh_limit_mem) prdh_limit_mem_save = TRUE;
    input.prdh_limit_mem = TRUE;

    for (nspect = 0; nspect < io.ray_nwave_sel; nspect++) {
      idx = io.ray_wave_idx[nspect];
      as  = &spectrum.as[idx];

      alloc_as(idx, crosscoupling=FALSE);
      Opacity(idx, 0, to_obs=TRUE, initialize=TRUE);

      if (input.backgr_in_mem)
        loadBackground(idx, 0, to_obs=TRUE);
      else
        readBackground(idx, 0, to_obs=TRUE);

      if (!input.limit_memory) J = spectrum.J[idx];

      for (k = 0; k < atmos.Nspace; k++) {
        l = k + mpi.zcut;
        c->chi[l * io.ray_nwave_sel + nspect] =
          (float)(as->chi[k] + as->chi_c[k]);
        c->sca[l * io.ray_nwave_sel + nspect] =
          (float)(as->sca_c[k] * J[k]);
        c->S_ray[l * io.ray_nwave_sel + nspect] =
          (float)((as->eta[k] + as->eta_c[k] + as->sca_c[k] * J[k]) /
                  (as->chi[k] + as->chi_c[k]));
        c->Jnu[l * io.ray_nwave_sel + nspect] = (float) J[k];
      }
      free_as(idx, crosscoupling=FALSE);
    }

    if (input.PRD_angle_dep == PRD_ANGLE_APPROX && atmos.NPRDactive > 0)
      input.prdh_limit_mem = prdh_limit_mem_save;
    if (input.limit_memory) free(J);
  }

}
/* ------- end   --------------------------   poolbuf_store_ray  -------- */

/* ------- begin --------------------------   poolbuf_free  ------------- */
void poolbuf_free(PoolOutputBuf *buf) {
  int i, nact;

  for (i = 0; i < buf->ncols; i++) {
    PoolColumnBuf *c = &buf->cols[i];
    free(c->dpopsmax_hist);
    free(c->atmos_T);
    free(c->atmos_vz);
    free(c->atmos_z);
    free(c->atmos_ne);
    free(c->intensity);
    free(c->stokes_Q);
    free(c->stokes_U);
    free(c->stokes_V);
    free(c->tau_one);
    free(c->chi);
    free(c->S_ray);
    free(c->Jnu);
    free(c->sca);

    if (c->atom_n != NULL) {
      for (nact = 0; nact < atmos.Nactiveatom; nact++) {
        free(c->atom_n[nact]);
        free(c->atom_nstar[nact]);
        free(c->atom_RijL[nact]);
        free(c->atom_RjiL[nact]);
        free(c->atom_RijC[nact]);
        free(c->atom_RjiC[nact]);
        free(c->atom_Cij[nact]);
      }
      free(c->atom_n);
      free(c->atom_nstar);
      free(c->atom_RijL);
      free(c->atom_RjiL);
      free(c->atom_RijC);
      free(c->atom_RjiC);
      free(c->atom_Cij);
    }

    if (c->mol_nv != NULL) {
      for (nact = 0; nact < atmos.Nactivemol; nact++) {
        free(c->mol_nv[nact]);
        free(c->mol_nvstar[nact]);
      }
      free(c->mol_nv);
      free(c->mol_nvstar);
    }
  }

  free(buf->cols);
  buf->cols = NULL;
  buf->ncols = buf->capacity = 0;
}
/* ------- end   --------------------------   poolbuf_free  ------------- */

/* ------- begin --------------------------   poolbuf_reset ------------- */
void poolbuf_reset(PoolOutputBuf *buf) {
/* Frees all column data but keeps the PoolOutputBuf shell allocated
   so the buffer can be reused for the next batch of columns. */

  int i, nact;

  for (i = 0; i < buf->ncols; i++) {
    PoolColumnBuf *c = &buf->cols[i];
    free(c->dpopsmax_hist);
    free(c->atmos_T);
    free(c->atmos_vz);
    free(c->atmos_z);
    free(c->atmos_ne);
    free(c->intensity);
    free(c->stokes_Q);
    free(c->stokes_U);
    free(c->stokes_V);
    free(c->tau_one);
    free(c->chi);
    free(c->S_ray);
    free(c->Jnu);
    free(c->sca);

    if (c->atom_n != NULL) {
      for (nact = 0; nact < atmos.Nactiveatom; nact++) {
        free(c->atom_n[nact]);
        free(c->atom_nstar[nact]);
        free(c->atom_RijL[nact]);
        free(c->atom_RjiL[nact]);
        free(c->atom_RijC[nact]);
        free(c->atom_RjiC[nact]);
        free(c->atom_Cij[nact]);
      }
      free(c->atom_n);
      free(c->atom_nstar);
      free(c->atom_RijL);
      free(c->atom_RjiL);
      free(c->atom_RijC);
      free(c->atom_RjiC);
      free(c->atom_Cij);
    }

    if (c->mol_nv != NULL) {
      for (nact = 0; nact < atmos.Nactivemol; nact++) {
        free(c->mol_nv[nact]);
        free(c->mol_nvstar[nact]);
      }
      free(c->mol_nv);
      free(c->mol_nvstar);
    }
  }

  /* Zero-initialise all slots and reset counter, keep allocation */
  memset(buf->cols, 0, buf->capacity * sizeof(PoolColumnBuf));
  buf->ncols = 0;
}
/* ------- end   --------------------------   poolbuf_reset ------------- */


/* ======================================================================
   Collective write helpers
   ====================================================================== */

/* Compare function for sorting columns by (ix, iy) in row-major order,
   needed so that multi-block hyperslab selection matches file iteration
   order for the collective write. */
static int cmp_col_ixiy(const void *a, const void *b) {
  const PoolColumnBuf *ca = (const PoolColumnBuf *)a;
  const PoolColumnBuf *cb = (const PoolColumnBuf *)b;
  if (ca->ix != cb->ix) return ca->ix - cb->ix;
  return ca->iy - cb->iy;
}

/* Write a single dataset collectively.  Every rank must call this.
   - dset_id:    open HDF5 dataset handle
   - rank_ndim:  number of dimensions of the file dataset
   - col_offset: per-column base offset (ix, iy position) — remaining
                 offset components (level, z, etc.) added inside
   - membuf:     contiguous memory buffer holding all columns' data,
                 or NULL if this rank has nothing to write
   - mem_total:  total number of elements in membuf
   - plist_id:   collective transfer property list
*/
static void write_dataset_collective(
    hid_t dset_id, hid_t plist_id,
    PoolColumnBuf *cols, int ncols,
    int ndim,
    hsize_t *per_col_count,   /* count for one column's hyperslab */
    int ix_dim,               /* which dim holds ix */
    int iy_dim,               /* which dim holds iy */
    int z_dim,                /* which dim holds z (-1 if none) */
    hid_t mem_type,
    void *membuf,
    hsize_t mem_total)
{
  const char routineName[] = "write_dataset_collective";
  hid_t file_dspace, mem_dspace;
  hsize_t offset[4] = {0, 0, 0, 0};
  int first = 1;

  if (( file_dspace = H5Dget_space(dset_id) ) < 0) HERR(routineName);

  if (ncols > 0 && membuf != NULL) {
    /* Build multi-block file hyperslab selection */
    for (int i = 0; i < ncols; i++) {
      memset(offset, 0, sizeof(offset));
      offset[ix_dim] = cols[i].ix;
      offset[iy_dim] = cols[i].iy;
      if (z_dim >= 0) offset[z_dim] = cols[i].zcut;

      if (( H5Sselect_hyperslab(file_dspace,
              first ? H5S_SELECT_SET : H5S_SELECT_OR,
              offset, NULL, per_col_count, NULL) ) < 0)
        HERR(routineName);
      first = 0;
    }

    /* Memory dataspace: contiguous 1D */
    if (( mem_dspace = H5Screate_simple(1, &mem_total, NULL) ) < 0)
      HERR(routineName);
  } else {
    /* This rank has nothing to write — participate with empty selection */
    if (( H5Sselect_none(file_dspace) ) < 0) HERR(routineName);
    hsize_t zero = 0;
    if (( mem_dspace = H5Screate_simple(1, &zero, NULL) ) < 0)
      HERR(routineName);
  }

  if (( H5Dwrite(dset_id, mem_type, mem_dspace, file_dspace,
                  plist_id, membuf) ) < 0) HERR(routineName);

  if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
  if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
}


/* Variant for datasets where each column may contribute a different
   z-extent (variable Nspace with zcut offset).  Builds per-column
   hyperslabs with individual counts.

   ONLY valid when (ix, iy) are the two leading dimensions, i.e. shape
   [nx, ny, nz].  Then the file's row-major iteration order over the
   OR-combined selection is exactly the (ix, iy)-sorted column order,
   which is how the caller packs membuf.  For shapes with a dimension
   in front of ix — populations [Nlevel, nx, ny, nz], collision rates
   [Nlevel, Nlevel, nx, ny, nz] — that no longer holds; use
   write_dataset_collective_lead instead. */
static void write_dataset_collective_varz(
    hid_t dset_id, hid_t plist_id,
    PoolColumnBuf *cols, int ncols,
    int ndim,
    hsize_t *base_count,      /* count template — z_dim entry overridden */
    int ix_dim, int iy_dim, int z_dim,
    hid_t mem_type,
    void *membuf,
    hsize_t mem_total)
{
  const char routineName[] = "write_dataset_collective_varz";
  hid_t file_dspace, mem_dspace;
  hsize_t offset[4], count[4];
  int first = 1;

  if (( file_dspace = H5Dget_space(dset_id) ) < 0) HERR(routineName);

  if (ncols > 0 && membuf != NULL) {
    for (int i = 0; i < ncols; i++) {
      memcpy(count, base_count, ndim * sizeof(hsize_t));
      memset(offset, 0, sizeof(offset));
      offset[ix_dim] = cols[i].ix;
      offset[iy_dim] = cols[i].iy;
      if (z_dim >= 0) {
        offset[z_dim] = cols[i].zcut;
        count[z_dim]  = cols[i].Nspace;
      }

      if (( H5Sselect_hyperslab(file_dspace,
              first ? H5S_SELECT_SET : H5S_SELECT_OR,
              offset, NULL, count, NULL) ) < 0)
        HERR(routineName);
      first = 0;
    }

    if (( mem_dspace = H5Screate_simple(1, &mem_total, NULL) ) < 0)
      HERR(routineName);
  } else {
    if (( H5Sselect_none(file_dspace) ) < 0) HERR(routineName);
    hsize_t zero = 0;
    if (( mem_dspace = H5Screate_simple(1, &zero, NULL) ) < 0)
      HERR(routineName);
  }

  if (( H5Dwrite(dset_id, mem_type, mem_dspace, file_dspace,
                  plist_id, membuf) ) < 0) HERR(routineName);

  if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
  if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
}


/* Write a dataset that carries one or two leading dimensions in front
   of (nx, ny, nz) in a single collective H5Dwrite.  Shapes handled:

     [Nkr,    nx, ny, nz]        rates      (nlead = 1, lead_dims = Nkr)
     [Nlevel, nx, ny, nz]        populations(nlead = 1)
     [Nlevel, Nlevel, nx, ny, nz] collision rates (nlead = 2)

   For each (leading index, column) pair we select the hyperslab
   [.., ix, iy, zcut:zcut+Nspace].

   ORDERING — the reason this function exists.  HDF5 fills an
   OR-combined selection in FILE row-major order, not in the order the
   blocks were OR'd in.  With the leading dimension outermost that
   means: every column at leading index 0, then every column at index
   1, and so on.  The memory buffer must therefore be lead-major:

     for il in 0..Nlead-1, for col in 0..ncols-1, Nspace doubles

   The per-column source arrays are lead-major *within* a column
   (col_ptrs[c] holds Nlead * Nspace_c contiguous doubles), so packing
   here is a transpose.  Handing HDF5 a column-major buffer instead
   silently scrambles values across columns whenever ncols > 1 and
   Nlead > 1 — it does not error, the numbers just land in the wrong
   place.  Columns must already be sorted by (ix, iy).             */
static void write_dataset_collective_lead(
    hid_t dset_id, hid_t plist_id,
    PoolColumnBuf *cols, int ncols,
    const hsize_t *lead_dims, /* extents of the leading dims */
    int nlead,                /* 1 or 2 */
    double **col_ptrs)        /* col_ptrs[c] → Nlead*Nspace contiguous */
{
  const char routineName[] = "write_dataset_collective_lead";
  hid_t file_dspace, mem_dspace;
  int first = 1;
  long Nlead = 1;
  int  c, d, il;

  for (d = 0; d < nlead; d++) Nlead *= (long) lead_dims[d];

  if (( file_dspace = H5Dget_space(dset_id) ) < 0) HERR(routineName);

  if (ncols > 0 && Nlead > 0 && col_ptrs != NULL) {
    /* Compute total memory size */
    hsize_t mem_total = 0;
    for (c = 0; c < ncols; c++)
      mem_total += (hsize_t) Nlead * cols[c].Nspace;

    /* Build lead-major memory buffer (see ORDERING above) */
    double *membuf = (double *) malloc(mem_total * sizeof(double));
    long pos = 0;
    for (il = 0; il < Nlead; il++) {
      for (c = 0; c < ncols; c++) {
        memcpy(membuf + pos,
               col_ptrs[c] + (long) il * cols[c].Nspace,
               cols[c].Nspace * sizeof(double));
        pos += cols[c].Nspace;
      }
    }

    /* Build multi-block file hyperslab: one block per (lead index, col) */
    for (il = 0; il < Nlead; il++) {
      for (c = 0; c < ncols; c++) {
        hsize_t offset[5], count[5];
        int k = 0;
        if (nlead == 2) {
          offset[k] = (hsize_t)(il / (long) lead_dims[1]);  count[k++] = 1;
          offset[k] = (hsize_t)(il % (long) lead_dims[1]);  count[k++] = 1;
        } else {
          offset[k] = (hsize_t) il;                         count[k++] = 1;
        }
        offset[k] = cols[c].ix;    count[k++] = 1;
        offset[k] = cols[c].iy;    count[k++] = 1;
        offset[k] = cols[c].zcut;  count[k++] = cols[c].Nspace;
        if (( H5Sselect_hyperslab(file_dspace,
                first ? H5S_SELECT_SET : H5S_SELECT_OR,
                offset, NULL, count, NULL) ) < 0) HERR(routineName);
        first = 0;
      }
    }

    if (( mem_dspace = H5Screate_simple(1, &mem_total, NULL) ) < 0)
      HERR(routineName);
    if (( H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_dspace,
                    file_dspace, plist_id, membuf) ) < 0) HERR(routineName);
    if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
    free(membuf);
  } else {
    /* Empty participation for collective call */
    if (( H5Sselect_none(file_dspace) ) < 0) HERR(routineName);
    hsize_t zero = 0;
    if (( mem_dspace = H5Screate_simple(1, &zero, NULL) ) < 0)
      HERR(routineName);
    if (( H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_dspace,
                    file_dspace, plist_id, NULL) ) < 0) HERR(routineName);
    if (( H5Sclose(mem_dspace) ) < 0) HERR(routineName);
  }
  if (( H5Sclose(file_dspace) ) < 0) HERR(routineName);
}


/* ======================================================================
   Phase 2: Collective writer
   ====================================================================== */

/* ------- begin --------------------------   writeCollective_pool  ----- */
void writeCollective_pool(PoolOutputBuf *buf, bool_t flush) {
  const char routineName[] = "writeCollective_pool";
  int    c, nact;
  hid_t  plist_id;
  bool_t write_xtra;

  sprintf(messageStr, "Process %4d: --- START collective write (%d columns)\n",
          mpi.rank, buf->ncols);
  Error(MESSAGE, "main", messageStr);

  /* Sort buffered columns by (ix, iy) so multi-block hyperslab selection
     matches file row-major iteration order */
  if (buf->ncols > 1)
    qsort(buf->cols, buf->ncols, sizeof(PoolColumnBuf), cmp_col_ixiy);

  /* --- Create transfer property list (collective or independent per
         keyword 15D_POOL_COLLECTIVE_WRITE) --- */
  if (( plist_id = H5Pcreate(H5P_DATASET_XFER) ) < 0) HERR(routineName);
  {
    H5FD_mpio_xfer_t xfer_mode = input.pool_collective_write
                                 ? H5FD_MPIO_COLLECTIVE
                                 : H5FD_MPIO_INDEPENDENT;
    if (( H5Pset_dxpl_mpio(plist_id, xfer_mode) ) < 0)
      HERR(routineName);
  }

  write_xtra = (io.ray_nwave_sel > 0);

  /* --- Build index arrays for converged vs all columns --- */
  /* "conv" = converged columns (have ray + aux data)
     "all"  = all columns including crashed/non-converged (have MPI metadata) */
  int nconv = 0;
  for (c = 0; c < buf->ncols; c++)
    if (buf->cols[c].convergence == 1) nconv++;

  PoolColumnBuf *conv_cols = NULL;
  if (nconv > 0) {
    conv_cols = (PoolColumnBuf *) malloc(nconv * sizeof(PoolColumnBuf));
    int j = 0;
    for (c = 0; c < buf->ncols; c++)
      if (buf->cols[c].convergence == 1) conv_cols[j++] = buf->cols[c];
  }


  /* =================================================================
     RAY FILE
     ================================================================= */
  {
    /* --- Intensity [nx, ny, Nspect] --- */
    hsize_t count3[3] = {1, 1, spectrum.Nspect};
    hsize_t mem_total = (hsize_t)nconv * spectrum.Nspect;
    double *membuf = NULL;

    if (nconv > 0) {
      membuf = (double *) malloc(mem_total * sizeof(double));
      for (c = 0; c < nconv; c++)
        memcpy(membuf + (long)c * spectrum.Nspect,
               conv_cols[c].intensity,
               spectrum.Nspect * sizeof(double));
    }

    write_dataset_collective(io.ray_int_var, plist_id,
        conv_cols, nconv, 3, count3, 0, 1, -1,
        H5T_NATIVE_DOUBLE, membuf, mem_total);
    free(membuf);

    /* --- Stokes Q, U, V --- */
    if (atmos.Stokes || input.backgr_pol) {
      for (int s = 0; s < 3; s++) {
        hid_t dset = (s == 0) ? io.ray_stokes_q_var :
                     (s == 1) ? io.ray_stokes_u_var :
                                io.ray_stokes_v_var;
        membuf = NULL;
        if (nconv > 0) {
          membuf = (double *) malloc(mem_total * sizeof(double));
          for (c = 0; c < nconv; c++) {
            double *src = (s == 0) ? conv_cols[c].stokes_Q :
                          (s == 1) ? conv_cols[c].stokes_U :
                                     conv_cols[c].stokes_V;
            memcpy(membuf + (long)c * spectrum.Nspect, src,
                   spectrum.Nspect * sizeof(double));
          }
        }
        write_dataset_collective(dset, plist_id,
            conv_cols, nconv, 3, count3, 0, 1, -1,
            H5T_NATIVE_DOUBLE, membuf, mem_total);
        free(membuf);
      }
    }

    /* --- tau_one [nx, ny, Nspect] --- */
    if (input.p15d_wtau) {
      hsize_t fmem_total = (hsize_t)nconv * spectrum.Nspect;
      float *fmembuf = NULL;
      if (nconv > 0) {
        fmembuf = (float *) malloc(fmem_total * sizeof(float));
        for (c = 0; c < nconv; c++)
          memcpy(fmembuf + (long)c * spectrum.Nspect,
                 conv_cols[c].tau_one,
                 spectrum.Nspect * sizeof(float));
      }
      write_dataset_collective(io.ray_tau1_var, plist_id,
          conv_cols, nconv, 3, count3, 0, 1, -1,
          H5T_NATIVE_FLOAT, fmembuf, fmem_total);
      free(fmembuf);
    }

    /* --- Extra 4D datasets: chi, S, Jnu, sca [nx, ny, nz, nwave_sel] --- */
    if (write_xtra) {
      long xtra_per_col = (long)infile.nz * io.ray_nwave_sel;
      hsize_t count4[4] = {1, 1, infile.nz, io.ray_nwave_sel};
      hsize_t fmem_total = (hsize_t)nconv * xtra_per_col;

      hid_t xdsets[4] = {io.ray_chi_var, io.ray_S_var,
                          io.ray_j_var,   io.ray_sca_c_var};

      for (int d = 0; d < 4; d++) {
        float *fmembuf = NULL;
        if (nconv > 0) {
          fmembuf = (float *) malloc(fmem_total * sizeof(float));
          for (c = 0; c < nconv; c++) {
            float *src = (d == 0) ? conv_cols[c].chi :
                         (d == 1) ? conv_cols[c].S_ray :
                         (d == 2) ? conv_cols[c].Jnu :
                                    conv_cols[c].sca;
            memcpy(fmembuf + (long)c * xtra_per_col, src,
                   xtra_per_col * sizeof(float));
          }
        }
        write_dataset_collective(xdsets[d], plist_id,
            conv_cols, nconv, 4, count4, 0, 1, -1,
            H5T_NATIVE_FLOAT, fmembuf, fmem_total);
        free(fmembuf);
      }
    }
  }


  /* =================================================================
     INDATA FILE — atmosphere [nx, ny, nz]
     ================================================================= */
  {
    /* Atmosphere datasets have variable z-extent per column */
    hsize_t base_count[3] = {1, 1, 0};  /* z overridden per column */

    hid_t atmos_dsets[4] = {io.in_atmos_T, io.in_atmos_vz,
                            io.in_atmos_z, io.in_atmos_ne};

    for (int d = 0; d < 4; d++) {
      hsize_t mem_total = 0;
      for (c = 0; c < nconv; c++) mem_total += conv_cols[c].Nspace;

      double *membuf = NULL;
      if (nconv > 0 && mem_total > 0) {
        membuf = (double *) malloc(mem_total * sizeof(double));
        long pos = 0;
        for (c = 0; c < nconv; c++) {
          double *src = (d == 0) ? conv_cols[c].atmos_T :
                        (d == 1) ? conv_cols[c].atmos_vz :
                        (d == 2) ? conv_cols[c].atmos_z :
                                   conv_cols[c].atmos_ne;
          memcpy(membuf + pos, src, conv_cols[c].Nspace * sizeof(double));
          pos += conv_cols[c].Nspace;
        }
      }

      write_dataset_collective_varz(atmos_dsets[d], plist_id,
          conv_cols, nconv, 3, base_count, 0, 1, 2,
          H5T_NATIVE_DOUBLE, membuf, mem_total);
      free(membuf);
    }
  }


  /* =================================================================
     INDATA FILE — MPI metadata [nx, ny] and [nx, ny, NmaxIter]
     ================================================================= */
  {
    /* These are written for ALL columns (including crashed/non-converged) */
    int nall = buf->ncols;
    PoolColumnBuf *all_cols = buf->cols;
    hsize_t count2[2] = {1, 1};

    /* task_map (rank), task_number, iterations, convergence, z_cut, delta_max
       are all scalar per column at [ix, iy] */
    hid_t mpi_dsets[6] = {io.in_mpi_tm, io.in_mpi_tn, io.in_mpi_it,
                          io.in_mpi_conv, io.in_mpi_zc, io.in_mpi_dm};

    for (int d = 0; d < 6; d++) {
      hsize_t mem_total = nall;
      void *membuf = NULL;
      hid_t mem_type;

      if (nall > 0) {
        if (d == 5) {  /* delta_max is double */
          double *dbuf = (double *) malloc(nall * sizeof(double));
          for (c = 0; c < nall; c++) dbuf[c] = all_cols[c].dpopsmax;
          membuf = dbuf;
          mem_type = H5T_NATIVE_DOUBLE;
        } else {
          int *ibuf = (int *) malloc(nall * sizeof(int));
          for (c = 0; c < nall; c++) {
            switch (d) {
              case 0: ibuf[c] = all_cols[c].rank;         break;
              case 1: ibuf[c] = (int)all_cols[c].task_id; break;
              case 2: ibuf[c] = all_cols[c].niter;        break;
              case 3: ibuf[c] = all_cols[c].convergence;  break;
              case 4: ibuf[c] = all_cols[c].zcut;         break;
            }
          }
          membuf = ibuf;
          mem_type = H5T_NATIVE_INT;
        }
      } else {
        mem_type = (d == 5) ? H5T_NATIVE_DOUBLE : H5T_NATIVE_INT;
      }

      write_dataset_collective(mpi_dsets[d], plist_id,
          all_cols, nall, 2, count2, 0, 1, -1,
          mem_type, membuf, mem_total);
      free(membuf);
    }

    /* delta_max_history [nx, ny, NmaxIter] */
    {
      hsize_t count3[3] = {1, 1, input.NmaxIter};
      hsize_t mem_total = (hsize_t)nall * input.NmaxIter;
      double *membuf = NULL;

      if (nall > 0) {
        membuf = (double *) malloc(mem_total * sizeof(double));
        for (c = 0; c < nall; c++)
          memcpy(membuf + (long)c * input.NmaxIter,
                 all_cols[c].dpopsmax_hist,
                 input.NmaxIter * sizeof(double));
      }

      write_dataset_collective(io.in_mpi_dmh, plist_id,
          all_cols, nall, 3, count3, 0, 1, -1,
          H5T_NATIVE_DOUBLE, membuf, mem_total);
      free(membuf);
    }
  }


  /* =================================================================
     AUX FILE — populations and rates
     ================================================================= */
  {
    Atom *atom;
    Molecule *molecule;

    for (nact = 0; nact < atmos.Nactiveatom; nact++) {
      atom = atmos.activeatoms[nact];

      /* --- Populations [Nlevel, nx, ny, nz] --- */
      if (input.p15d_wpop) {
        hsize_t lead[1] = {(hsize_t) atom->Nlevel};
        double **col_ptrs = NULL;

        if (nconv > 0)
          col_ptrs = (double **) malloc(nconv * sizeof(double *));

        /* populations n */
        for (c = 0; c < nconv; c++) col_ptrs[c] = conv_cols[c].atom_n[nact];
        write_dataset_collective_lead(io.aux_atom_pop[nact], plist_id,
            conv_cols, nconv, lead, 1, col_ptrs);

        /* populations nstar */
        for (c = 0; c < nconv; c++) col_ptrs[c] = conv_cols[c].atom_nstar[nact];
        write_dataset_collective_lead(io.aux_atom_poplte[nact], plist_id,
            conv_cols, nconv, lead, 1, col_ptrs);

        free(col_ptrs);
      }

      /* --- Rates [Nline/Ncont, nx, ny, nz] --- */
      if (input.p15d_wrates) {
        /* Build per-column pointer arrays, then write all kr slabs for
           each rate dataset in a single collective H5Dwrite call. */
        hsize_t lead_line[1] = {(hsize_t) atom->Nline};
        hsize_t lead_cont[1] = {(hsize_t) atom->Ncont};
        double **col_ptrs = NULL;

        if (nconv > 0)
          col_ptrs = (double **) malloc(nconv * sizeof(double *));

        /* Line rates: 2 datasets × 1 H5Dwrite each (all kr combined) */
        hid_t line_dsets[2] = {io.aux_atom_RijL[nact],
                               io.aux_atom_RjiL[nact]};
        for (int r = 0; r < 2; r++) {
          for (c = 0; c < nconv; c++) {
            col_ptrs[c] = (r == 0) ? conv_cols[c].atom_RijL[nact]
                                   : conv_cols[c].atom_RjiL[nact];
          }
          write_dataset_collective_lead(line_dsets[r], plist_id,
              conv_cols, nconv, lead_line, 1, col_ptrs);
        }

        /* Continuum rates: 2 datasets × 1 H5Dwrite each */
        hid_t cont_dsets[2] = {io.aux_atom_RijC[nact],
                               io.aux_atom_RjiC[nact]};
        for (int r = 0; r < 2; r++) {
          for (c = 0; c < nconv; c++) {
            col_ptrs[c] = (r == 0) ? conv_cols[c].atom_RijC[nact]
                                   : conv_cols[c].atom_RjiC[nact];
          }
          write_dataset_collective_lead(cont_dsets[r], plist_id,
              conv_cols, nconv, lead_cont, 1, col_ptrs);
        }

        free(col_ptrs);
      }

      /* --- Collision rates [Nlevel, Nlevel, nx, ny, nz] --- */
      if (input.p15d_wcrates) {
        hsize_t lead[2] = {(hsize_t) atom->Nlevel, (hsize_t) atom->Nlevel};
        double **col_ptrs = NULL;

        if (nconv > 0) {
          col_ptrs = (double **) malloc(nconv * sizeof(double *));
          for (c = 0; c < nconv; c++) col_ptrs[c] = conv_cols[c].atom_Cij[nact];
        }
        write_dataset_collective_lead(io.aux_atom_Cij[nact], plist_id,
            conv_cols, nconv, lead, 2, col_ptrs);
        free(col_ptrs);
      }
    }

    /* --- Molecule populations [Nv, nx, ny, nz] --- */
    for (nact = 0; nact < atmos.Nactivemol; nact++) {
      molecule = atmos.activemols[nact];

      if (input.p15d_wpop) {
        hsize_t lead[1] = {(hsize_t) molecule->Nv};
        double **col_ptrs = NULL;

        if (nconv > 0)
          col_ptrs = (double **) malloc(nconv * sizeof(double *));

        /* nv */
        for (c = 0; c < nconv; c++) col_ptrs[c] = conv_cols[c].mol_nv[nact];
        write_dataset_collective_lead(io.aux_mol_pop[nact], plist_id,
            conv_cols, nconv, lead, 1, col_ptrs);

        /* nvstar */
        for (c = 0; c < nconv; c++) col_ptrs[c] = conv_cols[c].mol_nvstar[nact];
        write_dataset_collective_lead(io.aux_mol_poplte[nact], plist_id,
            conv_cols, nconv, lead, 1, col_ptrs);

        free(col_ptrs);
      }
    }
  }

  /* --- Cleanup --- */
  if (( H5Pclose(plist_id) ) < 0) HERR(routineName);
  free(conv_cols);

  /* Flush HDF5 files only on last batch (files are flushed at close).

     Order matters for restartability: the convergence flags live in the
     indata file and are what a rerun trusts to decide a column is done.
     Flush the payload (ray, aux) first so that a job killed mid-flush
     can only ever lose the flag of a column whose data is already on
     disk — that column is simply recomputed.  The reverse order can
     leave a column flagged converged whose populations were never
     written, and the rerun then skips it forever. */
  if (flush) {
    if (( H5Fflush(io.aux_ncid, H5F_SCOPE_LOCAL) ) < 0) HERR(routineName);
    if (( H5Fflush(io.ray_ncid, H5F_SCOPE_LOCAL) ) < 0) HERR(routineName);
    if (( H5Fflush(io.in_ncid,  H5F_SCOPE_LOCAL) ) < 0) HERR(routineName);
  }

  sprintf(messageStr, "Process %4d: *** END   collective write\n", mpi.rank);
  Error(MESSAGE, "main", messageStr);
}
/* ------- end   --------------------------   writeCollective_pool  ----- */
