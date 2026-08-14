/* ------- file: -------------------------- rh15d_iobench.c ----------

       I/O benchmark for RH 1.5D pool mode.

       Performs the same setup, atmosphere reads, and collective HDF5
       writes as rh15d_ray_pool, but skips ALL radiative transfer
       computations (Background_p, getProfiles, initSolution_p,
       Iterate_p, calculate_ray).  Output buffers are zero-filled at
       the correct shapes so that writeCollective_pool() exercises the
       full real I/O path.

       Use this binary to measure the read/write cost of pool mode in
       isolation from the compute cost.

       --------------------------                      ----------RH-- */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"
#include "io.h"

#define WORKTAG  1
#define DIETAG   2
#define BATCHTAG 3

#ifndef REV_ID
#define REV_ID "UNKNOWN"
#endif

/* --- Global variables -- */
enum Topology topology = ONE_D_PLANE;
Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
Input_Atmos_file infile;
CommandLine commandline;
char messageStr[MAX_MESSAGE_LENGTH];
BackgroundData bgdat;
MPI_data  mpi;
IO_data   io;
IO_buffer iobuf;
const float FILLVALUE = FILL;

/* --- Local helpers ------------------------------------------------ */

/* Fill a PoolColumnBuf with zero data of the right shapes for the
   current column.  Mirrors poolbuf_store_mpi + poolbuf_store_aux_atmos
   + poolbuf_store_ray, but does not depend on Iterate_p / calculate_ray
   having populated atom->n, spectrum.I, etc.  All scientific values
   are set to zero; only sizes and metadata matter for an I/O test. */
static void iobench_fill_fake(PoolOutputBuf *buf, long task_id) {
  int nact, kr;
  Atom *atom;
  Molecule *molecule;

  /* --- Grow buffer if needed (mirrors poolbuf_store_mpi) --- */
  if (buf->ncols >= buf->capacity) {
    buf->capacity *= 2;
    buf->cols = (PoolColumnBuf *) realloc(buf->cols,
                    buf->capacity * sizeof(PoolColumnBuf));
    memset(&buf->cols[buf->ncols], 0,
           (buf->capacity - buf->ncols) * sizeof(PoolColumnBuf));
  }

  PoolColumnBuf *c = &buf->cols[buf->ncols];
  memset(c, 0, sizeof(PoolColumnBuf));

  /* --- MPI metadata --- */
  c->ix          = mpi.ix;
  c->iy          = mpi.iy;
  c->task_id     = task_id;
  c->zcut        = mpi.zcut;
  c->Nspace      = atmos.Nspace;
  c->niter       = 0;
  c->convergence = 1;            /* always pretend converged */
  c->dpopsmax    = 0.0;
  c->rank        = mpi.rank;
  c->dpopsmax_hist = (double *) calloc(input.NmaxIter, sizeof(double));

  buf->ncols++;

  /* --- Atmosphere arrays: copy real values just read by readAtmos --- */
  c->atmos_T  = (double *) malloc(c->Nspace * sizeof(double));
  c->atmos_vz = (double *) malloc(c->Nspace * sizeof(double));
  c->atmos_z  = (double *) malloc(c->Nspace * sizeof(double));
  c->atmos_ne = (double *) malloc(c->Nspace * sizeof(double));
  memcpy(c->atmos_T,  atmos.T,         c->Nspace * sizeof(double));
  memcpy(c->atmos_vz, geometry.vel,    c->Nspace * sizeof(double));
  memcpy(c->atmos_z,  geometry.height, c->Nspace * sizeof(double));
  memcpy(c->atmos_ne, atmos.ne,        c->Nspace * sizeof(double));

  /* --- Aux: atoms (zero-filled, sized like the real path) --- */
  if (atmos.Nactiveatom > 0) {
    c->atom_n     = (double **) calloc(atmos.Nactiveatom, sizeof(double *));
    c->atom_nstar = (double **) calloc(atmos.Nactiveatom, sizeof(double *));
    c->atom_RijL  = (double **) calloc(atmos.Nactiveatom, sizeof(double *));
    c->atom_RjiL  = (double **) calloc(atmos.Nactiveatom, sizeof(double *));
    c->atom_RijC  = (double **) calloc(atmos.Nactiveatom, sizeof(double *));
    c->atom_RjiC  = (double **) calloc(atmos.Nactiveatom, sizeof(double *));
    c->atom_Cij   = (double **) calloc(atmos.Nactiveatom, sizeof(double *));

    for (nact = 0; nact < atmos.Nactiveatom; nact++) {
      atom = atmos.activeatoms[nact];

      if (input.p15d_wpop) {
        long psize = (long)atom->Nlevel * c->Nspace;
        c->atom_n[nact]     = (double *) calloc(psize, sizeof(double));
        c->atom_nstar[nact] = (double *) calloc(psize, sizeof(double));
      }
      if (input.p15d_wrates) {
        long Lsize = (long)atom->Nline * c->Nspace;
        long Csize = (long)atom->Ncont * c->Nspace;
        c->atom_RijL[nact] = (double *) calloc(Lsize, sizeof(double));
        c->atom_RjiL[nact] = (double *) calloc(Lsize, sizeof(double));
        c->atom_RijC[nact] = (double *) calloc(Csize, sizeof(double));
        c->atom_RjiC[nact] = (double *) calloc(Csize, sizeof(double));
      }
      if (input.p15d_wcrates) {
        long Ksize = (long)atom->Nlevel * atom->Nlevel * c->Nspace;
        c->atom_Cij[nact] = (double *) calloc(Ksize, sizeof(double));
      }
    }
  }

  /* --- Aux: molecules --- */
  if (atmos.Nactivemol > 0) {
    c->mol_nv     = (double **) calloc(atmos.Nactivemol, sizeof(double *));
    c->mol_nvstar = (double **) calloc(atmos.Nactivemol, sizeof(double *));
    for (nact = 0; nact < atmos.Nactivemol; nact++) {
      molecule = atmos.activemols[nact];
      if (input.p15d_wpop) {
        long msize = (long)molecule->Nv * c->Nspace;
        c->mol_nv[nact]     = (double *) calloc(msize, sizeof(double));
        c->mol_nvstar[nact] = (double *) calloc(msize, sizeof(double));
      }
    }
  }

  /* --- Ray data: intensity, Stokes, tau, extras --- */
  c->intensity = (double *) calloc(spectrum.Nspect, sizeof(double));
  if (atmos.Stokes || input.backgr_pol) {
    c->stokes_Q = (double *) calloc(spectrum.Nspect, sizeof(double));
    c->stokes_U = (double *) calloc(spectrum.Nspect, sizeof(double));
    c->stokes_V = (double *) calloc(spectrum.Nspect, sizeof(double));
  }
  if (input.p15d_wtau) {
    c->tau_one = (float *) calloc(spectrum.Nspect, sizeof(float));
  }
  if (io.ray_nwave_sel > 0) {
    long xsize = (long)infile.nz * io.ray_nwave_sel;
    c->chi   = (float *) calloc(xsize, sizeof(float));
    c->S_ray = (float *) calloc(xsize, sizeof(float));
    c->Jnu   = (float *) calloc(xsize, sizeof(float));
    c->sca   = (float *) calloc(xsize, sizeof(float));
  }
}


/* Compute total bytes that one column contributes to the READ phase.
   Mirrors the datasets fetched by readAtmos_hdf5 (always at full nz). */
static double iobench_read_bytes_per_column(void) {
  double bytes = 0.0;
  /* T, z, ne, vz, vturb */
  bytes += 5.0 * (double)infile.nz * sizeof(double);
  /* nH */
  bytes += (double)atmos.NHydr * (double)infile.nz * sizeof(double);
  /* B-field components if Stokes */
  if (atmos.Stokes)
    bytes += 3.0 * (double)infile.nz * sizeof(double);
  return bytes;
}


/* Compute total bytes that one filled column contributes to disk.
   Mirrors the dataset list in writeCollective_pool. */
static double iobench_bytes_per_column(void) {
  double bytes = 0.0;
  int nact;

  /* Ray file */
  bytes += spectrum.Nspect * sizeof(double);                    /* intensity */
  if (atmos.Stokes || input.backgr_pol)
    bytes += 3.0 * spectrum.Nspect * sizeof(double);            /* Q, U, V  */
  if (input.p15d_wtau)
    bytes += spectrum.Nspect * sizeof(float);                   /* tau_one  */
  if (io.ray_nwave_sel > 0)
    bytes += 4.0 * (double)infile.nz * io.ray_nwave_sel * sizeof(float);

  /* Indata atmosphere (T, vz, z, ne) — uses Nspace */
  bytes += 4.0 * (double)atmos.Nspace * sizeof(double);

  /* Indata MPI metadata: 5 ints + 1 double + dpopsmax_hist */
  bytes += 5.0 * sizeof(int) + sizeof(double);
  bytes += (double)input.NmaxIter * sizeof(double);

  /* Aux: atom populations + rates */
  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    Atom *atom = atmos.activeatoms[nact];
    if (input.p15d_wpop)
      bytes += 2.0 * (double)atom->Nlevel * atmos.Nspace * sizeof(double);
    if (input.p15d_wrates) {
      bytes += 4.0 * (double)atom->Nline * atmos.Nspace * sizeof(double);
      bytes += 4.0 * (double)atom->Ncont * atmos.Nspace * sizeof(double);
    }
  }
  for (nact = 0; nact < atmos.Nactivemol; nact++) {
    Molecule *mol = atmos.activemols[nact];
    if (input.p15d_wpop)
      bytes += 2.0 * (double)mol->Nv * atmos.Nspace * sizeof(double);
  }
  return bytes;
}


/* ------- begin -------------------------- main ------------------- */
int main(int argc, char *argv[])
{
  bool_t run_ray = FALSE, writej = FALSE;
  long   task, total_tasks_local = 0;
  double t_setup0, t_setup1;
  double t_read_total = 0.0, t_fill_total = 0.0, t_write_total = 0.0;
  double t_run0, t_run1;
  double bytes_written_local = 0.0, bytes_written_global = 0.0;
  double bytes_read_local    = 0.0, bytes_read_global    = 0.0;
  double read_bytes_per_col  = 0.0, write_bytes_per_col  = 0.0;

  /* --- Set up MPI --- */
  mpi.main_logfile = stderr;
  initParallel(&argc, &argv, run_ray=FALSE);
  memset(&spectrum, 0, sizeof(spectrum));

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();
  mpi.main_logfile     = commandline.logfile;
  commandline.logfile  = mpi.logfile;
  strcpy(mpi.rev_id, REV_ID);

  t_setup0 = MPI_Wtime();

  /* --- Read input data and initialize (same path as ray_pool) --- */
  readInput(NULL);
  if (input.p15d_rerun) readSavedKeywords();

  /* --- I/O benchmark: force-disable depth_refine.
     readAtmos() otherwise calls depth_refine() per column, which does
     spline interpolation + Hminus_bf opacity — real compute work that
     would pollute the I/O measurement. */
  if (input.p15d_refine) {
    if (mpi.rank == 0) {
      fprintf(mpi.main_logfile,
              "  [iobench] forcing 15D_DEPTH_REFINE = FALSE for pure "
              "I/O measurement\n");
    }
    input.p15d_refine = FALSE;
  }

  spectrum.updateJ = TRUE;
  getCPU(1, TIME_START, NULL);
  init_atmos(&atmos, &geometry, &infile);
  distribute_jobs();
  if (input.p15d_rerun) {
      readSavedInput();
  } else {
      readRayInput();
  }
  mpi.Ntasks = 1;
  atmos.moving = TRUE;

  /* --- I/O benchmark: select read mode --------------------------------
     `readAtmos_hdf5` issues collective H5Dread when mpi.coll_atmos_read
     is TRUE and independent reads when FALSE.  In the real rh15d_ray_pool
     the overlord is excluded so 262144 / 2047 ≈ unbalanced → independent
     reads.  Our static round-robin makes the iobench balanced by
     accident, which would measure a code path the real run does not
     take.  Default to independent (real pool semantics); allow override
     via env var for A/B testing.

     RH_IOBENCH_READ_MODE = "independent" (default) | "collective"
  */
  {
    const char *mode = getenv("RH_IOBENCH_READ_MODE");
    if (mode != NULL && strcmp(mode, "collective") == 0) {
      mpi.coll_atmos_read = TRUE;
      if (mpi.rank == 0) {
        fprintf(mpi.main_logfile,
                "  [iobench] read mode: COLLECTIVE "
                "(forced via RH_IOBENCH_READ_MODE)\n");
      }
    } else {
      mpi.coll_atmos_read = FALSE;
      if (mpi.rank == 0) {
        fprintf(mpi.main_logfile,
                "  [iobench] read mode: INDEPENDENT "
                "(matches real rh15d_ray_pool)\n");
      }
    }
  }

  /* --- Populate the node-shared atmosphere cache.  Must come AFTER
     distribute_jobs (which sets node_ix0/ix1) and BEFORE the first
     readAtmos call so the dimension-discovery read also goes through
     the cache. */
  init_atmos_node_cache(&atmos, &infile);

  /* First read just to get dimensions.  Use a node-local column so
     the cache lookup hits — column (node_ix0, 0) is always owned by
     this node. */
  readAtmos(mpi.xnum[mpi.node_ix0], mpi.ynum[0], &atmos, &geometry, &infile);
  if (atmos.Stokes) Bproject();
  readAtomicModels();
  readMolecularModels();
  SortLambda();
  checkValuesRayInput();
  initParallelIO(run_ray=FALSE, writej=FALSE);

  t_setup1 = MPI_Wtime();

  /* --- Benchmark loop ------------------------------------------------
     Every rank (including rank 0) iterates over its own slice of the
     task map.  No master/drone scheme: this is purely an I/O test, so
     a static decomposition keeps the benchmark deterministic.        */
  PoolOutputBuf poolbuf;
  poolbuf_init(&poolbuf, 16);

  if (mpi.rank == 0) {
    sprintf(messageStr,
            "I/O benchmark: %ld total tasks, %d MPI ranks\n",
            mpi.total_tasks, mpi.size);
    fprintf(mpi.main_logfile, "%s", messageStr);
    Error(MESSAGE, "main", messageStr);
  }

  t_run0 = MPI_Wtime();

  /* Node-local static decomposition: each rank processes a stride
     of its node's task slice (mpi.node_task_start ..
     mpi.node_task_start + mpi.node_task_count).  Rank R on the node
     starts at offset node_rank, stride node_size.  Every column
     processed lies within this node's owned row range, so the cache
     lookup always hits and no rank ever falls back to file I/O.  */
  long node_start = mpi.node_task_start;
  long node_end   = mpi.node_task_start + mpi.node_task_count;

  /* Progress reporting: rank 0 prints every ~5% of its share */
  long my_share    = (mpi.node_task_count + mpi.node_size - 1 - mpi.node_rank)
                     / mpi.node_size;
  long progress_iv = my_share / 20;
  if (progress_iv < 1) progress_iv = 1;
  long progress_done = 0;

  for (task = node_start + mpi.node_rank; task < node_end; task += mpi.node_size) {
    double t0, t1;

    mpi.task = task;
    mpi.ix   = mpi.taskmap[task][0];
    mpi.iy   = mpi.taskmap[task][1];
    mpi.task = 0;  /* index into Ntasks=1 arrays */

    /* Real atmosphere read */
    t0 = MPI_Wtime();
    readAtmos(mpi.xnum[mpi.ix], mpi.ynum[mpi.iy], &atmos, &geometry, &infile);
    t1 = MPI_Wtime();
    t_read_total += (t1 - t0);
    /* Capture per-column sizes once we have a real atmos to measure. */
    if (read_bytes_per_col == 0.0) {
      read_bytes_per_col  = iobench_read_bytes_per_column();
      write_bytes_per_col = iobench_bytes_per_column();
    }
    bytes_read_local += read_bytes_per_col;

    /* Fake-fill output buffers (no compute) */
    t0 = MPI_Wtime();
    iobench_fill_fake(&poolbuf, task);
    bytes_written_local += iobench_bytes_per_column();
    t1 = MPI_Wtime();
    t_fill_total += (t1 - t0);

    total_tasks_local++;
    progress_done++;
    if (mpi.rank == 0 && (progress_done % progress_iv == 0)) {
      double elapsed = MPI_Wtime() - t_run0;
      double frac    = (double)progress_done / (double)my_share;
      double eta     = (frac > 0.0) ? elapsed * (1.0 / frac - 1.0) : 0.0;
      fprintf(mpi.main_logfile,
              "  [iobench] rank 0: %ld / %ld cols read+filled "
              "(%5.1f%%, elapsed %.1fs, eta %.1fs)\n",
              progress_done, my_share, 100.0 * frac, elapsed, eta);
      fflush(mpi.main_logfile);
    }
  }

  if (mpi.rank == 0) {
    fprintf(mpi.main_logfile,
            "  [iobench] rank 0: read+fill loop done, entering "
            "collective write...\n");
    fflush(mpi.main_logfile);
  }

  /* Single collective write of everything we buffered.  All ranks
     participate; ranks with fewer columns just contribute smaller
     hyperslab selections (empty selection if zero).                    */
  {
    double t0 = MPI_Wtime();
    writeCollective_pool(&poolbuf, TRUE);
    double t1 = MPI_Wtime();
    t_write_total += (t1 - t0);
  }
  poolbuf_reset(&poolbuf);

  t_run1 = MPI_Wtime();

  poolbuf_free(&poolbuf);
  closeParallelIO(run_ray=FALSE, writej=FALSE);
  finish_jobs();

  /* --- Aggregate timings across ranks --- */
  double t_read_max, t_fill_max, t_write_max, t_setup_max, t_run_max;
  double t_setup_local = t_setup1 - t_setup0;
  double t_run_local   = t_run1   - t_run0;

  MPI_Reduce(&t_setup_local, &t_setup_max, 1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&t_read_total,  &t_read_max,  1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&t_fill_total,  &t_fill_max,  1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&t_write_total, &t_write_max, 1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&t_run_local,   &t_run_max,   1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&bytes_written_local, &bytes_written_global, 1, MPI_DOUBLE,
             MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&bytes_read_local, &bytes_read_global, 1, MPI_DOUBLE,
             MPI_SUM, 0, MPI_COMM_WORLD);

  if (mpi.rank == 0) {
    double GB_w = bytes_written_global / (1024.0 * 1024.0 * 1024.0);
    double GB_r = bytes_read_global    / (1024.0 * 1024.0 * 1024.0);
    double KiB_per_col_r = read_bytes_per_col  / 1024.0;
    double KiB_per_col_w = write_bytes_per_col / 1024.0;
    double write_throughput = (t_write_max > 0) ? (GB_w / t_write_max) : 0.0;
    double read_throughput  = (t_read_max  > 0) ? (GB_r / t_read_max ) : 0.0;
    double end_to_end_throughput =
        (t_run_max > 0) ? ((GB_r + GB_w) / t_run_max) : 0.0;

    fprintf(mpi.main_logfile,
      "\n"
      "================================================================\n"
      "  RH 1.5D I/O BENCHMARK SUMMARY\n"
      "================================================================\n"
      "  MPI ranks                     : %d\n"
      "  Total atmosphere columns      : %ld\n"
      "  Atmosphere grid (nx,ny,nz)    : %ld x %ld x %ld\n"
      "  Spectrum Nspect               : %d\n"
      "  Active atoms / molecules      : %d / %d\n"
      "  Write populations (p15d_wpop) : %s\n"
      "  Write rates      (p15d_wrates): %s\n"
      "  Write tau_one    (p15d_wtau)  : %s\n"
      "  Extra ray wavelengths         : %d\n"
      "  Read mode                     : %s\n"
      "----------------------------------------------------------------\n"
      "  Per-column read  size         : %10.2f KiB\n"
      "  Per-column write size         : %10.2f KiB\n"
      "----------------------------------------------------------------\n"
      "  Setup time      (max)         : %10.3f s\n"
      "  Atmos read      (max)         : %10.3f s\n"
      "  Buffer fill     (max)         : %10.3f s\n"
      "  Collective write(max)         : %10.3f s\n"
      "  Total run loop  (max)         : %10.3f s\n"
      "----------------------------------------------------------------\n"
      "  Total bytes read              : %12.3f GiB\n"
      "  Total bytes written           : %12.3f GiB\n"
      "  Throughput (read only)        : %10.3f GiB/s\n"
      "  Throughput (write only)       : %10.3f GiB/s\n"
      "  Throughput (end-to-end loop)  : %10.3f GiB/s\n"
      "================================================================\n",
      mpi.size,
      mpi.total_tasks,
      (long)infile.nx, (long)infile.ny, (long)infile.nz,
      spectrum.Nspect,
      atmos.Nactiveatom, atmos.Nactivemol,
      input.p15d_wpop   ? "yes" : "no",
      input.p15d_wrates ? "yes" : "no",
      input.p15d_wtau   ? "yes" : "no",
      io.ray_nwave_sel,
      mpi.coll_atmos_read ? "collective" : "independent",
      KiB_per_col_r, KiB_per_col_w,
      t_setup_max, t_read_max, t_fill_max, t_write_max, t_run_max,
      GB_r, GB_w, read_throughput, write_throughput, end_to_end_throughput);
  }

  sprintf(messageStr, "*** I/O benchmark finished. Rank %d processed %ld "
          "columns, wrote %.3f MiB.\n", mpi.rank, total_tasks_local,
          bytes_written_local / (1024.0 * 1024.0));
  Error(MESSAGE, "main", messageStr);

  printTotalCPU();
  MPI_Info_free(&mpi.info);
  MPI_Finalize();
  return 0;
}
/* ------- end ---------------------------- main ------------------- */
