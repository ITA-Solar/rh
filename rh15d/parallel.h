/* ------- file: -------------------------- parallel.h --------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Tue Nov 24 10:59:58 2010 --

       --------------------------                      ----------RH-- */

#ifndef __PARALLEL_H__
#define __PARALLEL_H__

#include <mpi.h>
#include <hdf5.h>

typedef struct {
  char     name[MPI_MAX_PROCESSOR_NAME], rev_id[MAX_LINE_SIZE];
  bool_t   single_log, stop, isfirst,isbalanced;
  int      nx, ny;
  int      size, rank, namelen, ix, iy, *xnum, *ynum, *niter, zcut, ndims_z;
  int     *zcut_hist, **rh_converged, StokesMode_save, *convergence, snap_number;
  long     nconv, nnoconv, ncrash, my_start, backgrrecno;
  long   **taskmap, task, Ntasks, total_tasks;
  double  *dpopsmax, **dpopsmax_hist;
  FILE    *logfile, *main_logfile;
  MPI_Comm comm;
  MPI_Info info;
  /* --- Node-level communicator (ranks sharing physical memory) ----
     Populated by initParallel via MPI_Comm_split_type(MPI_COMM_TYPE_SHARED).
     Used by node-partitioned pool mode and the node-shared atmosphere
     cache.  All ranks on a node share the same node_comm; node_id is
     a unique [0, n_nodes) index per node, the same for every rank on
     that node. */
  MPI_Comm node_comm;
  int      node_rank;        /* this rank's position in node_comm        */
  int      node_size;        /* number of ranks on this node             */
  int      n_nodes;          /* total number of nodes in MPI_COMM_WORLD  */
  int      node_id;          /* unique [0, n_nodes) id for this node    */
  /* --- Node-owned slice of the global taskmap (set by distribute_jobs) -
     Each node owns a contiguous range of reduced rows ix in [ix0, ix1).
     node_task_start / node_task_count index into mpi.taskmap and select
     the columns assigned to this node.  Used by the node-local pool
     dispatcher.  Cyclic decomposition is a planned follow-up. */
  int      node_ix0, node_ix1;
  long     node_task_start, node_task_count;
} MPI_data;

void init_Background();
void Background_p(bool_t analyzeoutput, bool_t equilibria_only);
void close_Background();

void writeJlambda_single(int nspect, double *J);
void writeJ20_single(int nspect, double *J);
void readJlambda_single(int nspect, double *J);
void readJ20_single(int nspect, double *J);

void initSolution_p(void);

void distribute_jobs(void);
void finish_jobs(void);
void readConvergence(void);

void init_hdf5_indata(void);
void init_hdf5_indata_new(void);
void init_hdf5_indata_existing(void);
void close_hdf5_indata(void);
void writeAtmos_all(void);
void writeAtmos_p(void);
void writeMPI_all(void);
void writeMPI_p(int task);

void init_hdf5_aux(void);
void init_aux_new(void);
void init_aux_existing(void);
void close_hdf5_aux(void);
void writeAux_all(void);
void writeAux_p(void);
void writeOpacity_p(void);

void initParallel(int *argc, char **argv[], bool_t run_ray);
void initParallelIO(bool_t run_ray, bool_t writej);
void closeParallelIO(bool_t run_ray, bool_t writej);
hid_t create_hdf5_fapl(void);
hid_t create_hdf5_fapl_indep(void);
void UpdateAtmosDep(void);
void RequestStop_p(void);
bool_t StopRequested_p(void);
void ERR(int ierror, const char *rname);
void HERR(const char *rname);
void copyBufVars(bool_t writej);
void writeOutput(bool_t writej);

void Iterate_p(int NmaxIter, double iterLimit);
double solveSpectrum_p(bool_t eval_operator, bool_t redistribute);

void SolveLinearEq_p(int N, double **A, double *b, bool_t improve);


#define MPILOG_TEMPLATE     "scratch/rh_p%d.log"
#define RAY_MPILOG_TEMPLATE "scratch/solveray_p%d.log"
#define PRD_FILE_TEMPLATE   "scratch/PRD_%s_%d-%d_p%d.dat"
#define PRD_FILE_TEMPLATE1  "scratch/PRD_%.1s_%d-%d_p%d.dat"

#endif /* !__PARALLEL_H__ */

/* ---------------------------------------- parallel.h -------------- */
