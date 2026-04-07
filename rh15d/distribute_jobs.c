/* ------- file: -------------------------- distribute_jobs.c ------------

       Version:       rh2.0, 1.5-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Wed Dec 29 20:58:00 2010 --

       --------------------------                      ----------RH-- */

/* --- Divides all the 1D columns between the different MPI processes --- */

#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>


#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "inputs.h"
#include "parallel.h"
#include "error.h"
#include "io.h"


/* --- Function prototypes --                          -------------- */

int   *intrange(int start, int end, int step, int *N);
long  *get_tasks(long ntotal, int size);
long **get_taskmap(long remain_tasks, long *ntasks, long *my_start);
long **matrix_long(long Nrow, long Ncol);

/* --- Global variables --                             -------------- */

extern MPI_data mpi;
extern InputData input;
extern Input_Atmos_file infile;
extern Geometry geometry;
extern char messageStr[];


/* ------- begin --------------------------   distribute_jobs.c   --- */
void distribute_jobs(void)
/* Distributes the work between the available processes */
{
  const char routineName[] = "distribute_jobs";
  long *tasks, remain_tasks, i, j;
  long minNtasks, maxNtasks;
  bool_t abort = FALSE;

  mpi.backgrrecno = 0;

  /* Sanitise input */
  if ((input.p15d_x0 < 0)) input.p15d_x0 = 0;
  if ((input.p15d_y0 < 0)) input.p15d_y0 = 0;

  if ((input.p15d_x0 > input.p15d_x1)) {
      sprintf(messageStr, "\n X_START larger than X_END.");
      abort = TRUE;
  }
  if ((input.p15d_y0 > input.p15d_y1)) {
      sprintf(messageStr + strlen(messageStr),
              "\n Y_START larger than Y_END.");
      abort = TRUE;
  }

  if ((input.p15d_x0 >= infile.nx)) {
      sprintf(messageStr + strlen(messageStr),
              "\n X_START larger than NX points in atmosphere.");
      abort = TRUE;
  }
  if ((input.p15d_y0 >= infile.ny)) {
      sprintf(messageStr + strlen(messageStr),
              "\n Y_START larger than NY points in atmosphere.");
      abort = TRUE;
  }
  if (abort) Error(ERROR_LEVEL_2, routineName, messageStr);


  if ((input.p15d_x1 <= 0) || (input.p15d_x1 > infile.nx))
    input.p15d_x1 = infile.nx;
  if ((input.p15d_y1 <= 0) || (input.p15d_y1 > infile.ny))
    input.p15d_y1 = infile.ny;

  if ((input.p15d_xst < 1)) input.p15d_xst = 1;
  if ((input.p15d_yst < 1)) input.p15d_yst = 1;

  /* Calculate array maps of (mpi.ix/iy) > xi/yi */
  mpi.xnum = intrange(input.p15d_x0, input.p15d_x1, input.p15d_xst, &mpi.nx);
  mpi.ynum = intrange(input.p15d_y0, input.p15d_y1, input.p15d_yst, &mpi.ny);

  /* Populate x and y scales */
  geometry.xscale = (double *) malloc(mpi.nx * sizeof(double));
  geometry.yscale = (double *) malloc(mpi.ny * sizeof(double));
  for (i=0; i < mpi.nx; i++) geometry.xscale[i] = infile.x[mpi.xnum[i]];
  for (i=0; i < mpi.ny; i++) geometry.yscale[i] = infile.y[mpi.ynum[i]];

  /* Find how many tasks to perform */
  if (input.p15d_rerun) {
    /* If rerun, read convergence info, calculate only for non-converged columns */
    readConvergence();
    remain_tasks = 0;

    for (i = 0; i < mpi.nx; i++) {
      for (j = 0; j < mpi.ny; j++) {
       if (mpi.rh_converged[i][j] != 1) remain_tasks++;
      }
    }

  } else {
    /* If running first time, use all columns */
    remain_tasks = mpi.nx * mpi.ny;
    mpi.rh_converged = matrix_int(mpi.nx, mpi.ny);
  }

  mpi.total_tasks = remain_tasks;

  /* Calculate tasks and distribute */
  tasks       = get_tasks(remain_tasks, mpi.size);
  mpi.Ntasks  = tasks[mpi.rank];
  minNtasks   = mpi.Ntasks;
  maxNtasks   = mpi.Ntasks;
  mpi.taskmap = get_taskmap(remain_tasks, tasks, &mpi.my_start);

  /* Find max/min Ntask per processor */
  for (int i = 0; i < mpi.size; ++i) {
    if (minNtasks > tasks[i]) {
      minNtasks = tasks[i];
    }
    if (maxNtasks < tasks[i]) {
      maxNtasks = tasks[i];
    }
  }
  free(tasks);

  if (minNtasks == maxNtasks) {
    mpi.isbalanced = true;
  } else {
    mpi.isbalanced = false;
  }

  /* --- Compute this node's contiguous range of reduced rows ---------
     Each node owns rows ix in [node_ix0, node_ix1).  The split is as
     even as possible (the first `remainder` nodes get one extra row).
     This range is the basis for the node-shared atmosphere cache and
     the node-local pool dispatcher.

     For the rerun case the row range is the same; the actual taskmap
     section may be sparser because converged columns are filtered out,
     but row-major ordering of the taskmap means the entries belonging
     to a row range are still a contiguous slice in the taskmap. */
  {
    int base = mpi.nx / mpi.n_nodes;
    int rem  = mpi.nx % mpi.n_nodes;
    if (mpi.n_nodes > mpi.nx) {
      sprintf(messageStr,
              "\nNumber of nodes (%d) exceeds number of reduced rows (%d).\n"
              "  Each node must own at least one row.  Reduce node count\n"
              "  or increase grid size in X.\n",
              mpi.n_nodes, mpi.nx);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    mpi.node_ix0 = mpi.node_id * base + ((mpi.node_id < rem) ? mpi.node_id : rem);
    mpi.node_ix1 = mpi.node_ix0 + base + ((mpi.node_id < rem) ? 1 : 0);
  }

  /* Locate this node's slice of the global taskmap.  Because the
     taskmap is built row-major over (ix, iy), entries with ix in
     [node_ix0, node_ix1) are contiguous and we can find the slice
     with two linear scans (cheap, runs once at startup).            */
  mpi.node_task_start = 0;
  mpi.node_task_count = 0;
  if (mpi.taskmap != NULL && remain_tasks > 0) {
    long k;
    for (k = 0; k < remain_tasks; k++) {
      if (mpi.taskmap[k][0] >= mpi.node_ix0) {
        mpi.node_task_start = k;
        break;
      }
    }
    for (k = mpi.node_task_start; k < remain_tasks; k++) {
      if (mpi.taskmap[k][0] >= mpi.node_ix1) break;
      mpi.node_task_count++;
    }
  }

  return;
}
/* ------- end   --------------------------   distribute_jobs.c   --- */

/* ------- begin --------------------------   get_tasks.c ------- --- */
long *get_tasks(long ntotal, int size)
/* Divides the ntotal tasks by 'size' processes */
{
  long i, *tasks;

  tasks = (long *) calloc(size, sizeof(long));

  if (size > ntotal) {   /* More processes thank tasks, use ntotal processes */
    for (i = 0; i < ntotal; i++)
      tasks[i] = 1;
  } else {
    for (i=0; i < size; i++)
      tasks[i] = ntotal/size;
    /* Distribute remaining */
    if ((ntotal % size) > 0) {
      for (i=0; i < ntotal % size; i++) ++tasks[i];
    }
  }

  return tasks;
}
/* ------- end   --------------------------   get_tasks.c ------- --- */

/* ------- begin --------------------------   get_retaskmap.c --- --- */
long **get_taskmap(long remain_tasks, long *ntasks, long *my_start)
/* Distributes the tasks by the processes (defined by global taskmap and local
   my_start and mpi.Ntasks). Uses mpi.rh_converged to filter out already converged
   columns (if rerun is used).  */
{
  long i, j, k, *start, **taskmap;

  if (remain_tasks < 1)  /* For cases when there is no work */
    return NULL;
  taskmap = matrix_long(remain_tasks, (long) 2);
  start   = (long *) malloc(mpi.size * sizeof(long));

  /* Create map of tasks */
  k = 0;
  for (i=0; i < mpi.nx; i++) {
    for (j=0; j < mpi.ny; j++) {
      if (mpi.rh_converged[i][j] != 1) {
	       taskmap[k][0] = i;
	       taskmap[k][1] = j;
	       ++k;
      }
    }
  }

  /* distribute tasks */
  k = 0;
  for (i=0; i < mpi.size; i++) {
    start[i] = k;
    k += ntasks[i];
  }

  /* pointer to each process's starting task */
  *my_start = start[mpi.rank];

  free(start);

  return taskmap;
}
/* ------- end   --------------------------   get_retaskmap.c --- --- */

/* ------- begin --------------------------   intrange.c -------- --- */
int *intrange(int start, int end, int step, int *N)
/* Mimics Python's range function. Also gives a pointer
   to the number of elements. */
{
  int i, *arange;

  *N = (end - start) / step;
  if  ((end - start) % step > 0) ++*N;

  arange = (int *) malloc(*N * sizeof(int));

  arange[0] = start;
  for (i=1; i < *N; i++) {
    arange[i] = arange[i-1] + step;
  }

  return arange;
}
/* ------- end   --------------------------   intrange.c -------- --- */

/* ------- begin -------------------------- finish_jobs.c ------- --- */
void finish_jobs(void)
/* Frees from memory stuff used for job control */
{

  /* Get total number of tasks and convergence statuses */
  MPI_Allreduce(MPI_IN_PLACE, &mpi.Ntasks,  1, MPI_LONG, MPI_SUM, mpi.comm);
  MPI_Allreduce(MPI_IN_PLACE, &mpi.ncrash,  1, MPI_LONG, MPI_SUM, mpi.comm);
  MPI_Allreduce(MPI_IN_PLACE, &mpi.nconv,   1, MPI_LONG, MPI_SUM, mpi.comm);
  MPI_Allreduce(MPI_IN_PLACE, &mpi.nnoconv, 1, MPI_LONG, MPI_SUM, mpi.comm);

  free(mpi.xnum);
  free(mpi.ynum);
  if (mpi.taskmap != NULL)
    freeMatrix((void **) mpi.taskmap);
  if (mpi.rank_node != NULL) {
    free(mpi.rank_node);
    mpi.rank_node = NULL;
  }

}
/* ------- end   -------------------------- finish_jobs.c ------- --- */



/* ------- begin -------------------------- matrix_intl.c ----------- */

long **matrix_long(long Nrow, long Ncol)
{
  register long i;

  long *theMatrix, **Matrix;
  long typeSize = sizeof(long), pointerSize = sizeof(long *);

  theMatrix = (long *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (long **) malloc(Nrow * pointerSize);
  for (i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}
/* ------- end ---------------------------- matrix_intl.c ----------- */
