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
#include <sys/types.h>
#include <sys/stat.h>


#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "inputs.h"
#include "parallel.h"
#include "io.h"


/* --- Function prototypes --                          -------------- */

int  *intrange(int start, int end, int step, int *N);
int  *get_tasks(long ntotal, int size);
int **get_taskmap(int *ntasks, int *my_start);

/* --- Global variables --                             -------------- */

extern MPI_data mpi;
extern InputData input;
extern char messageStr[];
extern NCDF_Atmos_file infile;


/* ------- begin --------------------------   distribute_jobs.c   --- */
void distribute_jobs(void)
{
  /* The end product of this routine must be:

     For each process:
         * mpi.Ntasks --DONE
	 * An array of rank (mpi.Ntasks,2) that gives (mpi.ix, mpi.iy) 
           for each task index --DONE (mpi.my_tasks)
	 * A map of (mpi.ix, mpi.iy) -> (xi, yi) --DONE (two 1D arrays: mpi.xnum, mpi.ynum)

 */

  int  *tasks, my_start;

 
  mpi.backgrrecno = 0; 

  /* Sanitise input */
  if ((input.p15d_x0 < 0)) input.p15d_x0 = 0;
  if ((input.p15d_y0 < 0)) input.p15d_y0 = 0;

  if ((input.p15d_x0 > input.p15d_x1)) input.p15d_x0 = input.p15d_x1;
  if ((input.p15d_y0 > input.p15d_y1)) input.p15d_y0 = input.p15d_y1;

  if ((input.p15d_x1 <= 0) || (input.p15d_x1 > infile.nx))
    input.p15d_x1 = infile.nx;
  if ((input.p15d_y1 <= 0) || (input.p15d_y1 > infile.ny))
    input.p15d_y1 = infile.ny;

  if ((input.p15d_xst < 1)) input.p15d_xst = 1;
  if ((input.p15d_yst < 1)) input.p15d_yst = 1;

  /* Calculate array maps of (mpi.ix/iy) > xi/yi */
  mpi.xnum = intrange(input.p15d_x0, input.p15d_x1, input.p15d_xst, &mpi.nx);
  mpi.ynum = intrange(input.p15d_y0, input.p15d_y1, input.p15d_yst, &mpi.ny);

  /* Calculate tasks and distribute */ 
  tasks        = get_tasks(mpi.nx*mpi.ny, mpi.size);
  mpi.Ntasks   = tasks[mpi.rank];
  mpi.taskmap  = get_taskmap(tasks, &my_start);

  //printf("MMM = %d\n",mpi.my_tasks[10][1]);
 
  /*
  printf("Process %d: %d tasks, mystart = %d\n",mpi.rank,mpi.Ntasks, my_start);
  for (i=0; i < mpi.Ntasks; i++){
    printf("Process %d: task %d : (x,y) = (%d, %d)\n",mpi.rank,i, mpi.taskmap[i+my_start][0],
	   mpi.taskmap[i+my_start][1]);
    printf("Process %d: task %d : (x,y) = (%d, %d)\n",mpi.rank,i, mpi.my_tasks[i][0],
	   mpi.my_tasks[i][1]);
  }
  */

  free(tasks);

  return;
}
/* ------- end   --------------------------   distribute_jobs.c   --- */

/* ------- begin --------------------------   get_taskmap.c ----- --- */
int **get_taskmap(int *ntasks, int *my_start)
{
  int i, j, k, *start, **taskmap;

  taskmap = matrix_int(mpi.nx*mpi.ny, 2);
  start   = (int *) malloc(mpi.size * sizeof(int));

  /* Create map of tasks */
  k = 0;
  for (i=0; i < mpi.nx; i++) {
    for (j=0; j < mpi.ny; j++) {
      taskmap[k][0] = j;
      taskmap[k][1] = i;
      ++k;
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
/* ------- end   --------------------------   get_taskmap.c ----- --- */

/* ------- begin --------------------------   get_tasks.c ------- --- */
int *get_tasks(long ntotal, int size)
/* Divides the ntotal tasks by 'size' processes */
{
  int i, *tasks;

  tasks = (int *) malloc(size * sizeof(int));
  
  for (i=0; i < size; i++) tasks[i] = ntotal/size;

  /* Distribute remaining */
  if ((ntotal % size) > 0) {
    for (i=0; i < ntotal % size; i++) ++tasks[i];
  }
  
  return tasks;
}
/* ------- end   --------------------------   get_tasks.c ------- --- */

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

  free(mpi.xnum);
  free(mpi.ynum);
  free(mpi.taskmap);

}
/* ------- end   -------------------------- finish_jobs.c ------- --- */
