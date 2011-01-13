/* ------- file: -------------------------- parallel.h --------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Tue Nov 24 10:59:58 2010 --

       --------------------------                      ----------RH-- */

#ifndef __PARALLEL_H__
#define __PARALLEL_H__

#include <mpi.h>

typedef struct {
  int      size, rank, namelen, Ntasks, task, my_start, nx, ny, ix, iy;
  int     *xnum, *ynum, **taskmap, niter, nconv, nnoconv, ncrash;
  double   dpopsmax;
  bool_t   single_log, convergence, stop;
  FILE    *logfile, *main_logfile;
  long     backgrrecno; /* record number for background file */
  char     name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm comm;
  MPI_Info info;
} MPI_data;

void init_Background(void);
void Background_p(bool_t analyzeoutput, bool_t equilibria_only);
void close_Background();

void writeBRS_p(void);
void readBRS_p(void);
void init_ncdf_BRS(void);
void close_ncdf_BRS(void);
void writeBRS_ncdf(void);

void init_ncdf_J(void);
void close_ncdf_J(void);
void writeJlambda_ncdf(int nspect, double *J);
void writeJ20_ncdf(int nspect, double *J);
void readJlambda_ncdf(int nspect, double *J);
void readJ20_ncdf(int nspect, double *J);

void initSolution_p(void);

void distribute_jobs(void);
void finish_jobs(void);

void init_ncdf_spec(void);
void close_ncdf_spec(void);
void writeSpectrum_p(void);

void init_ncdf_indata(void);
void close_ncdf_indata(void);
void writeAtmos_p(void);
void writeMPI_p(void);

void init_ncdf_aux(void);
void init_aux_new(void);
void init_aux_old(void);
void close_ncdf_aux(void);
void writeAux_p(void);
void writeOpacity_p(void);
//void readPopulations_p(Atom *atom);

void initParallel(int *argc, char **argv[]);
void initParallelIO(void);
void closeParallelIO(void);
void UpdateAtmosDep(void);
void RequestStop_p(void);
bool_t StopRequested_p(void);

void Iterate_p(int NmaxIter, double iterLimit);
double solveSpectrum_p(bool_t eval_operator, bool_t redistribute);

void SolveLinearEq_p(int N, double **A, double *b, bool_t improve);


#define ERR(e,r) {printf("(EEE) %s: NetCDF: %s\n", r, nc_strerror(e)); exit(EXIT_FAILURE);}
#define MPILOG_TEMPLATE "scratch/rh_p%d.log"

#endif /* !__PARALLEL_H__ */

/* ---------------------------------------- parallel.h -------------- */
