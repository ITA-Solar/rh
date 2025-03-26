/* ------- file: -------------------------- writeAux_p.c ---------

       Version:       rh2.0, 1.5-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       --------------------------                      -----------RH-- */

/* --- Set of tools to deal with IO initialisation and closure ------- */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "accelerate.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "spectrum.h"
#include "parallel.h"
#include "io.h"


/* --- Function prototypes --                          -------------- */
void distribute_nH(void);
void  allocBufVars(bool_t writej);
void   freeBufVars(bool_t writej);

/* --- Global variables --                             -------------- */
extern Atmosphere atmos;
extern Geometry geometry;
extern Input_Atmos_file infile;
extern IO_data io;
extern IO_buffer iobuf;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];
extern MPI_data mpi;
extern Spectrum spectrum;

/* ------- begin --------------------------   initParallel.c --   --- */
void initParallel(int *argc, char **argv[], bool_t run_ray) {
  const char routineName[] = "initParallel";
  char   logfile[MAX_LINE_SIZE];

  /* Initialise MPI */
  MPI_Init(argc,argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi.size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
  MPI_Get_processor_name(mpi.name, &mpi.namelen);
  mpi.comm = MPI_COMM_WORLD;
  mpi.info = MPI_INFO_NULL;
  /* Open log files */
  sprintf(logfile, (run_ray) ? RAY_MPILOG_TEMPLATE : MPILOG_TEMPLATE, mpi.rank);
  if ((mpi.logfile = fopen(logfile, "w")) == NULL) {
    sprintf(messageStr, "Process %4d: Unable to open log file %s",
	    mpi.rank, logfile);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* _IOFBF for full buffering, _IOLBF for line buffering */
  setvbuf(mpi.logfile, NULL, _IOLBF, BUFSIZ_MPILOG);
  mpi.stop    = FALSE;
  mpi.nconv   = 0;
  mpi.nnoconv = 0;
  mpi.ncrash  = 0;
  mpi.taskmap = NULL;

  /* Initialise input */
  input.atoms_file_contents = NULL;
  input.atomic_file_contents = NULL;
  input.kurucz_file_contents = NULL;
  input.kurucz_line_file_contents = NULL;
  input.Nkurucz_files = 0;
  input.wavetable = NULL;
  return;
}
/* ------- end   --------------------------   initParallel.c --   --- */

/* ------- begin --------------------------   initParallelIO.c    --- */
void initParallelIO(bool_t run_ray, bool_t writej) {
  int i;

  init_hdf5_aux();
  init_Background();
  if (!run_ray) {
    init_hdf5_indata();
  }
  init_hdf5_ray();
  /* Save StokesMode (before adjustStokesMode changes it...) */
  mpi.StokesMode_save = input.StokesMode;
  mpi.single_log      = FALSE;
  /* Allocate some mpi. arrays */
  mpi.niter       = (int *)    calloc(mpi.Ntasks , sizeof(int));
  mpi.convergence = (int *)    calloc(mpi.Ntasks , sizeof(int));
  mpi.zcut_hist   = (int *)    calloc(mpi.Ntasks , sizeof(int));
  mpi.dpopsmax    = (double *) calloc(mpi.Ntasks , sizeof(double));
  /* max with 1 is used to make sure array is allocated even with Ntasks = 0 */
  mpi.dpopsmax_hist = matrix_double(MAX(mpi.Ntasks, 1), input.NmaxIter);
  mpi.zcut_hist[mpi.task] = mpi.zcut;

  /* Fill mpi.niter with ones, to avoid problems with crashes on 1st iteration */
  for (i=0; i < mpi.Ntasks; i++) mpi.niter[i] = 1;
  /* buffer quantities for final writes */
  allocBufVars(writej);
}
/* ------- end   --------------------------   initParallelIO.c    --- */


/* ------- begin --------------------------  closeParallelIO.c    --- */
void closeParallelIO(bool_t run_ray, bool_t writej) {
  int i;

  if (!run_ray) {
    close_hdf5_indata();
  }
  close_atmos(&atmos, &geometry, &infile);
  close_hdf5_aux();
  close_hdf5_ray();

  //free(io.atom_file_pos);
  free(mpi.niter);
  free(mpi.dpopsmax);
  free(mpi.convergence);
  freeMatrix((void **) mpi.dpopsmax_hist);

  /* Free buffer variables */
  freeBufVars(writej);

  /* Free RLK stuff */
  for (i = 0;  i < atmos.Nrlk;  i++) {
    if (atmos.rlk_lines[i].zm != NULL) freeZeeman(atmos.rlk_lines[i].zm);
  }
  if (atmos.rlk_lines != NULL) free(atmos.rlk_lines); 
}
/* ------- end   --------------------------  closeParallelIO.c    --- */

/* ------- begin --------------------------  updateAtmosDep.c     --- */
void UpdateAtmosDep(void) {
/* Updates the atmos-dependent factors for the atoms and molecules */
  const char routineName[] = "UpdateAtomsDep";
  int       nact, k, kr, la, Nlamu;
  double    vtherm;
  Atom     *atom;
  Molecule *molecule;
  AtomicLine    *line;
  MolecularLine *mrt;

  /* Put back initial Stokes mode */
  input.StokesMode = mpi.StokesMode_save;

  mpi.zcut_hist[mpi.task] = mpi.zcut;

  /* For single 1D atmosphere, this is not needed */
  //if (geometry.atmos_format == MULTI) return;

  /* Recalculate magnetic field projections */
  if (atmos.Stokes) {
    if (atmos.cos_gamma != NULL) {
      freeMatrix((void **) atmos.cos_gamma);
      atmos.cos_gamma = NULL;
    }
    if (atmos.cos_2chi != NULL) {
      freeMatrix((void **) atmos.cos_2chi);
      atmos.cos_2chi = NULL;
    }
    if (atmos.sin_2chi != NULL) {
      freeMatrix((void **) atmos.sin_2chi);
      atmos.sin_2chi = NULL;
    }
    Bproject();
  }

  /* Update atmos-dependent atomic  quantities --- --------------- */
  for (nact = 0; nact < atmos.Natom; nact++) {
    atom = &atmos.atoms[nact];

    /* Reallocate some stuff (because of varying Nspace) */
    atom->ntotal = (double *) realloc(atom->ntotal, atmos.Nspace * sizeof(double));
    atom->vbroad = (double *) realloc(atom->vbroad, atmos.Nspace * sizeof(double));

    if (atom->nstar != NULL) {
      freeMatrix((void **) atom->nstar);
      atom->nstar = matrix_double(atom->Nlevel, atmos.Nspace);
    }

    /* Only allocate n again for active atoms or atoms with read populations
       Hydrogen is a special case because populations can be read differently.
       When H_LTE is on, n is just a pointer to nstar and needs no freeing.  */
    if ((atom->active) || (atom->popsinFile) ||
        ((nact == 0) && (!atmos.H_LTE))) {
         if (atom->n != NULL) freeMatrix((void **) atom->n);
         atom->n = matrix_double(atom->Nlevel, atmos.Nspace);
    } else {
     /* alias to nstar again, just in case */
         atom->n = atom->nstar;
    }


    for (k = 0;  k < atmos.Nspace;  k++)
      atom->ntotal[k] = atom->abundance * atmos.nHtot[k];

    if (atom->Nline > 0) {
      vtherm = 2.0*KBOLTZMANN/(AMU * atom->weight);
      for (k = 0;  k < atmos.Nspace;  k++)
	atom->vbroad[k] = sqrt(vtherm*atmos.T[k] + SQ(atmos.vturb[k]));
    }
  }

  /* Now only for active atoms */
  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];
    /* Reallocate some stuff (because of varying Nspace) */

    /* Free collision rate array, will be reallocated by calls in Background_p */
    if (atom->C != NULL) {
      freeMatrix((void **) atom->C);
      atom->C = NULL;
    }

    /* Allocate Gamma, as iterate released the memory */
    atom->Gamma = matrix_double(SQ(atom->Nlevel), atmos.Nspace);

    /* Initialise some continuum quantities */
    for (kr = 0; kr < atom->Ncont; kr++) {
      atom->continuum[kr].Rij = (double *) realloc(atom->continuum[kr].Rij,
						   atmos.Nspace * sizeof(double));
      atom->continuum[kr].Rji = (double *) realloc(atom->continuum[kr].Rji,
						   atmos.Nspace * sizeof(double));
      for (k = 0;  k < atmos.Nspace;  k++) {
          atom->continuum[kr].Rij[k] = 0.0;
          atom->continuum[kr].Rji[k] = 0.0;
      }
    }


    /* Initialise some line quantities */
    for (kr = 0;  kr < atom->Nline;  kr++) {
      line = &atom->line[kr];

      if (line->phi  != NULL) {
          freeMatrix((void **) line->phi);
          line->phi = NULL;
      }
      if (line->wphi != NULL) {
          free(line->wphi);
          line->wphi = NULL;
      }

      if (atmos.moving && line->polarizable && (input.StokesMode >= FIELD_FREE)) {

	if (line->phi_Q != NULL) {
	  freeMatrix((void **) line->phi_Q);
	  line->phi_Q = NULL;
	}
	if (line->phi_U != NULL) {
	  freeMatrix((void **) line->phi_U);
	  line->phi_U = NULL;
	}
	if (line->phi_V != NULL) {
	  freeMatrix((void **) line->phi_V);
	  line->phi_V = NULL;
	}


	if (input.magneto_optical) {
	  if (line->psi_Q != NULL) {
	    freeMatrix((void **) line->psi_Q);
	    line->psi_Q = NULL;
	  }
	  if (line->psi_U != NULL) {
	    freeMatrix((void **) line->psi_U);
	    line->psi_U = NULL;
	  }
	  if (line->psi_V != NULL) {
	    freeMatrix((void **) line->psi_V);
	    line->psi_V = NULL;
	  }
	}
      }
      /* realloc because of varying Nspace */
      line->Rij = (double *) realloc(line->Rij, atmos.Nspace * sizeof(double));
      line->Rji = (double *) realloc(line->Rji, atmos.Nspace * sizeof(double));

      for (k = 0;  k < atmos.Nspace;  k++) {
          line->Rij[k] = 0.0;
          line->Rji[k] = 0.0;
      }

      if (line->PRD) {
	if (line->Ng_prd != NULL) {
	  NgFree(line->Ng_prd);
	  line->Ng_prd = NULL;
	}

	if (line->fp_GII != NULL) {
	  fclose(line->fp_GII);
	  line->fp_GII = NULL;
	}

	if (input.PRD_angle_dep == PRD_ANGLE_DEP)
	  Nlamu = 2*atmos.Nrays * line->Nlambda;
	else
	  Nlamu = line->Nlambda;

	/* Idea: instead of doing this, why not free and just set line->rho_prd = NULL,
	   (and also line->Qelast?), because profile.c will act on that and reallocate */
	if (line->rho_prd != NULL) freeMatrix((void **) line->rho_prd);
	line->rho_prd = matrix_double(Nlamu, atmos.Nspace);

	if (line->Qelast != NULL) {
	  line->Qelast = (double *) realloc(line->Qelast, atmos.Nspace * sizeof(double));
	} else {
	  line->Qelast = (double *) malloc(atmos.Nspace * sizeof(double));
	}

	/* Initialize the ratio of PRD to CRD profiles to 1.0 */
	for (la = 0;  la < Nlamu;  la++) {
	  for (k = 0;  k < atmos.Nspace;  k++)
	    line->rho_prd[la][k] = 1.0;
	}

	/* reset interpolation weights */
	if (input.PRD_angle_dep == PRD_ANGLE_APPROX) {
	  Nlamu = 2*atmos.Nrays * line->Nlambda;
	  if (line->frac != NULL) {
	    freeMatrix((void **) line->frac);
	    line->frac = NULL;
	  }
	  if (line->id0 != NULL){
	    freeMatrix((void **) line->id0);
	    line->id0 = NULL;
	  }
	  if (line->id1 != NULL) {
	    freeMatrix((void **) line->id1);
	    line->id1 = NULL;
	  }
	}

      }
    }
  }

  distribute_nH();

  /* Update atmos-dependent molecular  quantities --- --------------- */
  for (nact = 0; nact < atmos.Nmolecule; nact++) {
    molecule = &atmos.molecules[nact];

    /* Reallocate some stuff, because of varying Nspace */
    molecule->vbroad = (double *) realloc(molecule->vbroad, atmos.Nspace * sizeof(double));
    molecule->pf     = (double *) realloc(molecule->pf,     atmos.Nspace * sizeof(double));
    molecule->n      = (double *) realloc(molecule->n,      atmos.Nspace * sizeof(double));
    if (molecule->nv != NULL) {
      freeMatrix((void **) molecule->nv);
      molecule->nv = matrix_double(molecule->Nv, atmos.Nspace);
    }
    if (molecule->nvstar != NULL) {
      freeMatrix((void **) molecule->nvstar);
      molecule->nvstar = matrix_double(molecule->Nv, atmos.Nspace);
    }
    if (molecule->pfv != NULL) {
      freeMatrix((void **) molecule->pfv);
      molecule->pfv = matrix_double(molecule->Nv, atmos.Nspace);
    }

    vtherm = 2.0*KBOLTZMANN / (AMU * molecule->weight);
    for (k = 0;  k < atmos.Nspace;  k++)
      molecule->vbroad[k] = sqrt(vtherm*atmos.T[k] + SQ(atmos.vturb[k]));

    if (molecule->active) {
      /* Allocate Gamma, as iterate released the memory */
      molecule->Gamma = matrix_double(SQ(molecule->Nv), atmos.Nspace);

      LTEmolecule(molecule);

      /* Free CO collision rate array, will be reallocated in initSolution_p */
      if (strstr(molecule->ID, "CO")) {
	free(molecule->C_ul);
	molecule->C_ul = NULL;
      }

      /* Free some line quantities */
      for (kr = 0;  kr < molecule->Nrt;  kr++) {
	mrt = &molecule->mrt[kr];
	if (mrt->phi  != NULL) {
	  freeMatrix((void **) mrt->phi);
	  mrt->phi = NULL;
	}
	if (mrt->wphi != NULL) {
	  free(mrt->wphi);
	  mrt->wphi = NULL;
	}
      }

    } else {
      if (molecule->Npf > 0) {
	for (k = 0;  k < atmos.Nspace;  k++)
	  molecule->pf[k] = partfunction(molecule, atmos.T[k]);
      }
    }
  }
}
/* ------- end   --------------------------  updateAtmosDep.c     --- */

/* ------- begin --------------------------  HRRR.c -----         --- */
void HERR(const char *rname) {
  /* Processes NetCDF errors */

  printf("Process %4d: (EEE) %s: HDF5 error.\n", mpi.rank, rname);
  MPI_Abort(mpi.comm, 2);

}
/* ------- end   --------------------------  HERR.c -----         --- */

/* ------- begin --------------------------  Error_p.c --         --- */
void Error(enum errorlevel level, const char *routineName,
	   const char *messageStr)
{
  char errorStr[MAX_MESSAGE_LENGTH];
  enum errorlevel defaultLevel = ERROR_LEVEL_1;

  switch (level) {
  case MESSAGE:
    if ((mpi.single_log) && (mpi.rank != 0)) return;
    if (!commandline.quiet)
      fprintf(commandline.logfile, "%s", (messageStr) ? messageStr : "");
    return;
  case WARNING:
    if ((mpi.single_log) && (mpi.rank != 0))
      fprintf(mpi.logfile, "\nProcess %4d: -WARNING in routine %s\n %s\n",
	      mpi.rank, routineName, (messageStr) ? messageStr : " (Undocumented)\n");
    fprintf(commandline.logfile, "\nProcess %d: -WARNING in routine %s\n %s\n",
	    mpi.rank, routineName, (messageStr) ? messageStr : " (Undocumented)\n");
    return;
  default:
    if (level < defaultLevel) {
      fprintf(commandline.logfile, "\a\n-ERROR in routine %s\n %s \n %s\n",
	      routineName,(messageStr) ? messageStr : " (Undocumented)\n",
	      "Trying to continue.....");
      if (commandline.logfile == mpi.logfile)
	fprintf(mpi.main_logfile, "\a\n-Process %4d: ERROR in routine %s\n %s \n %s\n",
		mpi.rank, routineName,(messageStr) ? messageStr : " (Undocumented)\n",
		"Trying to continue.....");
      return;
    } else {
      sprintf(errorStr, "\a\n\n-FATAL_ERROR in routine %s\n %s \n %s\n",
	      routineName,(messageStr) ? messageStr : " (Undocumented)\n",
	      "Exiting.....");
      if (commandline.logfile == mpi.logfile)
	fprintf(mpi.main_logfile, "\nProcess %4d: %s", mpi.rank, errorStr);

      fprintf(commandline.logfile, "%s", errorStr);
      if (commandline.logfile != stderr) fprintf(stderr, "%s", errorStr);

      /* Make exception for Singular matrix error */
      if (!strstr(messageStr,"Singular matrix")) {
	if (errno) perror(routineName);
	MPI_Abort(mpi.comm, level);
      }
    }
  }
}
/* ------- end    --------------------------  Error_p.c --         --- */

/* ------- begin -------------------------- copyBufVars.c ------------ */

void copyBufVars(bool_t writej) {
/* Copies output variables to buffer arrays, to be written only at the end */
  const  char routineName[] = "copyBufVars";
  static long ind = 0;
  int         i, ndep, nspect, nact, kr;
  Atom       *atom;
  Molecule   *molecule;
  AtomicLine      *line;
  AtomicContinuum *continuum;

  if (writej) {
    for (nspect=0; nspect < spectrum.Nspect; nspect++) {
      i = 0;
      for (ndep=mpi.zcut; ndep < infile.nz; ndep++, i++) {
        iobuf.J[(mpi.task*spectrum.Nspect + nspect)*infile.nz + ndep] =
            (float) spectrum.J[nspect][i];
      }
    }

    if (input.backgr_pol) {
      memcpy((void *) &iobuf.J20[ind*spectrum.Nspect], (void *) spectrum.J20[0],
      	   spectrum.Nspect * atmos.Nspace * sizeof(double));
    }
  }

  /* --- ATOM loop --- */
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    /* n, nstar */
    memcpy((void *) &iobuf.n[nact][ind*atom->Nlevel],     (void *) atom->n[0],
	   atom->Nlevel * atmos.Nspace * sizeof(double));
    memcpy((void *) &iobuf.nstar[nact][ind*atom->Nlevel], (void *) atom->nstar[0],
	   atom->Nlevel * atmos.Nspace * sizeof(double));

    /* Rij, Rji for lines */
    for (kr=0; kr < atom->Nline; kr++) {
      line = &atom->line[kr];
      memcpy((void *) &iobuf.RijL[nact][ind*atom->Nline + kr*atmos.Nspace],
	     (void *) line->Rij, atmos.Nspace * sizeof(double));
      memcpy((void *) &iobuf.RjiL[nact][ind*atom->Nline + kr*atmos.Nspace],
	     (void *) line->Rji, atmos.Nspace * sizeof(double));
    }

    /* Rij, Rji for continua */
    for (kr=0; kr < atom->Ncont; kr++) {
      continuum = &atom->continuum[kr];
      memcpy((void *) &iobuf.RijC[nact][ind*atom->Ncont + kr*atmos.Nspace],
	     (void *) continuum->Rij, atmos.Nspace * sizeof(double));
      memcpy((void *) &iobuf.RjiC[nact][ind*atom->Ncont + kr*atmos.Nspace],
	     (void *) continuum->Rji, atmos.Nspace * sizeof(double));
    }

  }

  /* --- MOLECULE loop --- */
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];

    /* nv, nvstar */
    memcpy((void *) &iobuf.nv[nact][ind*molecule->Nv],
	   (void *) molecule->nv[0],
	   molecule->Nv * atmos.Nspace * sizeof(double));
    memcpy((void *) &iobuf.nvstar[nact][ind*molecule->Nv],
	   (void *) molecule->nvstar[0],
	   molecule->Nv * atmos.Nspace * sizeof(double));

  }
  ind += atmos.Nspace;

  return;
}

/* ------- end ---------------------------- copyBufVars.c ------------ */

/* ------- begin -------------------------- allocBufVars.c ----------- */

void allocBufVars(bool_t writej) {
/* Allocates buffer arrays, to be written only at the end */
  const char routineName[] = "allocBufVars";
  long jsize = mpi.Ntasks*spectrum.Nspect*infile.nz;
  long nsize, RLsize, RCsize;
  int  nact;
  Atom      *atom;
  Molecule  *molecule;

  /* J, J20 */
  if (writej) {
    iobuf.J = (float *) calloc(jsize, sizeof(float));
    if (iobuf.J == NULL) Error(ERROR_LEVEL_2, routineName, "Out of memory\n");

    if (input.backgr_pol) {
      iobuf.J20 = (float *) calloc(jsize, sizeof(float));
      if (iobuf.J20 == NULL) Error(ERROR_LEVEL_2, routineName, "Out of memory\n");
    }
  }


  if (atmos.Nactiveatom > 0) {
    iobuf.n     = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
    iobuf.nstar = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
    iobuf.RijL  = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
    iobuf.RjiL  = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
    iobuf.RijC  = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
    iobuf.RjiC  = (double **) malloc(atmos.Nactiveatom * sizeof(double *));
  }

  if (atmos.Nactivemol > 0) {
    iobuf.nv     = (double **) malloc(atmos.Nactivemol * sizeof(double *));
    iobuf.nvstar = (double **) malloc(atmos.Nactivemol * sizeof(double *));
  }

  /* --- Loop over active ATOMS --- */
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    nsize  = mpi.Ntasks * atom->Nlevel * infile.nz * sizeof(double);
    RLsize = mpi.Ntasks * atom->Nline  * infile.nz * sizeof(double);
    RCsize = mpi.Ntasks * atom->Ncont  * infile.nz * sizeof(double);

    /* n, nstar */
    iobuf.n[nact]     = (double *) malloc(nsize);
    if (iobuf.n[nact] == NULL)
      Error(ERROR_LEVEL_2, routineName, "Out of memory\n");

    iobuf.nstar[nact] = (double *) malloc(nsize);
    if (iobuf.nstar[nact] == NULL)
      Error(ERROR_LEVEL_2, routineName, "Out of memory\n");

    /* Rij, Rji for lines */
    iobuf.RijL[nact] = (double *) malloc(RLsize);
    if (iobuf.RijL[nact] == NULL)
      Error(ERROR_LEVEL_2, routineName, "Out of memory\n");

    iobuf.RjiL[nact] = (double *) malloc(RLsize);
    if (iobuf.RjiL[nact] == NULL)
      Error(ERROR_LEVEL_2, routineName, "Out of memory\n");

    /* Rij, Rji for continua */
    iobuf.RijC[nact] = (double *) malloc(RCsize);
    if (iobuf.RijC[nact] == NULL)
      Error(ERROR_LEVEL_2, routineName, "Out of memory\n");

    iobuf.RjiC[nact] = (double *) malloc(RCsize);
    if (iobuf.RjiC[nact] == NULL)
      Error(ERROR_LEVEL_2, routineName, "Out of memory\n");
  }

  /* --- Loop over active MOLECULES --- */
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    nsize  = mpi.Ntasks * molecule->Nv * infile.nz * sizeof(double);

    /* nv, nvstar */
    iobuf.nv[nact]     = (double *) malloc(nsize);
    if (iobuf.nv[nact] == NULL)
      Error(ERROR_LEVEL_2, routineName, "Out of memory\n");

    iobuf.nvstar[nact] = (double *) malloc(nsize);
    if (iobuf.nvstar[nact] == NULL)
      Error(ERROR_LEVEL_2, routineName, "Out of memory\n");
  }
  return;
}

/* ------- end ---------------------------- allocBufVars.c ----------- */

/* ------- begin -------------------------- freeBufVars.c ------------ */

void freeBufVars(bool_t writej) {
  int nact;

  if (input.backgr_pol) free(iobuf.J20);

  /* Loop over active ATOMS */
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    free(iobuf.n[nact]);
    free(iobuf.nstar[nact]);
    free(iobuf.RijL[nact]);
    free(iobuf.RjiL[nact]);
    free(iobuf.RijC[nact]);
    free(iobuf.RjiC[nact]);
  }

  free(iobuf.n);
  free(iobuf.nstar);
  free(iobuf.RijL);
  free(iobuf.RjiL);
  free(iobuf.RijC);
  free(iobuf.RjiC);

  /* Loop over active MOLECULES */
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    free(iobuf.nv[nact]);
    free(iobuf.nvstar[nact]);
  }

  free(iobuf.nv);
  free(iobuf.nvstar);


  return;
}

/* ------- end ---------------------------- freeBufVars.c ------------ */

/* ------- begin -------------------------- writeOutput.c ------------ */
void writeOutput(bool_t writej) {
/* Writes all output files, in the case where output all at once is active */

  /* Write output in order of rank. First 0, then send to next, until all
     processes have written the output.

     This can also be used with an integer, like if mpi.rank > ml,
     and then adding ml to mpi.rank in the send. But for now not using
  if (mpi.rank > 0)
    MPI_Recv(&msg, 0, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mpi.comm, &status);
  */

  if (mpi.Ntasks == 0) {
    sprintf(messageStr, "Process %4d: *** NO WORK (more processes than "
                        "tasks!)\n", mpi.rank);
    fprintf(mpi.main_logfile, "%s", messageStr);
    Error(MESSAGE, "main", messageStr);
  } else {
    sprintf(messageStr, "Process %4d: --- START output\n", mpi.rank);
    fprintf(mpi.main_logfile, "%s", messageStr);
    Error(MESSAGE, "main", messageStr);

    writeMPI_all();
    writeAux_all();

    sprintf(messageStr, "Process %4d: *** END output\n", mpi.rank);
    fprintf(mpi.main_logfile, "%s", messageStr);
    Error(MESSAGE, "main", messageStr);
  }



  return;
}
/* ------- end ---------------------------- writeOutput.c ------------ */
