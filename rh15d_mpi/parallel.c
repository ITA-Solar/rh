/* ------- file: -------------------------- writeAux_p.c ---------

       Version:       rh2.0, 1.5-D plane-parallel
       Author:        Tiago Pereira (tiago.pereira@nasa.gov)
       Last modified: Mon Jan 03 14:28:25 2011 --

       --------------------------                      -----------RH-- */

/* --- Set of tools to deal with IO initialisation and closure ------- */

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
#include "parallel.h"
#include "io.h"


/* --- Function prototypes --                          -------------- */
void  distribute_nH(void);

/* --- Global variables --                             -------------- */
extern Atmosphere atmos;
extern Geometry geometry;
extern NCDF_Atmos_file infile;
extern IO_data io;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];
extern MPI_data mpi;

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
    sprintf(messageStr, "Process %d: Unable to open log file %s", 
	    mpi.rank, logfile);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  setvbuf(mpi.logfile, NULL, _IOLBF, BUFSIZ);

  mpi.stop    = FALSE;
  mpi.nconv   = 0;
  mpi.nnoconv = 0;
  mpi.ncrash  = 0;


  return;
}
/* ------- end   --------------------------   initParallel.c --   --- */

/* ------- begin --------------------------   initParallelIO.c    --- */
void initParallelIO(bool_t run_ray) {
  int    nact;
  Atom  *atom;

  init_ncdf_aux();
  init_ncdf_J();
  init_Background();
  if (!run_ray) {
    if (input.p15d_wspec) init_ncdf_spec();
    init_ncdf_indata();
  }

  /* Get file position of atom files (to re-read collisions) */
  io.atom_file_pos = (long *) malloc(atmos.Nactiveatom * sizeof(long));

  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];
    io.atom_file_pos[nact] = ftell(atom->fp_input);
  }

  /* Save StokesMode (before adjustStokesMode changes it...) */
  mpi.StokesMode_save = input.StokesMode;
  
  mpi.single_log       = FALSE;  

  return;
}
/* ------- end   --------------------------   initParallelIO.c    --- */



/* ------- begin --------------------------  closeParallelIO.c    --- */
void closeParallelIO(bool_t run_ray) {

  close_ncdf_atmos(&atmos, &geometry, &infile);
  close_ncdf_aux();
  close_ncdf_J();
  close_Background();
  if (!run_ray) {
    if (input.p15d_wspec) close_ncdf_spec();
    close_ncdf_indata();
  }

  free(io.atom_file_pos);

  return;
}
/* ------- end   --------------------------  closeParallelIO.c    --- */

/* ------- begin --------------------------  updateAtmosDep.c     --- */
void UpdateAtmosDep(void) {
/* Updates the atmos-dependent factors for the atoms and molecules */
  const char routineName[] = "UpdateAtomsDep";
  int       ierror, nact, k, kr, la, Nlamu;
  double    vtherm;
  Atom     *atom;
  Molecule *molecule;
  AtomicLine    *line;
  MolecularLine *mrt;

  /* Put back initial Stokes mode */
  input.StokesMode = mpi.StokesMode_save;


  /* Update atmos-dependent atomic  quantities --- --------------- */
  for (nact = 0; nact < atmos.Natom; nact++) {
    atom = &atmos.atoms[nact];

    /* Reallocate some stuff (because of varying Nspace) */
    atom->ntotal = (double *) realloc(atom->ntotal, atmos.Nspace * sizeof(double));
    atom->vbroad = (double *) realloc(atom->vbroad, atmos.Nspace * sizeof(double));

    if (atom->n != NULL) {
      freeMatrix((void **) atom->n);
      atom->n = matrix_double(atom->Nlevel, atmos.Nspace);
    }
    
    if (atom->nstar != NULL) {
      freeMatrix((void **) atom->nstar);
      atom->nstar = matrix_double(atom->Nlevel, atmos.Nspace);
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

    
    /* Rewind atom files to point just before collisional data */
    if ((ierror = fseek(atom->fp_input, io.atom_file_pos[nact], SEEK_SET))) {
      sprintf(messageStr, "Unable to rewind atom file for %s", atom->ID);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    
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
      
      if (atmos.moving && line->polarizable && (input.StokesMode>FIELD_FREE)) {
	
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

	if (input.PRD_angle_dep)
	  Nlamu = 2*atmos.Nrays * line->Nlambda;
	else
	  Nlamu = line->Nlambda;
	

	// Idea: instead of doing this, why not free and just set line->rho_prd = NULL,
	// (and also line->Qelast?), because profile.c will act on that and reallocate
	freeMatrix((void **) line->rho_prd);
	line->rho_prd = matrix_double(Nlamu, atmos.Nspace);
	
	line->Qelast = (double *) realloc(line->Qelast, atmos.Nspace * sizeof(double));

	/* Initialize the ratio of PRD to CRD profiles to 1.0 */
	for (la = 0;  la < Nlamu;  la++) {
	  for (k = 0;  k < atmos.Nspace;  k++)
	    line->rho_prd[la][k] = 1.0;
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


  return;
}
/* ------- end   --------------------------  updateAtmosDep.c     --- */

/* ------- begin --------------------------  ERR.c ------         --- */
void ERR(int ierror, const char *rname) {
  /* Processes NetCDF errors */


  /* For conversion not representable (ie, NaNs, Inf, etc.), return so
     that process can go to next task */
  if (ierror == NC_ERANGE) {
    Error(WARNING, rname, "NetCDF: Numeric conversion not representable");
    return; 
  } 

  printf("Process %d: (EEE) %s: %d %s\n", mpi.rank, rname, ierror,
	 nc_strerror(ierror));
  MPI_Abort(mpi.comm, 2);

}
/* ------- end   --------------------------  ERR.c ------         --- */

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
      fprintf(mpi.logfile, "\nProcess %d: -WARNING in routine %s\n %s\n",
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
	fprintf(mpi.main_logfile, "\a\n-Process %d: ERROR in routine %s\n %s \n %s\n",
		mpi.rank, routineName,(messageStr) ? messageStr : " (Undocumented)\n",
		"Trying to continue.....");
      return;
    } else {
      sprintf(errorStr, "\a\n\n-FATAL_ERROR in routine %s\n %s \n %s\n",
	      routineName,(messageStr) ? messageStr : " (Undocumented)\n",
	      "Exiting.....");
      if (commandline.logfile == mpi.logfile) 
	fprintf(mpi.main_logfile, "\nProcess %d: %s", mpi.rank, errorStr);

      fprintf(commandline.logfile, "%s", errorStr);
      if (commandline.logfile != stderr) fprintf(stderr, "%s", errorStr);

      /* Make exception for Singular matrix error */
      if (!strstr(messageStr,"Singular matrix")) {
	if (errno) perror(routineName);
	exit(level);
      }
    }
  }
}
/* ------- end    --------------------------  Error_p.c --         --- */
