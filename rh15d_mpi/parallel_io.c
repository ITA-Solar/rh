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
extern IO_data io;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];
extern MPI_data mpi;

/* ------- begin --------------------------   initParallel.c --   --- */
void initParallel(int *argc, char **argv[]) {
  const char routineName[] = "initParallelIO";
  char   logfile[MAX_LINE_SIZE];

  /* Initialise MPI */
  MPI_Init(argc,argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi.size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
  MPI_Get_processor_name(mpi.name, &mpi.namelen);

  mpi.comm = MPI_COMM_WORLD;
  mpi.info = MPI_INFO_NULL;

  /* Open log files */ 
  sprintf(logfile, MPILOG_TEMPLATE, mpi.rank);
  if ((mpi.logfile = fopen(logfile, "w")) == NULL) {
    sprintf(messageStr, "Process %d: Unable to open log file %s", 
	    mpi.rank, logfile);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  setvbuf(mpi.logfile, NULL, _IOLBF, BUFSIZ);

  /* Keep log stream in main log */
  mpi.single_log = TRUE;

  mpi.stop    = FALSE;
  mpi.nconv   = 0;
  mpi.nnoconv = 0;
  mpi.ncrash  = 0;


  return;
}
/* ------- end   --------------------------   initParallel.c --   --- */

/* ------- begin --------------------------   initParallelIO.c    --- */
void initParallelIO(void) {
  int    nact;
  Atom  *atom;

  init_ncdf_aux();
  init_ncdf_J();
  init_Background();
  init_ncdf_spec();
  init_ncdf_indata();

  /* Get file position of atom files (to re-read collisions) */
  io.atom_file_pos = (long *) malloc(atmos.Nactiveatom * sizeof(long));

  for (nact = 0; nact < atmos.Nactiveatom; nact++) {
    atom = atmos.activeatoms[nact];
    io.atom_file_pos[nact] = ftell(atom->fp_input);
  }

  /* Direct log stream into MPI log files */
  commandline.logfile  = mpi.logfile;
  mpi.single_log       = FALSE;

  

  return;
}
/* ------- end   --------------------------   initParallelIO.c    --- */



/* ------- begin --------------------------  closeParallelIO.c    --- */
void closeParallelIO(void) {

  close_ncdf_aux();
  close_ncdf_J();
  close_Background();
  close_ncdf_spec();
  close_ncdf_indata();

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
  

  /* Update atmos-dependent atomic  quantities --- --------------- */
  for (nact = 0; nact < atmos.Natom; nact++) {
    atom = &atmos.atoms[nact];
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

      /* Free collision rate array, will be reallocated by calls in Background_p */
      freeMatrix((void **) atom->C);

      /* Allocate Gamma, as iterate released the memory */
      atom->Gamma = matrix_double(SQ(atom->Nlevel), atmos.Nspace);

      /* Initialise some line quantities */
      for (kr = 0;  kr < atom->Nline;  kr++) {
	line = &atom->line[kr];

	freeMatrix((void **) line->phi);
	free(line->wphi);

	if (atmos.moving && line->polarizable && (input.StokesMode>FIELD_FREE)) {
	  
	  freeMatrix((void **) line->phi_Q);
	  freeMatrix((void **) line->phi_U);
	  freeMatrix((void **) line->phi_V);

	  if (input.magneto_optical) {
	    freeMatrix((void **) line->psi_Q);
	    freeMatrix((void **) line->psi_U);
	    freeMatrix((void **) line->psi_V);
	  }
	}

	
	for (k = 0;  k < atmos.Nspace;  k++) {
	  line->Rij[k] = 0.0;
	  line->Rji[k] = 0.0;
	}
	

	if (line->PRD) {
	  NgFree(line->Ng_prd);
	  line->Ng_prd = NULL;

	  fclose(line->fp_GII);
	  line->fp_GII = NULL;

	  if (input.PRD_angle_dep)
	    Nlamu = 2*atmos.Nrays * line->Nlambda;
	  else
	    Nlamu = line->Nlambda;

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

    vtherm = 2.0*KBOLTZMANN / (AMU * molecule->weight);
    for (k = 0;  k < atmos.Nspace;  k++)
      molecule->vbroad[k] = sqrt(vtherm*atmos.T[k] + SQ(atmos.vturb[k]));
    
    if (molecule->active) {
      LTEmolecule(molecule);

      /* Free CO collision rate array, will be reallocated in initSolution_p */
      if (strstr(molecule->ID, "CO")) free(molecule->C_ul);

      /* Free some line quantities */
      for (kr = 0;  kr < molecule->Nrt;  kr++) {
	mrt = &molecule->mrt[kr];
	freeMatrix((void **) mrt->phi);
	free(mrt->wphi);
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
      fprintf(mpi.logfile, "\n-WARNING in routine %s\n %s\n",
	      routineName, (messageStr) ? messageStr : " (Undocumented)\n");
    fprintf(commandline.logfile, "\n-WARNING in routine %s\n %s\n",
	    routineName, (messageStr) ? messageStr : " (Undocumented)\n");
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
