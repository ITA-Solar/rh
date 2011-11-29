/* ------- file: -------------------------- iterate.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek  (huitenbroek@nso.edu)
       Last modified: Wed Apr  1 14:08:31 2009 --

       --------------------------                      ----------RH-- */

/* --- Main iteration routine --                       -------------- */

 
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "accelerate.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"


typedef struct {
  bool_t eval_operator, redistribute;
  int    nspect;
  double dJ;
} threadinfo;

/* --- Function prototypes --                          -------------- */

void *Formal_pthread(void *argument);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- Iterate.c --------------- */

void Iterate(int NmaxIter, double iterLimit)
{
  const char routineName[] = "Iterate";
  register int niter, nact;
  double cswitch;

  bool_t eval_operator;
  double dpopsmax, PRDiterlimit;
  Atom *atom;
  Molecule *molecule;
  AtomicLine *line;        // Tiago: DELETE
  int i, mu, to_obs, lamu; // Tiago: DELETE

  if (NmaxIter <= 0) return;
  getCPU(1, TIME_START, NULL);

  /* --- Initialize structures for Ng acceleration of population
         convergence --                                  ------------ */

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    atom->Ng_n = NgInit(atom->Nlevel*atmos.Nspace, input.Ngdelay,
			input.Ngorder, input.Ngperiod, atom->n[0]);
  }
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    molecule->Ng_nv = NgInit(molecule->Nv*atmos.Nspace, input.Ngdelay,
			       input.Ngorder, input.Ngperiod,
			       molecule->nv[0]);
  }
  /* --- Start of the main iteration loop --             ------------ */

  niter = 1;
  
  /* Collisional-radiative switching ? */
  if (input.crsw != 0.0)
    cswitch = input.crsw_ini;
  else
    cswitch = 1.0;
    
  /* PRD switching ? */
  if (input.prdsw > 0.0)
    input.prdswitch = 0.0;
  else
    input.prdswitch = 1.0;
  
  while (niter <= NmaxIter && !StopRequested()) {
    getCPU(2, TIME_START, NULL);


    for (nact = 0;  nact < atmos.Nactiveatom;  nact++)
      initGammaAtom(atmos.activeatoms[nact], cswitch);
    for (nact = 0;  nact < atmos.Nactivemol;  nact++)
      initGammaMolecule(atmos.activemols[nact]);

    /* --- Formal solution for all wavelengths --      -------------- */

    solveSpectrum(eval_operator=TRUE, FALSE);

    /* --- Solve statistical equilibrium equations --  -------------- */

    sprintf(messageStr, "\n -- Iteration %3d, switch = %.2f, prd switch = %.2f\n",
	    niter, cswitch, input.prdswitch);
    Error(MESSAGE, routineName, messageStr);
    dpopsmax = updatePopulations(niter);

    if (atmos.NPRDactive > 0) {
      
      /* --- Redistribute intensity in PRD lines if necessary -- ---- */

      if (input.PRDiterLimit < 0.0)
	PRDiterlimit = MAX(dpopsmax, -input.PRDiterLimit);
      else
	PRDiterlimit = input.PRDiterLimit;
      Redistribute(input.PRD_NmaxIter, PRDiterlimit);
    }

    sprintf(messageStr, "Total Iteration %3d", niter);
    getCPU(2, TIME_POLL, messageStr);

    if (dpopsmax < iterLimit && cswitch <= 1.0 && input.prdswitch == 1.0 ) break;
    niter++;
    
    /* Update collisional multiplier factor */
    if (input.crsw > 0)
      cswitch = MAX(1.0, cswitch * pow(0.1, 1./input.crsw));

    /* Update PRD switching */ 
    if (input.prdsw > 0.0)
      input.prdswitch = MIN(1.0, input.prdsw * (double) (niter * niter) ); // quadratic, for now



    if (atmos.hydrostatic) {
      if (!atmos.atoms[0].active) {
	sprintf(messageStr, "Can only perform hydrostatic equilibrium"
                            " for hydrogen active");
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      Hydrostatic(N_MAX_HSE_ITER, HSE_ITER_LIMIT);
    }
  }
  
      // Tiago: temporary printouts to get PRD rho after iteration
      /*
       atom = atmos.activeatoms[0];
       line = &atom->line[0];
       
       switch (input.PRD_angle_dep) {
	case PRD_ANGLE_INDEP:
	  printf("rho_prd = \n");
	  for (i = 0; i < line->Nlambda; i++) {
	    printf("%8.4f   %e   %e   %e   %e   %e\n", line->lambda[i], line->rho_prd[i][105], line->rho_prd[i][110], line->rho_prd[i][120], line->rho_prd[i][150], line->rho_prd[i][155]);
	  }
          //exit(1);
	  break;
	
	case PRD_ANGLE_DEP:
	    for (mu = 0; mu < atmos.Nrays; mu++) {
	      for (to_obs = 0; to_obs <= 1; to_obs++) {
	       for (i = 0; i < line->Nlambda; i++) {
		lamu = 2*(atmos.Nrays*i + mu) + to_obs;
		if ((to_obs == 1) && (mu == 4))
		printf("%8.4f  %e   %e   %e   %e   %e\n", line->lambda[i], line->rho_prd[lamu][105],line->rho_prd[lamu][110],line->rho_prd[lamu][120], line->rho_prd[lamu][150], line->rho_prd[lamu][155] );
	      }
	    }
	  }
	  //exit(1);
          break;
       }
      */
      
      
  

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    freeMatrix((void **) atom->Gamma);
    NgFree(atom->Ng_n);
  } 
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    freeMatrix((void **) molecule->Gamma);
    NgFree(molecule->Ng_nv);
  }

  getCPU(1, TIME_POLL, "Iteration Total");
}
/* ------- end ---------------------------- Iterate.c --------------- */

/* ------- begin -------------------------- solveSpectrum.c --------- */

double solveSpectrum(bool_t eval_operator, bool_t redistribute)
{
  register int nspect, n, nt;

  int         Nthreads, lambda_max;
  double      dJ, dJmax;
  pthread_t  *thread_id;
  threadinfo *ti;

  /* --- Administers the formal solution for each wavelength. When
         input.Nthreads > 1 the solutions are performed concurrently
         in Nthreads threads. These are POSIX style threads.

    See: - David R. Butenhof, Programming with POSIX threads,
           Addison & Wesley.

         - Multithreaded Programming Guide, http://sun.docs.com
           (search for POSIX threads).

         When solveSpectrum is called with redistribute == TRUE only
         wavelengths that have an active PRD line are solved. The
         redistribute key is passed to the addtoRates routine via
         Formal so that only the radiative rates of PRD lines are
         updated. These are needed for the emission profile ratio \rho.
         --                                            -------------- */

  getCPU(3, TIME_START, NULL);

  /* --- First zero the radiative rates --             -------------- */

  zeroRates(redistribute);
  lambda_max = 0;
  dJmax = 0.0;

  if (input.Nthreads > 1) {

    /* --- If input.Nthreads positive then solve Nthreads wavelengths
           concurrently in separate threads --         -------------- */

    ti = (threadinfo *) malloc(input.Nthreads * sizeof(threadinfo));
    for (nt = 0;  nt < input.Nthreads;  nt++) {
      ti[nt].eval_operator = eval_operator;
      ti[nt].redistribute  = redistribute;
    }
    thread_id = (pthread_t *) malloc(input.Nthreads * sizeof(pthread_t));

    /* --- Thread management is very simple. Submit a batch of as many
           as input.Nthreads at the same time, then wait till all of
           these have finished. There is no check on successful
           submission nor completion. --               -------------- */

    for (nspect = 0;  nspect < spectrum.Nspect;  nspect += input.Nthreads) {
      if (nspect + input.Nthreads <= spectrum.Nspect)
	Nthreads = input.Nthreads;
      else
	Nthreads = (spectrum.Nspect % input.Nthreads);

      /* --- Start batch of concurrent threads --      -------------- */

      for (nt = 0;  nt < Nthreads;  nt++) {
	ti[nt].nspect = nspect + nt;
	if (!redistribute ||
	    (redistribute && containsPRDline(&spectrum.as[nspect+nt]))) {
	  pthread_create(&thread_id[nt], &input.thread_attr,
			 Formal_pthread, &ti[nt]);
	} else
	  thread_id[nt] = 0;
      }
      /* --- Let the finished threads of the batch join again -- ---- */

      for (nt = 0;  nt < Nthreads;  nt++) {
	if (thread_id[nt]) {
	  pthread_join(thread_id[nt], NULL);
	  if (ti[nt].dJ > dJmax) {
	    dJmax = ti[nt].dJ;
	    lambda_max = nspect + nt;
	  }
	}
      }
    }
    free(thread_id);
    free(ti);
  } else {
      
    /* --- Else call the solution for wavelengths sequentially -- --- */
      
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
      if (!redistribute ||
	  (redistribute && containsPRDline(&spectrum.as[nspect]))) {
	dJ = Formal(nspect, eval_operator, redistribute);
	if (dJ > dJmax) {
	  dJmax = dJ;
	  lambda_max = nspect;
	}
      }
    }
  }

  sprintf(messageStr, " Spectrum max delta J = %6.4E (lambda#: %d)\n",
	  dJmax, lambda_max);
  Error(MESSAGE, NULL, messageStr);

  getCPU(3, TIME_POLL,
	 (eval_operator) ? "Spectrum & Operator" : "Solve Spectrum");

  return dJmax;
}
/* ------- end ---------------------------- solveSpectrum.c --------- */

/* ------- begin -------------------------- Formal_pthread.c -------- */

void *Formal_pthread(void *argument)
{
  threadinfo *ti = (threadinfo *) argument;

  /* --- Threads wrapper around Formal --              -------------- */ 

  ti->dJ = Formal(ti->nspect, ti->eval_operator, ti->redistribute);

  return (NULL);
}
/* ------- end ---------------------------- Formal_pthread.c -------- */
