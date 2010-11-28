/* ------- file: -------------------------- solvene.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Sep 30 16:40:26 2005 --

       --------------------------                      ----------RH-- */

/* --- Various routines to solve for the electron density, given:

       T     -- The kinetic electron temperature
       nHtot -- The total hydrogen density

       Assuming Saha - Boltzmann for the equilibrium between the
       ionization stages.

  See: D. Mihalas (1978), in "Stellar Atmospheres", pp. 114-119


       When keyword fromscratch is set the electron density is
       calculated from scratch, using pure Hydrogen ionization
       as initial guess. Otherwise, the values passed in ne are used.
       --                                              -------------- */
 
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "error.h"
#include "statistics.h"

#define MAX_ELECTRON_ERROR         1.0E-2
#define N_MAX_ELECTRON_ITERATIONS  10
#define N_MAX_ELEMENT              26


/* --- Function prototypes --                          -------------- */

double getKuruczpf(Element *element, int stage, int k);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- Solve_ne.c -------------- */

void Solve_ne(double *ne, bool_t fromscratch)
{
  const char routineName[] = "Solvene";
  register int k, n, j;

  int     Nmaxstage, niter;
  double *fjk, *dfjk, error, ne_old, akj, sum, PhiH, C1, Uk,
          dne, dnemax;

  getCPU(3, TIME_START, NULL);

  C1 = (HPLANCK/(2.0*PI*M_ELECTRON)) * (HPLANCK/KBOLTZMANN);

  /* --- Figure out the largest array size needed so that we do not
         have to allocate and free memory all the time -- ----------- */

  Nmaxstage = 0;
  for (n = 0;  n < N_MAX_ELEMENT;  n++)
    Nmaxstage = MAX(Nmaxstage, atmos.elements[n].Nstage);
  fjk  = (double *) malloc(Nmaxstage * sizeof(double));
  dfjk = (double *) malloc(Nmaxstage * sizeof(double));

  for (k = 0;  k < atmos.Nspace;  k++) {
    if (fromscratch) {

      /* --- Get the initial solution from ionization of H only -- -- */

      Uk = getKuruczpf(&atmos.elements[0], 0, k);
      PhiH = 0.5 * pow(C1/atmos.T[k], 1.5) *
	exp(Uk + atmos.elements[0].ionpot[0]/(KBOLTZMANN*atmos.T[k]));
      ne_old = (sqrt(1.0 + 4.0*atmos.nHtot[k]*PhiH) - 1.0) / (2.0*PhiH);

      /* --- Copy into ne as well to calculate first fij and dfij - - */

      ne[k] = ne_old;
    } else {
      /* --- Use original electron density as starting guess -- ----- */

      ne_old = ne[k]; 
    }

    niter = 0;
    while (niter < N_MAX_ELECTRON_ITERATIONS) {
      error = ne_old / atmos.nHtot[k];
      sum   = 0.0;

      for (n = 0;  n < N_MAX_ELEMENT;  n++) {
	getfjk(&atmos.elements[n], ne_old, k, fjk, dfjk);
	for (j = 1;  j < atmos.elements[n].Nstage;  j++) {
	  akj = atmos.elements[n].abund * j;
	  error -= akj * fjk[j];
	  sum   += akj * dfjk[j];
	}
      }
      ne[k] = ne_old -
	atmos.nHtot[k] * error / (1.0 - atmos.nHtot[k] * sum);
      dne = fabs((ne[k] - ne_old)/ne_old);
      ne_old = ne[k];
    
      if (dne <= MAX_ELECTRON_ERROR) break;
      niter++;
    }

    if (dne > MAX_ELECTRON_ERROR) {
      sprintf(messageStr, "Electron density iteration not converged:\n"
	      " spatial location: %d, temperature: %6.1f [K], \n"
	      " density: %9.3E [m^-3],\n dnemax: %9.3E\n",
	      k, atmos.T[k], atmos.nHtot[k], dne);
      Error(WARNING, routineName, messageStr);
    }
  }

  free(fjk);  free(dfjk);

  getCPU(3, TIME_POLL, "Electron density");
}
/* ------- end ---------------------------- Solve_ne.c -------------- */

/* ------- begin -------------------------- getfjk.c ---------------- */

void getfjk(Element *element, double ne, int k, double *fjk, double *dfjk)
{
  register int j;

  double C1, sum1, sum2, CT_ne, Uk, Ukp1;

  /* --- Get the fractional population f_j(ne, T) = N_j/N for element
         element and its partial derivative with ne. -- ------------- */

  C1 = (HPLANCK/(2.0*PI*M_ELECTRON)) * (HPLANCK/KBOLTZMANN);

  CT_ne   = 2.0 * pow(C1/atmos.T[k], -1.5) / ne;
  sum1    = 1.0;
  sum2    = 0.0;
  fjk[0]  = 1.0;
  dfjk[0] = 0.0;

  Uk = getKuruczpf(element, 0, k);
  for (j = 1;  j < element->Nstage;  j++) {
    Ukp1 = getKuruczpf(element, j, k);
    
    fjk[j]  = fjk[j-1] * CT_ne *
      exp(Ukp1 - Uk - element->ionpot[j-1]/(KBOLTZMANN*atmos.T[k]));
    dfjk[j] = -j * fjk[j] / ne;
    sum1   += fjk[j];
    sum2   += dfjk[j];
    Uk      = Ukp1;
  }

  for (j = 0;  j < element->Nstage;  j++) {
    fjk[j]  /= sum1;
    dfjk[j]  = (dfjk[j] - fjk[j] * sum2) / sum1;
  }
}
/* ------- end ---------------------------- getfjk.c ---------------- */

/* ------- begin -------------------------- getKuruczpf.c ----------- */

double getKuruczpf(Element *element, int stage, int k)
{
  bool_t hunt = TRUE;
  double Uk;

  Linear(atmos.Npf, atmos.Tpf, element->pf[stage], 
	 1, &atmos.T[k], &Uk, hunt);

  return Uk;
}
/* ------- end ---------------------------- getKuruczpf.c ----------- */
