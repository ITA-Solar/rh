/* ------- file: -------------------------- collision.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Sat Sep 19 15:54:07 2009 --

       --------------------------                      ----------RH-- */

/* --- Reads collisional data from atomic data file and computes
       collisional rates.

       Format is similar to MULTI's routine GENCOL:

         Input line can be either

           TEMP  Nitem  T[0]     ...   T[Nitem-1]

         to read in the temperature grid for collisional coefficients, or

           KEYWORD  i1  i2   coeff[0]    ...   coeff[Nitem-1]

         to read the coefficients for transitions between i1 and i2.
         Multiple entries of the temperature grid are allowed, but
         at least one entry with the correct number of grid points
         has to precede an entry with coefficients.

         Allowed keywords are:

         Keyword    Transition type
         ----------------------------------------------------

         TEMP  -->  Temperature grid

         OMEGA -->  Collisional de-excitation of ions by electrons
         CE    -->  Collisional de-excitation of neutrals by electrons
         CI    -->  Collisional ionization by electrons
         CP    -->  Collisional de-excitation by protons

         CH0   -->  Charge exchange of ion with neutral hydrogen
         CH+   -->  Charge exchange of neutral with protons

         END   -->  End of input data
         ----------------------------------------------------

 Note: Unit of number density is m^-3.

       Convention: C_ij = C[i][j] represents the
                   transition j --> i

       --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "statistics.h"

#define COMMENT_CHAR "#"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- CollisionRate.c --------- */

void CollisionRate(struct Atom *atom, FILE *fp_atom)
{
  const char routineName[] = "CollisionRate";
  register int k, n;

  char    inputLine[MAX_LINE_SIZE], *keyword, *pointer,
          labelStr[MAX_LINE_SIZE];
  bool_t  hunt, exit_on_EOF;
  int     nitem, i1, i2, i, j, ij, ji, Nlevel = atom->Nlevel, Nitem,
          status;
  long    Nspace = atmos.Nspace;
  double  dE, C0, *T, *coeff, *C, Cdown, Cup, gij, *np;

  getCPU(3, TIME_START, NULL);

  C0 = ((E_RYDBERG/sqrt(M_ELECTRON)) * PI*SQ(RBOHR)) *
    sqrt(8.0/(PI*KBOLTZMANN));

  atom->C = matrix_double(SQ(Nlevel), Nspace);
  for (ij = 0;  ij < SQ(Nlevel);  ij++) {
    for (k = 0;  k < Nspace;  k++) {
      atom->C[ij][k] = 0.0;
    }
  }
  C = (double *) malloc(Nspace * sizeof(double));

  T = coeff = NULL;
  while ((status = getLine(fp_atom, COMMENT_CHAR,
		  inputLine, exit_on_EOF=FALSE)) != EOF) {
    keyword = strtok(inputLine, " ");

    if (strstr(keyword, "TEMP")) {

      /* --- Read temperature grid --                  -------------- */

      Nitem = atoi(strtok(NULL, " "));
      T = (double *) realloc(T, Nitem*sizeof(double));
      for (n = 0, nitem = 0;  n < Nitem;  n++) {
        if ((pointer = strtok(NULL, " ")) == NULL) break;
	nitem += sscanf(pointer, "%lf", T+n);
      }
    } else if (strstr(keyword, "OMEGA") || strstr(keyword, "CE") ||
	       strstr(keyword, "CI")    || strstr(keyword, "CP") ||
	       strstr(keyword, "CH0")   || strstr(keyword, "CH+")) {

      /* --- Read level indices and collision coefficients -- ------- */

      i1 = atoi(strtok(NULL, " "));
      i2 = atoi(strtok(NULL, " "));
      coeff = (double *) realloc(coeff, Nitem*sizeof(double));
      for (n = 0, nitem = 0;  n < Nitem;  n++) {
        if ((pointer = strtok(NULL, " ")) == NULL) break;
	nitem += sscanf(pointer, "%lf", coeff+n);
      }
      /* --- Transitions i -> j are stored at index ji, transitions
	     j -> i are stored under ij. --            -------------- */

      i  = MIN(i1, i2);
      j  = MAX(i1, i2);
      ij = i*Nlevel + j;
      ji = j*Nlevel + i;

    } else if (strstr(keyword, "END")) {
      break;
    } else {
      sprintf(messageStr, "Unknown keyword: %s", keyword);
      Error(ERROR_LEVEL_1, routineName, messageStr);
    }

    if (nitem != Nitem) {
      sprintf(messageStr, "\n Read %d, not %d items (keyword = %s)\n",
	      nitem, Nitem, keyword);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    /* --- Spline interpolation in temperature T for all spatial
           locations. Linear if only 2 interpolation points given - - */

    if (!strstr(keyword, "TEMP")) {
      if (Nitem > 2) {
	splineCoef(Nitem, T, coeff);
	splineEval(Nspace, atmos.T, C, hunt=TRUE);
      } else
	Linear(Nitem, T, coeff, Nspace, atmos.T, C, hunt=TRUE);
    }

    if (strstr(keyword, "OMEGA")) {

      /* --- Collisional excitation of ions --         -------------- */ 

      for (k = 0;  k < Nspace;  k++) {
        Cdown = C0 * atmos.ne[k] * C[k] /
                                 (atom->g[j] * sqrt(atmos.T[k]));
	atom->C[ij][k] += Cdown;
	atom->C[ji][k] += Cdown * atom->nstar[j][k]/atom->nstar[i][k];
      }
    } else if (strstr(keyword, "CE")) {      

      /* --- Collisional excitation of neutrals --     -------------- */ 

      gij = atom->g[i] / atom->g[j];
      for (k = 0;  k < Nspace;  k++) {
        Cdown = C[k] * atmos.ne[k] * gij * sqrt(atmos.T[k]);
	atom->C[ij][k] += Cdown;
	atom->C[ji][k] += Cdown * atom->nstar[j][k]/atom->nstar[i][k];
      }
    } else if (strstr(keyword, "CI")) {      

      /* --- Collisional ionization --                 -------------- */

      dE = atom->E[j] - atom->E[i];
      for (k = 0;  k < Nspace;  k++) {
        Cup = C[k] * atmos.ne[k] *
	  exp(-dE/(KBOLTZMANN*atmos.T[k])) * sqrt(atmos.T[k]);
	atom->C[ji][k] += Cup;
	atom->C[ij][k] += Cup * atom->nstar[i][k]/atom->nstar[j][k];
      }
    } else if (strstr(keyword, "CP")) {

      /* --- Collisions with protons --                -------------- */

      np = atmos.H->n[atmos.H->Nlevel-1];
      for (k = 0;  k < Nspace;  k++) {
        Cdown = np[k] * C[k];
	atom->C[ij][k] += Cdown;
	atom->C[ji][k] += Cdown * atom->nstar[j][k]/atom->nstar[i][k];
      }
    } else if (strstr(keyword, "CH0")) {

      /* --- Charge exchange with neutral hydrogen --  -------------- */

      for (k = 0;  k < Nspace;  k++)
	atom->C[ij][k] += atmos.H->n[0][k] * C[k];

    } else if (strstr(keyword, "CH+")) {

      /* --- Charge exchange with protons --           -------------- */

      np = atmos.H->n[atmos.H->Nlevel-1];
      for (k = 0;  k < Nspace;  k++)
	atom->C[ji][k] += np[k] * C[k];
    }
  }
  if (status == EOF) {
    sprintf(messageStr, "Reached end of datafile before all data was read");
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }
  /* --- Clean up --                                   -------------- */

  free(C);
  free(T);
  free(coeff);

  sprintf(labelStr, "Collision Rate %2s", atom->ID);
  getCPU(3, TIME_POLL, labelStr);
}
/* ------- end ---------------------------- CollisionRate.c --------- */
