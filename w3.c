/* ------- file: -------------------------- w3.c --------------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Mar 31 14:44:25 1999 --

       --------------------------                      ----------RH-- */

/* --- Functions for piecewise quadratic integration -- ------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rh.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- w2.c -------------------- */

void w2(double dtau, double *w)
{
  double expdt;

  if (dtau < 5.0E-4) {
    w[0] = dtau*(1.0 - 0.5*dtau);
    w[1] = SQ(dtau) * (0.5 - 0.33333333*dtau);
  } else if (dtau > 50.0) {
    w[1] = w[0] = 1.0;
  } else {
    expdt = exp(-dtau);
    w[0]  = 1.0 - expdt;
    w[1]  = w[0] - dtau*expdt;
  }
}
/* ------- end ---------------------------- w2.c -------------------- */

/* ------- begin -------------------------- w3.c -------------------- */

void w3(double dtau, double *w)
{
  double expdt, delta;

  if (dtau < 5.0E-4) {
    w[0]   = dtau*(1.0 - 0.5*dtau);
    delta  = SQ(dtau);
    w[1]   = delta*(0.5 - 0.33333333*dtau);
    delta *= dtau;
    w[2]   = delta*(0.33333333 - 0.25*dtau);
  } else if (dtau > 50.0) {
    w[1] = w[0] = 1.0;
    w[2] = 2.0;
  } else {
    expdt = exp(-dtau);
    w[0]  = 1.0 - expdt;
    w[1]  = w[0] - dtau*expdt;
    w[2]  = 2.0*w[1] - SQ(dtau) * expdt;
  }
}
/* ------- end ---------------------------- w3.c -------------------- */
