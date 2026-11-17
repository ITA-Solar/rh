/* ------- file: -------------------------- fpehandler.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu May 20 11:29:23 2021 --

       --------------------------                      ----------RH-- */

/* --- Trap floating point exceptions on various machines:

       --                                              -------------- */

#include <math.h>
#include <stdio.h>

#include "rh.h"
#include "error.h"

extern char messageStr[];

void SetFPEtraps(void)
{
  /* --- Explicitly do not set traps --                -------------- */
  /*   Tiago: did not set this for 1.5D version (where it should not stop execution!) */
}

/* ------- end ------------------------------------------------------ */
