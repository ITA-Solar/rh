/* ------- file: -------------------------- voigt_cmplx.c ----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
                      C99 _Complex rewrite: M. Szydlarski
       Last modified: 2026 --

       --------------------------                      ----------RH-- */

/* --- Humlicek (1982) and Hui, Armstrong & Wray (1978) Voigt/Faraday-
       Voigt generators, expressed with C99 native complex arithmetic
       (double _Complex).  Mathematically identical to the region logic
       in humlicek.c and the rational approximation in voigt.c, but the
       complex operations are inlined by the compiler instead of going
       through the struct-by-value helpers in complex.c.  This recovers
       the speed of the optional FORTRAN path (hui_.f90/humlicek_.f90)
       with no FORTRAN dependency, so HAVE_F90 is no longer needed.

       Output (both routines):
         Voigt:            H(a, v) = Re[W(v + ia)]
         Faraday-Voigt:  2*F(a, v) = Im[W(v + ia)]

    See: Humlicek 1982, JQSRT 27, p. 437
         Hui, Armstrong & Wray 1978, JQSRT 19, pp. 509-516
         --                                            -------------- */

#include <math.h>
#include <complex.h>
#include <stddef.h>

/* --- This translation unit MUST see the C99 <complex.h>, not RH's own
       struct "complex.h".  The system header defines the macro I; RH's
       does not.  If a stray -I<rh-src-dir> ever shadows it, fail loudly
       here rather than miscompile. --                 -------------- */
#ifndef I
  #error "voigt_cmplx.c: system <complex.h> was shadowed (RH complex.h?). \
Do not add -I. / -I<rh source dir> to CFLAGS for this file."
#endif


/* ------- begin -------------------------- VoigtHumlicek.c --------- */

double VoigtHumlicek(double a, double v, double *F)
{
  double _Complex z = a - v * I;
  double _Complex W;
  double s = fabs(v) + a;

  if (s >= 15.0) {
    /* --- Approximation in region I --                -------------- */
    W = (0.5641896 * z) / (0.5 + z * z);
  } else if (s >= 5.5) {
    /* --- Approximation in region II --               -------------- */
    double _Complex u = z * z;
    W = (z * (1.410474 + u * 0.5641896)) / (0.75 + u * (3.0 + u));
  } else if (a >= 0.195 * fabs(v) - 0.176) {
    /* --- Approximation in region III --              -------------- */
    W = (16.4955 + z*(20.20933 + z*(11.96482 + z*(3.778987 + 0.5642236*z)))) /
        (16.4955 + z*(38.82363 + z*(39.27121 + z*(21.69274 +
                                                  z*(6.699398 + z)))));
  } else {
    /* --- Approximation in region IV --               -------------- */
    double _Complex u = z * z;
    W = cexp(u) -
        (z*(36183.31 - u*(3321.99 - u*(1540.787 - u*(219.031 -
             u*(35.7668 - u*(1.320522 - u*0.56419))))))) /
        (32066.6 - u*(24322.84 - u*(9022.228 - u*(2186.181 -
             u*(364.2191 - u*(61.57037 - u*(1.841439 - u)))))));
  }

  if (F != NULL) *F = cimag(W);
  return creal(W);
}
/* ------- end ---------------------------- VoigtHumlicek.c --------- */

/* ------- begin -------------------------- VoigtHui.c -------------- */

double VoigtHui(double a, double v, double *F)
{
  static const double ah[7] =
        {122.607931777104326, 214.382388694706425, 181.928533092181549,
          93.155580458138441,  30.180142196210589,   5.912626209773153,
           0.564189583562615};
  static const double bh[7] =
        {122.607931773875350, 352.730625110963558, 457.334478783897737,
         348.703917719495792, 170.354001821091472,  53.992906912940207,
          10.479857114260399};
  double _Complex z = a - v * I;
  double _Complex W1 = 0.0, W2 = z;
  double _Complex W;
  int n;

  /* --- Rational approximation of the complex error function W(a - iv),
         evaluated by Horner's scheme in the complex variable z. This
         routine is faster for large a (> 1.5). --      -------------- */

  for (n = 6; n >= 0; n--) {
    W1 = (W1 + ah[n]) * z;
    W2 = (W2 + bh[n]) * z;
  }
  W = W1 / W2;

  if (F != NULL) *F = cimag(W);
  return creal(W);
}
/* ------- end ---------------------------- VoigtHui.c -------------- */
