#include <stdio.h>
#include "idl_export.h"
#include <math.h>

#define SQRTPI 1.772453850905516
#define SQRT2 1.414213562373095145
#define cc 2.9979246e+10
#define cc2 8.9875517873681767e+20
#define K2OMP 1.8368511e-13
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

void Hunt(IDL_LONG n, double *array, double value, IDL_LONG *ilow)
{
  unsigned ascend;
  int    ihigh, index, increment;

  ascend = (array[n-1] > array[0]) ? TRUE : FALSE;
  if ((*ilow <= 0)  ||  (*ilow > n-1)) {

    /* --- Input guess not useful here, go to bisection --  --------- */

    *ilow = 0;
    ihigh = n;
  } else {

    /* --- Else hunt up or down to bracket value --    -------------- */

    increment = 1;
    if (((value >= array[*ilow]) ? TRUE : FALSE) == ascend) {
      ihigh = *ilow + increment;
      if (*ilow == n-1) return;

      /* --- Hunt up --                                -------------- */

      while (((value >= array[ihigh]) ? TRUE : FALSE) == ascend) {
        *ilow = ihigh;
        increment += increment;
        ihigh = *ilow + increment;
        if (ihigh >= n) { ihigh = n;  break; }
      }
    } else {
      ihigh = *ilow;
      if (*ilow == 0) return;

      /* --- Hunt down --                              -------------- */

      while (((value <= array[*ilow]) ? TRUE : FALSE) == ascend) {
        ihigh = *ilow;
        increment += increment;
        *ilow = ihigh - increment;
        if (*ilow <= 0) { *ilow = 0;  break; }
      }
    }
  }
  /* --- Bisection algorithm --                        -------------- */

  if (ascend) {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value >= array[index])
        *ilow = index;
      else
        ihigh = index;
    }
  } else {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value <= array[index])
        *ilow = index;
      else
        ihigh = index;
    }
  }
}

double BiLinear(IDL_LONG Na, IDL_LONG i, IDL_LONG j, double *fab,
                IDL_LONG ll, IDL_LONG nl, double *f)
{
  IDL_LONG idx00,idx01,idx10,idx11;


  /* --- Bi-linear interpolation of function f[][] given on 
 *          rectangular grid (a_table, b_table) --        -------------- */

  idx00 = Na*nl*j + Na*ll + i;
  idx10 = idx00+1;
  idx01 = Na*nl*(j+1) + Na*ll + i;
  idx11 = idx01+1;

  return  fab[0] * f[idx00] + fab[1]* f[idx01] + fab[2] * f[idx10] + fab[3] * f[idx11];
}

int uvlines_natural(double *tg1t, double *ne1t, double *d1t, double *vz1t, double *dzt, double *vbroad, double *lines_gofnt, double *gdens, 
                    double *gtemp, double *lines_wvl, double *lambda, double *spec, double *linesint, double *abund, IDL_LONG *iz, double *minintp, 
                    double *grphp, double *xlimp, IDL_LONG contribfunc, IDL_LONG ndep, IDL_LONG nt, IDL_LONG nlines, 
                    IDL_LONG ngdens, IDL_LONG ngtemp, IDL_LONG nlambda, IDL_LONG nabund, IDL_LONG nthreads)
{
   IDL_LONG i,j,k, i1,i2, timeidx, spidx, depthidx, iidx, ret, it, ine;

   double atomic_weights[] = {1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,58,63,64};
   double totbroad, broad, wvl, inten, thermalbroad2, convolwl, em, gofnt, binsize,grph, xlim, wvl1, wvl2, 
          minintensity, fab[4], fa, fb, lgt,lgne;
  
   grph = *grphp;
   xlim = *xlimp;
   minintensity = *minintp;
   
   binsize = lambda[1] - lambda[0];
   broad = 0;
   i1 = contribfunc ? ndep * nt * nlambda : nlambda * nt;
   for (i = 0; i< i1; i++) spec[i] = 0; // initialize spec to 0
   i1 = contribfunc ? ndep * nt * nlines : nlines * nt;
   for (i = 0; i< i1; i++) linesint[i] = 0; // initialize linesint to 0

   #pragma omp parallel for private(j,i, timeidx, depthidx, lgt, lgne, it, ine, fa, fb, fab, em, gofnt, inten, iidx, wvl, thermalbroad2, totbroad, wvl1, wvl2, i1, i2, spidx, k) num_threads(nthreads) schedule(dynamic) 
   for (j = 0; j<ndep*nt; j++)
   {
      if (tg1t[j] < 1e4) continue;
      if (ne1t[j] <= 0) continue;
      if (d1t[j] <= 0) continue;
      timeidx = j / ndep;
      depthidx = j % ndep;
      lgt = log10(tg1t[j]);
      lgne = log10(ne1t[j]);
      Hunt(ngtemp, gtemp, lgt, &it);
      Hunt(ngdens, gdens, lgne, &ine);
      it = MAX(0,MIN(it,ngtemp-2));
      ine = MAX(0,MIN(ine,ngdens-2));
      fa = (gtemp[it+1] - lgt) / (gtemp[it+1] - gtemp[it]);
      fb = (gdens[ine+1] - lgne) / (gdens[ine+1] - gdens[ine]);
      fab[0] = fa*fb; fab[1] = fa*(1.0 -fb); fab[2] = (1.0 -fa)*fb; fab[3] = (1-fa)*(1-fb);
      em = ne1t[j] * d1t[j]/grph;
      if (!contribfunc) em *=dzt[j];
      for (i= 0; i< nlines; i++)
      {
         if (lines_wvl[i] < lambda[0] || lines_wvl[i] > lambda[nlambda-1]) continue;

         gofnt = BiLinear(ngtemp, it, ine, fab, i, nlines, lines_gofnt); 
         iidx = contribfunc ? ndep * nt * i + ndep * timeidx + depthidx :  nt * i + timeidx; 
         inten = gofnt * abund[iz[i]] * em;
         #pragma omp atomic
         linesint[iidx] += inten;
         if (inten < minintensity) continue; // if below minintensity don't add to spectrum
         wvl = (1-vz1t[j]/cc) * lines_wvl[i];
         thermalbroad2 = K2OMP * tg1t[j] /atomic_weights[iz[i]];
         totbroad = wvl*sqrt(thermalbroad2 + vbroad[j]*vbroad[j]/cc2);
         wvl1 = -xlim * SQRT2 * totbroad + wvl;
         wvl2 = xlim * SQRT2 * totbroad + wvl;
         Hunt(nlambda, lambda,wvl1,&i1);
         Hunt(nlambda, lambda,wvl2,&i2);
         i1-=1; i1 = MAX(i1,0);
         i2+=1; i2 = MIN(i2,nlambda);
         spidx = (contribfunc) ? nlambda * ndep * timeidx + nlambda * depthidx  : nlambda * timeidx;
         for (k=i1; k<i2; k++)
         {
            wvl1 = lambda[k] - binsize/2; wvl2 = lambda[k] + binsize/2;
            #pragma omp atomic
            spec[spidx+k] += inten/binsize*.5*(erf((wvl2 - wvl)/totbroad) - erf((wvl1 - wvl)/totbroad));
         }
      }
   }
   ret = 1;
   return ret;
}
   
int uvlinesc(int argc,void* argv[])
/* 
 * Version with IDL portable calling convention.
 *
 * entry:
 *      argc - Must be 27.
 *      argv[0] - tg1t
 *      argv[1] - ne1t
 *      argv[2] - d1t
 *      argv[3] - vz1t
 *      argv[4] - dzt
 *      argv[5] - vbroad
 *      argv[6] - lines_gofnt
 *      argv[7] - gdens
 *      argv[8] - gtemp
 *      argv[9] - lines_wvl
 *      argv[10]- lambda
 *      argv[11]- spec
 *      argv[12]- linesint
 *      argv[13]- abund
 *      argv[14]- iz
 *      argv[15]- minintensity
 *      argv[16]- grph
 *      argv[17]- xlim
 *      argv[18]- contribfunc
 *      argv[19]- ndep
 *      argv[20]- nt
 *      argv[21]- nlines
 *      argv[22]- ngdens
 *      argv[23]- ngtemp
 *      argv[24]- nlamdba
 *      argv[25]- nabund 
 *      argv[26]- nthreads
 *
 *
 * exit:
 *      calculates a spectrum from gofnt and stores it in spectrum
 */
{
  IDL_LONG ret;
  ret = 0;
  if (argc != 27) return ret;

  return uvlines_natural((double *) argv[0], (double *) argv[1], (double *) argv[2], (double *) argv[3], (double *) argv[4], (double *) argv[5], 
                         (double *) argv[6], (double *) argv[7], (double *) argv[8], (double *) argv[9], (double *) argv[10], (double *) argv[11],
                         (double *) argv[12], (double *) argv[13], (IDL_LONG*) argv[14], (double*) argv[15], (double *) argv[16], (double *) argv[17],
                         (IDL_LONG) argv[18], (IDL_LONG) argv[19], (IDL_LONG) argv[20], (IDL_LONG) argv[21], (IDL_LONG) argv[22], (IDL_LONG) argv[23],
                         (IDL_LONG) argv[24], (IDL_LONG) argv[25],(IDL_LONG) argv[26]);
}

