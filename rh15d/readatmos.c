/*
    High level routines for selection of atmos reading functions.
*/
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "constant.h"
#include "background.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"
#include "io.h"

#define MULTI_COMMENT_CHAR  "*"

/* --- Global variables --                             -------------- */
extern MPI_data mpi;
extern InputData input;
extern char messageStr[];


void init_atmos(Atmosphere *atmos, Geometry *geometry, Input_Atmos_file *infile) {
    /* Finds out what kind of atmosphere is supplied, calls appropriate routines */
    const char routineName[] = "init_atmos";
    htri_t file_hdf5;
    bool_t exit_on_EOF;
    FILE *fp_atmos;
    int Nread = 0;
    char  scaleStr[20], inputLine[MAX_LINE_SIZE];

    file_hdf5 = H5Fis_hdf5(input.atmos_input);
    if (file_hdf5 > 0) {  /* File is HDF5 */
        geometry->atmos_format = HDF5;
    }
    else if (file_hdf5 < 0) {  /* Read error */
        Error(ERROR_LEVEL_2, routineName,
              "Could not read input file. Check if file exists.\n");
    }
    else {  /* File exists but is not HDF5 */
        if ((fp_atmos = fopen(input.atmos_input, "r")) == NULL) {
          sprintf(messageStr, "Unable to open inputfile %s\n", input.atmos_input);
          Error(ERROR_LEVEL_2, routineName, messageStr);
        } else {
            getLine(fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
            getLine(fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
            Nread += sscanf(inputLine, "%19s", scaleStr);
            fclose(fp_atmos);
            switch(toupper(scaleStr[0])) {
                case 'M':
                  geometry->atmos_format = MULTI;
                  break;
                case 'T':
                  geometry->atmos_format = MULTI;
                  break;
                case 'H':
                  geometry->atmos_format = MULTI;
                  break;
                default:
                  sprintf(messageStr,
                          "Failed to read %s as MULTI atmosphere failed.\n",
                          input.atmos_input);
                  Error(ERROR_LEVEL_2, routineName, messageStr);
            }
        }
    }
    /* Get abundances of background elements */
    readAbundance(atmos);
    switch (geometry->atmos_format) {
        case HDF5:
            init_hdf5_atmos(atmos, geometry, infile);
            break;
        case MULTI:
            /* For 1D MULTI no init necessary, read all in one go */
            readAtmos_multi(atmos, geometry, infile);
            break;
    }
}


void readAtmos(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
                Input_Atmos_file *infile) {
    /* Select read_atmos routine */
    switch (geometry->atmos_format) {
        case HDF5:
            readAtmos_hdf5(xi, yi, atmos, geometry, infile);
            break;
        case MULTI:
            /* Read again for consistency purposes */
            readAtmos_multi(atmos, geometry, infile);
            break;
        default:
            break;
    }
}


void close_atmos(Atmosphere *atmos, Geometry *geometry,
                 Input_Atmos_file *infile) {
    /* Select close_atmos routine */
    switch (geometry->atmos_format) {
        case HDF5:
            close_hdf5_atmos(atmos, geometry, infile);
            break;
        default:
            break;
    }
}


void setTcut(Atmosphere *atmos, Geometry *geometry, double Tmax) {
  /* Find the point where temperature (in TR region) gets below Tmax,
     set atmos.Nspace to remaining points, repoint all the atmospheric
     quantities to start at that point.  */
  const char routineName[] = "setTcut";
  int   i, cutpoint;

  mpi.zcut =  0;
  if (Tmax < 0) return; /* Do nothing when Tmax < 0 */
  cutpoint = -1;
  for (i=0; i < atmos->Nspace; i++) {
    if (atmos->T[i] <= Tmax) {
      cutpoint = i;
      break;
    }
  }
  if (cutpoint < 0) {
    sprintf(messageStr,
	    "\n-Could not find temperature cut point! Aborting.\n");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  mpi.zcut = cutpoint;
  /* subtract from total number of points */
  atmos->Nspace  -= mpi.zcut;
  geometry->Ndep -= mpi.zcut;
  /* Reallocate arrays of Nspace */
  realloc_ndep(atmos, geometry);
}


void realloc_ndep(Atmosphere *atmos, Geometry *geometry) {
  /* Reallocates the arrays of Nspace */
  atmos->T = (double *) realloc(atmos->T, atmos->Nspace * sizeof(double));
  atmos->ne = (double *) realloc(atmos->ne, atmos->Nspace * sizeof(double));
  atmos->nHtot = (double *) realloc(atmos->nHtot,
                                    atmos->Nspace * sizeof(double));
  geometry->vel = (double *) realloc(geometry->vel,
                                     atmos->Nspace * sizeof(double));
  geometry->height = (double *) realloc(geometry->height,
					                    atmos->Nspace * sizeof(double));
  if (atmos->nHmin != NULL)
      atmos->nHmin = (double *) realloc(atmos->nHmin,
                                        atmos->Nspace * sizeof(double));
  /* default zero */
  free(atmos->vturb);
  atmos->vturb  = (double *) calloc(atmos->Nspace , sizeof(double));
  if (atmos->Stokes) {
    atmos->B = (double *) realloc(atmos->B, atmos->Nspace * sizeof(double));
    atmos->gamma_B = (double *) realloc(atmos->gamma_B,
                                        atmos->Nspace * sizeof(double));
    atmos->chi_B = (double *) realloc(atmos->chi_B,
				                      atmos->Nspace * sizeof(double));
  }
}


void depth_refine(Atmosphere *atmos, Geometry *geometry, double Tmax) {
  /* Performs depth refinement to optimise for gradients in temperature,
     density and optical depth. In the same fashion as multi23's ipol_dscal
     (in fact, shamelessly copied).
  */
  bool_t nhm_flag, hunt;
  long    i, k, k0=0, k1=0;
  size_t  bufsize;
  double  CI, PhiHmin, *chi, *eta, *tau, tdiv, rdiv, taudiv, *aind, *xpt;
  double *new_height, *buf;
  const double taumax = 100.0, lg1 = log10(1.1);

  if (Tmax < 0) Tmax = 1.e20; /* No zcut if Tmax is negative */
  chi = (double *) malloc(atmos->Nspace * sizeof(double));
  eta = (double *) malloc(atmos->Nspace * sizeof(double));
  tau = (double *) calloc(atmos->Nspace , sizeof(double));
  xpt = (double *) calloc(atmos->Nspace , sizeof(double));
  buf = (double *) malloc(atmos->Nspace * sizeof(double));
  aind = (double *) calloc(atmos->Nspace , sizeof(double));
  new_height = (double *) calloc(atmos->Nspace , sizeof(double));
  bufsize = atmos->Nspace * sizeof(double);
  nhm_flag = FALSE;
  CI = (HPLANCK/(2.0*PI*M_ELECTRON)) * (HPLANCK/KBOLTZMANN);
  k1 = atmos->Nspace;
  /* --- Calculate tau500 scale from H-bf opacity alone. --- */
  /* First, build nHmin because ChemicalEquilibrium has not been called yet.
     Unless H level pops are given in hdf5 file, all H is assumed neutral.  */
  if (atmos->nHmin == NULL) {
    atmos->nHmin = (double *) malloc(atmos->Nspace * sizeof(double));
    nhm_flag = TRUE;
  }

  for (k = 0; k < atmos->Nspace; k++) {
    PhiHmin = 0.25*pow(CI/atmos->T[k], 1.5) *
        exp(0.754 * EV / (KBOLTZMANN * atmos->T[k]));
    atmos->nHmin[k] = atmos->ne[k] * atmos->nH[0][k] * PhiHmin;
  }

  Hminus_bf(500.0, chi, eta);

  /* integrate for optical depth */
  for (k = 1; k < atmos->Nspace; k++) {
    if ((atmos->T[k] > Tmax) && (k0 == k - 1)) k0 = k;
    tau[k] = tau[k-1] + 0.5 * (chi[k] + chi[k-1]) *
      (fabs(geometry->height[k-1] - geometry->height[k]));
    if (tau[k] < taumax) k1 = k;
    xpt[k] = (double) k;
  }

  tau[0] = tau[1];

  /* --- Compute log variations --- */
  for (k = k0+1; k <= k1; k++) {
    tdiv = fabs(log10(atmos->T[k]) - log10(atmos->T[k-1]))/lg1;
    /* rho is not available, so nH[0] used instead */
    rdiv = fabs(log10(atmos->nH[0][k]) - log10(atmos->nH[0][k-1]))/lg1;
    taudiv = fabs(log10(tau[k]) - log10(tau[k-1]))/0.1;
    aind[k] = aind[k-1] + MAX(MAX(tdiv,rdiv),taudiv);
  }

  for (k = 1; k <= k1; k++)
    aind[k] *= (atmos->Nspace-1)/aind[k1];

  /* --- Create new height scale --- */
  splineCoef(k1-k0+1, &aind[k0], &geometry->height[k0]);
  splineEval(atmos->Nspace, xpt, new_height, hunt=FALSE);

  /* --- Interpolate quantities to new scale --- */
  /* Take logs of densities to avoid interpolation to negatives */
  for (k = 0; k < atmos->Nspace; k++) {
    atmos->ne[k] = log(atmos->ne[k]);
    for (i = 0; i < atmos->NHydr; i++)
      atmos->nH[i][k] = log(atmos->nH[i][k]);
  }

  splineCoef(atmos->Nspace, geometry->height, atmos->T);
  splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
  memcpy((void *) atmos->T, (void *) buf, bufsize);
  /* condition for T < 1000 K ? */

  splineCoef(atmos->Nspace, geometry->height, atmos->ne);
  splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
  memcpy((void *) atmos->ne, (void *) buf, bufsize);

  for (i = 0; i < atmos->NHydr; i++) {
    splineCoef(atmos->Nspace, geometry->height, atmos->nH[i]);
    splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
    memcpy((void *) atmos->nH[i], (void *) buf, bufsize);
  }

  splineCoef(atmos->Nspace, geometry->height, geometry->vel);
  splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
  memcpy((void *) geometry->vel, (void *) buf, bufsize);

  splineCoef(atmos->Nspace, geometry->height, atmos->vturb);
  splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
  memcpy((void *) atmos->vturb, (void *) buf, bufsize);

  if (atmos->Stokes) {
    splineCoef(atmos->Nspace, geometry->height, atmos->B);
    splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
    memcpy((void *) atmos->B, (void *) buf, bufsize);

    splineCoef(atmos->Nspace, geometry->height, atmos->gamma_B);
    splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
    memcpy((void *) atmos->gamma_B, (void *) buf, bufsize);

    splineCoef(atmos->Nspace, geometry->height, atmos->chi_B);
    splineEval(atmos->Nspace, new_height, buf, hunt=FALSE);
    memcpy((void *) atmos->chi_B, (void *) buf, bufsize);
  }

  /* Take back logs */
  for (k = 0; k < atmos->Nspace; k++) {
    atmos->ne[k] = exp(atmos->ne[k]);
    for (i = 0; i < atmos->NHydr; i++)
      atmos->nH[i][k] = exp(atmos->nH[i][k]);
  }

  memcpy((void *) geometry->height, (void *) new_height, bufsize);

  if (nhm_flag) {
    free(atmos->nHmin);
    atmos->nHmin = NULL;
  }
  free(buf); free(new_height);
  free(chi); free(eta);
  free(tau); free(aind); free(xpt);
}
