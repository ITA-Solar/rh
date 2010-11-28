/* ------- file: -------------------------- shortchar.c -------------

       Version:       rh2.0, 2-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Sat Sep 19 15:00:02 2009 --

       --------------------------                      ----------RH-- */

/* --- Formal 2-D radiative transfer solver using the method of short
       characteristics. The transfer equation is solved in one direction
       along the ray only. Use keyword to_observer=TRUE to solve in the
       upward direction, and to_observer=FALSE to solve in opposite
       direction.

  See: P. B. Kunasz and  L. H. Auer 1988, JQSRT 39, 67-79

       For now it implements fixed and periodic boundary
       conditions in the horizontal direction.

  See: L. H. Auer, P. Fabiani Bendicho and J. Trujillo Bueno 1994,
       A&A, 292, 599-615

       At the top and bottom only a fixed boundary is allowed.


 Note: Rays with to_observer == TRUE go from BOTTOM to TOP.

       --                                              -------------- */

 
#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "inputs.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;


/* ------- begin -------------------------- ShortChar.c ------------- */

void ShortChar(Geometry *geometry, int nspect, int mu,
               bool_t to_observer, double *chi, double *S,
	       double *I, double *Psi)
{
  const char routineName[] = "ShortChar";
  register int l, k, n;

  bool_t  monotonic = TRUE;
  enum    boundval bvalue;
  enum    direction upwind, downwind;
  enum    sweep sweep;
  int     Nx = geometry->Nx, Nz = geometry->Nz, local, dl, lstart,
          Nrow, row_left, row_right, dk, kstart, kend;
  double  mux = geometry->mux[mu], I_uw, *Iboundary, dtau_uw, dtau_dw,
          c1, c2, dS_uw, dS_dw, S_uw, chi_uw, S_dw, chi_dw, w[3],
          Bnu_uw, Bnu, T_uw;
  Stencil  *st;
  LongChar *lc = NULL;

  /* --- In case of FIXED left and right boundary conditions treat
         the left and right columns --                 -------------- */

  if (geometry->hboundary == FIXED) {
    if ((to_observer && mux > 0.0) || (!to_observer && mux < 0.0)) {
      
      /* --- Rays that begin at the LEFT boundary and end at the RIGHT
	     boundary --                               -------------- */

      for (k = 0;  k < Nz;  k++) I[k*Nx] = geometry->Ileft[nspect][k];
      if (Psi) for (k = 0;  k < Nz;  k++) Psi[k*Nx] = 0.0;
    } else if ((to_observer && mux < 0.0) || (!to_observer && mux > 0.0)) {

      /* --- Rays that begin at the RIGHT boundary and end at the
             LEFT boundary --                          -------------- */ 

      for (k = 0;  k < Nz;  k++)
	I[(k+1)*Nx - 1] = geometry->Iright[nspect][k];
      if (Psi) for (k = 0;  k < Nz;  k++) Psi[(k+1)*Nx - 1] = 0.0;
    }
  }
  /* --- Boundary conditions in row at the top or bottom -- --------- */

  if (to_observer) {
    row_left = (Nz-1)*Nx;
    bvalue = geometry->bvalue[BOTTOM];
  } else {
    row_left = 0;
    bvalue = geometry->bvalue[TOP];
  }
  row_right = row_left + Nx;

  switch (bvalue) {
  case IRRADIATED:
    if (to_observer)
      Iboundary = geometry->Ibottom[nspect];
    else
      Iboundary = geometry->Itop[nspect];
    for (l = 0;  l < Nx;  l++)  I[row_left + l] = Iboundary[l];
    break;
  case ZERO:
    for (l = row_left;  l < row_right;  l++) I[l] = 0.0;
    break;
  case THERMALIZED:
    if (to_observer) {
      local = (Nz-1) * Nx;
      for (l = 0;  l < Nx;  l++) {
	st = &geometry->stencil[mu][local];

        T_uw = Quadr(atmos.T, st, UPWIND, monotonic);    
	Planck(1, &T_uw, spectrum.lambda[nspect], &Bnu_uw);
	Planck(1, &atmos.T[local], spectrum.lambda[nspect], &Bnu);

	chi_uw  = Quadr(chi, st, UPWIND, monotonic);
	dtau_uw = 0.5 * (chi_uw + chi[local]) * st->ds[UPWIND];

	I[local] = Bnu - (Bnu_uw - Bnu) / dtau_uw;
        local++;
      }
    } else {
      Error(ERROR_LEVEL_2, routineName, 
	    "Boundary condition THERMALIZED not implemented for TOP");
    }
    break;
  }

  if (Psi) for (l = row_left;  l < row_right;  l++) Psi[l] = 0.0;

  /* --- Finally, go through the grid to solve the non-local transfer
         equation --                                   -------------- */

  if (to_observer) {
    if (mux >= 0) dl = 1;  else  dl= -1;
    sweep  = UP;
    dk     = -1;
    kstart = Nz - 2;
    kend   = 0;
    upwind   = DOWNWIND;
    downwind = UPWIND;
  } else {
    if (mux <= 0) dl = 1;  else  dl= -1;
    sweep  = DOWN;
    dk     = 1;
    kstart = 1;
    kend   = Nz - 1;
    upwind   = UPWIND;
    downwind = DOWNWIND;
  }

  for (k = kstart;  k != kend + dk;  k += dk) {

    if (geometry->hboundary == PERIODIC) {
      lc = &geometry->longchar[mu][sweep*Nz + k];
      lstart = lc->lstart;
      Nrow = Nx;
    } else {
      if (mux == 0.0) {
        lstart = 0;
	Nrow = Nx;
      } else {
	if ((to_observer && mux > 0.0) || (!to_observer && mux < 0.0))
	  lstart = 1;
	else
	  lstart = Nx - 2;

        Nrow = Nx - 1;
      }
    }

    for (n = 0, l = lstart;  n < Nrow;  n++, l += dl) {
      if (geometry->hboundary == PERIODIC)
	local = k*Nx + MODULO(l, Nx);
      else
	local = k*Nx + l;
      st = &geometry->stencil[mu][local];

      chi_uw  = Quadr(chi, st, upwind, monotonic);
      dtau_uw = 0.5 * (chi_uw + chi[local]) * st->ds[upwind];
      S_uw    = Quadr(S, st, upwind, monotonic);
      dS_uw   = (S_uw - S[local]) / dtau_uw;

      if (geometry->hboundary == PERIODIC && l == lc->lstart && lc->Nlc > 0)
	I_uw = SolveLong(lc, local, chi, S, I);
      else
	I_uw = Quadr(I, st, upwind, monotonic);

      if (k == kend || (geometry->hboundary == FIXED && n == Nrow-1)) {
	w2(dtau_uw, w);

        /* --- Piecewise linear at end of ray --       -------------- */

	c1 = (S_uw - S[local]) / dtau_uw;
	I[local]   = (1.0 - w[0])*I_uw + w[0]*S[local] + w[1]*c1;
	if (Psi) Psi[local] = w[0] - w[1]/dtau_uw;

      } else {
	w3(dtau_uw, w);

        /* --- Piecewise quadratic elsewhere --        -------------- */

	chi_dw  = Quadr(chi, st, downwind, monotonic);
	dtau_dw = 0.5 * (chi[local] + chi_dw) * st->ds[downwind];
	S_dw    = Quadr(S, st, downwind, monotonic);
	dS_dw   = (S[local] - S_dw) / dtau_dw;

	c1 = (dS_uw*dtau_dw + dS_dw*dtau_uw);
	c2 = dS_uw - dS_dw;

	I[local] = (1.0 - w[0])*I_uw + w[0]*S[local] +
	  (w[1]*c1 + w[2]*c2) / (dtau_uw + dtau_dw);

        /* --- Try piecewise linear if quadratic gives negative
	       monochromatic intensity --              -------------- */ 

	if (I[local] < 0.0) {
	  c1 = dS_uw;
	  I[local] = (1.0 - w[0])*I_uw + w[0]*S[local] + w[1]*c1;

	  if (Psi) Psi[local] = w[0] - w[1]/dtau_uw;
	} else {
	  if (Psi) {
	    c1 = dtau_uw - dtau_dw;
	    Psi[local] = w[0] + (w[1]*c1 - w[2]) / (dtau_uw * dtau_dw);
	  }
	}
      }
    }
  }
}
/* ------- end ---------------------------- ShortChar.c ------------- */

/* ------- begin -------------------------- Quadr.c ----------------- */

/* --- Evaluates function value with polynomial interpolation. -- --- */

double Quadr(double *values, Stencil *stencil,
	     enum direction direction, bool_t monotonic)
{
  const char routineName[] = "Quadr";
  register int n;

  int    *index;
  double  result, v[3], v_max, v_min, fraction, *coeff; 

  /* --- If monotonic = TRUE the interpolated value is forced to be
         monotonic with the original values.

    See: Auer & Paletou, A&A, 285, 675-686 --          -------------- */

  index    = stencil->index[direction];
  coeff    = stencil->coeff[direction];
  fraction = stencil->fraction[direction];

  result = 0.0;
  switch (stencil->order[direction]) {
  case SC_LINEAR:

    result = (1.0 - fraction)*values[index[0]] + fraction*values[index[1]];
    break;

  case SC_QUADRATIC:

    if (monotonic) {
      for (n = 0;  n < 3;  n++) {
	v[n]    = values[index[n]];
	result += coeff[n]*v[n];
      }
      v_max = MAX(MAX(v[0], v[1]), v[2]);
      v_min = MIN(MIN(v[0], v[1]), v[2]);

      if (result > v_max || result < v_min) {
	if (stencil->triplet[direction] == UPPER_TRIPLET)
	  result =
	    (1.0 - fraction)*values[index[0]] + fraction*values[index[1]];
	else
	  result =
	    (1.0 - fraction)*values[index[1]] + fraction*values[index[2]];
      }
    } else {
      for (n = 0;  n < 3;  n++)
	result += coeff[n] * values[index[n]];
    }
    break;

  default:
    Error(ERROR_LEVEL_2, routineName, "Illegal interpolation order");
  }

  return result;
}
/* ------- end ---------------------------- Quadr.c ----------------- */

/* ------- begin -------------------------- PiecewiseStokes.c ------- */

void PiecewiseStokes(Geometry *geometry, int nspect, int mu,
		     bool_t to_observer,
		     double *chi, double **S, double **I, double *Psi)
{
  /* --- Piecewise integration of the coupled Stokes transfer equations
         in two dimensions. Method is quasi-parabolic DELO method.

    See: - D. E. Rees, G. A. Murphy and C. J. Durrant 1989, ApJ 339,
           1093-1106.

         - H. Socas Navarro, J. Trujillo Bueno and B. Ruiz Cobo 2000,
           "Non-LTE Inversion of Stokes Profiles", ApJ 530, 977.
      --                                               -------------- */

  const char routineName[] = "PiecewiseStokes";
  register int l, lp, k, n, m;

  bool_t  monotonic = TRUE;
  enum    boundval bvalue;
  enum    direction upwind, downwind;
  enum    sweep sweep;
  int     local, dl, lstart, row_left, row_right, dk, kstart, kend;
  double  I_uw[4], dtau_uw, dtau_dw, c1, c2, dS_uw[4], dS_dw[4],
          S_uw[4], chi_uw, S_dw[4], chi_dw, w[3], P[4], Q[4][4],
          **R, K[4][4], K_uw[4][4], Bnu_uw, Bnu, T_uw;
  Stencil  *st;
  LongChar *lc = NULL;

  R = matrix_double(4, 4);

  /* --- Boundary conditions in row at the top or bottom -- --------- */

  if (to_observer) {
    row_left = (geometry->Nz-1) * geometry->Nx;
    bvalue = geometry->bvalue[BOTTOM];
  } else {
    row_left = 0;
    bvalue = geometry->bvalue[TOP];
  }
  row_right = row_left + geometry->Nx;

  switch (bvalue) {
  case IRRADIATED:
    if (to_observer) {
      for (l = 0;  l < geometry->Nx;  l++)
	I[0][row_left + l] = geometry->Ibottom[nspect][l];
    } else {
      for (l = 0;  l < geometry->Nx;  l++)
	I[0][row_left + l] = geometry->Itop[nspect][l];
    }
    break;
  case ZERO:
    for (l = row_left;  l < row_right;  l++) I[0][l] = 0.0;
    break;
  case THERMALIZED:
    if (to_observer) {
      local = (geometry->Nz-1) * geometry->Nx;
      for (l = 0;  l < geometry->Nx;  l++) {
	st = &geometry->stencil[mu][local];

        T_uw = Quadr(atmos.T, st, UPWIND, monotonic);    
	Planck(1, &T_uw, spectrum.lambda[nspect], &Bnu_uw);
	Planck(1, &atmos.T[local], spectrum.lambda[nspect], &Bnu);

	chi_uw  = Quadr(chi, st, UPWIND, monotonic);
	dtau_uw = 0.5 * (chi_uw + chi[local]) * st->ds[UPWIND];

	I[0][local] = Bnu - (Bnu_uw - Bnu) / dtau_uw;
        local++;
      }
    } else {
      Error(ERROR_LEVEL_2, routineName, 
	    "Boundary condition THERMALIZED not implemented for TOP");
    }
    break;

  }
  if (Psi)
    for (l = row_left;  l < row_right;  l++) Psi[l] = 0.0;

  /* --- Assume irradiation is unpolarized --          -------------- */

  for (n = 1;  n < 4;  n++) {
    for (l = row_left;  l < row_right;  l++) I[n][l] = 0.0;
  } 

  if (to_observer) {
    if (geometry->mux[mu] >= 0) dl = 1;  else  dl = -1;
    sweep  = UP;
    dk     = -1;
    kstart = geometry->Nz - 2;
    kend   = 0;
    upwind   = DOWNWIND;
    downwind = UPWIND;
  } else {
    if (geometry->mux[mu] <= 0) dl = 1;  else  dl= -1;
    sweep  = DOWN;
    dk     = 1;
    kstart = 1;
    kend   = geometry->Nz - 1;
    upwind   = UPWIND;
    downwind = DOWNWIND;
  }
  /* --- Finally, go through the grid to solve the non-local transfer
         equation --                                   -------------- */

  for (k = kstart;  k != kend + dk;  k += dk) {
    lc = &geometry->longchar[mu][sweep*geometry->Nz + k];
    lstart = lc->lstart;
 
    for (lp = 0, l = lstart;  lp < geometry->Nx;  lp++, l += dl) {
      local = k*geometry->Nx + MODULO(l, geometry->Nx);
      st    = &geometry->stencil[mu][local];

      chi_uw  = Quadr(chi, st, upwind, monotonic);
      dtau_uw = 0.5 * (chi_uw + chi[local]) * st->ds[upwind];
      StokesK_2D(nspect, st, upwind, monotonic, chi_uw, K_uw);

      for (n = 0;  n < 4;  n++) {
	S_uw[n]  = Quadr(S[n], st, upwind, monotonic);
	dS_uw[n] = (S_uw[n] - S[n][local]) / dtau_uw;
      }
      if (l == lc->lstart  &&  lc->Nlc > 0)
	SolveLongStokes(nspect, lc, local, chi, S, I, I_uw);
      else {
	for (n = 0;  n < 4;  n++)
	  I_uw[n] = Quadr(I[n], st, upwind, monotonic);
      }
      StokesK(nspect, local, chi[local], K);

      if (k == kend) {
	w2(dtau_uw, w);

        /* --- Linear piecewise at end of ray --       -------------- */

        for (n = 0;  n < 4;  n++) {
	  c1 = (S_uw[n] - S[n][local]) / dtau_uw;
	  P[n] = w[0]*S[n][local] + w[1]*dS_uw[n];
	}
	if (Psi) Psi[local] = w[0] - w[1]/dtau_uw;
      } else {
	w3(dtau_uw, w);

        /* --- Quadratic piecewise elsewhere --        -------------- */

	chi_dw  = Quadr(chi, st, downwind, monotonic);
	dtau_dw = 0.5 * (chi[local] + chi_dw) * st->ds[downwind];

	for (n = 0;  n < 4;  n++) {
	  S_dw[n]  = Quadr(S[n], st, downwind, monotonic);
	  dS_dw[n] = (S[n][local] - S_dw[n]) / dtau_dw;

	  c1 = dS_uw[n]*dtau_dw + dS_dw[n]*dtau_uw;
	  c2 = dS_uw[n] - dS_dw[n];

          P[n] = w[0]*S[n][local] + (w[1]*c1 + w[2]*c2) /
	    (dtau_uw + dtau_dw);
	}
	if (Psi) {
	  c1 = dtau_uw - dtau_dw;
	  Psi[local] = w[0] + (w[1]*c1 - w[2]) / (dtau_uw * dtau_dw);
	}
      }

      for (n = 0;  n < 4;  n++) {
	for (m = 0;  m < 4;  m++) {
	  Q[n][m] = -w[1]/dtau_uw * K_uw[n][m];
	  R[n][m] = (w[0] - w[1]/dtau_uw) * K[n][m];
	}
	Q[n][n] = 1.0 - w[0];
	R[n][n] = 1.0;
      }
      for (n = 0;  n < 4;  n++) {
	for (m = 0;  m < 4;  m++) 
	  P[n] += Q[n][m] * I_uw[m];
      }
      /* --- Solve linear equations for I --           -------------- */
      
      SolveLinearEq(4, R, P, TRUE);
      
      /* --- Store results for Stokes vector --        -------------- */
      
      for (n = 0;  n < 4;  n++) I[n][local] = P[n];
    }
  }
  freeMatrix((void **) R);
}
/* ------- end ---------------------------- PiecewiseStokes.c ------- */

/* ------- begin -------------------------- StokesK_2D.c ------------ */

void StokesK_2D(int nspect, Stencil *st, enum direction direction,
		bool_t monotonic, double chi_I, double K[4][4])
{
  register int i, j;

  ActiveSet *as;

  /* --- Return the elements of the reduced 4x4 Stokes
         opacity matrix K', which is defined as:

           =           =         =              =
           K' = (chi_c*1 + chi_l*Phi) / chi_I - 1,

	   for wavelength# nspect, spatial point k, and ray mu.

         This is the 2-D short characteristics version.

    See: Rees, Murphy, & Durrant, 1989, ApJ 339, 1093-1106.

         --                                            -------------- */
  as = &spectrum.as[nspect];

  for (j = 0;  j < 4;  j++)
    for (i = 0;  i < 4;  i++) K[j][i] = 0.0;

  /* --- First the contribution from the active set -- -------------- */

  if (containsPolarized(as)) {
    K[0][1] = Quadr(as->chi + atmos.Nspace,
                    st, direction, monotonic);
    K[0][2] = Quadr(as->chi + 2*atmos.Nspace,
                    st, direction, monotonic);
    K[0][3] = Quadr(as->chi + 3*atmos.Nspace,
                    st, direction, monotonic);

    if (input.magneto_optical) {
      K[1][2] = Quadr(as->chip + 2*atmos.Nspace,
                      st, direction, monotonic);
      K[1][3] = Quadr(as->chip + atmos.Nspace,
                      st, direction, monotonic);
      K[2][3] = Quadr(as->chip,
                      st, direction, monotonic);
    }
  }
  /* --- Add possible contribution from the background -- ----------- */

  if (atmos.backgrflags[nspect].ispolarized) {
    K[0][1] += Quadr(as->chi_c + atmos.Nspace,
                     st, direction, monotonic);
    K[0][2] += Quadr(as->chi_c + 2*atmos.Nspace,
                     st, direction, monotonic);
    K[0][3] += Quadr(as->chi_c + 3*atmos.Nspace,
                     st, direction, monotonic);

    if (input.magneto_optical) {
      K[1][2] += Quadr(as->chip_c + 2*atmos.Nspace,
                      st, direction, monotonic);
      K[1][3] += Quadr(as->chip_c + atmos.Nspace,
                      st, direction, monotonic);
      K[2][3] += Quadr(as->chip_c,
                      st, direction, monotonic);
    }
  }
  /* --- Divide by Stokes I opacity and fill lower diagonal part -- - */

  for (j = 0;  j < 3;  j++) {
    for (i = j+1;  i < 4;  i++) {
      K[j][i] /= chi_I;
      K[i][j]  = K[j][i];
    }
  }
  /* --- Anti-symmetric magneto-optical elements --    -------------- */

  if (input.magneto_optical) {
    K[1][3] *= -1.0;
    K[2][1] *= -1.0;
    K[3][2] *= -1.0;
  }
}
/* ------- end ---------------------------- StokesK_2D.c ------------ */
