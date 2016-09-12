/* ------- file: -------------------------- iterate.c ---------------

       Version:       rh2.0 (1.5D)
       Author:        Tiago Pereira  (tiago.pereira@nasa.gov)
       Last modified: Tue Apr 19 13:28:31 2019 --

       --------------------------                      ----------RH-- */

/* --- Main iteration routine --                       -------------- */

 
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "background.h"
#include "spectrum.h"
#include "error.h"
#include "inputs.h"
#include "constant.h"
#include "parallel.h"

/* --- Function prototypes --                          -------------- */

void loadBackground(int la, int mu, bool_t to_obs);

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin ----------------------------- Pesc ----------------- */
double Pesc(double tau){
    /* Quick estimation of the escape probability */
    
    if (fabs(tau) < 1.01)
        return 1.0;            /* optically-thin lines, Pesc = 1 */
    else {
        if (tau < 0)
            return 1.0/(-tau); /* hack for negative opacities */
        else
            return 1.0/tau;
    }
}
/* ------- end   ----------------------------- Pesc --------------------- */

/* ------- begin ----------------------------- Escape ------------------- */
void Escape(Atom *atom) {
    /*
      
      Calculates radiative rates using the escape probability approximation

      
      Notes: Escape is inserted in an atom loop. Whatever takes place here should be at
             atom-only level. Why? Because different atoms can have different starting solutions.
             Calling Opacity directly is probably overkill, as it does its own loops, and
             works per wavelength, not per atom.
             
              The line opacity has to be per transition, and thus averaged in wavelength.
              For bound-bound the procedure is the following:
                 
              1. Loop over wavelength to calculate opacity at all wavelengths. Keep an
                 array for each transition, accumulate there the opacity * wla for each
                 transition. (this opacity is depth-dependent as well)
                    
              2. Loop over transitions (and not wavelength!), calculate tau by integrating
                 over depth for each transition
                    
              3. Once we have tau for each transition, add to the gamma of each transition
                 (in the same loop) the nrb recipe.
                    
              For bound-free the procedure is the following:
              
              1. Calculate intensity only for the relevant wavelenths.
              
              2. Using intensity, calculate proper radiative rates and add them to rate matrix.
              
              
    Notes for continuum:

          * To see how rates are added, look at fillgamma.addtoRates or fillgamma.addtoGamma
          * For a simplified version of how to calculate intensity and update rates, look
            at formal.c:233 (the angle-independent case). They could be calculated by a single
            call of Formal, but there is too much rubbish in the main routine, and it could
            call opacity and readbackground more than once -- duplicating many tasks.
      
      
    */
    
    
    const char routineName[] = "Escape";
    register int     n, k, kr, i, ij, ji, nspect;
    int              la, j, nt,  mu, nact;
    double         **opa, *chi, *I, *S, *Psi, *Jdag, *J, tau, wlambda, hc_4PI, twohc,
                     twohnu3_c2, wmu, wlamu, Ieff;
    bool_t           initialize, to_obs;
    ActiveSet       *as;
    AtomicLine      *line;
    AtomicContinuum *continuum;
    
    
    /* --- Some useful constants --                        -------------- */
    hc_4PI = HPLANCK * CLIGHT / (4.0 * PI);
    twohc = 2.0*HPLANCK*CLIGHT / CUBE(NM_TO_M);
    

    opa  = matrix_double(atom->Nline, atmos.Nspace);
    Jdag = (double *) malloc(atmos.Nspace * sizeof(double));
    Psi  = (double *) malloc(atmos.Nspace * sizeof(double));
    chi  = (double *) malloc(atmos.Nspace * sizeof(double));
    I    = (double *) malloc(atmos.Nspace * sizeof(double));
    S    = (double *) malloc(atmos.Nspace * sizeof(double));


    nact = atom->activeindex;
    
    /* Calculate opacities */
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) { 
        as = spectrum.as + nspect;
        alloc_as(nspect, FALSE);
        
        nt = nspect % input.Nthreads;

        /* Get line and background opacity */
        Opacity(nspect, 0, to_obs=TRUE, initialize=TRUE); 

	if (input.backgr_in_mem) {
	  loadBackground(nspect, 0, to_obs=TRUE);
	} else {
	  readBackground(nspect, 0, to_obs=TRUE);
	}
        
        /* --- For bound-bound: store opacity in array, per transition -- */
        for (kr = 0; kr < atom->Nline; kr++) {
            line = atom->line + kr;
            la = nspect - line->Nblue;
            
            if (la == 0)
                wlambda = getwlambda_line(line, la);
                
            if ((la >= 0) && (la < line->Nlambda)) {
                /* increment opacities with each wavelength, multiplying by the
                   integration weights */
                for (k = 0 ; k < atmos.Nspace; k++) 
                    opa[kr][k] += (as->chi[k] + as->chi_c[k]) * (wlambda * line->wphi[k] /
                                                                    hc_4PI);
            }
        }
        
        /* --- For bound-free: calculate intensity and update rates ----- */        
        for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {      
            if (as->art[nact][n].type == ATOMIC_CONTINUUM) {
                
                continuum = as->art[nact][n].ptype.continuum;
                la = nspect - continuum->Nblue;
                i = continuum->i;
                j = continuum->j;
                ij = i*atom->Nlevel + j;
                ji = j*atom->Nlevel + i;
                
                twohnu3_c2 = twohc / CUBE(spectrum.lambda[nspect]);
                
                /* Use old J and zero new array */
                if (input.limit_memory) {
                    J = (double *) malloc(atmos.Nspace *sizeof(double));
                    //readJlambda_single(nspect, Jdag);
                } else {
                    J = spectrum.J[nspect];
                    for (k = 0;  k < atmos.Nspace;  k++) Jdag[k] = J[k];
                }
                for (k = 0;  k < atmos.Nspace;  k++) J[k] = 0.0;
                
                for (mu = 0; mu < atmos.Nrays; mu++) {
                    wmu = 0.5 * geometry.wmu[mu];
                        
                    for (k = 0;  k < atmos.Nspace;  k++) {
                        chi[k] = as->chi[k] + as->chi_c[k];
                        S[k]   = (as->eta[k] +
                            as->eta_c[k] + as->sca_c[k]*Jdag[k]) / chi[k];
                    } 
                     
                    /* Get intensity */
                    Piecewise_1D(nspect, mu, to_obs=TRUE, chi, S, I, Psi);
                    
                    /* Update rates */
                    for (k = 0;  k < atmos.Nspace;  k++) {
                        
                        Ieff = I[k] - Psi[k] * atom->rhth[nt].eta[k];
                        wlamu = atom->rhth[nt].Vij[n][k] * atom->rhth[nt].wla[n][k] * wmu;
                        
                        atom->Gamma[ji][k] += Ieff * wlamu;
                        atom->Gamma[ij][k] += (twohnu3_c2 + Ieff) * atom->rhth[nt].gij[n][k] * wlamu;
                        
                        /*  Accumulate mean intensity */
                        J[k] += wmu * I[k];
                    }                    
                }
                if (input.limit_memory) free(J);
            }
        }
        
        free_as(nspect, FALSE);
    }  

    /* Bound-bound: calculate optical depth and add approximation to rates */  
    for (kr = 0; kr < atom->Nline; kr++) {
        line = atom->line + kr;
        i  = line->i;
        j  = line->j;
        ij = i*atom->Nlevel + j;

        tau = 0.0;
        
        for (k = 0; k < atmos.Nspace ; k++) {
            if (k > 0) {
                tau += 0.5*(opa[kr][k-1] + opa[kr][k]) *
                    (geometry.height[k-1] - geometry.height[k]);
            }
            
            /* add escape probability approximation to the rates matrix */
            atom->Gamma[ij][k] += line->Aji * Pesc(tau);   
        }

    }
    
    freeMatrix((void **) opa);
    free(Jdag);
    free(chi);
    free(Psi);
    free(I);
    free(S);
}

/* ------- end   ----------------------------- Escape --------------- */
