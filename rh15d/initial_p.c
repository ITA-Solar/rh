/* ------- file: -------------------------- initial_p.c -----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Jul 21 12:31:49 2009 --

       --------------------------                      ----------RH-- */

/* --- Reads and/or computes the initial solution (populations and/or
       mean intensity J).

       NetCDF version.

       Possible options:

         LTE_POPULATIONS    -- Assume LTE populations initially.
         ZERO_RADIATION     -- Solve statistical equilibrium with
                               zero radiation field
         OLD_POPULATIONS    -- Read old populations from file
         OLD_POPS_AND_J     -- Read both old populations and J from file
         ESCAPE_PROBABILITY -- Not yet implemented
         OLD_J              -- Use mean intensities from previous solution
                               (Only implemented for wavelength_table).

       --                                              -------------- */


#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "constant.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"

#define IMU_FILE_TEMPLATE "scratch/Imu_p%d.dat"


/* --- Function prototypes --                          -------------- */
void Escape(Atom *atom);

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];
extern MPI_data mpi;
extern enum Topology topology;


/* ------- begin -------------------------- initSolution_alloc.c ---- */
void initSolution_alloc(void) {
  const char routineName[] = "initSolution_p";
  register int nspect, nact,k,kr,mu;
  char    permission[3], file_imu[MAX_MESSAGE_LENGTH];
  int     Nsr, Nplane, index, oflag,la;
  int to_obs,lamuk,sign,ncoef,ilow,Nlamu,lamu;
  Molecule   *molecule;
  Atom       *atom;
  AtomicLine *line;
  long int idx, lc;
  double *lambda,fac,lambda_prv,lambda_gas,lambda_nxt,dl,lag;

  /* --- Allocate space for angle-averaged mean intensity -- -------- */

  if (!input.limit_memory) {
    if (spectrum.J != NULL) freeMatrix((void **) spectrum.J);
    spectrum.J = matrix_double(spectrum.Nspect, atmos.Nspace);

    /* --- If we do background polarization we need space for the
           anisotropy --                               -------------- */

    if (input.backgr_pol){
      if (spectrum.J20 != NULL) freeMatrix((void **) spectrum.J20);
      spectrum.J20 = matrix_double(spectrum.Nspect, atmos.Nspace);
    }
  }

 /* --- For the PRD angle approximation  we need to store J in
     the gas frame,                                   -------- */
  if (input.PRD_angle_dep == PRD_ANGLE_APPROX &&  atmos.NPRDactive > 0) {

    if (spectrum.Jgas != NULL) freeMatrix((void **) spectrum.Jgas);
    spectrum.Jgas  = matrix_double(spectrum.Nspect, atmos.Nspace);
    if (spectrum.v_los != NULL) freeMatrix((void **) spectrum.v_los);
    spectrum.v_los = matrix_double(    atmos.Nrays, atmos.Nspace);

    /* Calculate line of sight velocity */
    for (mu = 0;  mu < atmos.Nrays;  mu++) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	spectrum.v_los[mu][k] = vproject(k, mu); // / vbroad[k];
      }
    }

    /* --- allocate space for gII redistribution function in PRD
       lines                                                     --- */
    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {

      atom = atmos.activeatoms[nact];

      for (kr = 0;  kr < atom->Nline;  kr++) {

	line = &atom->line[kr];

	if (line->PRD) {

	  /*	  line->gII  = (double **) malloc(atmos.Nrays * sizeof(double *));

	  for (k = 0;  k < atmos.Nspace;  k++) {
	    for (la = 0 ;  la < line->Nlambda;  la++) {

	      // second index in line.gII array
	      lak= k * line->Nlambda + line->Nlambda;

	      q_emit = (line->lambda[la] - line->lambda0) * CLIGHT /
		(line->lambda0 * atom->vbroad[k]);


	      if (fabs(q_emit) < PRD_QCORE) {
		q0 = -PRD_QWING;
		qN =  PRD_QWING;
	      } else {
		if (fabs(q_emit) < PRD_QWING) {
		  if (q_emit > 0.0) {
		    q0 = -PRD_QWING;
		    qN = q_emit + PRD_QSPREAD;
		  } else {
		    q0 = q_emit - PRD_QSPREAD;
		    qN = PRD_QWING;
		  }
		} else {
		  q0 = q_emit - PRD_QSPREAD;
		  qN = q_emit + PRD_QSPREAD;
		}
	      }
	      Np = (int) ((qN - q0) / PRD_DQ) + 1;

	      line->gII[lak]= (double *) malloc(Np*sizeof(double));

	    }
	    }*/

	  Nsr=MAX(3*PRD_QWING,2*PRD_QSPREAD)/PRD_DQ+1;
	  if (line->gII != NULL) freeMatrix((void **) line->gII);
	  line->gII  = matrix_double(atmos.Nspace*line->Nlambda,Nsr);
	  line->gII[0][0] = -1.0;
	}
      }
    }

   /* precompute prd_rho interpolation coefficients if requested */
    if (!input.prdh_limit_mem) {

      for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {

	atom = atmos.activeatoms[nact];

	for (kr = 0;  kr < atom->Nline;  kr++) {

	  line = &atom->line[kr];

	  if (line->PRD) {

	    Nlamu = 2*atmos.Nrays * line->Nlambda;

	    if (line->frac != NULL) freeMatrix((void **) line->frac);
	    line->frac = matrix_double(Nlamu, atmos.Nspace);

	    if (line->id0 != NULL) freeMatrix((void **) line->id0);
	    line->id0  = matrix_int(Nlamu, atmos.Nspace);

	    if (line->id1 != NULL) freeMatrix((void **) line->id1);
	    line->id1  = matrix_int(Nlamu, atmos.Nspace);

	    for (la = 0;  la < line->Nlambda;  la++) {
	      for (mu = 0;  mu < atmos.Nrays;  mu++) {
		for (to_obs = 0;  to_obs <= 1;  to_obs++) {
		  sign = (to_obs) ? 1.0 : -1.0;
		  lamu = 2*(atmos.Nrays*la + mu) + to_obs;

		  for (k = 0;  k < atmos.Nspace;  k++) {

		    // wavelength in local rest frame
		    lag=line->lambda[la] * (1.+spectrum.v_los[mu][k]*sign/CLIGHT);

		    if (lag <= line->lambda[0]) {
		      // out of the lambda table, constant extrapolation
		      line->frac[lamu][k]=0.0;
		      line->id0[lamu][k]=0;
		      line->id1[lamu][k]=1;
		    } else if (lag >= line->lambda[line->Nlambda-1] ) {
		      // out of the lambda table, constant extrapolation
		      line->frac[lamu][k]=1.0;
		      line->id0[lamu][k]=line->Nlambda-2;
		      line->id1[lamu][k]=line->Nlambda-1;
		    } else {
		      // Locate index of line->lambda of point directly to the left of lag
		      Locate(line->Nlambda,line->lambda,lag,&ilow);
		      line->frac[lamu][k] = (lag-line->lambda[ilow])/ (line->lambda[ilow+1]-line->lambda[ilow]);
		      line->id0[lamu][k]=ilow;
		      line->id1[lamu][k]=ilow+1;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    /* precompute Jgas interpolation coefficients if requested */
    if (!input.prdh_limit_mem) {

      lambda = spectrum.lambda;

      /* --- keeps track of where to get indices and interpolation
             coefficients in spectrum.iprhh and spectrum.cprdh --- */
      if (spectrum.nc != NULL) free(spectrum.nc);
      spectrum.nc=  (int *) malloc( 2*atmos.Nrays*spectrum.Nspect*atmos.Nspace * sizeof(int));

      for (la = 0;  la < spectrum.Nspect;  la++) {
	for (mu = 0;  mu < atmos.Nrays;  mu++) {
	  for (to_obs = 0;  to_obs <= 1;  to_obs++) {

	    sign = (to_obs) ? 1.0 : -1.0;

	    for (k = 0;  k < atmos.Nspace;  k++) {

	      lamuk = la * (atmos.Nrays*2*atmos.Nspace)
		+ mu     * (2*atmos.Nspace)
		+ to_obs * (atmos.Nspace)
		+ k ;

	      ncoef=0;

	      // previous, current and next wavelength shifted to gas rest frame
	      fac = (1.+spectrum.v_los[mu][k]*sign/CLIGHT);
	      lambda_prv = lambda[ MAX(la-1,0)                 ]*fac;
	      lambda_gas = lambda[ la                          ]*fac;
	      lambda_nxt = lambda[ MIN(la+1,spectrum.Nspect-1) ]*fac;

	      // do lambda_prv and lambda_gas bracket lambda points?
	      if (lambda_prv !=  lambda_gas) {
		dl= lambda_gas - lambda_prv;
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] > lambda_prv && lambda[idx] <= lambda_gas) ncoef=ncoef+1;
		}
	      } else {
		// edge case, use constant extrapolation for lambda[idx]<lambda gas
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] <=  lambda_gas) ncoef=ncoef+1;
		}
	      }

	      // do lambda_gas and lambda_nxt bracket lambda points?
	      if (lambda_gas != lambda_nxt) {
		dl= lambda_nxt - lambda_gas;
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] > lambda_gas && lambda[idx] < lambda_nxt) ncoef=ncoef+1;
		}
	      } else {
		// edge case, use constant extrapolation for lambda[idx]>lambda gas
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] >=  lambda_gas) ncoef=ncoef+1;
		}
	      }

	      /* --- number of point this lambda contributes to is
	 	 computed as a difference --- */
	      if (lamuk == 0) {
		spectrum.nc[lamuk] = ncoef;
	      } else {
		spectrum.nc[lamuk]=spectrum.nc[lamuk-1]+ncoef;
	      }

	    } // k
	  } // to_obs
	} // mu
      } // la

      /* --- now we know the number of interpolation coefficients,
             it's stored in the last element of spectrum.nc,
	     so allocate space                                     --- */

      idx=spectrum.nc[2*atmos.Nrays*spectrum.Nspect*atmos.Nspace-1];

      if (spectrum.iprdh != NULL) free(spectrum.iprdh);
      spectrum.iprdh= (int *)    malloc( idx * sizeof(int   ));

      if (spectrum.cprdh != NULL) free(spectrum.cprdh);
      spectrum.cprdh= (double *) malloc( idx * sizeof(double));

     /* --- Run through all lamuk points again, and now store indices
            to lambda array and the corresponding interpolation
            coefficients                                          --- */
      for (la = 0;  la < spectrum.Nspect;  la++) {
	for (mu = 0;  mu < atmos.Nrays;  mu++) {
	  for (to_obs = 0;  to_obs <= 1;  to_obs++) {

	    sign = (to_obs) ? 1.0 : -1.0;

	    for (k = 0;  k < atmos.Nspace;  k++) {

	      lamuk = la * (atmos.Nrays*2*atmos.Nspace)
		+ mu     * (2*atmos.Nspace)
		+ to_obs * (atmos.Nspace)
		+ k ;

	      // starting index for storage for this lamuk point
	      lc = (lamuk==0) ? 0 : spectrum.nc[lamuk-1];

	      // previous, current and next wavelength shifted to gas rest frame
	      fac = (1.+spectrum.v_los[mu][k]*sign/CLIGHT);
	      lambda_prv = lambda[ MAX(la-1,0)                 ]*fac;
	      lambda_gas = lambda[ la                          ]*fac;
	      lambda_nxt = lambda[ MIN(la+1,spectrum.Nspect-1) ]*fac;

	      // do lambda_prv and lambda_gas bracket lambda points?
	      if (lambda_prv !=  lambda_gas) {
		dl= lambda_gas - lambda_prv;
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] > lambda_prv && lambda[idx] <= lambda_gas) {
		    // bracketed point found
		    spectrum.iprdh[lc]=idx;
		    spectrum.cprdh[lc]=(lambda[idx]-lambda_prv)/dl;
		    lc++;
		  }
		}
	      } else {
		// edge case, use constant extrapolation for lambda[idx]<lambda gas
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] <=  lambda_gas)  {
		    spectrum.iprdh[lc]=idx;
		    spectrum.cprdh[lc]=1.0;
		    lc++;
		  }
		}
	      }

	      // do lambda_gas and lambda_nxt bracket lambda points?
	      if (lambda_gas != lambda_nxt) {
		dl= lambda_nxt - lambda_gas;
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] > lambda_gas && lambda[idx] < lambda_nxt) {
		    // bracketed point found
		    spectrum.iprdh[lc]=idx;
		    spectrum.cprdh[lc]=1.0 - (lambda[idx]-lambda_gas)/dl;
		    lc++;
		  }
		}
	      } else {
		// edge case, use constant extrapolation for lambda[idx]>lambda gas
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] >=  lambda_gas)  {
		    spectrum.iprdh[lc]=idx;
		    spectrum.cprdh[lc]=1.0;
		    lc++;
		  }
		}
	      }

	    } // k
	  } // to_obs
	} // mu
      } // la

    }  //input.prdh_limit_mem if switch
  } // PRD_ANGLE_APPROX if switch


  /* --- Allocate space for the emergent intensity --  -------------- */

  if (atmos.Stokes || input.backgr_pol) {
    if (spectrum.Stokes_Q != NULL) {
      freeMatrix((void **) spectrum.Stokes_Q);
      spectrum.Stokes_Q = NULL;
    }
    if (spectrum.Stokes_U != NULL) {
      freeMatrix((void **) spectrum.Stokes_U);
      spectrum.Stokes_U = NULL;
    }
    if (spectrum.Stokes_V != NULL) {
      freeMatrix((void **) spectrum.Stokes_V);
      spectrum.Stokes_V = NULL;
    }
  }

  switch (topology) {
  case ONE_D_PLANE:
    if (spectrum.I != NULL) freeMatrix((void **) spectrum.I);
    spectrum.I = matrix_double(spectrum.Nspect, atmos.Nrays);
    if (atmos.Stokes || input.backgr_pol) {
      spectrum.Stokes_Q = matrix_double(spectrum.Nspect, atmos.Nrays);
      spectrum.Stokes_U = matrix_double(spectrum.Nspect, atmos.Nrays);
      spectrum.Stokes_V = matrix_double(spectrum.Nspect, atmos.Nrays);
    }
    break;
  case TWO_D_PLANE:
    Nsr = spectrum.Nspect * atmos.Nrays;
    spectrum.I = matrix_double(Nsr, atmos.N[0]);
    if (atmos.Stokes || input.backgr_pol) {
      spectrum.Stokes_Q = matrix_double(Nsr, atmos.N[0]);
      spectrum.Stokes_U = matrix_double(Nsr, atmos.N[0]);
      spectrum.Stokes_V = matrix_double(Nsr, atmos.N[0]);
    }
    break;
  case THREE_D_PLANE:
    spectrum.I = matrix_double(spectrum.Nspect * atmos.Nrays,
			       atmos.N[0] * atmos.N[1]);
    if (atmos.Stokes || input.backgr_pol) {
      Nsr    = spectrum.Nspect * atmos.Nrays;
      Nplane = atmos.N[0] * atmos.N[1];

      spectrum.I = matrix_double(Nsr, Nplane);
      if (atmos.Stokes || input.backgr_pol) {
	spectrum.Stokes_Q = matrix_double(Nsr, Nplane);
	spectrum.Stokes_U = matrix_double(Nsr, Nplane);
	spectrum.Stokes_V = matrix_double(Nsr, Nplane);
      }
    }
    break;
  case SPHERICAL_SYMMETRIC:
    spectrum.I = matrix_double(spectrum.Nspect, atmos.Nrays);
    if (atmos.Stokes) {
      Error(ERROR_LEVEL_2, routineName,
	    "Cannot do a full Stokes solution in spherical geometry");
    }
    break;
  default:
    sprintf(messageStr, "Unknown topology (%d)", topology);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }



  /* Things to be done only for the first task */
  if (mpi.isfirst) {
    /* --- Need storage for angle-dependent specific intensities for
       angle-dependent PRD --                        -------------- */

    if (atmos.NPRDactive > 0 && input.PRD_angle_dep == PRD_ANGLE_DEP) {
      oflag = 0;
      if (input.startJ == OLD_J) {
	if (spectrum.updateJ) {
	  strcpy(permission, "r+");
	  oflag |= O_RDWR;
	} else {
	  strcpy(permission, "r");
	  oflag |= O_RDONLY;
	}
      } else {
	strcpy(permission, "w+");
	oflag |= (O_RDWR | O_CREAT);
      }
      /* Imu file name, this may not work very well in the mpi version... */
      sprintf(file_imu, IMU_FILE_TEMPLATE, mpi.rank);

      if ((spectrum.fd_Imu = open(file_imu, oflag, PERMISSIONS)) == -1) {
	sprintf(messageStr, "Unable to open %s file %s with permission %s",
		(spectrum.updateJ) ? "update" : "input",
		file_imu, permission);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      /* --- Fill the index list that keeps track of the location
	 of intensity Imu in file spectrum.fd_Imu at wavelength
	 corresponding to nspect. --                 -------------- */

      spectrum.PRDindex = (int *) malloc(spectrum.Nspect * sizeof(int));
      index = 0;
      for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
	if (containsPRDline(&spectrum.as[nspect])) {
	  spectrum.PRDindex[nspect] = index++;
	}
      }
    }


    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      atom = atmos.activeatoms[nact];

      /* --- Allocate memory for the rate equation matrix -- ---------- */
      atom->Gamma = matrix_double(SQ(atom->Nlevel), atmos.Nspace);
    }


    for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
      molecule = atmos.activemols[nact];

      /* --- Allocate memory for the rate equation matrix -- ---------- */
      molecule->Gamma = matrix_double(SQ(molecule->Nv), atmos.Nspace);
    }

  } /* End of first task condition */



  return;
}
/* ------- end   -------------------------- initSolution_alloc.c ---- */

/* ------- begin -------------------------- initSolution.c ---------- */

void initSolution_p(void)
{
  const char routineName[] = "initSolution_p";
  register int k, i, ij, nspect, n, nact;
  int     la, j, status;
  double  gijk, wla, twohnu3_c2, hc_k, twoc, fourPI;
  ActiveSet  *as;
  Molecule   *molecule;
  Atom       *atom;
  AtomicLine *line;
  AtomicContinuum *continuum;

  getCPU(2, TIME_START, NULL);

  /* allocate memory always (because of dynamic Nspace) */
  initSolution_alloc();

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    /* --- Initialize the mutex lock for the operator Gamma if there
           are more than one threads --                -------------- */

    if (input.Nthreads > 0) {
      if ((status = pthread_mutex_init(&atom->Gamma_lock, NULL))) {
	sprintf(messageStr, "Unable to initialize mutex_lock, status = %d",
		status);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }


    switch(atom->initial_solution) {
    case LTE_POPULATIONS:
      for (i = 0;  i < atom->Nlevel;  i++) {
	for (k = 0;  k < atmos.Nspace;  k++)
	  atom->n[i][k] = atom->nstar[i][k];
      }
      break;

    case ZERO_RADIATION:
      hc_k   = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M);
      twoc   = 2.0*CLIGHT / CUBE(NM_TO_M);
      fourPI = 4.0 * PI;

      initGammaAtom(atom, 1.0);

      /* --- Then add radiative contributions of active transitions --  */

      for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
	as = spectrum.as + nspect;

	for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {
	  switch (as->art[nact][n].type) {
	  case ATOMIC_LINE:
	    line = as->art[nact][n].ptype.line;
	    la = nspect - line->Nblue;
	    i  = line->i;
	    j  = line->j;
	    ij = i*atom->Nlevel + j;

	    if (la == 0) {
	      for (k = 0;  k < atmos.Nspace;  k++)
		atom->Gamma[ij][k] += line->Aji;
	    }
	    break;

	  case ATOMIC_CONTINUUM:
	    continuum = as->art[nact][n].ptype.continuum;
	    la = nspect - continuum->Nblue;
	    i  = continuum->i;
	    j  = continuum->j;
	    ij = i*atom->Nlevel + j;

	    wla = fourPI * getwlambda_cont(continuum, la) /
	      continuum->lambda[la];
	    twohnu3_c2 = twoc / CUBE(continuum->lambda[la]);
	    for (k = 0;  k < atmos.Nspace;  k++) {
	      gijk = atom->nstar[i][k]/atom->nstar[j][k] *
		exp(-hc_k/(continuum->lambda[la] * atmos.T[k]));
	      atom->Gamma[ij][k] += gijk * twohnu3_c2 *
		continuum->alpha[la]*wla;
	    }
	    break;
	  default:
	    break;
	  }
	}
      }
      /* --- Solve statistical equilibrium equations --  ------------ */

      statEquil(atom, (input.isum == -1) ? 0 : input.isum);
      break;

    case ESCAPE_PROBABILITY:
      for (k=0; k < input.NpescIter; k++) {
	/* set the Gamma rate matrix equal to the collisional rates */
	initGammaAtom(atom, 1.0);

	/* Perform escape probability approximation */
	Escape(atom);

	/* Get populations from statistical equilibrium */
	statEquil(atom, (input.isum == -1) ? 0 : input.isum);
      }

      break;

    case OLD_POPULATIONS:
      readPopulations(atom);
      break;

    default:;
    break;
    }
  }

  /* --- Now the molecules that are active --          -------------- */

  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];

    /* --- Calculate the LTE vibration level populations here. They
           cannot be calculated yet in readMolecule since chemical
           equilibrium has to be established first --  -------------- */

    for (i = 0;  i < molecule->Nv;  i++) {
      for (k = 0;  k < atmos.Nspace;  k++)
	molecule->nvstar[i][k] = molecule->n[k] *
	  molecule->pfv[i][k] / molecule->pf[k];
    }

    /* --- Initialize the mutex lock for the operator Gamma if there
           are more than one thread --                 -------------- */

    if (input.Nthreads > 0) {
      if ((status = pthread_mutex_init(&molecule->Gamma_lock, NULL))) {
	sprintf(messageStr, "Unable to initialize mutex_lock, status = %d",
		status);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }

    switch(molecule->initial_solution) {

    case LTE_POPULATIONS:
      for (i = 0;  i < molecule->Nv;  i++) {
	for (k = 0;  k < atmos.Nspace;  k++)
	  molecule->nv[i][k] = molecule->nvstar[i][k];
      }
      break;

    case OLD_POPULATIONS:
      readMolPops(molecule);
      break;

    default:;
    }

    /* --- Calculate collisions for molecule (must be done here because
           rotation-vibration transitions are dominated by hydrogen and
           H2 collisions for which chemical equilibrium needs to be
           established first --                        -------------- */

    if (strstr(molecule->ID, "CO"))
      COcollisions(molecule);
    else if (strstr(molecule->ID, "H2"))
      H2collisions(molecule);
    else {
      sprintf(messageStr, "Collisions for molecule %s not implemented\n",
	      molecule->ID);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }
}
/* ------- end ---------------------------- initSolution.c ---------- */
