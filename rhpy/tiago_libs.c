/* 

Tiago's set of programs to interface RH and Python

get_nstar: get hydrogen levels in LTE using LTEpops(Atom). 

--Tiago, 20101117

*/

 
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../rh.h"
#include "../atom.h"
#include "../atmos.h"
#include "../rhsc2d/geometry.h"
#include "../spectrum.h"
#include "../statistics.h"
#include "../inputs.h"
#include "../constant.h"


/* --- Global variables --                             -------------- */

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_MESSAGE_LENGTH];
enum Topology topology = ONE_D_PLANE;


/* ------- begin ---------------------------------------------------- */


double **get_nstar(char *hatom_file, char *pf_file, char *abund_file,
		   int nspace, double *temp, double *nne, double *toth)
{
/*
  Calculates the Hydrogen LTE populations, for a given atom file (typically you want
  this to be 6-level). Expects location of H atom file, [Kurucz] partition functions,
  abundances, plus number of points, temperature [K], electron density [m^-3],
  and total number of hydrogen atoms [m^-3].

 */

  Atom *atom;
  int i;

  /* Allocate space for atom */
  atom = (Atom *) malloc(sizeof(Atom));

  /* initialise input  */
  input.Nthreads = 1;
  input.PRD_NmaxIter = 0;
  input.PRD_angle_dep = FALSE;
  input.limit_memory = FALSE;
  input.XRD = FALSE;
  input.isum = 0;
  strcpy(input.pfData,pf_file);
  strcpy(input.abund_input,abund_file);
  input.metallicity = 0.;
  commandline.logfile = stderr; 
  
  /* initialise atmos */
  atmos.Nspace = nspace;
  readAbundance(&atmos);

  atmos.Stokes = FALSE;
  atmos.moving = FALSE;
  atmos.vturb = (double *) calloc(atmos.Nspace, sizeof(double));  
  atmos.T     = temp;
  atmos.ne    = nne;
  atmos.nHtot = toth;

  for (i = 0; i < atmos.Nspace; i++){
    atmos.vturb[i] = 0.0;
  }

  readAtom(atom,hatom_file,(bool_t )FALSE);
  
  /* set LTE pops */
  LTEpops(atom,TRUE);

  /* Free stuff   */
  free(atmos.vturb);

  return atom->nstar;
    
}


/*----------------------------------------------------------------------------*/

void main(void)
{
  double *temp, *nne, *toth, **nstar;
  int i,nspace=5;
  char *hatom_file = "Atoms/H_6.atom";
  char *pf_file    = "Atoms/pf_Kurucz.input";
  char *abund_file = "Atoms/abundance.input";
  
  temp  = (double *) calloc(nspace, sizeof(double));
  nne   = (double *) calloc(nspace, sizeof(double));  
  toth  = (double *) calloc(nspace, sizeof(double));

  /* Manually assigning some test values
  // FALC higher 5 layers:
  atmos.T[0] = 1.00e5;
  atmos.T[1] = 9.56e4;
  atmos.T[2] = 9.08e4;
  atmos.T[3] = 8.39e4;
  atmos.T[4] = 7.59e4;
  atmos.ne[0] = 1.251891E10;
  atmos.ne[1] = 1.304293E10;
  atmos.ne[2] = 1.366348E10;
  atmos.ne[3] = 1.467464E10;
  atmos.ne[4] = 1.603707E10;
  atmos.nHtot[0] = 1.0457E10;
  atmos.nHtot[1] = 1.0940E10;
  atmos.nHtot[2] = 1.1518E10;
  atmos.nHtot[3] = 1.2472E10;
  atmos.nHtot[4] = 1.3782E10;  */

  // FALC lower 5 layers (in LTE):
  temp[0] = 8.220000E+03 ;
  temp[1] = 8.540000E+03 ;
  temp[2] = 8.860000E+03 ;
  temp[3] = 9.140000E+03 ;
  temp[4] = 9.400000E+03 ;
  nne[0]  = 1.041290E+15 ;
  nne[1]  = 1.531806E+15 ;
  nne[2]  = 2.194603E+15 ;
  nne[3]  = 2.952398E+15 ;
  nne[4]  = 3.831726E+15 ;
  toth[0] = 1.28402E17 ;
  toth[1] = 1.29180E17 ;
  toth[2] = 1.3002E17 ;
  toth[3] = 1.3128E17 ;
  toth[4] = 1.3266E17 ;  
  

  nstar = get_nstar(hatom_file,pf_file,abund_file,nspace,temp,nne,toth);


  for (i=0; i < nspace; i++){
    printf("%e  %e  %e  %e  %e  %e\n",
	   nstar[0][i]*CUBE(CM_TO_M),nstar[1][i]*CUBE(CM_TO_M),
	   nstar[2][i]*CUBE(CM_TO_M),nstar[3][i]*CUBE(CM_TO_M),
	   nstar[4][i]*CUBE(CM_TO_M),nstar[5][i]*CUBE(CM_TO_M));
  }

  printf("%d\n",sizeof(int));

}

