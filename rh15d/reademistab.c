/* ------- file: -------------------------- reademistab.c ---------------------
 *  
 *           Version:      rh_ita
 *           Author:       Graham S. Kerr 
 *                         (graham.s.kerr@nasa.gov ; grahamkerr.astro@gmail.com)
 *                                
 *                         First written: Fri May 17th 2019 --
 *        
 * ----------------------------------------------------------------------RH-- */

/* --- Reads the table of emissivities and performs a trilinear interpolation 
 *     to obtain the emissivity in each grid cell for each of the wavelengths 
 *     in the active set. Integrating through height then yeilds the intensity 
 *     that can be passed for a downward directed radiation source.
 *
 *     This was developed for use with flare simulations where the hot, denser, 
 *     corona can potentially photoionise chromospheric radiation, given the 
 *     increase in x-ray, and (extreme-) ultraviolet radiation. 
 *
 *     There are some hardcoded unit conversions (the grid is in cgs, but the
 *     RH variables are not). At the very last step, convert the intensity to 
 *     RH units. 
 *
 *     Can probably be rewritten to make more efficient!
 *
 *     The input table (emiss_grid.dat) in should contain the following: 
 *
 *     nlambda -- number of wavelength points
 *     n_temp  -- number of temperature bins
 *     n_dens  -- number of density bins
 *     lambda  -- the wavelengths for which the emissivities are defined
 *     temper  -- the temperatures for which the emissivities are defined
 *     edens   -- the electron densities for which the emissivities are defined
 *     emiss_grid -- the emissivities [n_temp x n_dens x nlambda]
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "statistics.h"
#include "spectrum.h"
#include "geometry.h"

#define COMMENT_CHAR "#"

/* --- Function prototypes --                          -------------- */
extern double trapezoidal(double *x, double *y, int n);

/* --- Global variables --                             -------------- */
extern InputData input;
extern char messageStr[];

/* ------- begin ---------------------------------------ReadEmisTab.c ----- */
void ReadEmisTab(Atmosphere *atmos, Spectrum *spectrum, Geometry *geometry)
{
  const char routineName[] = "ReadEmisTab";
  int nlambda, n_temp, n_dens;
  double *lambda, *emiss_intpl, *emiss_grid, *temper, *edens, *int_summed, *wavel;
  double p, dlambda;
  double leftedge, rightedge, res, c, tot_int, avg_int, *new_wave, *new_emiss;
  int Npoints;  
  int i, j, k, n;
  int tind, dind, wind; 
  int result;
  FILE *fp, *fptr;
  unsigned int Nspect = spectrum->Nspect;
  unsigned int Ndep = geometry->Ndep;
  unsigned int Nrays = geometry->Nrays;


  printf("\n>>>> The file in input.emistab_file is %s\n",input.emistab_file);

//if ((fp = fopen("../Atmos/emiss_grid.dat", "r")) == NULL) {
 //   sprintf(messageStr, "Unable to open input file %s", "emiss_grid.dat");
  //  Error(ERROR_LEVEL_2, routineName, messageStr);
  //}

if ((fp = fopen(input.emistab_file, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", "emiss_grid.dat");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  printf("\n\n>>> Reading Emistab file and calculating Irradiation\n");
  
  /* Read the file */
  result=fread(&(nlambda), sizeof(int), 1, fp);

  result=fread(&(n_temp), sizeof(int), 1, fp);

  result=fread(&(n_dens), sizeof(int), 1, fp);

  lambda = malloc(sizeof(double)*nlambda);
  result=fread(lambda, sizeof(double), nlambda, fp);
  
  temper = malloc(sizeof(double)*n_temp);
  result=fread(temper, sizeof(double), n_temp, fp);
  
  edens = malloc(sizeof(double)*n_dens);
  result=fread(edens, sizeof(double), n_dens, fp);
 
  emiss_grid = malloc(sizeof(double)*n_temp*n_dens*nlambda);
  result=fread(emiss_grid, sizeof(double), n_temp*n_dens*nlambda, fp);
  
  fclose(fp);

  /* Create the array to hold the intensity to inject */
  emiss_intpl = malloc(sizeof(double)*Nspect*Ndep);
  int_summed = malloc(sizeof(double)*Nspect);
  wavel = malloc(sizeof(double)*Nspect);

  dlambda = lambda[2]-lambda[1];

/* Loop through the height grid, and if the temperature or density is outwith
 * the bounds of temper or edens, then dont do anything.
 * This will hunt for the index in temperature, density, and wavelength space
 * that the atmosphere requires. Loop through each depth point to get the local 
 * temperature and density, then loop through wavelength space. Is it faster to 
 * loop through wavelength first?? Hunt returns the index j, such that the desired
 * value lies between j and j + 1. This can give us the relevant values to be 
 * passed into a trilinear interpolation routine.
 */

/* We want to make sure flux is conserved within the wavelength bins. So, we
 * will:
 *   -- Loop through each wavelength bin i = 0->Nlambda-1
 *   -- Define the cell edges to be the average between the cells
 *         left -> (lambda[i]+lambda[i-1])/2
 *         rifght -> (lambda[i] + lambda[i+1])/2
 *   -- Make a new finer grid that is linearly spaced between the edges. It should have as 
 *      resolution of {delta Wave/4}, where delta Wave is the wavelengh spacing of the
 *      grid that we have defined the emissivities on. That gives N = (edge1-edge2)/res. If 
 *      RH natively has very fine spacing then N = 0 (integer math) and we can skip the 
 *      recasting step
 *   -- Perform a trilinear interpolation for
 *            N > 0, the recasted grid
 *            N < 0, the RH wavelength point
 *      That gives the specific intensity at each point
 *   -- Integrate over the grid to get the total integrated intensity in that bin
 *   -- Divide by the binsize to get the average specific intensity in the RH bin, 
 *      which is injeced at lambda[i]         
*/
  res = dlambda/10.0;
  for (k = 0; k< Ndep; k++){
      if (atmos->T[k] < 50000.00 || atmos->ne[k]/1.0E06 < edens[0] || atmos->ne[k]/1.0E06 > edens[n_dens-1])  {
         for (i = 0; i<Nspect; i++){
             emiss_intpl[i+ k*Nspect]= 0.0;
         }
     } else {
         Hunt(n_temp, temper, atmos->T[k], &tind);
         Hunt(n_dens, edens, atmos->ne[k]/1.0E06, &dind);
         for (i = 0; i<Nspect; i++){
             if (spectrum->lambda[i]*10.0 < lambda[0] || spectrum->lambda[i]*10.0 > lambda[nlambda-1]) {
                 emiss_intpl[i+ k*Nspect]= 0.0;
             } else {
                     //Find the cell edges
                     if (i == 0){
                        rightedge = (spectrum->lambda[i]*10.0 + spectrum->lambda[i+1]*10.0)/2.0;
                        leftedge =  spectrum->lambda[i]*10.0 - (rightedge - spectrum->lambda[i]*10.0); 
                     } else if (i == Nspect-1){
                        leftedge  = (spectrum->lambda[i-1]*10.0 + spectrum->lambda[i]*10.0)/2.0;
                        rightedge = spectrum->lambda[i]*10.0 + (spectrum->lambda[i]*10.0-leftedge);
                     } else {
                        leftedge  = (spectrum->lambda[i-1]*10.0 + spectrum->lambda[i]*10.0)/2.0;
                        rightedge = (spectrum->lambda[i]*10.0 + spectrum->lambda[i+1]*10.0)/2.0;
                     }
                     //Number of points
                     Npoints = (rightedge-leftedge)/res;
                     if (Npoints > 2){
                        //Make a new grid
                        new_wave =  malloc(sizeof(double)*Npoints);
                        new_emiss = malloc(sizeof(double)*Npoints);
                        c = (rightedge - leftedge)/(Npoints - 1);
                        for(j = 0; j < Npoints - 1; ++j)
                           new_wave[j] = leftedge + j*c;
                        new_wave[Npoints-1] = rightedge;
                         
                  
                        for (j = 0; j < Npoints; j++){
                        //Peform the interpolation of wavelength, density and temperature
                            Hunt(nlambda-1, lambda, new_wave[j], &wind);
                            TrilinearInterp(n_temp, n_dens, nlambda,
                                            emiss_grid, temper, edens, lambda, 
                                            tind, dind, wind,
                                            atmos->T[k], atmos->ne[k]/1.0E06, new_wave[j],
                                            &p);
                           // new_emiss[j] = p*ERG_TO_JOULE / CUBE(CM_TO_M) * SQ(new_wave[j]) / CLIGHT / 1.0E10 / 2.0/ PI * fabs(geometry->height[k]-geometry->height[k+1]);
                              new_emiss[j] = p;
                         }
                        //Perform an integration over wavelength to get total intensity in the bin
                        //then divide by wavelength to get avg int
                        tot_int = 0.0;
                        avg_int = 0.0;
                        tot_int = trapezoidal(new_wave, new_emiss, Npoints);
                        avg_int = tot_int/(rightedge-leftedge);
                        emiss_intpl[i+ k*Nspect] = avg_int * ERG_TO_JOULE / CUBE(CM_TO_M) * SQ(spectrum->lambda[i]*10.0) / CLIGHT / 1.0E10 / 2.0/ PI * fabs(geometry->height[k]-geometry->height[k+1]);   
                        free(new_emiss);
                        free(new_wave);
                     } else {
                          Hunt(nlambda-1, lambda, spectrum->lambda[i]*10.0, &wind);
                          TrilinearInterp(n_temp, n_dens, nlambda, 
                                       emiss_grid, temper, edens, lambda, 
                                       tind, dind, wind, 
                                       atmos->T[k], atmos->ne[k]/1.0E06, spectrum->lambda[i]*10.0,
                                       &p);
                  // Add to the array and convert to SI units, and turn into per steradian 
                  // erg/s/cm^3/ang -> J/s/m^2/sr/Hz */ 
                         emiss_intpl[i+ k*Nspect] = p * ERG_TO_JOULE / CUBE(CM_TO_M) * SQ(spectrum->lambda[i]*10.0) / CLIGHT / 1.0E10 / 2.0/ PI * fabs(geometry->height[k]-geometry->height[k+1]) ;
                    }
         }
     }
  }
}


/* Sum the intensities through height */ 
 for (i = 0; i<Nspect; i++){
    int_summed[i] = 0.0;
    }
 for (i = 0; i<Nspect; i++){ 
     wavel[i] = spectrum->lambda[i]*10.0;
     for (k = 0; k< Ndep; k++){
         int_summed[i] += emiss_intpl[i+ k*Nspect];
         //int_summed[i] = 1e10*ERG_TO_JOULE/CUBE(CM_TO_M) * SQ(spectrum->lambda[i]*10.0)/ CLIGHT / 1.0E10 / 2.0/ PI * fabs(geometry->height[k]-geometry->height[k+1]) ;
     }
 } 
 
/* Place the intensity into the Itop variable, to be used in the rest of RH */
  geometry->Itop = matrix_double(Nspect, Nrays);
  for (i = 0; i<Nspect; i++){
      for (j = 0; j<Nrays; j++){
           geometry->Itop[i][j] = int_summed[i] / geometry->muz[j];
      }
  }


/* Write out the intensity and wavelength to a file */
/*  fptr = fopen("int_summed.dat","wb");
  fwrite(int_summed, sizeof(double), Nspect, fptr);
  fwrite(wavel, sizeof(double),Nspect,fptr);
  fclose(fptr);
*/
  /* Free up the allocated arrays */
  /*printf("\n");*/
  free(lambda);
//  free(dlam);
  free(temper);
  free(edens);
  free(emiss_grid);
  free(emiss_intpl);
  free(int_summed);
}
/* ------- end ------------------------------------ ReadEmisTab.c --------- */
