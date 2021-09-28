/*
    Functions for HDF5 reading of input Non-thermal Collisional Rates
*/

#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

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

#define FAIL -1
#define MULTI_COMMENT_CHAR  "*"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */
extern Atmosphere atmos;
extern char messageStr[];
extern MPI_data mpi;
extern InputData input;
Geometry geometry;
Input_Atmos_file infile;

/* ------- begin -------------------------- readCnt_hdf5_alt  ------- */

// void readCnt_hdf5(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
//                 Input_Atmos_file *infile, Atom *atom) {
void readCnt_hdf5_alt(int xi, int yi, Atom *atom) {

	 /* Reads the variables pops for a given (xi,yi) pair */
    const char routineName[] = "readCnt_hdf5_alt";
    hsize_t     start_rate[] = {0, 0, 0, 0, 0};
    hsize_t     count_rate[] = {1, 1, 1, 1, 1};
    hsize_t    dims_memory[2];
    hid_t      ncid, dataspace_id, memspace_id;
  
    int        ierror;
    register int        ji, i, j, k;
    register int        n, nact;
    double     **Cnt;
  
    ncid = infile.ncid;

    // Atom *atom;
  
    Cnt = matrix_double(SQ(atom->Nlevel), atmos.Nspace);

    printf("\n\n>>> (alt) Including Non-Thermal Collisional Rates for %2s\n",atom->ID);
   
   // /* Memory dataspace */
   //  dims_memory[0] = infile.nz;
   //  if ((memspace_id = H5Screate_simple(1, dims_memory, NULL)) < 0)
   //    HERR(routineName);
    /* read nH, all at once */
    start_rate[0] = input.p15d_nt; count_rate[0] = 1;
    start_rate[1] = 0;             count_rate[1] = SQ(atom->Nlevel);
    start_rate[2] = (size_t) xi;   count_rate[2] = 1;
    start_rate[3] = (size_t) yi;   count_rate[3] = 1;
    start_rate[4] = mpi.zcut;      count_rate[4] = atmos.Nspace;

    if (strcmp(atom->ID,"HE") == 0){
        if (H5LTfind_dataset(ncid, He_CNT_NAME)) {
          if ((infile.He_Cnt_varid = H5Dopen2(ncid, He_CNT_NAME, H5P_DEFAULT)) < 0)
             HERR(routineName);
          } else {
           infile.He_Cnt_varid = -1;
          }
          
        if (infile.He_Cnt_varid != -1){
          
            dataspace_id = H5Dget_space(infile.He_Cnt_varid);

            ierror = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_rate,
                                         NULL, count_rate, NULL);
            
            dims_memory[0] = SQ(atom->Nlevel);
            dims_memory[1] = atmos.Nspace;
            
            memspace_id = H5Screate_simple(2, dims_memory, NULL);
            ierror = H5Dread(infile.He_Cnt_varid, H5T_NATIVE_DOUBLE,
                 memspace_id, dataspace_id, H5P_DEFAULT, Cnt[0]);
            if (( H5Sclose(dataspace_id) ) < 0) HERR(routineName);
            if (( H5Sclose(memspace_id) ) < 0) HERR(routineName);
            // ierror = H5Dclose(infile->He_Cnt_varid)
             
          // k = 45;
          j = 21;
          i = 0;
          ji = j*atom->Nlevel + i;
          // printf("\n\n>>> CIJ therm = %.17f",atom->C[ji][k]);
          // printf("\n\n>>> CIJ non-therm = %.17f",Cnt[ji][k]);
          for (k = 0; k< atmos.Nspace; k++){
              atom->C[ji][k] += Cnt[ji][k];
              }
          // k = 45;
          // printf("\n\n>>> CIJ tot = %.17f\n\n",atom->C[ji][k]);


          // k = 45;
          j = 32;
          i = 21;
          ji = j*atom->Nlevel + i;
          // printf("\n\n>>> CIJ therm = %.17f",atom->C[ji][k]);
          // printf("\n\n>>> CIJ non-therm = %.17f",Cnt[ji][k]);
          for (k = 0; k< atmos.Nspace; k++){
              atom->C[ji][k] += Cnt[ji][k];
              }
          // k = 45;
          // printf("\n\n>>> CIJ tot = %.17f\n\n",atom->C[ji][k]);
         
          // free(Cnt);
      

      }

    }

freeMatrix((void **) Cnt);


}

/* ------- end ---------------------------- readCntin_hdf5_alt  --------- */

