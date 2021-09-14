/*
    High level routines for selection of input population reading functions.
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

void init_popsin(Atmosphere *atmos, Geometry *geometry, Input_Pops_file *infilepop, Atom *atom) {
    const char routineName[] = "init_popsin";
    htri_t file_hdf5_pops;
    bool_t exit_on_EOF;
    
    int Nread = 0;
    char  scaleStr[20], inputLine[MAX_LINE_SIZE];
    
    printf("\n>>> The filename is %s",input.popsin_file);

    file_hdf5_pops = H5Fis_hdf5(input.popsin_file);
    if (file_hdf5_pops < 0) {  /* Read error */
        Error(ERROR_LEVEL_2, routineName,
              "Could not read input file. Check if file exists.\n");
    }
   
    
       init_popsin_hdf5(atmos, geometry, infilepop, atom);
       
}

void readPopsin(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
                Input_Pops_file *infilepop, Atom *atom) {
      
      readPopsin_hdf5(xi, yi, atmos, geometry, infilepop, atom);

}

void closePopsin(Atmosphere *atmos, Geometry *geometry,
                 Input_Pops_file *infilepop, Atom *atom) {
    
      closePopsin_hdf5(atmos, geometry, infilepop, atom);
 
}