/*
    High level routines for selection of input nonthermal collisional rates reading functions.
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
// extern InputData input;
extern char messageStr[];



// void readCntin(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
//                 Input_Atmos_file *infile) {
void readCntin(int xi, int yi) {
      // readCnt_hdf5_alt(xi, yi, atmos, geometry, infile);
      readCnt_hdf5_alt(xi, yi);

}

