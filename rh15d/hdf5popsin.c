/*
    Functions for HDF5 reading of input populations
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
#include "atom.h"

#define FAIL -1
#define MULTI_COMMENT_CHAR  "*"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */
extern MPI_data mpi;
extern InputData input;
extern char messageStr[];

/* ------- begin --------------------------   init_popsin_hdf5   ----- */
void init_popsin_hdf5(Atmosphere *atmos, Geometry *geometry,
    Input_Pops_file *infilepop, Atom *atom) {
  
  const char routineName[] = "init_popsin_hdf5";
  struct  stat statBuffer;
  hid_t plist_id2, p_ncid;
  hsize_t dims[5];
  size_t nn = 0;
  size_t start[] = {0, 0};
  size_t count[] = {1, 1};
  char *filename;

  /* --- Open input file for model pops ---       -------------- */
  if ((plist_id2 = H5Pcreate(H5P_FILE_ACCESS)) < 0) HERR(routineName);
  // if ((H5Pset_fapl_mpio(plist_id2, mpi.comm, mpi.info)) < 0) HERR(routineName);
    printf("\n\n>>>> MADE IT TO HERE");
    printf("\n\n>>> filename = %s",input.popsin_file);
  if ((p_ncid = H5Fopen(input.popsin_file, H5F_ACC_RDONLY, plist_id2)) < 0)
    HERR(routineName);
  infilepop->p_ncid = p_ncid;
  if ((H5Pclose(plist_id2)) < 0) HERR(routineName); /* plist no longer needed */

  /* Get dimensions from hydrogen_populations (if present) */
  if (H5LTfind_dataset(p_ncid, "pops")) {
    if ((H5LTget_dataset_info(p_ncid, "pops", dims,
			      NULL, NULL)) < 0)
      HERR(routineName);
    nn = dims[1];
    infilepop->pnx = dims[2];
    infilepop->pny = dims[3];
    infilepop->pnz = dims[4];
  } else {
    sprintf(messageStr, "Popsfile missing pops"
	    " array, aborting.\n");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  if ((infilepop->popsin_varid = H5Dopen2(p_ncid, POPSIN_NAME, H5P_DEFAULT)) < 0)
    HERR(routineName);


  /* read things that don't depend on x, y */
  start[0] = input.p15d_nt;
  count[0] = 1;
  start[1] = 0;
  count[1] = infilepop->pnz;

}
/* ------- end ---------------------------- init_popsin_hdf5  -------- */


/* ------- begin -------------------------- readPopsin_hdf5  --------- */
void readPopsin_hdf5(int xi, int yi, Atmosphere *atmos, Geometry *geometry,
                Input_Pops_file *infilepop, Atom *atom) {
  /* Reads the variables pops for a given (xi,yi) pair */
  const char routineName[] = "readPopsin_hdf5";
  hsize_t     start_pop[] = {0, 0, 0, 0, 0};
  hsize_t     count_pop[] = {1, 1, 1, 1, 1};
  hsize_t    dims_memory[2];
  hid_t      p_ncid, dataspace_id, memspace_id;
  int        ierror, i, j;
  bool_t     old_moving;

  p_ncid = infilepop->p_ncid;
 
  /* Memory dataspace */
  dims_memory[0] = infilepop->pnz;
  if ((memspace_id = H5Screate_simple(1, dims_memory, NULL)) < 0)
    HERR(routineName);
  
  /* read nH, all at once */
  start_pop[0] = input.p15d_nt; count_pop[0] = 1;
  start_pop[1] = 0;             count_pop[1] = atom->Nlevel;
  start_pop[2] = (size_t) xi;   count_pop[2] = 1;
  start_pop[3] = (size_t) yi;   count_pop[3] = 1;
  start_pop[4] = mpi.zcut;      count_pop[4] = atmos->Nspace;
  dataspace_id = H5Dget_space(infilepop->popsin_varid);
  ierror = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_pop,
                               NULL, count_pop, NULL);
  dims_memory[0] = atom->Nlevel;
  dims_memory[1] = atmos->Nspace;
  memspace_id = H5Screate_simple(2, dims_memory, NULL);
  ierror = H5Dread(infilepop->popsin_varid, H5T_NATIVE_DOUBLE,
		   memspace_id, dataspace_id, H5P_DEFAULT, atom->n[0]);
  if (( H5Sclose(dataspace_id) ) < 0) HERR(routineName);
  if (( H5Sclose(memspace_id) ) < 0) HERR(routineName);
  /* Depth grid refinement */
  // if (input.p15d_refine)
  //   depth_refine(atmos, geometry, input.p15d_tmax);


}
/* ------- end ---------------------------- readPopsin_hdf5  --------- */


/* ------- begin -------------------------- closePopsin_hdf5  ------------- */
void closePopsin_hdf5(Atmosphere *atmos, Geometry *geometry,
                      Input_Pops_file *infilepop, Atom *atom) {
  /* Closes the HDF5 file and frees memory */
  int ierror;
  /* Close the file. */
  ierror = H5Dclose(infilepop->popsin_varid);

  ierror = H5Fclose(infilepop->p_ncid);
  /* Free stuff */
  free(atom->n);
}