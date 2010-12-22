/* ------- file: -------------------------- readj.c -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jul 23 15:52:59 2009 --

       --------------------------                      ----------RH-- */

/* --- Routines to read and write angle-averaged mean intensity and
       background opacities and emissivities. --       -------------- */
 

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "geometry.h"
#include "inputs.h"
#include "error.h"
#include "parallel.h"

#define J_FILE_TEMPLATE "%s.ncdf"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];
extern MPI_data mpi;
extern BackgroundData bgdat;

/* ------- begin -------------------------- init__ncdf_J.c ----------- */
void init_ncdf_J(void)
/* Creates the netCDF file for J 
   Must also put option to create or not, depending on OLD_J
   This will be the first file that will be asynchronously read and written
   in parallel
*/
{
  const char routineName[] = "init_netcdf_J";
  int     ierror, jncid, nx_id, ny_id, nspect_id, nspace_id,
          Jlambda_var, J20_var, dimids[4];
  FILE   *test;
  char    file_J[MAX_MESSAGE_LENGTH];


  /* add the .ncdf extension to J.dat */
  sprintf(file_J, J_FILE_TEMPLATE, input.JFile);
 
  /* Check if we can open the file */
  if ((test = fopen(file_J, "a")) == NULL) {
    sprintf(messageStr, "Unable to open J output file %s", input.JFile);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(test);
  }

  /* Create the file */
  if ((ierror = nc_create_par(file_J, NC_NETCDF4 | NC_CLOBBER | NC_MPIIO, 
			      mpi.comm, mpi.info, &jncid))) ERR(ierror,routineName);

  /* Write atmos.ID as global attribute */
  if ((ierror = nc_put_att_text(jncid, NC_GLOBAL, "atmosID", strlen(atmos.ID),
				atmos.ID ))) ERR(ierror,routineName);

  /* Create dimensions */ 
  if ((ierror = nc_def_dim(jncid, "nx",    mpi.nx,          &nx_id     ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(jncid, "ny",    mpi.ny,          &ny_id     ))) 
    ERR(ierror,routineName);
  // in the future, this nz should be a maximum value of atmos.Nspace to account
  // for the cases where more points are added (interpolated) to different columns
  if ((ierror = nc_def_dim(jncid, "nz",    atmos.Nspace,    &nspace_id ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(jncid, "nwave", spectrum.Nspect, &nspect_id ))) 
    ERR(ierror,routineName);

  /* Create variables */
  dimids[0] = nx_id;
  dimids[1] = ny_id;
  dimids[2] = nspace_id;
  dimids[3] = nspect_id;

  /* Jlambda */
  if ((ierror = nc_def_var(jncid, "Jlambda", NC_FLOAT, 4, dimids, &Jlambda_var)))
    ERR(ierror,routineName);

  /* J20 */
  if (input.backgr_pol)
      if ((ierror = nc_def_var(jncid, "J20", NC_FLOAT, 4, dimids, &J20_var    )))
	ERR(ierror,routineName);


  /* End define mode */
  if ((ierror = nc_enddef(jncid))) ERR(ierror,routineName);

  /* Copy stuff to Background data struct */
  printf("j_ncid = %d\n", jncid);
  bgdat.j_ncid      = jncid;
  bgdat.jlambda_var = Jlambda_var;
  bgdat.j20_var     = J20_var;

  return;
}

/* ------- end ---------------------------- init_ncdf_J.c ----------- */


/* ------- begin -------------------------- close_ncdf_J.c ---------- */
void close_ncdf_J(void)
/* Closes the J netCDF file */
{
  const char routineName[] = "close_ncdf_J";
  int        ierror;
  if ((ierror = nc_close(bgdat.j_ncid))) ERR(ierror,routineName);
  return;
}
/* ------- end ---------------------------- close_ncdf_J.c ---------- */


/* ------- begin -------------------------- writeJlambda_ncdf.c ----- */
void writeJlambda_ncdf(int nspect, double *J)
{
  const char routineName[] = "writeJlambda_ncdf";
  int   ierror;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};

  
  /* Get variable */
  start[0] = mpi.task; // temporary, to be changed to mpi.ix 
  start[1] = mpi.task;
  start[3] = nspect;
  count[2] = atmos.Nspace;

  //printf("start = %d %d %d %d\n", start[0],start[1],start[2],start[3]);
  //printf("count = %d %d %d %d\n", count[0],count[1],count[2],count[3]);

  if ((ierror = nc_put_vara_double(bgdat.j_ncid, bgdat.jlambda_var, start, count,
				   J ))) ERR(ierror,routineName);

  return;
}
/* ------- end ---------------------------- writeJlambda_ncdf.c ----- */


/* ------- begin -------------------------- writeJ20_ncdf.c --------- */
void writeJ20_ncdf(int nspect, double *J)
{
  const char routineName[] = "writeJ20_ncdf";
  int   ierror;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};

  
  /* Get variable */
  start[0] = mpi.task; // temporary, to be changed to mpi.ix 
  start[1] = mpi.task;
  start[3] = nspect;
  count[2] = atmos.Nspace;

  if ((ierror = nc_put_vara_double(bgdat.j_ncid, bgdat.j20_var, start, count,
				   J ))) ERR(ierror,routineName);

  return;
}
/* ------- end ---------------------------- writeJ20_ncdf.c --------- */



/* ------- begin -------------------------- readJlambda_ncdf.c ------ */
void readJlambda_ncdf(int nspect, double *J)
{
  const char routineName[] = "readJlambda_ncdf";
  int   ierror;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};

  
  /* Get variable */
  start[0] = mpi.task; // temporary, to be changed to mpi.ix 
  start[1] = mpi.task;
  start[3] = nspect;
  count[2] = atmos.Nspace;

  /* read as double, although it is written as float */
  if ((ierror = nc_get_vara_double(bgdat.j_ncid, bgdat.jlambda_var, start, count, J )))
    ERR(ierror,routineName);


  return;
}
/* ------- end ---------------------------- readJlambda_ncdf.c ------ */


/* ------- begin -------------------------- readJ20_ncdf.c ---------- */
void readJ20_ncdf(int nspect, double *J)
{
  const char routineName[] = "readJ20_ncdf";
  int   ierror;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};

  
  /* Get variable */
  start[0] = mpi.task; // temporary, to be changed to mpi.ix 
  start[1] = mpi.task;
  start[3] = nspect;
  count[2] = atmos.Nspace;

  /* read as double, although it is written as float */
  if ((ierror = nc_get_vara_double(bgdat.j_ncid, bgdat.j20_var, start, count, J )))
    ERR(ierror,routineName);

  return;
}
/* ------- end ---------------------------- readJ20_ncdf.c ---------- */



/* ------- begin -------------------------- readJlambda.c ----------- */

void readJlambda(int nspect, double *J)
{
  const char routineName[] = "readJlambda";

  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);
  offset     = recordsize * nspect;

  result &= (pread(spectrum.fd_J, J, recordsize, offset) == recordsize);

  if (!result) {
    sprintf(messageStr,
	    "Error reading file: offset = %lu, recordsize = %zu",
	    offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- readJlambda.c ----------- */

/* ------- begin -------------------------- writeJlambda.c ---------- */

void writeJlambda(int nspect, double *J)
{
  const char routineName[] = "writeJlambda";

  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);
  offset     = recordsize * nspect;

  result &= (pwrite(spectrum.fd_J, J, recordsize, offset) == recordsize);

  if (!result) {
    sprintf(messageStr,
	    "Error writing file: offset = %lu, recordsize = %zu",
	    offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- writeJlambda.c ---------- */

/* ------- begin -------------------------- readJ20lambda.c --------- */

void readJ20lambda(int nspect, double *J20)
{
  const char routineName[] = "readJ20lambda";

  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);
  offset     = recordsize * nspect;

  result &= (pread(spectrum.fd_J20, J20,
		   recordsize, offset) == recordsize);

  if (!result) {
    sprintf(messageStr,
	    "Error reading file: offset = %lu, recordsize = %zu",
	    offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- readJ20lambda.c --------- */

/* ------- begin -------------------------- writeJ20lambda.c -------- */

void writeJ20lambda(int nspect, double *J20)
{
  const char routineName[] = "writeJ20lambda";

  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);
  offset     = recordsize * nspect;

  result &= (pwrite(spectrum.fd_J20, J20,
		    recordsize, offset) == recordsize);

  if (!result) {
    sprintf(messageStr,
	    "Error writing file: offset = %lu, recordsize = %zu",
	    offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- writeJ20lambda.c -------- */

/* ------- begin -------------------------- readImu.c --------------- */

void readImu(int nspect, int mu, bool_t to_obs, double *I)
{
  const char routineName[] = "readImu";

  bool_t result = TRUE;
  int    index;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);

  index  = spectrum.PRDindex[nspect];
  offset = 2*(index*atmos.Nrays + mu) * recordsize;
  if (to_obs)  offset =+ recordsize;

  result &= (pread(spectrum.fd_Imu, I, 
		    recordsize, offset) == recordsize);

 if (!result) {
    sprintf(messageStr,
	    "Error reading file: offset = %lu, recordsize = %zu",
	    offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- readImu.c --------------- */

/* ------- begin -------------------------- writeImu.c -------------- */

void writeImu(int nspect, int mu, bool_t to_obs, double *I)
{
  const char routineName[] = "writeImu";

  bool_t result = TRUE;
  int    index;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);

  index  = spectrum.PRDindex[nspect];
  offset = 2*(index*atmos.Nrays + mu) * recordsize;
  if (to_obs)  offset =+ recordsize;

  result &= (pwrite(spectrum.fd_Imu, I, 
		    recordsize, offset) == recordsize);

 if (!result) {
    sprintf(messageStr,
	    "Error writing file: offset = %lu, recordsize = %zu",
	    offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- writeImu.c -------------- */

/* ------- begin -------------------------- readBackground.c -------- */

void readBackground(int nspect, int mu, bool_t to_obs)
{
  const char routineName[] = "readBackground";

  int    NrecStokes, NskipStokes;
  long   recordno;
  size_t recordsize;
  off_t  offset;
  bool_t result = TRUE;
  ActiveSet *as;

  /* --- Read background opacity, emissivity and scattering opacity
         into ActiveSet structure as. This cannot be done sequentially
         with relative offsets because the file may be read by
         different threads at the same time. Therefore, we use absolute
         offsets here and in writeBackground. --       -------------- */

  as = &spectrum.as[nspect];

  recordsize = atmos.Nspace * sizeof(double);

  if (atmos.moving || atmos.Stokes)
    recordno = atmos.backgrrecno[2*(nspect*atmos.Nrays + mu) + to_obs];
  else
    recordno = atmos.backgrrecno[nspect];
  offset = recordno * recordsize;

  /* --- Read emissivity and opacity for Q, U, and V only if we are
         solving explicitly for all four Stokes quantities, but always
         skip the proper amount of records --           ------------- */

  NrecStokes = 1;
  if (atmos.backgrflags[nspect].ispolarized) {
    NskipStokes = 4;
    if (input.StokesMode == FULL_STOKES) NrecStokes = 4;
  } else
    NskipStokes = 1;

  /* --- Read background opacity --                    -------------- */

  result &= (pread(atmos.fd_background, as->chi_c,
                   NrecStokes * recordsize,
                   offset) == NrecStokes * recordsize);
  offset += NskipStokes * recordsize;

  /* --- Read off-diagonal elements propagation matrix K -- --------- */

  if (atmos.backgrflags[nspect].ispolarized && input.magneto_optical) {
    if (input.StokesMode == FULL_STOKES)
      result &= (pread(atmos.fd_background, as->chip_c, 3*recordsize,
                       offset) == 3*recordsize);
    offset += 3 * recordsize;
  }
  /* --- Read background emissivity --                 -------------- */

  result &= (pread(atmos.fd_background, as->eta_c,
                   NrecStokes * recordsize,
                   offset) == NrecStokes * recordsize);
  offset += NskipStokes * recordsize;

  /* --- Read background scattering opacity --         -------------- */

  result &= (pread(atmos.fd_background,
                   as->sca_c, recordsize, offset) == recordsize);

  /* --- Exit if reading is unsuccessful --            -------------- */

  if (!result) Error(ERROR_LEVEL_2, routineName, "Error reading file");
}
/* ------- end ---------------------------- readBackground.c -------- */

/* ------- begin -------------------------- writeBackground.c ------- */

int writeBackground(int nspect, int mu, bool_t to_obs,
                    double *chi_c, double *eta_c, double *sca_c,
                    double *chip_c)
{
  const char routineName[] = "writeBackground";

  long   recordno, Nwrite, NrecStokes;
  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  /* --- Writes background opacity, emissivity, and scattering 
         opacity. Returns Nwrite, the number of records written,
         with a record length of atmos.Nspace * sizeof(double) -- --- */

  recordsize = atmos.Nspace * sizeof(double);

  if (atmos.moving || atmos.Stokes)
    recordno = atmos.backgrrecno[2*(nspect*atmos.Nrays + mu) + to_obs];
  else
    recordno = atmos.backgrrecno[nspect];
  offset = recordno * recordsize;

  if (atmos.backgrflags[nspect].ispolarized)
    NrecStokes = 4;
  else
    NrecStokes = 1;

  Nwrite = 0;

  /* --- Write background opacity --                   -------------- */

  result &= (pwrite(atmos.fd_background, chi_c, NrecStokes * recordsize,
                    offset) == NrecStokes * recordsize);
  Nwrite += NrecStokes;
  offset += NrecStokes * recordsize;

  /* --- Write off-diagonal elements of propagation matrix K -- ----- */

  if (atmos.backgrflags[nspect].ispolarized && input.magneto_optical) {
    result &= (pwrite(atmos.fd_background, chip_c, 3*recordsize,
                      offset) == 3*recordsize);
    Nwrite += 3;
    offset += 3 * recordsize;
  }
  /* --- Write background emissivity --                -------------- */

  result &= (pwrite(atmos.fd_background, eta_c, NrecStokes * recordsize,
                    offset) == NrecStokes * recordsize);
  Nwrite += NrecStokes;
  offset += NrecStokes * recordsize;

  /* --- Write background scattering opacity --        -------------- */

  result &= (pwrite(atmos.fd_background, sca_c, recordsize,
                    offset) == recordsize);
  Nwrite += 1;

  if (!result) Error(ERROR_LEVEL_2, routineName, "Error writing file");

  /* --- Return the number of written records --       -------------- */

  return Nwrite;
}
/* ------- end ---------------------------- writeBackground.c ------- */

/* ------- begin -------------------------- readProfile.c ----------- */

void readProfile(AtomicLine *line, int lamu, double *phi)
{
  const char routineName[] = "readProfile";

  int    Nrecphi, NrecSkip;
  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  if (line->polarizable && (input.StokesMode > FIELD_FREE)) {
    if (input.magneto_optical)
      NrecSkip = 7;
    else
      NrecSkip = 4;
  } else
    NrecSkip = 1;

  if (line->polarizable && input.StokesMode == FULL_STOKES) {
    if (input.magneto_optical)
      Nrecphi = 7;
    else
      Nrecphi = 4;
  } else
    Nrecphi = 1;
  
  recordsize = Nrecphi * atmos.Nspace * sizeof(double);
  offset     = NrecSkip * atmos.Nspace * sizeof(double) * lamu;

  result &= (pread(line->fd_profile, phi, recordsize, offset) ==
	     recordsize);

  if (!result) Error(ERROR_LEVEL_2, routineName, "Error reading file");
}
/* ------- end ---------------------------- readProfile.c ----------- */

/* ------- begin -------------------------- writeProfile.c ---------- */

void writeProfile(AtomicLine *line, int lamu, double *phi)
{
  const char routineName[] = "writeProfile";

  int    Nrecphi;
  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  if (line->polarizable && (input.StokesMode > FIELD_FREE)) {
    if (input.magneto_optical)
      Nrecphi = 7;
    else
      Nrecphi = 4;
  } else
    Nrecphi = 1;

  recordsize = Nrecphi * atmos.Nspace * sizeof(double);
  offset     = recordsize * lamu;

  result &= (pwrite(line->fd_profile, phi, recordsize, offset) ==
	     recordsize);

  if (!result) Error(ERROR_LEVEL_2, routineName, "Error writing file");
}
/* ------- end ---------------------------- writeProfile.c ---------- */
