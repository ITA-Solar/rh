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

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "geometry.h"
#include "inputs.h"
#include "error.h"
#include "parallel.h"
#include "io.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];
extern Input_Atmos_file infile;
extern IO_data   io;
extern IO_buffer iobuf;
extern MPI_data  mpi;

/* ------- begin -------------------------- init_ncdf_J.c ----------- */
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
  size_t  len_id;
  char    file_J[MAX_MESSAGE_LENGTH],  *atmosID;


  /* add the .ncdf extension to J.dat */
  sprintf(file_J, J_FILE_TEMPLATE, input.JFile);

  
  if (input.startJ == OLD_J) {

    /* Read existing file */
    if ((ierror = nc_open_par(file_J, NC_WRITE | NC_MPIPOSIX, mpi.comm, mpi.info,
    //if ((ierror = nc_open_par(file_J, NC_WRITE | NC_MPIIO, mpi.comm, mpi.info,
			      &jncid))) ERR(ierror,routineName); 

    /* Get the variable ids */
    if ((ierror = nc_inq_varid(jncid, "Jlambda", &Jlambda_var))) 
      ERR(ierror,routineName);
    if (input.backgr_pol)
      if ((ierror = nc_inq_varid(jncid, "J20", &J20_var))) 
	ERR(ierror,routineName);    

    /* consistency checks */
    if ((ierror = nc_inq_attlen(jncid, NC_GLOBAL, "atmosID", &len_id ))) 
      ERR(ierror,routineName);

    atmosID = (char *) malloc(len_id+1);

    if ((ierror = nc_get_att_text( jncid, NC_GLOBAL, "atmosID", atmosID ))) 
      ERR(ierror,routineName);

    if (!strstr(atmosID, atmos.ID)) {
      sprintf(messageStr,
	      "J data was derived from different atmosphere (%s) than current",
	      atmosID);
      Error(WARNING, routineName, messageStr);
    }
    free(atmosID);

  } else {  

    /* Create the file, when NEW_J is used */
    if ((ierror = nc_create_par(file_J, NC_NETCDF4 | NC_CLOBBER | NC_MPIPOSIX, 
				mpi.comm, mpi.info, &jncid))) ERR(ierror,routineName);

    /* Write atmos.ID as global attribute */
    if ((ierror = nc_put_att_text(jncid, NC_GLOBAL, "atmosID", strlen(atmos.ID),
				  atmos.ID ))) ERR(ierror,routineName);

    /* Create dimensions */ 
    if ((ierror = nc_def_dim(jncid, "nwave", spectrum.Nspect, &nspect_id ))) 
      ERR(ierror,routineName);
    if ((ierror = nc_def_dim(jncid, "nx",    mpi.nx,          &nx_id     ))) 
      ERR(ierror,routineName);
    if ((ierror = nc_def_dim(jncid, "ny",    mpi.ny,          &ny_id     ))) 
      ERR(ierror,routineName);
    // in the future, this nz should be a maximum value of atmos.Nspace to account
    // for the cases where more points are added (interpolated) to different columns
    if ((ierror = nc_def_dim(jncid, "nz",    infile.nz,       &nspace_id ))) 
      ERR(ierror,routineName);


    /* Create variables */
    dimids[0] = nx_id;
    dimids[1] = ny_id;
    dimids[2] = nspect_id;
    dimids[3] = nspace_id;

    /* Jlambda */
    if ((ierror = nc_def_var(jncid, "Jlambda", NC_FLOAT, 4, dimids, &Jlambda_var)))
      ERR(ierror,routineName);

    /* J20 */
    if (input.backgr_pol)
      if ((ierror = nc_def_var(jncid, "J20", NC_FLOAT, 4, dimids, &J20_var    )))
	ERR(ierror,routineName);

    /* End define mode */
    if ((ierror = nc_enddef(jncid))) ERR(ierror,routineName);
  }

  /* Copy stuff to Background data struct */
  io.j_ncid        = jncid;
  io.j_jlambda_var = Jlambda_var;
  io.j_j20_var     = J20_var;

  return;
}

/* ------- end ---------------------------- init_ncdf_J.c ----------- */

/* ------- begin -------------------------- writeJ_all.c ------------ */
void writeJ_all(void) {
/* Writes J for all the tasks */

  const char routineName[] = "writeJ_p";
  int       ierror, task;
  long      ind;
  size_t    start[] = {0, 0, 0, 0};
  size_t    count[] = {1, 1, 1, 1};

  
  start[2] = 0;   count[2] = spectrum.Nspect;

  //ind = 0;

  /* Set collective access for variables
  if ((ierror = nc_var_par_access(io.j_ncid, io.j_jlambda_var, NC_COLLECTIVE)))
    ERR(ierror,routineName);
  */

  for (task = 0; task < mpi.Ntasks; task++) {
    
    /* If there was a crash, no data were written into buffer variables */
    if (mpi.convergence[task] < 0) continue;

    start[0] = mpi.taskmap[task + mpi.my_start][0];  count[0] = 1;
    start[1] = mpi.taskmap[task + mpi.my_start][1];  count[1] = 1;
    //start[3] = mpi.zcut_hist[task];   count[3] = infile.nz - start[3];
    start[3] = 0;   count[3] = infile.nz;
    
    ind = task*spectrum.Nspect*infile.nz;

    //if ((ierror = nc_put_vara_double(io.j_ncid, io.j_jlambda_var, start, count,
    //				     &iobuf.J[ind] ))) ERR(ierror,routineName);
    if ((ierror = nc_put_vara_float(io.j_ncid, io.j_jlambda_var, start, count,
				     &iobuf.J[ind] ))) ERR(ierror,routineName);

    //ind += spectrum.Nspect*count[3];
  }


  return;
}
/* ------- end ---------------------------- writeJ_all.c ------------ */

/* ------- begin -------------------------- writeJ_p.c -------------- */
void writeJ_p(void) {
/* Writes J, J20 for all wavelengths */
  const char routineName[] = "writeJ_p";
  int ierror;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};

  
  /* Get variable */

  start[0] = mpi.ix;  
  start[1] = mpi.iy;
  start[2] = 0;         count[2] = spectrum.Nspect;	
  start[3] = mpi.zcut;  count[3] = atmos.Nspace;


  if ((ierror = nc_put_vara_double(io.j_ncid, io.j_jlambda_var, start, count,
				   spectrum.J[0] ))) ERR(ierror,routineName);

  if (input.backgr_pol) 
    if ((ierror = nc_put_vara_double(io.j_ncid, io.j_j20_var, start, count,
        		         spectrum.J20[0] ))) ERR(ierror,routineName);

  return;
}
/* ------- end   -------------------------- writeJ_p.c -------------- */


/* ------- begin -------------------------- readJ_p.c -------------- */
void readJ_p(void) {
/* Reads J, J20 for all wavelengths */
  const char routineName[] = "readJ_p";
  int ierror;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};

  
  /* Get variable */
  start[0] = mpi.ix;  
  start[1] = mpi.iy;
  start[2] = 0;         count[2] = spectrum.Nspect;
  start[3] = mpi.zcut;  count[3] = atmos.Nspace;


  if ((ierror = nc_get_vara_double(io.j_ncid, io.j_jlambda_var, start, count,
				   spectrum.J[0] ))) ERR(ierror,routineName);

  if (input.backgr_pol) 
    if ((ierror = nc_get_vara_double(io.j_ncid, io.j_j20_var, start, count,
        		         spectrum.J20[0] ))) ERR(ierror,routineName);

  return;
}
/* ------- end   -------------------------- readJ_p.c --------------- */


/* ------- begin -------------------------- close_ncdf_J.c ---------- */
void close_ncdf_J(void)
/* Closes the J netCDF file */
{
  const char routineName[] = "close_ncdf_J";
  int        ierror;
  if ((ierror = nc_close(io.j_ncid))) ERR(ierror,routineName);
  return;
}
/* ------- end ---------------------------- close_ncdf_J.c ---------- */


/* ------- begin -------------------------- writeJlambda_single.c --- */
void writeJlambda_single(int nspect, double *J)
{
  const char routineName[] = "writeJlambda_ncdf";
  int   ierror;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};

  
  /* Get variable */
  start[1] = mpi.ix;  
  start[2] = mpi.iy;
  start[0] = nspect;
  start[3] = mpi.zcut;
  count[3] = atmos.Nspace;

  if ((ierror = nc_put_vara_double(io.j_ncid, io.j_jlambda_var, start, count,
				   J ))) ERR(ierror,routineName);

  return;
}
/* ------- end ---------------------------- writeJlambda_single.c --- */


/* ------- begin -------------------------- writeJ20_single.c ------- */
void writeJ20_single(int nspect, double *J)
{
  const char routineName[] = "writeJ20_ncdf";
  int   ierror;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};

  
  /* Get variable */
  start[1] = mpi.ix;  
  start[2] = mpi.iy;
  start[0] = nspect;
  start[3] = mpi.zcut;
  count[3] = atmos.Nspace;

  if ((ierror = nc_put_vara_double(io.j_ncid, io.j_j20_var, start, count,
				   J ))) ERR(ierror,routineName);

  return;
}
/* ------- end ---------------------------- writeJ20_single.c ------- */



/* ------- begin -------------------------- readJlambda_single.c ---- */
void readJlambda_single(int nspect, double *J)
{
  const char routineName[] = "readJlambda_ncdf";
  int   ierror;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};

  
  /* Get variable */
  start[1] = mpi.ix;
  start[2] = mpi.iy;
  start[0] = nspect;
  start[3] = mpi.zcut;
  count[3] = atmos.Nspace;

  /* read as double, although it is written as float */
  if ((ierror = nc_get_vara_double(io.j_ncid, io.j_jlambda_var, start, count, J )))
    ERR(ierror,routineName);


  return;
}
/* ------- end ---------------------------- readJlambda_single.c ---- */


/* ------- begin -------------------------- readJ20_single.c -------- */
void readJ20_single(int nspect, double *J)
{
  const char routineName[] = "readJ20_ncdf";
  int   ierror;
  size_t start[] = {0, 0, 0, 0};
  size_t count[] = {1, 1, 1, 1};

  
  /* Get variable */
  start[1] = mpi.ix; 
  start[2] = mpi.iy;
  start[0] = nspect;
  start[3] = mpi.zcut;
  count[3] = atmos.Nspace;

  /* read as double, although it is written as float */
  if ((ierror = nc_get_vara_double(io.j_ncid, io.j_j20_var, start, count, J )))
    ERR(ierror,routineName);

  return;
}
/* ------- end ---------------------------- readJ20_single.c ------- */

