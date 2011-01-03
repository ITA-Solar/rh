/* ------- file: -------------------------- brs_p.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jan  3 14:35:46 2008 --

       --------------------------                      ----------RH-- */

/* --- Routines to write and read background record structure to file.

       XDR (external data representation) version. --  -------------- */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "error.h"
#include "xdr.h"
#include "geometry.h"
#include "parallel.h"
#include "io.h"

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern char messageStr[];
extern IO_data io; 
extern MPI_data mpi;

/* ------- begin -------------------------- init_ncdf_BRS.c --------- */
void init_ncdf_BRS(void)
/* Initialises the BRS netCDF file */
{
  const char routineName[] = "init_ncdf_BRS";
  int     ierror, ncid, task_id, nspect_id, nrec_id, nspace_id,
          hl_var, ip_var, nrec_var, Nrecno, dimids[2];
  FILE   *test;
  char   *fname = BRS_FILE, *fext  = BRS_EXT;
  char    file_brs[MAX_MESSAGE_LENGTH];

  /* get file name, with the MPI rank */
  sprintf(file_brs,"%s%d%s", fname, mpi.rank, fext);

  if (atmos.moving || atmos.Stokes)
    Nrecno = 2 * spectrum.Nspect * atmos.Nrays;
  else
    Nrecno = spectrum.Nspect;

  /* Check if we can open the file */
  if ((test = fopen((char *)file_brs, "a")) == NULL) {
    sprintf(messageStr, "Unable to open BRS output file %s", file_brs);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(test);
  }

  printf("%s\n",file_brs);

  /* Create the file */
  if ((ierror = nc_create(file_brs, NC_NETCDF4 | NC_CLOBBER, &ncid)))
    ERR(ierror,routineName);

  /* Write atmos.ID as global attribute */
  if ((ierror = nc_put_att_text(ncid, NC_GLOBAL, "atmosID", strlen(atmos.ID),
				atmos.ID ))) ERR(ierror,routineName);

  /* Create dimensions */ 
  if ((ierror = nc_def_dim(ncid, "task", NC_UNLIMITED,    &task_id   ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid, "nwave", spectrum.Nspect, &nspect_id ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid, "nrec", Nrecno,          &nrec_id   ))) 
    ERR(ierror,routineName);
  if ((ierror = nc_def_dim(ncid, "nspace", atmos.Nspace,  &nspace_id ))) 
    ERR(ierror,routineName);
  /* this last dimension is not used, but kept to save Nspace */

  /* Create variables */
  dimids[0] = task_id;
  dimids[1] = nspect_id;

  /* presence of background line */
  if ((ierror = nc_def_var(ncid, HASLINE_VAR, NC_UBYTE, 2, dimids, &hl_var   )))
    ERR(ierror,routineName);
  /* presence of polarised background */
  if ((ierror = nc_def_var(ncid, ISPOL_VAR,   NC_UBYTE, 2, dimids, &ip_var   )))
    ERR(ierror,routineName);
  /* record numbers */
  dimids[1] = nrec_id;
  if ((ierror = nc_def_var(ncid, BGREC_VAR,   NC_LONG,  2, dimids, &nrec_var )))
    ERR(ierror,routineName);
 
  /* End define mode */
  if ((ierror = nc_enddef(ncid))) ERR(ierror,routineName);
 
  /* Copy stuff to the IO data struct */
  io.brs_ncid     = ncid;
  io.brs_hl_var   = hl_var;
  io.brs_ip_var   = ip_var;
  io.brs_nrec_var = nrec_var;

  return;
}
/* ------- end ---------------------------- init_ncdf_BRS.c --------- */


/* ------- begin -------------------------- close_ncdf_BRS.c --------- */
void close_ncdf_BRS(void)
/* Closes the BRS netCDF file */
{
  const char routineName[] = "close_ncdf_BRS";
  int        ierror;

  if ((ierror = nc_close(io.brs_ncid))) ERR(ierror,routineName);
  return;
}
/* ------- end ---------------------------- close_ncdf_BRS.c --------- */



/* ------- begin -------------------------- writeBRS_ncdf.c -------------- */

void writeBRS_ncdf(void)
{
  const char routineName[] = "writeBRS_ncdf";
  unsigned char *hasline, *ispolarised;
  int     Nrecno, nspect, ierror, ncid;
  size_t  start[] = {0, 0};
  size_t  count[] = {1, 1};
  
  ncid = io.brs_ncid;

  /* For unlimited dimension */
  start[0] = (size_t) mpi.task; 
  
  if (atmos.moving || atmos.Stokes)
    Nrecno = 2 * spectrum.Nspect * atmos.Nrays;
  else
    Nrecno = spectrum.Nspect;

  /* create some arrays for data */
  hasline     =
    (unsigned char *) malloc(spectrum.Nspect * sizeof(unsigned char));
  ispolarised =
    (unsigned char *) malloc(spectrum.Nspect * sizeof(unsigned char));
  
  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++){
    ispolarised[nspect] = atmos.backgrflags[nspect].ispolarized;
    hasline[nspect]     = atmos.backgrflags[nspect].hasline;
  }

  
  /* write data to file */
  count[1] = (size_t) spectrum.Nspect;
  if ((ierror = nc_put_vara_ubyte(io.brs_ncid, io.brs_hl_var,   start, count, 
				  hasline )))           ERR(ierror,routineName);
  if ((ierror = nc_put_vara_ubyte(io.brs_ncid, io.brs_ip_var,   start, count,
				  ispolarised )))       ERR(ierror,routineName);
  count[1] = (size_t) Nrecno;
  if ((ierror = nc_put_vara_long( io.brs_ncid, io.brs_nrec_var, start, count,
				  atmos.backgrrecno ))) ERR(ierror,routineName);

  
  free(hasline);
  free(ispolarised);
  return;

}
/* ------- end ---------------------------- writeBRS_ncdf.c -------------- */


/* ------- begin -------------------------- readBRS_ncdf.c --------------- */
bool_t readBRS_ncdf(void)
{
  /****************************************************************
   **** WARNING: this routine has NOT been tested, must check  ****
   ****          with solveray!!!                              ****
   ****************************************************************/

  const char routineName[] = "readBRS_ncdf";
  unsigned char *hasline, *ispolarised;
  int     ierror, ncid, task_id, nrec_id, nspace_id, wave_id,
          hl_var, ip_var, nrec_var, Nrecno, nspect;
  size_t  task_len, wave_len, nrec_len, nspace_len;
  size_t  start[] = {0, 0};
  size_t  count[] = {1, 1};
  FILE   *test;
  char   *fname = BRS_FILE, *fext  = BRS_EXT, *atmosID;
  char    file_brs[MAX_MESSAGE_LENGTH];

 
  if (atmos.moving || atmos.Stokes)
    Nrecno = 2 * spectrum.Nspect * atmos.Nrays;
  else
    Nrecno = spectrum.Nspect;

  /* get file name, with the MPI rank */
  sprintf(file_brs,"scratch/%s%d%s", fname, mpi.rank, fext);

  /* Check if we can open the file */
  if ((test = fopen((const char *)file_brs, "r")) == NULL) {
    sprintf(messageStr, "Unable to open BRS output file %s", file_brs);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    fclose(test);
  }

  /* Open the file */
  if ((ierror = nc_open(file_brs, NC_NOWRITE, &ncid))) ERR(ierror,routineName);

  /* Get the dimids and values */
  if ((ierror = nc_inq_dimid( ncid, "task",     &task_id  ))) ERR(ierror,routineName);
  if ((ierror = nc_inq_dimlen(ncid, task_id,    &task_len ))) ERR(ierror,routineName);

  if ((ierror = nc_inq_dimid( ncid, "nwave",     &wave_id  ))) ERR(ierror,routineName);
  if ((ierror = nc_inq_dimlen(ncid, wave_id,    &wave_len ))) ERR(ierror,routineName);

  if ((ierror = nc_inq_dimid( ncid, "nrec",     &nrec_id  ))) ERR(ierror,routineName);
  if ((ierror = nc_inq_dimlen(ncid, nrec_id,    &nrec_len ))) ERR(ierror,routineName);

  if ((ierror = nc_inq_dimid( ncid, "nspace", &nspace_id  ))) ERR(ierror,routineName);
  if ((ierror = nc_inq_dimlen(ncid, nspace_id,&nspace_len ))) ERR(ierror,routineName);

  /* Get the varids */
  if ((ierror = nc_inq_varid(ncid, HASLINE_VAR, &hl_var   )))  ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(ncid, ISPOL_VAR,   &ip_var   )))  ERR(ierror,routineName);
  if ((ierror = nc_inq_varid(ncid, BGREC_VAR,   &nrec_var )))  ERR(ierror,routineName);

  /* consistency checks */
  if (nspace_len != atmos.Nspace || wave_len != spectrum.Nspect 
      || Nrecno != nrec_len ) return FALSE;
 
  if ((ierror = nc_get_att_text( ncid, NC_GLOBAL, "atmosID", atmosID ))) 
    ERR(ierror,routineName);

  if (!strstr(atmosID, atmos.ID)) {
    sprintf(messageStr,
	    "Input was derived from different atmosphere (%s) than current",
	    atmosID);
    Error(WARNING, routineName, messageStr);
  }
  free(atmosID);

  /* Get variables */
  atmos.backgrflags = (flags *) malloc(spectrum.Nspect * sizeof(flags));
  atmos.backgrrecno = (long  *) malloc(Nrecno          * sizeof(long) );

  hasline     = (unsigned char *) malloc(wave_len * sizeof(unsigned char));
  ispolarised = (unsigned char *) malloc(wave_len * sizeof(unsigned char));

  /* Background opacity record numbers  */
  start[0] = mpi.task;
  count[1] = nrec_len;
  if ((ierror = nc_get_vara_long(ncid, nrec_var, start, count, atmos.backgrrecno )))
    ERR(ierror,routineName);
  
  count[1] = wave_len;
  /* --- Flags for presence of background line --    -------------- */
  if ((ierror = nc_get_vara_ubyte(ncid, hl_var,  start, count, hasline     ))) 
    ERR(ierror,routineName);
  /* --- Flags for presence of polarized background -- ------------ */
  if ((ierror = nc_get_vara_ubyte(ncid, hl_var,  start, count, ispolarised ))) 
    ERR(ierror,routineName);

  /* Put back in structure */
  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++){
    atmos.backgrflags[nspect].ispolarized  =  ispolarised[nspect];
    atmos.backgrflags[nspect].hasline      =  hasline[nspect];
  }
  
  /* Close the file. */
  if ((ierror = nc_close(ncid))) ERR(ierror,routineName);
  
  return TRUE;
}

/* ------- end ---------------------------- readBRS_ncdf.c --------------- */

/* ------- begin -------------------------- writeBRS.c -------------- */

void writeBRS_p(void)
{
  const char routineName[] = "writeBRS";

  FILE *fp_out;
  XDR   xdrs;


  if ((fp_out = fopen(BRS_DOT_OUT, "a")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", BRS_DOT_OUT);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  if (!xdr_BRS(&xdrs)) {
    sprintf(messageStr, "Unable to write to output file %s", BRS_DOT_OUT);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp_out);
}
/* ------- end ---------------------------- writeBRS.c -------------- */

/* ------- begin -------------------------- readBRS.c --------------- */

void readBRS(void)
{
  const char routineName[] = "readBRS";

  FILE *fp_in;
  long  ftask_size;
  XDR   xdrs;

  if ((fp_in = fopen(BRS_DOT_OUT, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", BRS_DOT_OUT);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* Skip the number of bytes relative to previous tasks */
  if (mpi.task > 0){
    ftask_size = strlen(atmos.ID) + sizeof(long) + sizeof(int);
  }
  
  
  xdrstdio_create(&xdrs, fp_in, XDR_DECODE);

  if (!xdr_BRS(&xdrs)) {
    sprintf(messageStr, "Unable to read from input file %s", BRS_DOT_OUT);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp_in);
}
/* ------- end ---------------------------- readBRS.c --------------- */

/* ------- begin -------------------------- xdr_BRS.c --------------- */

bool_t xdr_BRS(XDR *xdrs)
{
  const char routineName[] = "xdr_BRS";
  register int nspect;

  char  *atmosID;
  bool_t result = TRUE, *hasline, *ispolarized;
  int    Nsp;
  long   Nspace, Nrecno;
  

  hasline     =
    (bool_t *) malloc(spectrum.Nspect * sizeof(bool_t));
  ispolarized =
    (bool_t *) malloc(spectrum.Nspect * sizeof(bool_t));

  if (atmos.moving || atmos.Stokes)
    Nrecno = 2 * spectrum.Nspect * atmos.Nrays;
  else
    Nrecno = spectrum.Nspect;

  if (xdrs->x_op == XDR_ENCODE) {
    atmosID = atmos.ID;
    result &= xdr_counted_string(xdrs, &atmosID);
    result &= xdr_long(xdrs, (int *) &atmos.Nspace);
    result &= xdr_int(xdrs, (int *) &spectrum.Nspect);

    /* --- Flags for presence of background line --    -------------- */
 
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
      hasline[nspect] = atmos.backgrflags[nspect].hasline;
    result &= xdr_vector(xdrs, (char *) hasline, spectrum.Nspect,
			 sizeof(bool_t), (xdrproc_t) xdr_bool);

    /* --- Flags for presence of polarized background -- ------------ */
 
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
      ispolarized[nspect] = atmos.backgrflags[nspect].ispolarized;
    result &= xdr_vector(xdrs, (char *) ispolarized, spectrum.Nspect,
			 sizeof(bool_t), (xdrproc_t) xdr_bool);
  } else {
    atmos.backgrflags = (flags *) malloc(spectrum.Nspect * sizeof(flags));
    atmos.backgrrecno = (long *) malloc(Nrecno * sizeof(long));

    result &= xdr_counted_string(xdrs, &atmosID);
    if (!strstr(atmosID, atmos.ID)) {
      sprintf(messageStr,
	      "Input was derived from different atmosphere (%s) than current",
              atmosID);
      Error(WARNING, routineName, messageStr);
    }
    free(atmosID);

    result &= xdr_long(xdrs, (int *) &Nspace);
    result &= xdr_int(xdrs, &Nsp);
    if (Nspace != atmos.Nspace || Nsp != spectrum.Nspect) {
      free(hasline);
      return FALSE;
    }
    /* --- Flags for presence of background line --    -------------- */
 
    result &= xdr_vector(xdrs, (char *) hasline, spectrum.Nspect,
			 sizeof(bool_t), (xdrproc_t) xdr_bool);
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
      atmos.backgrflags[nspect].hasline = hasline[nspect];

    /* --- Flags for presence of polarized background -- ------------ */
 
    result &= xdr_vector(xdrs, (char *) ispolarized, spectrum.Nspect,
			 sizeof(bool_t), (xdrproc_t) xdr_bool);
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
      atmos.backgrflags[nspect].ispolarized = ispolarized[nspect];
  }

  result &= xdr_vector(xdrs, (char *) atmos.backgrrecno, Nrecno,
		       sizeof(long), (xdrproc_t) xdr_int);

  free(hasline);
  free(ispolarized);

  return result;
}
/* ------- end ---------------------------- xdr_BRS.c --------------- */
