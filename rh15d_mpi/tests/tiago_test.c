/* tests to read atmosphere in netCDF format */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <netcdf.h>
#include <mpi.h>

#define FILE_NAME "/Users/tiago/data/bifrost/bifrost_s20_slice.ncdf"

#define TEMP_NAME "temperature"
#define VZ_NAME "velocity_z"
#define NE_NAME "electron_density"
#define NH_NAME "hydrogen_populations"

#define ERR(e) {printf("(EEE) NetCDF: %s\n", nc_strerror(e)); return;}

void main(int argc, char **argv)
{
  int ncid, temp_varid, x_varid, nx_id,ny_id,nz_id, i, retval;
  size_t nx, ny, nz;
  size_t start[] = {0, 0, 0}; /* starting values */
  size_t count[] = {1, 1, 1};

  double *temp, *x;
  char x_units[100];

  /* MPI stuff. */
  int mpi_namelen;		
  char mpi_name[MPI_MAX_PROCESSOR_NAME];
  int mpi_size, mpi_rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Get_processor_name(mpi_name, &mpi_namelen);

  /* Open the file. */
  if ((retval = nc_open_par(FILE_NAME, NC_NOWRITE, comm, info, &ncid))) ERR(retval);

  /* Get the varids and dimids */
  if (retval = nc_inq_varid(ncid, TEMP_NAME, &temp_varid)) ERR(retval);
  if (retval = nc_inq_varid(ncid, "x", &x_varid)) ERR(retval);
  if (retval = nc_inq_dimid(ncid, "nx", &nx_id)) ERR(retval);
  if (retval = nc_inq_dimid(ncid, "ny", &ny_id)) ERR(retval);
  if (retval = nc_inq_dimid(ncid, "nz", &nz_id)) ERR(retval);

  /* get the dimension values */
  if (retval = nc_inq_dimlen(ncid, nx_id, &nx)) ERR(retval);
  if (retval = nc_inq_dimlen(ncid, ny_id, &ny)) ERR(retval);
  if (retval = nc_inq_dimlen(ncid, nz_id, &nz)) ERR(retval);

  if (mpi_rank == 0){
    printf("nx is %d\n", nx);
    printf("ny is %d\n", ny);
    printf("nz is %d\n", nz);
  }

  /* read some stuff */
  temp  = (double *) calloc(nz, sizeof(double));
  x     = (double *) calloc(nx, sizeof(double));

  // read the whole variable, in this case
  if (retval = nc_get_var_double(ncid, x_varid, x)) ERR(retval);

  // read only one column, adjust size
  count[2] = nz;
  start[0] = mpi_rank;
  start[1] = mpi_rank;

  if (retval = nc_get_vara_double(ncid, temp_varid, start, count, temp)) ERR(retval);
  
  printf("%d\n", temp_varid);

  
  for (i=0; i < nz; i++){
    printf("Process %d: temp[%d] = %e\n",mpi_rank, i, temp[i]);
  }
  

   


  /* Close the file. */
  if (retval = nc_close(ncid)) ERR(retval);


  /* free stuff */
  free(temp);
  free(x);

  printf("-- End\n");
  MPI_Finalize();

  return;
}
