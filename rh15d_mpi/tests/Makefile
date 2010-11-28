readatmos: readatmos_ncdf.c 
	mpicc -Wall readatmos_ncdf.c -o readatmos_ncdf -lnetcdf -L/opt/local/lib/ ../../librh.a 
pread: sfc_pres_temp_rd.c
	mpicc sfc_pres_temp_rd.c -o pread -lnetcdf -L/opt/local/lib/
pwrite: sfc_pres_temp_wr.c
	mpicc sfc_pres_temp_wr.c -o pread -lnetcdf -L/opt/local/lib/
