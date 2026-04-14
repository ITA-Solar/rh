# Reading output in IDL

There are no specific IDL routines for reading the output from RH 1.5D.
However, there is a utility function that can be used to variables from
HDF5/netCDF4 files, under the `idl/` directory in a file named
`read_ncdf_var.pro`. The function `read_ncdf_var()` can be used to read
variables from an HDF5 or netCDF4 file, e.g.:

``` fortran
IDL> data = read_ncdf_var("output_ray.hdf5", "intensity")
IDL> help, data
DATA            FLOAT     = Array[902, 512, 512]
IDL> pops = read_ncdf_var("output_aux.hdf5", "populations", groupname="atom_CA")
IDL> help, pops
POPS            FLOAT     = Array[400, 512, 512, 5]
```

!!! note

    The IDL analysis suite of RH **does not work with RH 1.5D**.
