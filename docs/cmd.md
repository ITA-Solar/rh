# Command line tools

Two useful command line tools that come with HDF5 are
[h5dump](https://support.hdfgroup.org/HDF5/doc/RM/Tools.html#Tools-Dump)
and
[h5repack](https://support.hdfgroup.org/HDF5/doc/RM/Tools.html#Tools-Repack).

`h5dump` can be used with the `-H` option to look at the header of a
file: see the dimensions, variables, groups. It can also be used to
print a text version of any variable in an HDF5 file (e.g. this can be
redirected to a text file). When printing a variable (dataset in HDF5)
one uses the option `-d variable`, and the resulting output is the same
as in the `-H` mode, with the variable printed at the end. The NetCDF
[ncdump](https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf/ncdump.html)
program offers an even clearer look into the file (e.g. used with the
`-h` option to print out the header).

The `h5repack` program can be used to copy and modify the parameters of
HDF5 files. It can convert the files between different format versions,
compress variables, etc. Of particular importance is the option for
rechunking a file. Chunking in HDF5 files can be used to improve
performance by changing the disk structures to improve different read
patterns. It is analogous to fully or partially transposing the
variables along certain dimensions.

!!! note "See also"

    [h5dump guide](https://support.hdfgroup.org/HDF5/doc/RM/Tools.html#Tools-Dump): Detailed information about `h5dump`.

    [h5repack guide](https://support.hdfgroup.org/HDF5/doc/RM/Tools.html#Tools-Repack): Detailed information about `h5repack`.

    [Chunking in HDF5](http://www.hdfgroup.org/HDF5/doc/Advanced/Chunking/): Description on the advantages of chunking.
