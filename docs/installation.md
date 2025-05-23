# Installation

## Getting the code

The code is available on a git repository, hosted on github:
<https://github.com/ita-solar/rh>. If you don't have
[git](http://git-scm.com/) installed and just want to get started, the
easiest way is to download a zip file with the latest revision:
<https://github.com/ita-solar/rh/archive/master.zip>. If you have git
installed and would like to be up-to-date with the repository, you can
do a git clone:

    git clone https://github.com/ita-solar/rh.git

or using SSH:

    git clone git@github.com:ita-solar/rh.git

Whether you unpack the zip file or do one of the above it will create a
directory called `rh` in your current path. This directory will have the
following subdirectories:

|  Directory        |  Contents    |
|:------------------|:-----------------------------------------------------------------|
|  rh               |  Main RH source                               |
|  rh/Atmos         |  Used to keep atmosphere files  |
|  rh/Atoms_example |  Used to keep atom files, line and wavelength lists (example)  |
|  rh/idl           |  Old RH IDL routines, not used  |
|  rh/Molecules     |  Used to keep molecule files  |
|  rh/python        |  Utility Python programs  |
|  rh/rh15d.        |  Source files for RH 1.5D  |
|  rh/rhf1d         |  Source files for 1D geometry :material-close: **deprecated, do not use**  |
|  rh/rhsphere      | Source files for spherical geometry :material-close: **deprecated, do not use**  |
|  rh/tools         |  Associate C programs for RH, not tested. |


!!! warning

    The source code directories for other geometries (rhf1d, rhsphere) are
    still in the code tree, but they are deprecated and will be removed
    soon. With the latest changes related to rh15d, they are not guaranteed
    to work or even run. Do not use.

## Dependencies

### HDF5

RH 1.5D makes use of the [HDF5](http://www.hdfgroup.org/HDF5/) library
to read the atmosphere files and write the output. It is not possible to
run the code without this library. RH 1.5D requires HDF5 version 1.8.1
or newer (including versions 1.10.x).

!!! info

    RH 1.5D previously made use of the netCDF4 library for its output (which
    in turn also required HDF5). The latest changes mean RH 1.5D needs only
    HDF5. Because netCDF4 files are also HDF5 files, the output is still
    readable in the same way as before and input files in netCDF version 4
    format can still be read in the same way by RH 1.5D. If you used input
    atmospheres in netCDF version 3 format, then these will have to be
    converted to HDF5. It is recommended that new atmosphere files be
    created in HDF5 only.


Because HDF5 is commonly used in high-performance computing, many
supercomputers already have them available. Here are a few setups for
different supercomputers:

Betzy:

    module load iompi/2022a HDF5/1.12.2-iompi-2022a

Fram:

    module load HDF5/1.8.19-intel-2018a intel/2018a

Pleiades:

    module load hdf5/1.8.18_mpt

Vilje:

    module load intelcomp/18.0.1 mpt/2.14 hdf5/1.8.19

Hexagon:

    module load cray-hdf5-parallel

and at ITA's Linux system:

    module load intel/oneapi compiler/latest mpi/latest hdf5/Intel/1.14.3

### MPI

You need MPI to run RH 1.5D. In supercomputers and clusters these are
provided, but for your workstation or laptop you may need to install
manually both MPI and HDF5 (with parallel support). If you have an
existing conda or mamba installation, it should be possible to install
compilers, MPI libraries, and HDF5 parallel. For example, the following
will install the latest version of HDF5 with OpenMPI as the MPI library:

    conda install -c conda-forge 'hdf5=*=*openmpi*'

This will install HDF5 to the same directory of your python environment
(files under `lib/` and `bin/`). This setup has been tested on macOS
with Apple Silicon, and works for RH 1.5D. However, packages may change
and you may need to specify another way. Installing HDF5 via the system
package manager (e.g. apt-get in Linux) is not recommended, since those
HDF5 builds will probably not have parallel support.

If installing binaries fails, the safest bet is to download and compile
HDF5 from the source, enabling parallel builds in the `./configure`
script, e.g.:

    ./configure (...) --enable-parallel

## Compilation

Compilation of RH 1.5D consists of two steps:

1.  Compilation of the geometry-independent main libraries (`librh.a`
    and `librh_f90.a`)
2.  Compilation of the `rh15d_mpi` tree and main binaries

RH 1.5D has been compiled in a variety of architectures and compilers,
including gcc, the Intel compilers, and clang. As for MPI
implementations, it has been tested with SGI's mpt, OpenMPI, mpich,
mvapich, and Intel's MPI.

### Makefile configuration

RH 1.5D does not automatically look for the compilers and libraries. You
need to tell RH which compilers to use and where to find the HDF5
library by editing the file `rh/Makefile.config`. This file is also used
to set up any additional compiler or linker flags, if appropriate.
**Changes to any other Makefiles are not necessary.** It is also no
longer necessary to set the environment variables `OS` and `CPU`, as in
previous versions.

For `HDF5_DIR`, please enter the base directory for the library (not the
directory with the lib* files), so that both library and include files
are used. In Fram and Hexagon this is already stored in the `HDF5_DIR`
environment variable, so you can comment that line in `Makefile.config`.
If your version of HDF5 was not built as a shared binary, you need to
link HDF5 and other used libraries directly (you will need to set at
least `-lz` in `LDFLAGS`).

The following compiler flags are recommended for Betzy:

    CFLAGS = -O3 -DHAVE_F90 -qopt-prefetch -use-intel-optimized-headers -march=core-avx2 -fp-model source  
    F90FLAGS = -O3 -qopt-prefetch -use-intel-optimized-headers -march=core-avx2 -fp-model source  

And for the ITA linux system (RHEL 9.x):

    CFLAGS = -O3 -DHAVE_F90 -Wformat  -I/usr/include/tirpc/  -std=gnu89
    F90FLAGS = -O3 -Wformat  -I/usr/include/tirpc/  -std=gnu89
    LDFLAGS = -ltirpc
    HDF5_DIR = /astro/local/hdf5/rhel9/1.14.3/intel/

There are two steps in the compilation: main libraries and rh15d
binaries. To speed up compilation, you can use parallel builds (e.g.
`make -j8`) in all steps of the compilation.

### Main libraries

The common RH files are put in a library under the base directory. After
editing `Makefile.config`, build the main libraries with `make` on the
`rh` directory. If successful, the compilation will produce the two
library files `librh.a` and `librh_f90.a`.

### Program binaries

The `rh15d` contains the source files for the 1.5D version. After
compiling the main library, go to that directory and compile the
binaries with `make`. The following executables will be created:

|  <div style="width:100px">File</div>         |   Description |
:---------------|:---------------
  `rh15d_ray_pool` |  Main RH 1.5D binary, uses a job pool
  `rh15d_ray`      |  Alternative RH 1.5D binary. **Deprecated.** This program runs much slower than `rh15d_ray_pool` and is kept for backwards compatibility only. Will be removed in a future revision.
  `rh15d_lteray`   |  Special binary for running in LTE

## Run directory

Once compiled, you can copy or link the binaries to a run directory.
This directory will contain all the necessary input files, and it should
contain two subdirectories called `output` and `scratch`.

!!! warning

    If the subdirectories `output` and `scratch` do not exist in the
    directory where the code is run, the code will crash with an obscure
    error message.

