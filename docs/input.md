# Input files

## Configuration files

The configuration of an RH 1.5D is made primarily through several text
files that reside in the run directory. The main file is
`keyword.input`. All the other files and their locations are specified
in `keyword.input`. The source tree contains a sample `rh/rh15d/run/`
directory with the following typically used configuration files:

|  File            |  Description  |
|:-----------------|:---------------|
|  `atoms.input`     |  Lists the atom files to be used  |
|  `keyword.input`   |  Main configuration file  |
|  `kurucz.input`    |  Contains list of line lists to be used  |
|  `molecules.input` |  Lists the molecule files to be used  |
|  `ray.input`       |  Selects output mu and wavelengths for detailed output  |

The `kurucz.input` and the `molecules.input` files are identical under
RH, so we refer to the RH manual for more information about them. Most
of the other files behave very similarly in RH and RH 1.5D, with a few
differences.

The `atoms.input` file is identical in RH, but it can also have a new
starting solution, `ESCAPE_PROBABILITY`.

The `keyword.input` file functions in very much the same manner under RH
and RH 1.5D. The main difference is that there are new options for the
1.5D version, and some options should not be used.

The new `keyword.input` options for the 1.5D version are:


|  <div style="width:190px">Name</div>                 |    <div style="width:120px">Default value</div>  |    Description |
|:-----------------------|:-----------------|:-------------------------------------------
|  `SNAPSHOT`              |   `0`             |   Snapshot index from the atmosphere file.  |
| `X_START`                |  `0`              |  Starting column in the x direction. If < 0, will be set to 0. |
| `X_END`                   | `-1`            |   Ending column in the x direction. If <= 0, will be set to `NX` in the atmosphere file. |
|  `X_STEP`                 |  `1`              |  How many columns to sample in the y direction. If < 1, will be set to 1. |
|  `Y_START`                 | `0`             |   Starting column in the y direction. If < 0, will be set to 0. |
|  `Y_END`                   | `-1`           |    Ending column in the y direction. If <= 0, will be set to `NY` in the atmosphere file.  |
|  `Y_STEP`                  | `1`           |     How many columns to sample in the y direction. If < 1, will be set to 1. |  
|  `15D_WRITE_POPS`          | `FALSE`          |  If `TRUE`, will write the level populations (including LTE) for each active atom to `output_aux.hdf5`.  |
|  `15D_WRITE_RRATES`        | `FALSE`          |  If `TRUE`, will write the radiative rates (lines and continua) for each active atom to `output_aux.hdf5`.  |
|  `15D_WRITE_CRATES`        | `FALSE`          |  If `TRUE`, will write the collisional rates for each active atom to `output_aux.hdf5`.  |
|  `15D_WRITE_TAU1`          | `FALSE`           | If `TRUE`, will write the height of tau=1 to `output_ray.hdf5`, for all wavelengths (this takes up as much space as the intensity).  |
|  `15D_RERUN`               | `FALSE`           | If `TRUE`, will rerun for non-converged columns.  |
|  `15D_DEPTH_ZCUT`          | `TRUE`            | If `TRUE`, will perform a cut in z for  points above a threshold temperature  |
|  `15D_TMAX_CUT`            | `-1`             |  Threshold temperature (in K) over which the above depth cut ids made. If < 0, no temperature cut will be made.  |
|  `15D_DEPTH_REFINE`        | `FALSE`          |  If `TRUE`, will perform an optimisation of the depth scale, based on optical depth, density and temperature gradients.  |
|  `BACKGR_IN_MEM`           | `FALSE`          |  If `TRUE`, will keep background opacity coefficients in memory instead of scratch files on disk.  |
|  `BARKLEM_DATA_DIR`        | `../../Atoms`    |  Directory where the `Barklem_*data.dat` data files are saved.  |
|  `COLLRAD_SWITCH`          | `0.0`            |  Defines if collisional radiative switching is on. If < 0, switching parameter is constant (and equal to `COLLRAD_SWITCH_INI`). If = 0, no collisional radiative switching. If > 0, collisional radiative switching decreases by `COLLRAD_SWITCH` per log decade,  starting with `COLLRAD_SWITCH_INI`.  |
|  `COLLRAD_SWITCH_INIT`     | `1.0`            |  Initial increment for collisional-radiative  switching  |
|  `LIMIT_MEMORY`            | `FALSE`          |  If `TRUE`, will not keep several large  arrays in memory but rather save them to scratch files. Not recommended unless memory usage is critical.  |
|  `N_PESC_ITER`             | `3`              |  Number of escape probability iterations, if any atoms have it as initial solution.  |
|  `PRD_SWITCH`              | `0.0`            |  If > 0, the PRD effects will be added gradually, converging to the full PRD solution in `1/sqrt(PRD_SWITCH)` iterations.  |
|  `PRDH_LIMIT_MEM`          | `FALSE`          |  If `TRUE` and using `PRD_ANGLE_APPROX`, will not keep in memory quantities necessary to calculate the current PRD weights, but rather calculate them again. Will affect the performance, so should be used only when necessary.  |
|  `S_INTERPOLATION`         | `LINEAR`         |  Type of source function interpolation to use in formal solver. Can be `LINEAR`, `BEZIER`, `BEZIER3` or `CUBIC_HERMITE`.  |
|  `S_INTERPOLATION_STOKES` |  `DELO_PARABOLIC` |  Type of source function interpolation to  use in formal solver for polarised cases. Can be `DELO_PARABOLIC` or `DELO_BEZIER3`, see[^1] .  |
|  `VTURB_MULTIPLIER`       |  `1.0`           |   Atmospheric `vturb` will be multiplied by this value  |
|  `VTURB_ADD`             |   `0.0`            |  Value to be added to atmospheric `vturb`  |


The `X_START`, `X_END`, and `X_STEP` keywords (and the equivalent for
the y direction) define which columns of the atmosphere file are going
to be run. They can be used to calculate only a specific region. RH 1.5D
chooses the columns to calculate using the `(start, end, step)`
parameters as in the `range()` function in
[Python](http://docs.python.org/2/library/functions.html#range): the
result is `[start, start + step, start + 2 * step, ...]`. The last
element is the largest `start + i * step` less than `end`. This means
that the numbers given by `X_END` and `Y_END` are **not inclusive**
(e.g. if `nx = 50` and `X_END = 49`, the column with the index 49 will
not be calculated). One must set `X_END = nx` to calculate all the
columns.

The following options have a different meaning under RH 1.5D:

|  <div style="width:120px">Name</div>                 |    <div style="width:120px">Default value</div>    |  Description   |
|:-------------------|:------------------|:--------------------------------------------|
|  `PRD_ANGLE_DEP`   |  `PRD_ANGLE_INDEP` |  This keyword is no longer boolean. To accommodate for new options, it now takes the values `PRD_ANGLE_INDEP` for  angle-independent PRD, `PRD_ANGLE_DEP` for angle-dependent PRD, and `PRD_ANGLE_APPROX` for the approximate angle-dependent scheme of Leenaarts et al. (2012)[^2] .  |
|  `BACKGROUND_FILE`  |                  |   This keyword is no longer the name of the  background file, but the prefix of the background files. There will be one file per process, and the filenames are this prefix plus `_i.dat`, where `i` is the process number.  |
|  `STOKES_INPUT`    |                 |     This option is not used in RH 1.5D because  the magnetic fields are now written to the atmosphere file. However, it **must** be set to any string if one is using any `STOKES_MODE` other than `NO_STOKES` (RH  won't read B otherwise).  |


And the following options are valid for RH but may not work with RH
1.5D:

|  <div style="width:120px">Name</div>                 |    <div style="width:120px">Default value</div>    |  Description   |
|:-------------------|:------------------|:--------------------------------------------|
|  `LIMIT_MEMORY` |  `FALSE` |  This option has not been tested and may not work well with RH 1.5D. |
|  `PRINT_CPU`  |    `FALSE` |  This option does not work with RH 1.5D and **should always be** `FALSE`.
|  `N_THREADS` |     `0`     |  Thread parallelism will not work with RH 1.5D. This option should always be `0` or`1`. |


The `ray.input` has the same structure in RH1D and RH 1.5D. In RH it is
used as input for the `solveray` program, but in RH 1.5D it is used for
the main program. It should contain the following:

    1.00
    Nsource

The first line is the `mu` angle for the output ray, and it should
always be 1.00. The second line is `Nsource`, the number of wavelengths
for which detailed output (typically source function, opacity, and
emissivities) will be written. If `Nsource > 0`, it should be followed
in the same line by the indices of the wavelengths (e.g. `0 2 10 20`).

## Atom and molecule files

The atom and molecule files have the same format as in RH. In the
`rh/Atoms` and `rh/Molecules` directories there are a few sample files.
They are read by the procedures in `readatom.c` and `readmolecule.c`.
The atom files have the following basic structure:


|  <div style="width:190px">Input</div>         |  Format   |
|:-----------------------------|:------------------------------------------|
|  `ID`                        |  **(A2)**. Two-character atom identifier.
|  `Nlevel Nline Ncont Nfixed` |  **(4I)**. Number of levels, lines,  continua, and fixed radiation temperature transitions. |
|  `level_entries`   |            Nlevel * **(2F, A20, I)** |
|  `line entries`     |           Nline * **(2I, F, A, I, A, 2F, A, 6F)**| 
|  `continuum_entries` |          Ncont * **(I, I, F, I, A, F)**  |
|  `fixed_entries`      |         Ncont * **(2I, 2F, A)**  |


## Atmosphere files

The atmosphere files for RH 1.5D are a significant departure from RH.
They are written in the flexible and self-describing
[HDF5](https://www.hdfgroup.org/HDF5/) format. They can be written with
any version, except the 1.10.x development branch.

The atmosphere files contain all the atmospheric variables necessary for
RH 1.5D, and they may contain one or more simulation snapshots. The
basic dimensions of the file are:

|  |  |
|:--------|:----------------------------|
|  `nt`    |  Number of snapshots. |
|  `nx`    |  Number of x points  |
|  `ny`    |  Number of y points.  |
|  `nz`    |  Number of depth points.  |
|  `nhydr` |  Number of hydrogen levels.  |

While strictly 3D atmosphere files, 2D and 1D snapshots can also be used
provided that one or both of `nx` and `ny` are equal to 1.

!!! note 
   
    The atmosphere variables must be written to the file in a particular
    way. They should be written in a *height grid* (meaning the top of the
    atmosphere has a larger value of `z`), and must *start from the top*
    (meaning that the first height index of the arrays must be the **TOP**
    of the atmosphere). Failure to follow these two rules can either lead to
    RH aborting a run, or worse, getting wrong results without a clear error
    message.

The atmosphere file can contain the following variables:


|  <div style="width:140px">Name</div>  |    Dimensions       |    Units   |   Notes |
:--------------------|:---------------------|:----------|:-------------------------
  `B_x`               |     `(nt, nx, ny, nz)`    |      T   |       Magnetic field x component. **Optional**
  `B_y`                |    `(nt, nx, ny, nz)`      |    T       |   Magnetic field y component. **Optional**
  `B_z`                 |   `(nt, nx, ny, nz)`       |   T       |   Magnetic field z component. **Optional**
  `electron_density`    |   `(nt, nx, ny, nz)`       |   m^-3^   |   **Optional**.
  `hydrogen_populations` |  `(nt, nhydr, nx, ny, nz)`|   m^-3^    |  `nhydr` must correspond to the number of levels  in the hydrogen atom used. If `nhydr=1`, this  variable should contain the total number of hydrogen atoms (in all levels), and LTE populations will be calculated.
  `snapshot_number`       | `(nt)`                   |   None    |   The snapshot number is an array of integers to identify each snapshot in the output files.
  `temperature`          |  `(nt, nx, ny, nz)`      |    K       |
  `velocity_z`           |  `(nt, nx, ny, nz)`      |    m s^-1^  |  Vertical component of velocity. **Positive velocity is upflow.**
  `velocity_turbulent`   |  `(nt, nx, ny, nz)`      |    m s^-1^  |  Turbulent velocity (microturbulence).
  `z`                    |  `(nt, nx, ny, nz)` or  `(nt, nz)` |     m     |     Height grid. Can vary with column and snapshot. First `nz` index is **top of the atmosphere** (closest to observer).


Any other variable in the file will not be used. In addition, the
atmosphere file **must** have a global attribute called `has_B`. This
attribute should be 1 when the magnetic field variables are present, and
0 otherwise. Also recommended, but optional, is a global attribute
called `description` with a brief description of the atmosphere file
(e.g. how and from they were generated).

!!! note

    Variables in the atmosphere files can be compressed (zlib or szip), but compression is not recommended for performance reasons.

As HDF5 files, the contents of the atmosphere files can be examined with
the `h5dump` utility. To see a summary of what's inside a given file,
one can do:

    h5dump -H atmosfile

Here is the output of the above for a sample file:

    HDF5 "example.hdf5" {
    GROUP "/" {
       ATTRIBUTE "boundary_bottom" {
          DATATYPE  H5T_STD_I64LE
          DATASPACE  SCALAR
       }
       ATTRIBUTE "boundary_top" {
          DATATYPE  H5T_STD_I64LE
          DATASPACE  SCALAR
       }
       ATTRIBUTE "description" {
          DATATYPE  H5T_STRING {
             STRSIZE H5T_VARIABLE;
             STRPAD H5T_STR_NULLTERM;
             CSET H5T_CSET_UTF8;
             CTYPE H5T_C_S1;
          }
          DATASPACE  SCALAR
       }
       ATTRIBUTE "has_B" {
          DATATYPE  H5T_STD_I64LE
          DATASPACE  SCALAR
       }
       ATTRIBUTE "nhydr" {
          DATATYPE  H5T_STD_I64LE
          DATASPACE  SCALAR
       }
       ATTRIBUTE "nx" {
          DATATYPE  H5T_STD_I64LE
          DATASPACE  SCALAR
       }
       ATTRIBUTE "ny" {
          DATATYPE  H5T_STD_I64LE
          DATASPACE  SCALAR
       }
       ATTRIBUTE "nz" {
          DATATYPE  H5T_STD_I64LE
          DATASPACE  SCALAR
       }
       DATASET "electron_density" {
          DATATYPE  H5T_IEEE_F64LE
          DATASPACE  SIMPLE { ( 1, 512, 512, 425 ) / ( H5S_UNLIMITED, 512, 512, 425 ) }
       }
       DATASET "hydrogen_populations" {
          DATATYPE  H5T_IEEE_F32LE
          DATASPACE  SIMPLE { ( 1, 6, 512, 512, 425 ) / ( H5S_UNLIMITED, 6, 512, 512, 425 ) }
       }
       DATASET "snapshot_number" {
          DATATYPE  H5T_STD_I32LE
          DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
       }
       DATASET "temperature" {
          DATATYPE  H5T_IEEE_F32LE
          DATASPACE  SIMPLE { ( 1, 512, 512, 425 ) / ( H5S_UNLIMITED, 512, 512, 425 ) }
       }
       DATASET "velocity_z" {
          DATATYPE  H5T_IEEE_F32LE
          DATASPACE  SIMPLE { ( 1, 512, 512, 425 ) / ( H5S_UNLIMITED, 512, 512, 425 ) }
       }
       DATASET "x" {
          DATATYPE  H5T_IEEE_F32LE
          DATASPACE  SIMPLE { ( 512 ) / ( 512 ) }
       }
       DATASET "y" {
          DATATYPE  H5T_IEEE_F32LE
          DATASPACE  SIMPLE { ( 512 ) / ( 512 ) }
       }
       DATASET "z" {
          DATATYPE  H5T_IEEE_F32LE
          DATASPACE  SIMPLE { ( 1, 425 ) / ( H5S_UNLIMITED, 425 ) }
       }
    }
    }

All the floating point variables can be either double or single
precision.

## Line lists and wavelength files

Other auxiliary files that can be used are line lists files and
wavelength files.

The line list files are used to include additional lines not included in
the different atoms. These lines will be treated in LTE. The line lists
are specified in the `kurucz.input` file (one per line), and have the
Kurucz line list format
([link](http://kurucz.harvard.edu/linelists.html)).

Just adding new transitions doesn't mean that they will be included in
the synthetic spectra. The extra lines will only be included in the
existing wavelength grid, which depends on the active atoms used. The
calculation of additional wavelengths can be forced by using a
wavelength file. This file is specified in `keyword.input` using the
keyword `WAVETABLE`. The format is a binary XDR file. Its contents are,
in order: the number of new wavelengths (1 XDR int), vacuum wavelength
values (XDR doubles).

[^1]: de la Cruz Rodríguez, J.; Piskunov, N. 2013, ApJ, 764, 33, [ADS
    link](http://adsabs.harvard.edu/abs/2013ApJ...764...33D).

[^2]: Leenaarts, J., Pereira, T. M. D., & Uitenbroek, H. 2012, A&A, 543,
    A109, [ADS
    link](http://adsabs.harvard.edu/abs/2012A%26A...543A.109L).
