# Output file structure 


## Overview

Output is written into three files: `output_aux.hdf5`,
`output_indata.hdf5`, and `output_ray.hdf5`. This is a big departure
from RH, which contained several more output files. In particular, RH
1.5D will not write all the information that was written by RH, due to
the sheer size it would take for large 3D simulations. The files are
written in the machine-independent, self-describing HDF5 format. The
contents of the files are organised in groups, variables, and
attributes. Groups and variables can be imagined as directories and
files in a filesystem. Inside groups, different variables and dimensions
can be organised. The content of the output files can vary: some runs
will have more detailed information and therefore more variables written
to the files.

HDF5 is an open, platform-independent format, and therefore interfaces
to many programming languages are available. The [main interface
libraries](https://www.hdfgroup.org/HDF5/) are available in C, C++,
Fortran, and Java. But there are also interfaces for Python
([h5py](http://www.h5py.org/)),
[Julia](https://github.com/JuliaIO/HDF5.jl), IDL (from version 6.2),
MATLAB , Octave,
[Perl](http://search.cpan.org/~chm/PDL-IO-HDF5-0.6501/hdf5.pd), and
[R](https://cran.r-project.org/package=h5).

The RH 1.5D output format is standard HDF5 but it is also compatible
with NetCDF 4 readers: in most cases one needs to specify only the
variable or group name to read the data. The HDF5 and NetCDF libraries
provide useful command line tools, which can be used to gather
information about the RH 1.5D files or extract data. Additionally, there
is a more complete set of tools written in Python to read and analyse
these files.

!!! warning

    Because of the limitations of different languages, not all interfaces
    support all HDF5 features. Some libraries (e.g. IDL or plain h5py in
    Python) will not detect missing data in arrays (written with the fill
    value). In such cases, reading variables with missing data, the data are read with no
    warning or indication of those that have special fill values.


The structure of the three output files is given below.

!!! info "Missing data"

    When a column fails to converge, output for that column is not written.
    This means that the variables that depend on `(nx, ny)` will have some
    values missing. HDF5 marks these values as missing data
    and uses a fill value (of 9.9692e+36). When the `15D_DEPTH_ZCUT` option
    is used, not all heights will be used in the calculation. The code does
    not read the skipped parts of the atmosphere. When writing such
    variables of `nz`, only the points that were used are written to the
    file, and the rest will be marked as missing data (typically the z cut
    height varies with the column).


## `output_aux.hdf5`

This file contains the level populations and radiative rates. For each
active atom or molecule, it contains different groups called `atom_XX`
or `molecule_XX`, where `XX` is the identifier for the species (e.g.
`MG`, `CO`).

!!! note

    The atmosphere dimensions on many of the output files are not
    necessarily the same as in the atmosphere file. They depend on the
    number of columns calculated, which are a function of
    `X/Y_START/END/STEP`.


It has the following global attributes:

| | |
-----------|:--------------------------------------
  `atmosID`|   Identifier for the atmosphere file.
  `rev_id` |   Revision identifier.
  `nx`     |   Number of points in x dimension
  `ny`     |   Number of points in y dimension
  `nz`     |   Number of points in height dimension


Inside each of the atom/molecule groups, the following dimensions can
exist:

  | Name |                 Description|
|:-------------------|:--------------------------------------
  `x`               |   Horizontal x dimension.
  `y`               |   Horizontal y dimension.
  `height`          |   Vertical dimension.
  `level`           |   Number of atomic levels.
  `line`            |   Number of atomic transitions
  `continuum`       |   Number of bound-free transitions.
  `vibration_level` |   Number of molecule vibration levels.
  `molecular_line`  |   Number of molecular lines.
  `rotational_state` |  Number of rotational states.

The atom groups can contain the following optional variables:


|   <div style="width:120px">Name</div>             |   Dimensions   |         Description |
:------------------|:-------------------------------|:----------------------------
  `populations`     |  `(level, x, y, height)`      |    Atomic populations.
  `populations_LTE` |  `(level, x, y, height)`      |    Atomic LTE populations.
  `Rij_line`         | `(line, x, y, height)`        |   Radiative rates out of the line.
  `Rji_line`         | `(line, x, y, height)`        |   Radiative rates into the line.
  `Rij_continuum`   |  `(continuum, x, y, height)`   |   Radiative rates out of the bf transition.
  `Rji_continuum`   |  `(continuum, x, y, height)`   |   Radiative rates into the bf transition.
  `collision_rates` |  `(level, level, x, y, height)` |  Collisional rates.  Convention: first index is upper level.

The molecule groups can contain the following optional variables:



|   <div style="width:120px">Name</div>             |   Dimensions   |         Description |
:------------------|:-------------------------------|:----------------------------
  `populations`    |   `(vibration_level, x, y, height)` |   Molecular populations.
  `populations_LTE` |  `(vibration_level, x, y, height)` |  Molecular LTE populations.


All units are SI.

!!! note

    In older versions it was possible to specify the keyword
    `15D_WRITE_EXTRA` and get additional output written to `output_aux.hdf5`
    (e.g. a new `opacity` group and more rates). While the procedures are
    still in `writeAux_p.c`, the functionality is deprecated because other
    changes in the code were not compatible with this way of writing the
    output. It is possible that this functionality will return at a later
    version.


## `output_indata.hdf5`

This file contains data and metadata related to the run. It contains
three groups: `input` (mostly settings from `keyword.input`), `atmos`
(atmospheric variables), and `mpi` (several variables relating to the
run).

It has the following global attributes:

| | |
-----------|:--------------------------------------
  `atmosID`|   Identifier for the atmosphere file.
  `rev_id` |   Revision identifier.
  `nx`     |   Number of points in x dimension
  `ny`     |   Number of points in y dimension
  `nz`     |   Number of points in height dimension


The `input` group contains all the input files (except atmosphere and
molecular data), and a few attributes that are options from
`keyword.input`. It contains the following string variables:


|   <div style="width:170px">Variable</div>            |         Description |
:------------------|:-------------------------------
  `atom_groups`     |        Array with names of atom groups. Other is the same as atom order in `atoms.input`
  `atoms_file_contents` |    Contents of `atoms.input` saved into a string
  `keyword_file_contents` |  Contents of `keyword.input` saved into a string
  `kurucz_file contents`  |  Contents of `kurucz.input` saved into a string, only if used


The `input` group also has other groups inside. If Kurucz line lists are
used, it contains groups called `Kurucz_line_file0`, ...,
`Kurucz_line_fileN`, where N-1 is the total number of line list files.
The other groups are all atom files (`PASSIVE` and `ACTIVE`), and they
take the names of `atom_XX`, where `XX` is the element name (for a list
of these, see the variable `atom_groups` above). Inside all of these
groups (Kurucz and atom) there is one variable, called `file_contents`,
which contains the file saved intro a string and an attribute, called
`file_name`, which contains the file name and path. These input options
and files are read instead of the original files when doing a rerun.

The `atmos` groups contains the dimensions `x`, `y`, `height`, `element`
and `ray`. It also contains the following variables:


|  <div style="width:140px">Name</div>  |    Dimensions       |    Units   |   Description |
:--------------------|:---------------------|:----------|:-------------------------
  `temperature`      |   `(x, y, height)`  |  K        |   Temperatures
  `velocity_z`       |   `(x, y, height)`  |  m s^-1^  |   Vertical velocities
  `electron_density` |   `(x, y, height)`  |  m^-3^  | Electron   densities         
  `height_scale`     |   `(x, y, height)`  |  m      |     Height scale used. Can be different for every column when depth refine is used.
  `element_weight`    |  `(element)`       |  a.m.u.  |    Atomic weights
  `element_abundance` |  `(element)`       |          |    Element abundances relative to hydrogen.
  `muz`               |  `(ray)`           |          |    mu values for each ray.
  `muz`               |  `(ray)`           |          |    mu weights for each ray.
  `x`                 |  `(x)`             |  m       |    Spatial coordinates along x axis.
  `y`                 |  `(y)`             |  m       |    Spatial coordinates along y axis.


!!! note

    When `15D_DEPTH_REFINE` is used, each column will have a different
    (optimised) height scale, but they all have the same number of depth
    points (`nz`). In these cases, it is very important to save the `height`
    variable because otherwise one does not know how to relate the height
    relations of quantities from different columns.


The `atmos` group also contains the following attributes:

| | |
-----------|:--------------------------------------
  `moving` |  Unsigned int, 1 if velocity fields present.
  `stokes` |  Unsigned int, 1 if stokes output present.


The `mpi` group contains the dimensions `x`, `y`, and `iteration`
(maximum number of iterations).

!!! warning

    `iteration` is needs to be hardcoded for reruns. The default maximum value is 1500. If your `N_MAX_ITER` is larger than this, `iteration` will be adjusted. Reruns can never have a value larger than the maximum of 1500 or `N_MAX_ITER` (of the first run).


The `mpi` group also contains several variables:


|   <div style="width:120px">Name</div>             |   Dimensions   |         Description |
:------------------|:-------------------------------|:----------------------------
  `xnum`           |     `(x)`              |   Indices of x positions calculated.
  `xnum`            |    `(x)`              |   Indices of x positions calculated.
  `task_map`         |   `(x, y)`           |   Maps which process ran which column.
  `task_map_number`  |   `(x, y)`           |   Maps which task number each column was.
  `iterations`       |   `(x, y)`           |   Number of iterations used for each column.
  `convergence`       |  `(x, y)`           |   Indicates if each column converged or not. Possible values are `1` (converged), `0` (non converged),  or `-1` (crashed).
  `delta_max`         |  `(x, y)`            |  Final value for `delta_max` when iteration finished.
  `delta_max_history` |  `(x, y, iteration)` |  Evolution of `delta_max`
  `z_cut`             |  `(x, y)`        |      Height index of the temperature cut.


The `mpi` group also contains the following attributes: `x_start`,
`x_end`, `x_step`, `y_start`, `y_end`, and `y_step`, all of which are
options from `keyword.input`.

## `output_ray.hdf5`

This file contains the synthetic spectra and can also contain extra
information such as opacities and the source function. It contains only
the root group. Its dimensions are `x`, `y`, `wavelength`, and
eventually `wavelength_selected` and `height`. The latter two are only
present when `ray.input` specifies more than `0` wavelengths for
detailed output, and it matches `Nsource`, the number of those
wavelengths entered in `ray.input`.

It can contain the following variables:

|  <div style="width:140px">Name</div>  |    <div style="width:100px">Dimensions</div>       |    Units   |   Description |
|:--------------------|:---------------------|:----------|:-------------------------|
| `wavelength` | `(    wavelength)`      | nm                    |  Wavelength     scale.    |
| `intensity`  | `(x, y,wavelength)`      | W m^-2^ Hz^-1^ sr^-1^ | Synthetic  disk-centre    intensity (Stokes I)    |
| `stokes_Q`   | `(x, y, wavelength)`     | W m^-2^ Hz^-1^ sr^-1^ | Stokes Q.  **Optional.**       |
| `stokes_U`   | `(x, y, wavelength)`     | W m^-2^ Hz^-1^ sr^-1^ | Stokes U.  **Optional.**       |
| `stokes_V`   | `(x, y, wavelength)`     | W m^-2^ Hz^-1^ sr^-1^ | Stokes V.  **Optional.**       |
| `tau_one_height`   | `(x, y,wavelength)`  | m                     | Height where  optical depth reaches unity, for each column. **Optional.**      |
| `wavelength_selected`    | `(wavelength_selected)`  |     | Wavelength scale for the detailed  output variables   below.  **Optional.**      |
| `wavelength_indices`    | `(wavelength_selected)`  |     | Indices of wavelengths selected for variables below. **Optional.**     
| `chi`  | `(x, y, height, wavelength_selected)`           | m^-1^                 | Total opacity (line and continuum) **Optional.**    |
| `source_function`  | `(x, y, height, wavelength_selected)`           |W m^-2^ Hz^-1^ sr^-1^                 | Total source function (line and continuum) **Optional.**    |
| `Jlambda`  | `(x, y, height, wavelength_selected)`           |W m^-2^ Hz^-1^ sr^-1^                 | Angle-averaged radiation field **Optional.**    |
| `scattering`  | `(x, y, height, wavelength_selected)`           |                | Scattering term multiplied by Jlambda **Optional.**    |


The `wavelength` is in nm, air or vacuum units, depending if
`VACUUM_TO_AIR` is `TRUE` or `FALSE` (in `keyword.input`). `chi` is in
m^-1^ and `tau_one_height` in m.

Despite internally being calculated in double precision, all the output
(except the wavelength scale) is written in single precision to save
disk space.

The full Stokes vector is only written when in `keyword.input`
`STOKES_MODE` is not `NO_STOKES` and the `STOKES_INPUT` is set.

The `chi`, `source_function`, and `Jlambda` variables depend on the 3D
grid and on wavelength. Therefore, for even moderate grid sizes they can
take huge amounts of space. If `nx = ny = nz = 512` and
`wavelength_selected = 200`, each of these variables will need 100Gb of
disk space. For a simulation with a cubic grid of 1024^3^ points and
saving the full output for 1000 wavelength points, `output_ray.hdf5`
will occupy a whopping 12Tb per snapshot of disk space. To avoid such
problems, these large arrays are only written when `ray.input` contains
`Nsource > 0`, and for the wavelengths selected.

The `output_ray.hdf5` file contains the following global attributes:

| | |
-----------|:--------------------------------------
  `atmosID`             |  Identifier for the atmosphere file
  `snapshot_number`     |  Number of simulation snapshot (from atmosphere file)
  `rev_id`              |  Revision identifier
  `nx`                  |  Number of points in x dimension
  `ny`                  |  Number of points in y dimension
  `nz`                  |  Number of points in height dimension
  `nwave`               |  Number of wavelength points
  `wavelength_selected` |  Number of wavelength points selected for detailed output
  `creation_time`       |  Local time when file was created

