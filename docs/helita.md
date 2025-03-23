# `helita` interface

The [helita](https://github.com/ITA-Solar/helita) Python package
contains several routines to interface with RH 1.5D. Installation
instructions are [available in its
website](https://helita.readthedocs.io/en/latest/installation.html).

## Reading and writing input files

### Writing atmosphere files

The `rh15d` module in `helita.sim` contains a function to write an input
atmosphere in RH 1.5D format, assuming the user already has the required
data to write at hand. Its function definition is:

``` python
def make_xarray_atmos(outfile, T, vz, z, nH=None, x=None, y=None, Bz=None, By=None,
                      Bx=None, rho=None, ne=None, vx=None, vy=None, vturb=None,
                      desc=None, snap=None, boundary=None, append=False):
    """
    Creates HDF5 input file for RH 1.5D using xarray.

    Parameters
    ----------
    outfile : string
        Name of destination. If file exists it will be wiped.
    T : n-D array
        Temperature in K. Its shape will determine the output
        dimensions. Shape is generally (nt, nx, ny, nz), but any
        dimensions except nz can be omitted. Therefore the array can
        be 1D, 2D, or 3D, 4D but ultimately will always be saved as 4D.
    vz : n-D array
        Line of sight velocity in m/s. Same shape as T.
    z : n-D array
        Height in m. Can have same shape as T (different height scale
        for each column) or be only 1D (same height for all columns).
    nH : n-D array, optional
        Hydrogen populations in m^-3. Shape is (nt, nhydr, nx, ny, nz),
        where nt, nx, ny can be omitted but must be consistent with
        the shape of T. nhydr can be 1 (total number of protons) or
        more (level populations). If nH is not given, rho must be given!
    ne : n-D array, optional
        Electron density in m^-3. Same shape as T.
    rho : n-D array, optional
        Density in kg m^-3. Same shape as T. Only used if nH is not given.
    vx : n-D array, optional
        x velocity in m/s. Same shape as T. Not in use by RH 1.5D.
    vy : n-D array, optional
        y velocity in m/s. Same shape as T. Not in use by RH 1.5D.
    vturb : n-D array, optional
        Turbulent velocity (Microturbulence) in km/s. Not usually needed
        for MHD models, and should only be used when a depth dependent
        microturbulence is needed (constant microturbulence can be added
        in RH).
    Bx : n-D array, optional
        Magnetic field in x dimension, in Tesla. Same shape as T.
    By : n-D array, optional
        Magnetic field in y dimension, in Tesla. Same shape as T.
    Bz : n-D array, optional
        Magnetic field in z dimension, in Tesla. Same shape as T.
    x : 1-D array, optional
        Grid distances in m. Same shape as first index of T.
    y : 1-D array, optional
        Grid distances in m. Same shape as second index of T.
    x : 1-D array, optional
        Grid distances in m. Same shape as first index of T.
    snap : array-like, optional
        Snapshot number(s).
    desc : string, optional
        Description of file
    boundary : Tuple, optional
        Tuple with [bottom, top] boundary conditions. Options are:
        0: Zero, 1: Thermalised, 2: Reflective.
    append : boolean, optional
        If True, will append to existing file (if any).
    """
```

Note that while in this routine the writing of the hydrogen populations
is optional (they can be derived from the mass density, if available),
RH 1.5D does not support this yet.

!!! note

    The variables passed to `make_xarray_atmos` must be consistent with the
    height scale. The first height index must be the top of the atmosphere
    (closest to observer), and the height scale must be strictly decreasing.


### Reading atmosphere files

Once written, the input atmosphere files can be read in Python with
`xarray`, and do not require `helita`. For example:

``` python
>>> import xarray
>>> atmos = xarray.open_dataset('my_atmos.hdf5')
>>> atmos
<xarray.Dataset>
Dimensions:               (depth: 82, nhydr: 6, snapshot_number: 1, x: 5, y: 5)
Coordinates:
  * x                     (x) int64 0 1 2 3 4
  * y                     (y) int64 0 1 2 3 4
    z                     (snapshot_number, depth) float32 ...
  * snapshot_number       (snapshot_number) int32 0
Dimensions without coordinates: depth, nhydr
Data variables:
    temperature           (snapshot_number, x, y, depth) float32 ...
    velocity_z            (snapshot_number, x, y, depth) float32 ...
    electron_density      (snapshot_number, x, y, depth) float64 ...
    hydrogen_populations  (snapshot_number, nhydr, x, y, depth) float32 ...
    velocity_turbulent    (snapshot_number, x, y, depth) float32 ...
Attributes:
    comment:          Created with make_xarray_atmos on 2018-01-25 15:28:10.4...
    boundary_top:     0
    boundary_bottom:  1
    has_B:            0
    description:      FAL C model with 82 depth points replicated to 5x5 colu...
    nx:               5
    ny:               5
    nz:               82
    nt:               1
```

The amount of detail loaded by `xarray` will depend how the atmosphere
was written. Older atmosphere files may not have as much verbose
attributes or labeled coordinates (especially if written by plain HDF5
with no attaching of dimension scales), but they are still valid. Older
netCDF atmospheres should work fine with `xarray`.

It is also possible to modify the data with `xarray`, and saving and
updated atmosphere is done via the `to_netcdf()` method:

``` python
>>> atmos.to_netcdf("newfile.hdf5", format='NETCDF4')
```

Be sure to use `format='NETCDF4'` so that the file is internally HDF5!

### Writing wavelength files

Another utility function in `rh15d.py` is `make_wave_file`. This creates
an RH wavelength file (to be used with the option `WAVETABLE` in
`keyword.input`) that contains additional wavelengths to be calculated.
The function's usage is documented in its function call:

``` python
def make_wave_file(outfile, start=None, end=None, step=None, new_wave=None,
                   ewave=None, air=True):
   """
   Writes RH wave file (in xdr format). All wavelengths should be in nm.

   Parameters
   ----------
   start: number
       Starting wavelength.
   end: number
       Ending wavelength (non-inclusive)
   step: number
       Wavelength separation
   outfile: string
       Name of file to write.
   ewave: 1-D array, optional
       Array of existing wavelengths. Program will make discard points
       to make sure no step is enforced using these points too.
   air: boolean, optional
       If true, will at the end convert the wavelengths into vacuum
       wavelengths.
   """
```

You can either supply an array with the wavelengths, or give a range of
wavelengths and a fixed spacing, e.g.:

``` python
>>> from helita.sim import rh15d
>>> rr = rh15d.Rh15dout()
# this will write wavelenghts from 650 to 650 nm, 0.01 nm spacing:
>>> rh15d.make_wave_file('my.wave', 650, 660, 0.01)
# this will write an existing array "my_waves", if it exists
>>> rh15d.make_wave_file('my.wave', ewave=my_waves)
```

## Reading output files

The main class to read the output is called `Rh15dout`. It uses `xarray`
under the hood and populates an object with all the different datasets.
It can be initiated in the following way:

``` python
>>> from helita.sim import rh15d
>>> rr = rh15d.Rh15dout()
--- Read ./output_aux.hdf5 file.
--- Read ./output_indata.hdf5 file.
--- Read ./output_ray.hdf5 file.
```

By default, it will look for the three files in the directory specified
as main argument (defaults to current directory). Additionally, the
method `read_group(infile)` can be used to manually load the
`output_aux.hdf5` or `output_indata.hdf5` and the method and
`read_ray(infile)` can be used to manually load the `output_ray.hdf5`
file. The variables themselves are not read into memory, but are rather
a [memmap
object](http://docs.scipy.org/doc/numpy/reference/generated/numpy.memmap.html)
(file pointer; only read when needed) that `xarray` opens.

After loading the files, the `Rh15dout` instance loads each file as an
`xarray` dataset with the base name of each group (e.g. `ray`, `atmos`,
`atom_CA`, `mpi`).The `ray` attribute contains the same dataset as shown
in the `xarray` example above.

The attributes of each file are still accessible under the attributes of
each object, e.g.:

``` python
>>> rr.ray.creation_time
'2018-01-10T16:16:42+0100'
>>> rr.atmos.nrays
5
>>> rr.mpi.nprocesses
2048
```

With `xarray` it is easy to quickly inspect and plot different
quantities. For example, to plot the intensity at `(x, y) = (0, 0)`:

``` python
>>> rr.ray.intensity[0, 0].plot()
```

Or the intensity at a fixed wavelength:

``` py
>>> rr.ray.intensity.sel(wavelength=279.55, method='nearest').plot()
```

(This only shows a 2D image if you calculated the intensity from a 3D
model, otherwise an histogram or line plot is shown.)

## Visualisation and notebooks

`helita` includes a visualisation module, `helita.sim.rh15d_vis`, with
widgets that are meant to be used inside the [Jupyter
notebook](https://jupyter.org/). To use these, you will need to install
not only `helita` but also the [Matplotlib Jupyter
Extension](https://github.com/matplotlib/jupyter-matplotlib) and the
[IPython widgets for
Jupyter](https://github.com/jupyter-widgets/ipywidgets). If you have
Anaconda, both can be installed with conda:

    conda install -c conda-forge ipywidgets ipympl widgetsnbextension
    jupyter nbextension enable --py widgetsnbextension

You can also install them with `pip` (check their pages for details).

Currently we have the following Jupyter notebooks for visualisation of
RH 1.5D output:

-   [Basic
    output](https://github.com/ITA-Solar/rh/blob/master/doc/notebooks/BasicOutput.ipynb)
-   [Visualisation
    widgets](https://github.com/ITA-Solar/rh/blob/master/doc/notebooks/Visualisation.ipynb)

To use the above notebooks, you need to have run RH 1.5D and have the
output files ready!

You can also explore the input atmosphere files with the Jupyter widget
`rh15d_vis.InputAtmosphere` in `helita`:

``` py
>>> from helita.sim import rh15d_vis
>>> rh15d_vis.InputAtmosphere('my_atmos.hdf5')
```
