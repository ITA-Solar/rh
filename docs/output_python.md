# Reading output in Python

The [helita](https://github.com/ITA-Solar/helita) package has a complete
python interface to read the output, input, and visualise files from RH
1.5D. The `helita` tools are described in detail in section
[helita](helita.md).

If `helita` is not available, the easiest and fastest way to read
the RH 1.5D output (or input) files in Python is via the
[xarray](http://xarray.pydata.org) package. `xarray` can load the output
files as a dataset directly, but in the case of the `output_aux.hdf5`
and `output_indata.hdf5` one needs to specify which group to read (see
above).

Here is a quick example on how to read some output from RH 1.5D with
`xarray`:

``` python
>>> import xarray
>>> ray = xarray.open_dataset("output_ray.hdf5")
>>> ray
<xarray.Dataset>
Dimensions:              (height: 82, wavelength: 902, wavelength_selected: 10, x: 1, y: 1)
Coordinates:
  * wavelength           (wavelength) float64 28.0 31.4 32.8 33.7 34.3 35.3 ...
  * wavelength_selected  (wavelength_selected) float64 85.1 276.4 278.5
  * x                    (x) float64 0.0
  * y                    (y) float64 0.0
Dimensions without coordinates: height
Data variables:
    Jlambda              (x, y, height, wavelength_selected) float64 ...
    chi                  (x, y, height, wavelength_selected) float64 ...
    intensity            (x, y, wavelength) float64 ...
    scattering           (x, y, height, wavelength_selected) float64 ...
    source_function      (x, y, height, wavelength_selected) float64 ...
    wavelength_indices   (wavelength_selected) int32 ...
Attributes:
    atmosID:              FALC_82_5x5.hdf5 (Wed Jan 10 15:29:28 2018)
    snapshot_number:      0
    rev_id:               001d537  Tiago Pereira  2018-01-10 12:34:07 +0100
    nx:                   1
    ny:                   1
    nz:                   82
    nwave:                902
    wavelength_selected:  3
    creation_time:        2018-01-10T16:16:42+0100
>>> aux = xarray.open_dataset("output_aux.hdf5", group="atom_MG")
>>> aux
<xarray.Dataset>
Dimensions:          (continuum: 10, height: 82, level: 11, line: 15, x: 1, y: 1)
Coordinates:
  * x                (x) float64 0.0
  * y                (y) float64 0.0
Dimensions without coordinates: continuum, height, level, line
Data variables:
    Rij_continuum    (continuum, x, y, height) float64 ...
    Rij_line         (line, x, y, height) float64 ...
    Rji_continuum    (continuum, x, y, height) float64 ...
    Rji_line         (line, x, y, height) float64 ...
    populations      (level, x, y, height) float64 ...
    populations_LTE  (level, x, y, height) float64 ...
Attributes:
    nlevel:      11
    nline:       15
    ncontinuum:  10
```