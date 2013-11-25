# Set of programs and tools to read the outputs from RH 1.5D
import os
import numpy as np
import netCDF4 as nc


class Rh15dout:
    """
    Main class for reading RH 1.5D output.

    Parameters
    ----------
    fdir : string
        Directory where to look for files.
    verbose: bool
        Will set self.verbose to allow extra output.

    Attributes
    ----------
    files : list
        Files that are currently open.
    params : dictionary
        Parameters read from output_indata.ncdf attributes
    fdir : string
        File directory
    verbose : bool
        If True, will display more information.
    atmos : DataHolder instance
        sub-class with data for output_atmos.ncdf (if loaded)
    ray : DataHolder instance
        sub-class with data for output_ray.ncdf (if loaded)
    spectrum : DataHolder instance
        sub-class with data for output_spectrum.ncdf (if loaded)

    Methods
    -------
    read_aux(infile=None)
        Reads Aux file.
    read_indata(infile=None)
        Reads indata file.
    read_spectrum(infile=None)  **DEPRECATED**
        Reads spectrum file.
    read_ray(infile=None)
        Reads ray file.
    read_J(infile='scratch/J.dat.ncdf')  **DEPRECATED**
        Reads angle averaged intensity file
    close()
        Closes the open NetCDF files.
    """
    def __init__(self, fdir='.', verbose=True):
        self.files = []
        self.params = {}
        self.verbose = verbose
        self.fdir = fdir
        if os.path.isfile('%s/output_aux.ncdf' % self.fdir):
            self.read_aux()
        if os.path.isfile('%s/output_indata.ncdf' % self.fdir):
            self.read_indata()
        if os.path.isfile('%s/output_spectrum.ncdf' % self.fdir):
            self.read_spectrum()
        if os.path.isfile('%s/output_ray.ncdf' % self.fdir):
            self.read_ray()

    def read_aux(self, infile=None):
        ''' Reads Aux file. '''
        if infile is None:
            infile = '%s/output_aux.ncdf' % self.fdir
        self.files.append(read_ncdf(self, infile))
        if self.verbose:
            print('--- Read %s file.' % infile)
        return

    def read_indata(self, infile=None):
        ''' Reads indata file. '''
        if infile is None:
            infile = '%s/output_indata.ncdf' % self.fdir
        self.files.append(read_ncdf(self, infile))
        if self.verbose:
            print('--- Read %s file.' % infile)
        return

    def read_spectrum(self, infile=None):
        ''' Reads spectrum file. '''
        if infile is None:
            infile = '%s/output_spectrum.ncdf' % self.fdir
        self.spectrum = DataHolder()
        self.files.append(read_ncdf(self.spectrum, infile))
        if self.verbose:
            print('--- Read %s file.' % infile)
        return

    def read_ray(self, infile=None):
        ''' Reads ray file. '''
        if infile is None:
            infile = '%s/output_ray.ncdf' % self.fdir
        self.ray = DataHolder()
        self.files.append(read_ncdf(self.ray, infile))
        if self.verbose:
            print('--- Read %s file.' % infile)
        return

    def read_J(self, infile='scratch/J.dat.ncdf'):
        ''' Reads angle averaged intensity file '''
        self.files.append(read_ncdf(self, infile))
        if self.verbose:
            print('--- Read %s file.' % infile)
        return

    def close(self):
        ''' Closes the open NetCDF files '''
        for f in self.files:
            f.close()


class NcdfAtmos:
    """
    Class to hold information from an RH 1.5D atmosphere file. Uses
    DataHolder and read_ncdf()

    Attributes
    ----------
    file : string
        Filename read
    closed : bool
        False if file is open, True if closed.

    Methods
    -------
    close()
        Closes the netcdf file.
    read()
        Reads the netcdf file (ran on __init__())
    """
    def __init__(self, infile):
        self.file = read_ncdf(self, infile)
        self.closed = False

    def close(self):
        try:
            self.file.close()
            self.closed = True
        except RuntimeError:
            print('(WWW) NcdfAtmos: input file already closed.')

    def read(self, infile):
        if not self.closed:
            self.close()
        self.file = read_ncdf(self, infile)
        return

#############################################################################
###   TOOLS                                                               ###
#############################################################################


class DataHolder:
    def __init__(self):
        pass


def read_ncdf(inclass, infile):
    """
    Reads NetCDF file into inclass, instance of any class.
    Variables are read into class attributes, dimensions and attributes
    are read into params dictionary.
    """
    # internal attributes of NetCDF groups
    ncdf_internals = ['__class__', '__delattr__', '__doc__', '__format__',
        '__getattr__', '__getattribute__', '__hash__', '__init__', '__new__',
        '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__',
        '__str__', '__subclasshook__', '_enddef', '_grpid', '_redef', 'close',
        'cmptypes', 'createCompoundType', 'createDimension', 'createGroup',
        'createVLType', 'createVariable', 'delncattr', 'dimensions',
        'file_format', 'getncattr', 'groups', 'maskanscale', 'ncattrs',
        'parent', 'path', 'renameDimension', 'renameVariable', 'set_fill_off',
        'set_fill_on', 'setncattr', 'sync', 'variables', 'vltypes']
    if not os.path.isfile(infile):
        raise IOError('read_ncdf: File %s not found' % infile)
    f = nc.Dataset(infile, mode='r')
    if 'params' not in dir(inclass):
        inclass.params = {}
    # add dimensions as attributes
    for d in f.dimensions.keys():
        inclass.params[d] = len(f.dimensions[d])
    # add attributes
    attrs = [a for a in dir(f) if a not in ncdf_internals]
    for att in attrs:
        inclass.params[att] = getattr(f, att)
    # add variables
    for v in f.variables.keys():
        vname = v.replace(' ', '_')     # sanitise string for spaces
        setattr(inclass, vname, f.variables[v])
    # Now do the same for all groups
    for group in f.groups.keys():
        gname = group.replace(' ', '_')  # sanitise string for spaces
        setattr(inclass, gname, DataHolder())
        cur_group = f.groups[group]
        cur_class = getattr(inclass, gname)
        # add variables
        for v in cur_group.variables.keys():
            vname = v.replace(' ', '_')  # sanitise string for spaces
            setattr(cur_class, vname, cur_group.variables[v])
        # add dimensions as attributes
        for d in cur_group.dimensions.keys():
            inclass.params[d] = len(cur_group.dimensions[d])
        # add attributes
        attrs = [a for a in dir(cur_group) if a not in ncdf_internals]
        for att in attrs:
            inclass.params[att] = getattr(cur_group, att)
    return f


def make_ncdf_atmos(outfile, T, vz, ne, nH, z, x=None, y=None, Bz=None,
                    By=None, Bx=None, desc=None, comp=False, append=True,
                    snap=None, complev=2):
    """
    Creates NetCDF input file for rh15d.

    Parameters
    ----------
    outfile :  string
        Name of destination. If file exists it will be wiped.
    T : nD array
        Temperature. Its shape will determine the dimensionality written
    vz : nD array
        Vertical velocity. Same shape as T. In m/s.
    ne : nD array
        Electron density. Same shape as T. In m-3.
    nH : nD array
        Hydrogen populations. Shape [6, shape.T]. In m-3.
    z : 1D array
        Height scale. Same shape as last index of T. In m.
    x : 1D array
        Same shape as first index of T. In m.
    y : 1D array
        Same shape as second index of T. In m.
    snap : array-like
        Snapshot number(s)
    Bx, By, Bz : nD arrays
        Magnetic field components. Same shape as T. In T.
    append : bool
        If True, will append to existing file (if any).
    comp : bool
        If True, compress file.
    complev : integer
        Compression level.
    """
    import os
    mode = ['w', 'a']
    if (append and not os.path.isfile(outfile)):
        append = False
    rootgrp = nc.Dataset(outfile, mode[append], format='NETCDF4')
    complev = 2
    nt = 1
    if nH.shape == T.shape:
        nhydr = 1
    else:
        nhydr = nH.shape[0]
    if len(T.shape) == 1:
        x = [0.]
        y = [0.]
        # for these, only single snapshot is supported
        T = T[np.newaxis, np.newaxis, np.newaxis, :]
        ne = ne[np.newaxis, np.newaxis, np.newaxis, :]
        nH = nH[np.newaxis, np.newaxis, np.newaxis, :]
        vz = vz[np.newaxis, np.newaxis, np.newaxis, :]
        z = z[np.newaxis, np.newaxis, np.newaxis, :]
        if Bz is not None:
            Bx = Bx[np.newaxis, np.newaxis, np.newaxis, :]
            By = By[np.newaxis, np.newaxis, np.newaxis, :]
            Bz = Bz[np.newaxis, np.newaxis, np.newaxis, :]
    if len(T.shape) == 3:  # single snapshot
        T = T[np.newaxis, :]
        ne = ne[np.newaxis, :]
        nH = nH[np.newaxis, :]
        vz = vz[np.newaxis, :]
        z = z[np.newaxis, :]
        if Bz is not None:
            Bx = Bx[np.newaxis, :]
            By = By[np.newaxis, :]
            Bz = Bz[np.newaxis, :]
    elif len(T.shape) != 4:
        raise ValueError('Invalid shape for T')
    else:
        nt = T.shape[0]
        nx = T.shape[1]
        ny = T.shape[2]
    if snap is None:
        snap = np.arange(nt, dtype='i4')

    # for a new file, create dimensions and variables
    if not append:
        rootgrp.createDimension('nt', None)  # create unlimited dimension
        rootgrp.createDimension('nx', T.shape[-3])
        rootgrp.createDimension('ny', T.shape[-2])
        rootgrp.createDimension('nz', T.shape[-1])
        rootgrp.createDimension('nhydr', nhydr)
        T_var = rootgrp.createVariable('temperature', 'f4',
                                       ('nt', 'nx', 'ny', 'nz'), zlib=comp,
                                        least_significant_digit=1,
                                        complevel=complev)
        vz_var = rootgrp.createVariable('velocity_z', 'f4',
                                        ('nt', 'nx', 'ny', 'nz'), zlib=comp,
                                        least_significant_digit=1,
                                        complevel=complev)
        ne_var = rootgrp.createVariable('electron_density', 'f8',
                                        ('nt', 'nx', 'ny', 'nz'), zlib=comp,
                                        complevel=complev)
        nh_var = rootgrp.createVariable('hydrogen_populations', 'f4',
                                        ('nt', 'nhydr', 'nx', 'ny', 'nz'),
                                        zlib=comp, complevel=complev)
        x_var = rootgrp.createVariable('x', 'f4', ('nx',))
        y_var = rootgrp.createVariable('y', 'f4', ('ny',))
        z_var = rootgrp.createVariable('z', 'f4', ('nt', 'nz'))
        nt_var = rootgrp.createVariable('snapshot_number', 'i4', ('nt',))
        if Bz is not None:
            bx_var = rootgrp.createVariable('B_x', 'f4',
                                            ('nt', 'nx', 'ny', 'nz'),
                                            zlib=comp, complevel=complev,
                                            least_significant_digit=5)
            by_var = rootgrp.createVariable('B_y', 'f4',
                                            ('nt', 'nx', 'ny', 'nz'),
                                            zlib=comp, complevel=complev,
                                            least_significant_digit=5)
            bz_var = rootgrp.createVariable('B_z', 'f4',
                                            ('nt', 'nx', 'ny', 'nz'),
                                            zlib=comp, complevel=complev,
                                            least_significant_digit=5)
        if desc is None:
            rootgrp.description = "N/A."
        else:
            rootgrp.description = desc
        if Bz is None:
            rootgrp.has_B = 0
        else:
            rootgrp.has_B = 1
        nt = [0, nt]
    else:
        # get variables
        T_var = rootgrp.variables['temperature']
        vz_var = rootgrp.variables['velocity_z']
        ne_var = rootgrp.variables['electron_density']
        nh_var = rootgrp.variables['hydrogen_populations']
        nt_var = rootgrp.variables['snapshot_number']
        x_var = rootgrp.variables['x']
        y_var = rootgrp.variables['y']
        z_var = rootgrp.variables['z']
        if Bz is not None:
            bx_var = rootgrp.variables['B_x']
            by_var = rootgrp.variables['B_y']
            bz_var = rootgrp.variables['B_z']
        nti = len(rootgrp.dimensions['nt'])
        nt = [nti, nti + nt]

    T_var[nt[0]:nt[1]] = T
    vz_var[nt[0]:nt[1]] = vz
    ne_var[nt[0]:nt[1]] = ne
    nh_var[nt[0]:nt[1], :nhydr] = nH
    if Bz is not None:
        bx_var[nt[0]:nt[1]] = Bx
        by_var[nt[0]:nt[1]] = By
        bz_var[nt[0]:nt[1]] = Bz
    x_var[:] = x
    y_var[:] = y
    z_var[nt[0]:nt[1]] = z
    nt_var[nt[0]:nt[1]] = snap
    rootgrp.close()
    return


def waveconv(wave, mode='vac2air', gravred=False, verbose=False):
    """
    Converts given wavelength from air to vacuum or vacuum to air.
    Can also account for the Solar gravitational redshift.
    Uses the formula from P. Ciddor, Applied Optics vol 35, no 9, 1566 (1996)
    20070306: Coded --Tiago

    Parameters
    ----------
    wave : array-like
        Wavelength values in nm.
    mode : string
        Which option to use, either 'vac2air' or 'air2vac'.
    gravred : bool
        If True, will compensate for solar gravitational redshift.
    verbose: bool
        If True, prints more information.

    Returns
    -------
    result : ndarray
        Converted wavelengths.
    """
    # All these constants are in um^-2
    k = np.array([238.0185, 5792105., 57.362, 167917.], dtype='d')
    wn = np.array(1 / (1.e-3 * wave), dtype='d')  # to wave number in um^-1
    # Index of refraction of dry air (15deg C, 101325 Pa, 450 ppm CO2)
    n = 1. + 1.e-8 * (k[1] / (k[0] - wn ** 2) + k[3] / (k[2] - wn ** 2))
    if verbose:
        print '*** Index of refraction is %s' % n
    if mode == 'air2vac':
        result = wave * n
    elif mode == 'vac2air':
        result = wave / n
    else:
        print '(EEE) waveconv: Mode not valid, exiting...'
        return 0
    # Account for solar gravitational redshift?
    if gravred:
        # Multiplying by the 636m/s velocity shift / vac. speed of light
        result += result * 636. / 299792458.
    return result


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
    import xdrlib
    if new_wave is None:
        new_wave = np.arange(start, end, step)
        if None in [start, end, step]:
            raise ValueError(('Must specify either new_wave, or start, end, '
                              'step. Stopping.'))
    if step is None:
        step = np.median(np.diff(new_wave))
    if ewave is not None:  # ensure step is kept at most times
        keepers = []
        for w in new_wave:
            if np.min(np.abs(w - ewave)) > step * 0.375:
                keepers.append(w)
        new_wave = np.array(keepers)
    if air:
        new_wave = waveconv(new_wave, mode='air2vac')
    # write file
    p = xdrlib.Packer()
    nw = len(new_wave)
    p.pack_int(nw)
    p.pack_farray(nw, new_wave.astype('d'), p.pack_double)
    f = open(outfile, 'wb')
    f.write(p.get_buffer())
    f.close()
    print("Wrote %i wavelengths to file." % nw)
    return
