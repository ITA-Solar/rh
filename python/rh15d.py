"""
Set of programs and tools to read the outputs from RH, 1.5D version
"""
import os
import numpy as np



class Rh15dout:
    def __init__(self, fdir='.', verbose=True):
        self.files = []
        self.params = {}
        self.verbose = verbose
        self.fdir = fdir
        OUTFILE_FUNC = {"output_aux": self.read_aux,
                        "output_indata": self.read_indata,
                        "output_ray": self.read_ray,
                        "output_spectrum": self.read_spectrum}
        for outfile, func in OUTFILE_FUNC.items():
            if (os.path.isfile("%s/%s.ncdf" % (self.fdir, outfile)) or
                os.path.isfile("%s/%s.hdf5" % (self.fdir, outfile))):
                func()

    def read_aux(self, infile=None):
        ''' Reads Aux file. '''
        if infile is None:
            infile = '%s/output_aux.hdf5' % self.fdir
        self.files.append(read_hdf5(self, infile))
        if self.verbose:
            print(('--- Read %s file.' % infile))

    def read_indata(self, infile=None):
        ''' Reads indata file. '''
        if infile is None:
            infile = '%s/output_indata.hdf5' % self.fdir
        self.files.append(read_hdf5(self, infile))
        if self.verbose:
            print(('--- Read %s file.' % infile))

    def read_spectrum(self, infile=None):
        ''' Reads spectrum file. '''
        if infile is None:
            infile = '%s/output_spectrum.ncdf' % self.fdir
        self.spectrum = DataHolder()
        self.files.append(read_hdf5(self.spectrum, infile))
        if self.verbose:
            print(('--- Read %s file.' % infile))

    def read_ray(self, infile=None):
        ''' Reads ray file. '''
        if infile is None:
            infile = '%s/output_ray.hdf5' % self.fdir
        self.ray = DataHolder()
        self.files.append(read_hdf5(self.ray, infile))
        if self.verbose:
            print(('--- Read %s file.' % infile))

    def read_J(self, infile='scratch/J.dat.hdf5'):
        ''' Reads angle averaged intensity file '''
        self.files.append(read_hdf5(self, infile))
        if self.verbose:
            print(('--- Read %s file.' % infile))

    def close(self):
        ''' Closes the open NetCDF files '''
        for f in self.files:
            f.close()

class HDF5Atmos:
    def __init__(self, infile):
        self.file = read_hdf5(self, infile)
        self.closed = False

    def close(self):
        try:
            self.file.close()
            self.closed = True
        except RuntimeError:
            print('(WWW) HDF5Atmos: input file already closed.')

    def read(self, infile):
        if not self.closed:
            self.close()
        self.file = read_hdf5(self, infile)

    def write_multi(self, outfile, xi, yi, nti=0, writeB=False,
                    write_dscale=False, zcut=0, depth_optimise=False):
        ''' Writes MULTI atmosphere file from a column of the 3D model,
        in RH 1.5D HDF5 format. Also writes the binary XDR file with magnetic
        fields, if writeB is true.
        '''
        from helita.sim.multi import watmos_multi
        from helita.sim.rh import write_B
        writeB = writeB and self.params['has_B']
        # if only total H available, will have to use rhpy (which is sometimes
        # risky...)
        if self.params['nhydr'] == 1:
            try:
                import rhpy
            except ImportError:
                raise ValueError("This function depents on rhpy, which is not"
                                 " installed in this system.")
            nh = rhpy.nh_lte(self.temperature[nti, xi, yi, zcut:].astype('Float64'),
                             self.electron_density[
                                   nti, xi, yi, zcut:].astype('Float64'),
                             self.hydrogen_populations[
                                   nti, 0, xi, yi, zcut:].astype('Float64'))
        elif self.params['nhydr'] == 6:
            nh = self.hydrogen_populations[nti, :, xi, yi, zcut:]
        else:
            raise ValueError("(EEE) write_multi: found %i hydrogen levels."
                             " For multi, need 6 or 1 " % self.params['nhydr'])
        M_TO_CM3 = (100.)**3
        M_TO_KM = 0.001
        temp = self.temperature[nti, xi, yi, zcut:].copy()
        ne = self.electron_density[nti, xi, yi, zcut:].copy() / M_TO_CM3
        if len(self.z.shape) > 2:
            self.z = self.z[:, xi, yi]
        z = self.z[nti, zcut:].copy() * M_TO_KM * 1.e5    # in cm
        vz = self.velocity_z[nti, xi, yi, zcut:].copy() * M_TO_KM
        nh = nh / M_TO_CM3
        if writeB:
            bx = self.B_x[nti, xi, yi, zcut:].copy()
            by = self.B_y[nti, xi, yi, zcut:].copy()
            bz = self.B_z[nti, xi, yi, zcut:].copy()
        else:
            bx = by = bz = None
        if depth_optimise:
            rho = self.hydrogen_populations[
                nti, 0, xi, yi, zcut:] * 2.380491e-24 / M_TO_CM3
            res = depth_optim(z, temp, ne, vz, rho, nh=nh, bx=bx, by=by, bz=bz)
            z, temp, ne, vz, rho, nh = res[:6]
            if writeB:
                bx, by, bz = res[6:]
        watmos_multi(outfile, temp, ne, z * 1e-5, vz=vz, nh=nh,
                     write_dscale=write_dscale,
                     id='%s txy-slice: (t,x,y) = (%i,%i,%i)' %
                     (self.params['description'], nti, xi, yi))
        if writeB:
            write_B('%s.B' % outfile, bx, by, bz)
            print(('--- Wrote magnetic field to %s.B' % outfile))

    def write_multi_3d(self, outfile, nti=0, sx=None, sy=None, sz=None,
                       big_endian=False):
        ''' Writes atmosphere in multi_3d format (the same as the
            pre-Jorrit multi3d) '''
        from helita.sim import multi
        ul = 1e2  # m to cm
        uv = 1e-3  # m/s to km/s
        # slicing and unit conversion
        if sx is None:
            sx = [0, self.nx, 1]
        if sy is None:
            sy = [0, self.ny, 1]
        if sz is None:
            sz = [0, self.nz, 1]
        if self.params['nhydr'] > 1:
            nh = np.mean(self.hydrogen_populations[nti, :, sx[0]:sx[1]:sx[2],
                                                   sy[0]:sy[1]:sy[2],
                                                   sz[0]:sz[1]:sz[2]], axis=1) / (ul**3)
        else:
            nh = self.hydrogen_populations[nti, 0, sx[0]:sx[1]:sx[2],
                                           sy[0]:sy[1]:sy[2],
                                           sz[0]:sz[1]:sz[2]] / (ul**3)
        rho = nh * 2.380491e-24  # nH to rho [g cm-3]
        x = self.x[sx[0]:sx[1]:sx[2]] * ul
        y = self.y[sy[0]:sy[1]:sy[2]] * ul
        z = self.z[nti, sz[0]:sz[1]:sz[2]] * ul
        ne = self.electron_density[nti, sx[0]:sx[1]:sx[2], sy[0]:sy[1]:sy[2],
                                   sz[0]:sz[1]:sz[2]] / (ul**3)
        temp = self.temperature[nti, sx[0]:sx[1]:sx[2], sy[0]:sy[1]:sy[2],
                                sz[0]:sz[1]:sz[2]]
        vz = self.velocity_z[nti, sx[0]:sx[1]:sx[2], sy[0]:sy[1]:sy[2],
                             sz[0]:sz[1]:sz[2]] * uv
        # write to file
        multi.write_atmos3d(outfile, x, y, z, ne, temp, vz, rho=rho,
                            big_endian=big_endian)


class NcdfAtmos:
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

    def write_multi(self, outfile, xi, yi, nti=0, writeB=False, write_dscale=False,
                    zcut=0, depth_optimise=False):
        ''' Writes MULTI atmosphere file from a column of the 3D model,
        in RH 1.5D ncdf format. Also writes the binary XDR file with magnetic
        fields, if writeB is true.
        '''
        from helita.sim.multi import watmos_multi
        from helita.sim.rh import write_B
        writeB = writeB and self.params['has_B']
        # if only total H available, will have to use rhpy (which is sometimes
        # risky...)
        if self.params['nhydr'] == 1:
            try:
                import rhpy
            except ImportError:
                raise ValueError("This function depents on rhpy, which is not"
                                 " installed in this system.")
            nh = rhpy.nh_lte(self.temperature[nti, xi, yi, zcut:].astype('Float64'),
                             self.electron_density[
                                   nti, xi, yi, zcut:].astype('Float64'),
                             self.hydrogen_populations[
                                   nti, 0, xi, yi, zcut:].astype('Float64'))
        elif self.params['nhydr'] == 6:
            nh = self.hydrogen_populations[nti, :, xi, yi, zcut:]
        else:
            raise ValueError("(EEE) write_multi: found %i hydrogen levels."
                             " For multi, need 6 or 1 " % self.params['nhydr'])
        M_TO_CM3 = (100.)**3
        M_TO_KM = 0.001
        temp = self.temperature[nti, xi, yi, zcut:].copy()
        ne = self.electron_density[nti, xi, yi, zcut:].copy() / M_TO_CM3
        if len(self.z.shape) > 2:
            self.z = self.z[:, xi, yi]
        z = self.z[nti, zcut:].copy() * M_TO_KM * 1.e5    # in cm
        vz = self.velocity_z[nti, xi, yi, zcut:].copy() * M_TO_KM
        nh = nh / M_TO_CM3
        if writeB:
            bx = self.B_x[nti, xi, yi, zcut:].copy()
            by = self.B_y[nti, xi, yi, zcut:].copy()
            bz = self.B_z[nti, xi, yi, zcut:].copy()
        else:
            bx = by = bz = None
        if depth_optimise:
            rho = self.hydrogen_populations[
                nti, 0, xi, yi, zcut:] * 2.380491e-24 / M_TO_CM3
            res = depth_optim(z, temp, ne, vz, rho, nh=nh, bx=bx, by=by, bz=bz)
            z, temp, ne, vz, rho, nh = res[:6]
            if writeB:
                bx, by, bz = res[6:]
        watmos_multi(outfile, temp, ne, z * 1e-5, vz=vz, nh=nh,
                     write_dscale=write_dscale,
                     id='%s txy-slice: (t,x,y) = (%i,%i,%i)' %
                     (self.params['description'], nti, xi, yi))
        if writeB:
            write_B('%s.B' % outfile, bx, by, bz)
            print(('--- Wrote magnetic field to %s.B' % outfile))

    def write_multi_3d(self, outfile, nti=0, sx=None, sy=None, sz=None,
                       big_endian=False):
        ''' Writes atmosphere in multi_3d format (the same as the
            pre-Jorrit multi3d) '''
        from helita.sim import multi
        ul = 1e2  # m to cm
        uv = 1e-3  # m/s to km/s
        # slicing and unit conversion
        if sx is None:
            sx = [0, self.nx, 1]
        if sy is None:
            sy = [0, self.ny, 1]
        if sz is None:
            sz = [0, self.nz, 1]
        if self.params['nhydr'] > 1:
            nh = np.mean(self.hydrogen_populations[nti, :, sx[0]:sx[1]:sx[2],
                                                   sy[0]:sy[1]:sy[2],
                                                   sz[0]:sz[1]:sz[2]], axis=1) / (ul**3)
        else:
            nh = self.hydrogen_populations[nti, 0, sx[0]:sx[1]:sx[2],
                                           sy[0]:sy[1]:sy[2],
                                           sz[0]:sz[1]:sz[2]] / (ul**3)
        rho = nh * 2.380491e-24  # nH to rho [g cm-3]
        x = self.x[sx[0]:sx[1]:sx[2]] * ul
        y = self.y[sy[0]:sy[1]:sy[2]] * ul
        z = self.z[nti, sz[0]:sz[1]:sz[2]] * ul
        ne = self.electron_density[nti, sx[0]:sx[1]:sx[2], sy[0]:sy[1]:sy[2],
                                   sz[0]:sz[1]:sz[2]] / (ul**3)
        temp = self.temperature[nti, sx[0]:sx[1]:sx[2], sy[0]:sy[1]:sy[2],
                                sz[0]:sz[1]:sz[2]]
        vz = self.velocity_z[nti, sx[0]:sx[1]:sx[2], sy[0]:sy[1]:sy[2],
                             sz[0]:sz[1]:sz[2]] * uv
        # write to file
        multi.write_atmos3d(outfile, x, y, z, ne, temp, vz, rho=rho,
                            big_endian=big_endian)


#############################################################################
###   TOOLS                                                               ###
#############################################################################
class DataHolder:
    def __init__(self):
        pass


def read_ncdf(inclass, infile):
    ''' DEPRECATED. Use read_hdf5 instead.
        Reads NetCDF file into inclass, instance of any class.
        Variables are read into class attributes, dimensions and attributes
        are read into params dictionary. '''
    from warnings import warn
    import netCDF4 as nc
    warn("Please use read_hdf5 instead", DeprecationWarning)
    # internal attributes of NetCDF groups
    ncdf_internals = dir(nc.Dataset)
    if not os.path.isfile(infile):
        raise IOError('read_ncdf: File %s not found' % infile)
    f = nc.Dataset(infile, mode='r')
    if 'params' not in dir(inclass):
        inclass.params = {}
    # add dimensions as attributes
    for d in list(f.dimensions.keys()):
        inclass.params[d] = len(f.dimensions[d])
    # add attributes
    attrs = [a for a in dir(f) if a not in ncdf_internals]
    for att in attrs:
        inclass.params[att] = getattr(f, att)
    # add variables
    for v in list(f.variables.keys()):
        vname = v.replace(' ', '_')    # sanitise string for spaces
        setattr(inclass, vname, f.variables[v])
    # Now do the same for all groups
    for group in list(f.groups.keys()):
        gname = group.replace(' ', '_')  # sanitise string for spaces
        setattr(inclass, gname, DataHolder())
        cur_group = f.groups[group]
        cur_class = getattr(inclass, gname)
        # add variables
        for v in list(cur_group.variables.keys()):
            vname = v.replace(' ', '_')  # sanitise string for spaces
            setattr(cur_class, vname, cur_group.variables[v])
        # add dimensions as attributes
        for d in list(cur_group.dimensions.keys()):
            inclass.params[d] = len(cur_group.dimensions[d])
        # add attributes
        attrs = [a for a in dir(cur_group) if a not in ncdf_internals]
        for att in attrs:
            inclass.params[att] = getattr(cur_group, att)
    return f


def read_hdf5(inclass, infile):
    ''' Reads HDF5/netCDF4 file into inclass, instance of any class.
        Variables are read into class attributes, dimensions and attributes
        are read into params dictionary. '''
    import h5py
    if not os.path.isfile(infile):
        raise IOError('read_hdf5: File %s not found' % infile)
    f = h5py.File(infile, mode='r')
    if 'params' not in dir(inclass):
        inclass.params = {}
    # add attributes
    attrs = [a for a in f.attrs]
    for att in f.attrs:
        inclass.params[att] = f.attrs[att]
    # add variables and groups
    for element in f:
        name = element.replace(' ', '_')    # sanitise string for spaces
        if type(f[element]) == h5py._hl.dataset.Dataset:
            setattr(inclass, name, f[element])
        if type(f[element]) == h5py._hl.group.Group:
            setattr(inclass, name, DataHolder())
            cur_class = getattr(inclass, name)
            cur_class.params = {}
            for variable in f[element]:   # add group variables
                vname = variable.replace(' ', '_')
                setattr(cur_class, vname, f[element][variable])
            for att in f[element].attrs:  # add group attributes
                cur_class.params[att] = f[element].attrs[att]
    return f


def make_ncdf_atmos(outfile, T, vz, nH, z, x=None, y=None, Bz=None, By=None,
                    Bx=None, rho=None, ne=None, vx=None, vy=None, desc=None,
                    snap=None, boundary=[1, 0], comp=False, complev=2,
                    append=False):
    ''' Creates NetCDF input file for rh15d.

        IN:
           outfile:  string name of destination. If file exists it will be wiped.
           T:        Temperature. Its shape will determine the output dimensions
           vz:       Same shape as T. In m/s.
           ne:       Same shape as T. In m-3. Optional.
           nH:       Shape [6, shape.T]. In m-3.
           z:        Same shape as last index of T. In m.
           x:        Same shape as first index of T. In m.
           y:        Same shape as second index of T. In m.
           snap:     Snapshot number(s)
           Bx, By, Bz: same shape as T. In T.
           rho, vx, vy: same shape as T. Optional.
           desc:     Description string
           boundary: Tuple with [bottom, top] boundary conditions. Key is:
                     0: Zero, 1: Thermalised, 2: Reflective.
           append:   if True, will append to existing file (if any).
           comp:     if false, compress file.
           complev:  compression level.
    '''
    import os
    import netCDF4 as nc

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
    idx = [None] * (4 - len(T.shape)) + [Ellipsis]  # empty axes for 1D/2D/3D
    T = T[idx]
    if ne is not None:
        ne = ne[idx]
    nH = nH[idx]
    vz = vz[idx]
    z = z[idx]
    if Bz is not None:
        Bx = Bx[idx]
        By = By[idx]
        Bz = Bz[idx]
    if rho is not None:
        rho = rho[idx]
    if vx is not None:
        vx = vx[idx]
    if vy is not None:
        vy = vy[idx]
    if len(T.shape) != 4:
        raise ValueError('Invalid shape for T')
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
        if ne is not None:
            ne_var = rootgrp.createVariable('electron_density', 'f8',
                                            ('nt', 'nx', 'ny', 'nz'),
                                            zlib=comp, complevel=complev)
            ne_var.units = 'm^-3'
        nh_var = rootgrp.createVariable('hydrogen_populations', 'f4',
                                        ('nt', 'nhydr', 'nx', 'ny', 'nz'),
                                        zlib=comp, complevel=complev)
        x_var = rootgrp.createVariable('x', 'f4', ('nx',))
        y_var = rootgrp.createVariable('y', 'f4', ('ny',))
        z_var = rootgrp.createVariable('z', 'f4', ('nt', 'nz'))
        nt_var = rootgrp.createVariable('snapshot_number', 'i4', ('nt',))
        T_var.units = 'K'
        vz_var.units = 'm s^-1'
        nh_var.units = 'm^-3'
        z_var.units = 'm'
        x_var.units = 'm'
        y_var.units = 'm'
        if Bz is not None:
            bx_var = rootgrp.createVariable('B_x', 'f4',
                                            ('nt', 'nx', 'ny', 'nz'), zlib=comp,
                                            least_significant_digit=5,
                                            complevel=complev)
            by_var = rootgrp.createVariable('B_y', 'f4',
                                            ('nt', 'nx', 'ny', 'nz'), zlib=comp,
                                            least_significant_digit=5,
                                            complevel=complev)
            bz_var = rootgrp.createVariable('B_z', 'f4',
                                            ('nt', 'nx', 'ny', 'nz'), zlib=comp,
                                            least_significant_digit=5,
                                            complevel=complev)
            bx_var.units = 'T'
            by_var.units = 'T'
            bz_var.units = 'T'
        if rho is not None:
            rho_var = rootgrp.createVariable('density', 'f4',
                                             ('nt', 'nx', 'ny', 'nz'), zlib=comp,
                                             least_significant_digit=5,
                                             complevel=complev)
            rho_var.units = 'kg m^-3'
        if vx is not None:
            vx_var = rootgrp.createVariable('velocity_x', 'f4',
                                            ('nt', 'nx', 'ny', 'nz'), zlib=comp,
                                            least_significant_digit=5,
                                            complevel=complev)
            vx_var.units = 'm s^-1'
        if vy is not None:
            vy_var = rootgrp.createVariable('velocity_y', 'f4',
                                            ('nt', 'nx', 'ny', 'nz'), zlib=comp,
                                            least_significant_digit=5,
                                            complevel=complev)
            vy_var.units = 'm s^-1'
        if desc is None:
            rootgrp.description = \
                "BIFROST snapshot"
        else:
            rootgrp.description = desc
        if boundary is None:
            rootgrp.boundary_top = 0
            rootgrp.boundary_bottom = 1
        else:
            rootgrp.boundary_top = boundary[1]
            rootgrp.boundary_bottom = boundary[0]
        if Bz is None:
            rootgrp.has_B = 0
        else:
            rootgrp.has_B = 1
        nt = [0, nt]
    else:
        # get variables
        T_var = rootgrp.variables['temperature']
        vz_var = rootgrp.variables['velocity_z']
        nh_var = rootgrp.variables['hydrogen_populations']
        nt_var = rootgrp.variables['snapshot_number']
        x_var = rootgrp.variables['x']
        y_var = rootgrp.variables['y']
        z_var = rootgrp.variables['z']
        if ne is not None:
            ne_var = rootgrp.variables['electron_density']
        if Bz is not None:
            bx_var = rootgrp.variables['B_x']
            by_var = rootgrp.variables['B_y']
            bz_var = rootgrp.variables['B_z']
        if rho is not None:
            rho_var = rootgrp.variables['density']
        if vx is not None:
            vx_var = rootgrp.variables['velocity_x']
        if vy is not None:
            vy_var = rootgrp.variables['velocity_y']
        nti = len(rootgrp.dimensions['nt'])
        nt = [nti, nti + nt]
    T_var[nt[0]:nt[1]] = T
    vz_var[nt[0]:nt[1]] = vz
    nh_var[nt[0]:nt[1], :nhydr] = nH
    if ne is not None:
        ne_var[nt[0]:nt[1]] = ne
    if Bz is not None:
        bx_var[nt[0]:nt[1]] = Bx
        by_var[nt[0]:nt[1]] = By
        bz_var[nt[0]:nt[1]] = Bz
    if rho is not None:
        rho_var[nt[0]:nt[1]] = rho
    if vx is not None:
        vx_var[nt[0]:nt[1]] = vx
    if vy is not None:
        vy_var[nt[0]:nt[1]] = vy
    x_var[:] = x
    y_var[:] = y
    z_var[nt[0]:nt[1]] = z
    nt_var[nt[0]:nt[1]] = snap
    rootgrp.close()
    return


def make_hdf5_atmos(outfile, T, vz, nH, z, x=None, y=None, Bz=None, By=None,
                    Bx=None, rho=None, ne=None, vx=None, vy=None, desc=None,
                    snap=None, boundary=[1, 0], comp=None, complev=None,
                    append=False):
    """
    Creates HDF5 input file for RH 1.5D.

    Parameters
    ----------
    outfile : string
        Name of destination. If file exists it will be wiped.
    T : n-D array
        Temperature in K. Its shape will determine the output
        dimensions (can be 1D, 2D, or 3D).
    vz : n-D array
        Line of sight velocity in m/s. Same shape as T.
    nH : n-D array
        Hydrogen populations in m^-3. Shape is [nhydr, shape.T] where
        nydr can be 1 (total number of protons) or more (level populations).
    z : n-D array
        Height in m. Can have same shape as T (different height scale
        for each column) or be only 1D (same height for all columns).
    ne : n-D array, optional
        Electron density in m^-3. Same shape as T.
    rho : n-D array, optional
        Density in kg / m^-3. Same shape as T.
    vx : n-D array, optional
        x velocity in m/s. Same shape as T. Not in use by RH 1.5D.
    vy : n-D array, optional
        y velocity in m/s. Same shape as T. Not in use by RH 1.5D.
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
    comp : string, optional
        Options are: None (default), 'gzip', 'szip', 'lzf'.
    complev : integer or tuple, optional
        Compression level. Integer for 'gzip', 2-tuple for szip.
    """
    import os
    import datetime
    import h5py
    mode = ['w', 'a']
    if (append and not os.path.isfile(outfile)):
        append = False
    rootgrp = h5py.File(outfile, mode=mode[append])
    if nH.shape == T.shape:
        nhydr = 1
    else:
        nhydr = nH.shape[0]
    idx = [None] * (4 - len(T.shape)) + [Ellipsis]  # empty axes for 1D/2D/3D
    T = T[idx]
    if ne is not None:
        ne = ne[idx]
    nH = nH[idx]
    vz = vz[idx]
    z = z[idx]
    if Bz is not None:
        Bx = Bx[idx]
        By = By[idx]
        Bz = Bz[idx]
    if rho is not None:
        rho = rho[idx]
    if vx is not None:
        vx = vx[idx]
    if vy is not None:
        vy = vy[idx]
    if len(T.shape) != 4:
        raise ValueError('Invalid shape for T')
    nt = T.shape[0]
    if snap is None:
        snap = np.arange(nt, dtype='i4')
    if not append:
        # for a new file, create datasets
        max_dims = (None,) + T.shape[1:] # time is unlimited dimension
        rootgrp.attrs["nx"] = T.shape[-3]
        rootgrp.attrs["ny"] = T.shape[-2]
        rootgrp.attrs["nz"] = T.shape[-1]
        rootgrp.attrs["nhydr"] = nhydr
        T_var = rootgrp.create_dataset("temperature", dtype="f4",
                                       shape=T.shape, maxshape=max_dims,
                                       fletcher32=True, compression=comp,
                                       compression_opts=complev)
        vz_var = rootgrp.create_dataset("velocity_z", dtype="f4",
                                        shape=T.shape, maxshape=max_dims,
                                        fletcher32=True, compression=comp,
                                        compression_opts=complev)
        nh_var = rootgrp.create_dataset("hydrogen_populations", dtype="f4",
                                shape=(nt, nhydr,) + T.shape[1:],
                                maxshape=(None, nhydr) + T.shape[1:],
                                fletcher32=True, compression=comp,
                                compression_opts=complev)
        if ne is not None:
            ne_var = rootgrp.create_dataset("electron_density", dtype="f8",
                                            shape=T.shape, maxshape=max_dims,
                                            fletcher32=True, compression=comp,
                                            compression_opts=complev)
            ne_var.attrs["units"] = 'm^-3'
        x_var = rootgrp.create_dataset("x", dtype="f4", shape=(T.shape[1],))
        y_var = rootgrp.create_dataset("y", dtype="f4", shape=(T.shape[1],))
        z_var = rootgrp.create_dataset("z", dtype="f4", shape=(nt, T.shape[3]),
                                       maxshape=(None, T.shape[3]))
        nt_var = rootgrp.create_dataset("snapshot_number", dtype="i4",
                                        shape=(nt,))
        T_var.attrs["units"] = 'K'
        vz_var.attrs["units"] = 'm s^-1'
        nh_var.attrs["units"] = 'm^-3'
        z_var.attrs["units"] = 'm'
        x_var.attrs["units"] = 'm'
        y_var.attrs["units"] = 'm'
        if Bz is not None:
            bx_var = rootgrp.create_dataset("B_x", dtype="f4",
                                            shape=T.shape, maxshape=max_dims,
                                            fletcher32=True, compression=comp,
                                            compression_opts=complev)
            by_var = rootgrp.create_dataset("B_y", dtype="f4",
                                            shape=T.shape, maxshape=max_dims,
                                            fletcher32=True, compression=comp,
                                            compression_opts=complev)
            bz_var = rootgrp.create_dataset("B_z", dtype="f4",
                                            shape=T.shape, maxshape=max_dims,
                                            fletcher32=True, compression=comp,
                                            compression_opts=complev)
            bx_var.attrs["units"] = 'T'
            by_var.attrs["units"] = 'T'
            bz_var.attrs["units"] = 'T'
        if rho is not None:
            rho_var = rootgrp.create_dataset("density", dtype="f4",
                                             shape=T.shape, maxshape=max_dims,
                                             fletcher32=True, compression=comp,
                                             compression_opts=complev)
            rho_var.attrs["units"] = 'kg m^-3'
        if vx is not None:
            vx_var = rootgrp.create_dataset("velocity_x", dtype="f4",
                                            shape=T.shape, maxshape=max_dims,
                                            fletcher32=True, compression=comp,
                                            compression_opts=complev)
            vx_var.attrs["units"] = 'm s^-1'
        if vy is not None:
            vy_var = rootgrp.create_dataset("velocity_y", dtype="f4",
                                            shape=T.shape, maxshape=max_dims,
                                            fletcher32=True, compression=comp,
                                            compression_opts=complev)
            vy_var.attrs["units"] = 'm s^-1'
        if desc is None:
            rootgrp.attrs["description"] = ("Created with make_hdf5_atmos."
                                            "on %s" % datetime.datetime.now())
        else:
            rootgrp.attrs["description"] = desc
        if boundary is None:
            rootgrp.attrs["boundary_top"] = 0
            rootgrp.attrs["boundary_bottom"] = 1
        else:
            rootgrp.attrs["boundary_top"] = boundary[1]
            rootgrp.attrs["boundary_bottom"] = boundary[0]
        if Bz is None:
            rootgrp.attrs["has_B"] = 0
        else:
            rootgrp.attrs["has_B"] = 1
        nt = [0, nt]
    else:
        # get variables
        T_var = rootgrp['temperature']
        vz_var = rootgrp['velocity_z']
        nh_var = rootgrp['hydrogen_populations']
        nt_var = rootgrp['snapshot_number']
        x_var = rootgrp['x']
        y_var = rootgrp['y']
        z_var = rootgrp['z']
        if ne is not None:
            ne_var = rootgrp['electron_density']
        if Bz is not None:
            bx_var = rootgrp['B_x']
            by_var = rootgrp['B_y']
            bz_var = rootgrp['B_z']
        if rho is not None:
            rho_var = rootgrp['density']
        if vx is not None:
            vx_var = rootgrp['velocity_x']
        if vy is not None:
            vy_var = rootgrp['velocity_y']
        nti = int(rootgrp.attrs['nt'])
        nt = [nti, nti + nt]
    T_var[nt[0]:nt[1]] = T
    vz_var[nt[0]:nt[1]] = vz
    nh_var[nt[0]:nt[1], :nhydr] = nH
    if ne is not None:
        ne_var[nt[0]:nt[1]] = ne
    if Bz is not None:
        bx_var[nt[0]:nt[1]] = Bx
        by_var[nt[0]:nt[1]] = By
        bz_var[nt[0]:nt[1]] = Bz
    if rho is not None:
        rho_var[nt[0]:nt[1]] = rho
    if vx is not None:
        vx_var[nt[0]:nt[1]] = vx
    if vy is not None:
        vy_var[nt[0]:nt[1]] = vy
    if x is not None:
        x_var[:] = x
    if y is not None:
        y_var[:] = y
    z_var[nt[0]:nt[1]] = z
    nt_var[nt[0]:nt[1]] = snap
    rootgrp.attrs['nt'] = z_var.shape[0]
    rootgrp.close()
    return


def depth_optim(height, temp, ne, vz, rho, nh=None, bx=None, by=None, bz=None,
                tmax=5e4):
    """
    Performs depth optimisation of one single column (as per multi_3d).

        IN:
            height   [cm]
            temp     [K]
            ne       [cm-3]
            vz       [any]
            rho      [g cm-3]
            nh       [any] (optional)
            bx,by,bz [any] (optional)
            tmax     [K] maximum temperature of the first point

    """
    from scipy.integrate import cumtrapz
    import scipy.interpolate as interp
    ndep = len(height)
    # calculate optical depth from H-bf only
    taumax = 100
    grph = 2.26e-24
    crhmbf = 2.9256e-17
    ee = 1.602189E-12
    bk = 1.380662E-16
    xhbf = 1.03526e-16 * ne * crhmbf / temp**1.5 * \
        np.exp(0.754 * ee / bk / temp) * rho / grph
    tau = np.concatenate(([0.], cumtrapz(xhbf, -height)))
    idx = (tau < taumax) & (temp < tmax)
    # find maximum variance of T, rho, and tau for each depth
    tt = temp[idx]
    rr = rho[idx]
    ta = tau[idx]
    tdiv = np.abs(np.log10(tt[1:]) - np.log10(tt[:-1])) / np.log10(1.1)
    rdiv = np.abs(np.log10(rr[1:]) - np.log10(rr[:-1])) / np.log10(1.1)
    taudiv = np.abs(np.log10(ta[1:]) - np.log10(ta[:-1])) / 0.1
    taudiv[0] = 0.
    aind = np.concatenate(
        ([0.], np.cumsum(np.max(np.array([tdiv, rdiv, taudiv]), axis=0))))
    aind *= (ndep - 1) / aind[-1]
    # interpolate new height so it is constant in aind2
    nheight = interp.splev(np.arange(ndep), interp.splrep(
        aind, height[idx], k=3, s=0), der=0)
    # interpolate quantities for new depth scale
    ntemp = np.exp(interp.splev(nheight, interp.splrep(height[::-1], np.log(temp[::-1]),
                                                       k=3, s=0), der=0))
    nne = np.exp(interp.splev(nheight, interp.splrep(height[::-1], np.log(ne[::-1]),
                                                     k=3, s=0), der=0))
    nrho = np.exp(interp.splev(nheight, interp.splrep(height[::-1], np.log(rho[::-1]),
                                                      k=3, s=0), der=0))
    nvz = interp.splev(nheight, interp.splrep(height[::-1], vz[::-1],
                                              k=3, s=0), der=0)
    result = [nheight, ntemp, nne, nvz, nrho]
    if nh is not None:
        for k in range(nh.shape[0]):
            nh[k] = np.exp(interp.splev(nheight, interp.splrep(height[::-1],
                                                               np.log(nh[k, ::-1]), k=3, s=0), der=0))
        result += [nh]
    if bx is not None:
        nbx = interp.splev(nheight, interp.splrep(
            height[::-1], bx[::-1], k=3, s=0), der=0)
        nby = interp.splev(nheight, interp.splrep(
            height[::-1], by[::-1], k=3, s=0), der=0)
        nbz = interp.splev(nheight, interp.splrep(
            height[::-1], bz[::-1], k=3, s=0), der=0)
        result += [nbx, nby, nbz]
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
    new_wave: 1D array
        Alternatively to start/end, one can specify an array of
        wavelengths here.
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
    from helita.utils.waveconv import waveconv
    if new_wave is None:
        new_wave = np.arange(start, end, step)
        if None in [start, end, step]:
            raise ValueError('Must specify either new_wave, or start, end, '
                             'step. Stopping.')
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
    print(("Wrote %i wavelengths to file." % nw))
    return


def read_wave_file(infile):
    """
    Reads RH wavelength file.

    Parameters
    ----------
    infile - string
        Name of wavelength file to read.
    """
    import xdrlib
    import io
    from .rh import read_xdr_var
    f = io.open(infile, 'rb')
    buf = xdrlib.Unpacker(f.read())
    f.close()
    nw = read_xdr_var(buf, 'i')
    return read_xdr_var(buf, ('d', (nw,)))


def clean_var(data, only_positive=True):
    """
    Cleans a 2D or 3D variable filled with NaNs and other irregularities.
    """
    from ..utils import utilsfast
    data = np.ma.masked_invalid(data, copy=False)
    if only_positive:
        data = np.ma.masked_less(data, 0., copy=False)
    tmp = np.abs(data)
    thres = tmp.mean() + tmp.std() * 4  # points more than 4 std away
    data = np.ma.masked_where(tmp > thres, data, copy=False)
    if data.ndim == 2:
        data = data[..., np.newaxis]
    for k in range(data.shape[-1]):
        tmp = data[..., k].astype("d")
        tmp[data[..., k].mask] = np.nan
        data[..., k] = utilsfast.replace_nans(tmp, 15, 0.1, 3, "localmean")
    return np.squeeze(data)
