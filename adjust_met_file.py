#!/usr/bin/env python

"""
Create netcdf files to force CABLE-CNP for FACE MIP type experiment

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (17.10.2017)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import netCDF4 as nc
import datetime as dt
import numpy as np

class AdjustCableMetFile(object):

    def __init__(self, site,  met_fname, co2_ndep_fname):

        self.site = site
        self.met_fname = met_fname
        self.co2_ndep_fname = co2_ndep_fname

    def write_new_met_file(self, ofname):

        out = nc.Dataset(ofname, 'w', format='NETCDF4')
        ds = nc.Dataset(self.met_fname)
        time = nc.num2date(ds.variables['time'][:], ds.variables['time'].units)

        nc_attrs = ds.ncattrs()
        nc_dims = [dim for dim in ds.dimensions]
        nc_vars = [var for var in ds.variables]

        # Write dimensions (x, y, z, time)
        data = {}
        for dim in nc_dims:
            out.createDimension(dim, ds.variables[dim].size)
            data[dim] = out.createVariable(dim, ds.variables[dim].dtype,\
                                          (dim,))
            for ncattr in ds.variables[dim].ncattrs():
                data[dim].setncattr(ncattr, ds.variables[dim].getncattr(ncattr))

        out.variables['time'][:] = ds.variables['time'][:]
        out.variables['x'][:] = ds.variables['x'][:]
        out.variables['y'][:] = ds.variables['y'][:]
        out.variables['z'][:] = ds.variables['z'][:]
        [nc_vars.remove(i) for i in ["time", "x", "y", "z"]]

        # Write vars
        vars_without_time = []
        vars_with_time = []
        vars_with_time_and_zdim = []
        for v in nc_vars:
            if len(ds.variables[v].dimensions) == 2:
                vars_without_time.append(v)
            elif len(ds.variables[v].dimensions) == 3:
                vars_with_time.append(v)
            elif len(ds.variables[v].dimensions) == 4:
                vars_with_time_and_zdim.append(v)

        for v in vars_without_time:

            out.createVariable(v, ds.variables[v].dtype, ('y', 'x'))
            out.variables[v][:,:] = ds.variables[v][:,:]
            ncvar = ds.variables[v]
            out = self.write_attributes(v, ncvar, out)

        for v in vars_with_time:
            out.createVariable(v, ds.variables[v].dtype, ('time', 'y', 'x'))
            out.variables[v][:,:,:] = ds.variables[v][:,:,:]
            ncvar = ds.variables[v]
            out = self.write_attributes(v, ncvar, out)

        for v in vars_with_time_and_zdim:
            out.createVariable(v, ds.variables[v].dtype,
                               ('time', 'z', 'y', 'x'))
            out.variables[v][:,:,:,:] = ds.variables[v][:,:,:,:]
            ncvar = ds.variables[v]
            out = self.write_attributes(v, ncvar, out)

        # write global attributes
        for ncattr in ds.ncattrs():
            out.setncattr(ncattr, getattr(ds, ncattr))

        ds.close()

    def check_differences(self, ofname):
        ds1 = nc.Dataset(self.met_fname)
        ds2 = nc.Dataset(ofname)

        nc_vars = [var for var in ds1.variables]
        for v in nc_vars:
            if not np.array_equal(ds1.variables[v], ds2.variables[v]):
                print(v, "Differ")

        #for v in nc_vars:
        #    print(v,  ds1.variables[v].dtype, ds2.variables[v].dtype)

    def write_attributes(self, v, ncvar, out):

        if hasattr(ncvar, 'missing_value'):
            mval = ncvar.missing_value
            out.variables[v].setncatts({'missing_value': mval})
        if hasattr(ncvar, 'units'):
            uval = ncvar.units
            out.variables[v].setncatts({'units': uval})
        if hasattr(ncvar, 'long_name'):
            ln = ncvar.long_name
            out.variables[v].setncatts({'long_name': ln})
        if hasattr(ncvar, 'values'):
            vx = ncvar.values
            out.variables[v].setncatts({'values': vx})
        if hasattr(ncvar, 'CF_name'):
            cf = ncvar.CF_name
            out.variables[v].setncatts({'CF_name': cf})

        return out

    def setup_ini_spin_met_file(self, start_yr_spin, end_yr_spin, start_met_yr,
                                end_met_yr):

        out = nc.Dataset('test.nc', 'w', format='NETCDF4')
        ds = nc.Dataset(self.met_fname)

        total_yrs = end_yr_spin - start_yr_spin

        # Account for offset, first day stars 00:30 to 1.
        start = dt.datetime(start_met_yr,1,1,0,30)
        stop = dt.datetime(end_met_yr+1,1,1,0,0)
        tindex0 = nc.date2index(start, ds.variables['time'], select='nearest')
        tindex1 = nc.date2index(stop, ds.variables['time'], select='nearest')
        start_idx = nc.num2date(ds.variables['time'][tindex0],
                                ds.variables['time'].units)
        end_idx = nc.num2date(ds.variables['time'][tindex1],
                              ds.variables['time'].units)



        """
        out.createDimension(dim, nc_fid.variables[dim].size)


        print(time)
        out.variables['time'][:] = time
        #out.variables['lat'][:] = lats
        #out.variables['lon'][:] = lons


        #for ncattr in nc.variables['time'].ncattrs():
        #    print ( w_nc_dim.setncattr(ncattr, nc.variables['time'].getncattr(ncattr)))
        ds.close()
        sys.exit()

        vars = ds.variables.keys()
        for v in vars:
            var = ds.variables[v][tindex0:tindex1]
            print(var)
            print(ds.variables[v])
            print("\n")
            print(ds.variables[v].units)



            sys.exit()
        sys.exit()
        #for attr in ds.ncattrs():
        #    print(attr, '=', getattr(ds, attr))

        #temp = rootgrp.createVariable("temp","f4",("time","level","lat","lon",))


        #print(time[0], time[-1])
        total_yrs = end_yr_spin - start_yr_spin

        # Account for offset, first day stars 00:30 to 1.
        start = dt.datetime(start_met_yr,1,1,0,30)
        stop = dt.datetime(end_met_yr+1,1,1,0,0)
        tindex0 = nc.date2index(start, ds.variables['time'], select='nearest')
        tindex1 = nc.date2index(stop, ds.variables['time'], select='nearest')
        start_idx = nc.num2date(ds.variables['time'][tindex0],
                                ds.variables['time'].units)
        end_idx = nc.num2date(ds.variables['time'][tindex1],
                              ds.variables['time'].units)

        vars = ds.variables.keys()
        for v in vars:
            print(v)
            print(ds.variables[v])
            print("\n")
            print(ds.variables[v].units)

            p = ds.variables[v][tindex0:tindex1]
            print(p)
            sys.exit()
        sys.exit()

        print(start)
        print(stop)
        print(tindex0)
        print(tindex1)


        print (nc.num2date(ds.variables['time'][tindex0], ds.variables['time'].units), nc.num2date(ds.variables['time'][tindex1], ds.variables['time'].units))

        #print (nc.num2date(ds.variables['time'][tindex1], ds.variables['time'].units))



        vars = ds.variables.keys()

        print(vars)
        time = nc.num2date(ds.variables['time'][:], ds.variables['time'].units)
        print(time[0], time[-1])
        sys.exit()


        for v in ds.variables.keys():
            var = ds.variables[v][tindex0:tindex1]
            print(var)
            print(ds.variables[v])
            print("\n")
            print(ds.variables[v].units)

        for attr in dset.ncattrs():
            print(attr, '=', getattr(ds, attr))
        """




if __name__ == "__main__":

    site = "TumbaFluxnet"

    met_dir = "../../met_data/plumber_met/"
    co2_ndep_dir = "../../met_data/co2_ndep"

    met_fname = os.path.join(met_dir, '%s.1.4_met.nc' % (site))
    co2_ndep_fname = os.path.join(co2_ndep_dir,
                                  "AmaFACE_co2ndepforcing_1850_2015_AMB.csv")

    M = AdjustCableMetFile(site, met_fname, co2_ndep_fname)

    # Figured out in run_cable_site_CNP
    start_yr_spin = 1834
    end_yr_spin = 1853
    start_met_yr = 2002
    end_met_yr = 2005
    ofname = "test.nc"
    M.write_new_met_file(ofname)
    M.check_differences(ofname)
