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

        for v in nc_vars:
            out.createVariable(v, ds.variables[v].dtype,
                               ds.variables[v].dimensions)
            if len(ds.variables[v].dimensions) == 2:
                out.variables[v][:,:] = ds.variables[v][:,:]
            elif len(ds.variables[v].dimensions) == 3:
                out.variables[v][:,:,:] = ds.variables[v][:,:,:]
            elif len(ds.variables[v].dimensions) == 4:
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

        import matplotlib.pyplot as plt
        plt.plot(ds2.variables["Tair"][0:250,0,0,0]-273.15)
        plt.plot(ds1.variables["Tair"][0:250,0,0,0]-273.15, "r.")
        plt.show()

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

        #var = ds.variables[v][tindex0:tindex1]





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
