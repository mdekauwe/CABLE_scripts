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
import math

class GenerateMetFiles(object):

    def __init__(self, site,  met_fname, co2_ndep_fname):

        self.site = site
        self.met_fname = met_fname
        self.co2_ndep_fname = co2_ndep_fname

    def create_met_spin_file(self):

        ds = nc.Dataset(self.met_fname)
        out = nc.Dataset(ofname, 'w', format='NETCDF4')

        nc_attrs = ds.ncattrs()
        nc_dims = [dim for dim in ds.dimensions]
        #nc_vars = [var for var in ds.variables]
        nc_vars = ds.variables.keys()

        # Write vars
        for v in nc_vars:
            if len(ds.variables[v].dimensions) == 1:
                out.createDimension(v, ds.variables[v].size)
                out.createVariable(v, ds.variables[v].dtype, (v,))
                ncvar = ds.variables[v]
                if hasattr(ncvar, 'units'):
                    mval = ncvar.units
                    out.variables[v].setncatts({'units': mval})
                out.variables[v][:] = ds.variables[v][:]
            else:
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

        if hasattr(ncvar, 'units'):
            val = ncvar.units
            out.variables[v].setncatts({'units': val})
        if hasattr(ncvar, 'missing_value'):
            val = ncvar.missing_value
            out.variables[v].setncatts({'missing_value': val})
        if hasattr(ncvar, 'long_name'):
            val = ncvar.long_name
            out.variables[v].setncatts({'long_name': val})
        if hasattr(ncvar, 'values'):
            val = ncvar.values
            out.variables[v].setncatts({'values': val})
        if hasattr(ncvar, 'CF_name'):
            val = ncvar.CF_name
            out.variables[v].setncatts({'CF_name': val})

        return out






if __name__ == "__main__":

    site = "TumbaFluxnet"

    met_dir = "../../met_data/plumber_met/"
    co2_ndep_dir = "../../met_data/co2_ndep"

    met_fname = os.path.join(met_dir, '%s.1.4_met.nc' % (site))
    co2_ndep_fname = os.path.join(co2_ndep_dir,
                                  "AmaFACE_co2ndepforcing_1850_2015_AMB.csv")

    G = GenerateMetFiles(site, met_fname, co2_ndep_fname)

    # Figured out in run_cable_site_CNP
    start_yr_spin = 1834
    end_yr_spin = 1853
    start_met_yr = 2002
    end_met_yr = 2005
    ofname = "test.nc"

    #G.check_differences(ofname)
    G.create_met_spin_file()
