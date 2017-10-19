#!/usr/bin/env python

"""
Create netcdf files to force CABLE-CNP for FACE MIP type experiment.

Typically this process involves creating:
1. a file for the standard 30 year spin to equillibrium at 1850 (fixed CO2,
   Ndep, Pdep)
2. A transient met file, i.e. 1850-20XX, with time varying CO2, Ndep & Pdep.
3. An expermental run file, with whatever treatment added (i.e. eCO2,
   PPT mainpulation, etc).

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
        self.KG_2_G = 1000.0
        self.HA_2_M2 = 10000.0
        self.YR_2_DAY = 365.0

    def create_spin_file(self, ofname, co2_fixed, ndep_fixed, pdep_fixed):

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

        # Add CO2, NDEP & PDEP
        out_length = len(ds.variables["Rainf"][:,:])
        for v in ["CO2air", "Ndep", "Pdep"]:
            out.createVariable(v, 'float32', ('time', 'y', 'x'))

        out.variables["CO2air"][:,:,:] = np.ones((out_length,1,1)) * co2_fixed
        out.variables["CO2air"].setncatts({'units': "ppm"})
        out.variables["CO2air"].setncatts({'missing_value': "-9999"})
        out.variables["CO2air"].setncatts({'long_name':
                                           "Atmosphereic CO2 concentration"})

        # kg ha-1 y-1 -> gN m-2 d
        conv = self.KG_2_G / self.HA_2_M2 / self.YR_2_DAY
        out.variables["Ndep"][:,:,:] = np.ones((out_length,1,1)) * \
                                               (ndep_fixed * conv)
        out.variables["Ndep"].setncatts({'units': "gN m-2 d-1"})
        out.variables["Ndep"].setncatts({'missing_value': "-9999"})
        out.variables["Ndep"].setncatts({'long_name': "N deposition"})

        out.variables["Pdep"][:,:,:] = np.ones((out_length,1,1)) *  \
                                               (pdep_fixed * conv)
        out.variables["Pdep"].setncatts({'units': "gP m-2 d-1"})
        out.variables["Pdep"].setncatts({'missing_value': "-9999"})
        out.variables["Pdep"].setncatts({'long_name': "P deposition"})

        # write global attributes
        for ncattr in ds.ncattrs():
            out.setncattr(ncattr, getattr(ds, ncattr))

        ds.close()

    def create_transient_file(self, ofname):
        ds = nc.Dataset(self.met_fname)
        out = nc.Dataset(ofname, 'w', format='NETCDF4')

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
    local_met_dir = "met_files"
    ofname = os.path.join(local_met_dir, "%s_met_spin.nc" % (site))
    G = GenerateMetFiles(site, met_fname, co2_ndep_fname)

    co2_fixed = 284.7  # umol mol-1
    ndep_fixed = 0.79  # kg N ha-1 yr-1
    pdep_fixed = 0.144 # kg N ha-1 yr-1
    G.create_spin_file(ofname, co2_fixed, ndep_fixed, pdep_fixed)
    #G.create_transient_file(ofname)
    #G.check_differences(ofname)
