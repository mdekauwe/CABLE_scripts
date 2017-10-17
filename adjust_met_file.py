#!/usr/bin/env python

"""
Site run with full CNP (and Pop)
================================

- Model spin-up: using K34 tower info, CO2=284.7; NDEP-0.79 kg N ha-1 yr-1;
                 PDEP=0.144 kg P ha-1 yr-1
- Transient: 1851-1998, varying CO2 and NDEP, but just recycling the met data
- Historical: actual met (dates correspond to simulated dates) and actual CO2
- CNP + POP switched on.

During the spinup, we are recycling in 30 year chunks.

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (19.09.2017)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import shutil
import tempfile
import netCDF4 as nc
import math
import datetime as dt

class AdjustCableMetFile(object):

    def __init__(self, site,  met_fname, co2_ndep_fname):

        self.site = site
        self.met_fname = met_fname
        self.co2_ndep_fname = co2_ndep_fname

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

        time = nc.num2date(ds.variables['time'][:], ds.variables['time'].units)

        nc_attrs = ds.ncattrs()
        nc_dims = [dim for dim in ds.dimensions]
        nc_vars = [var for var in ds.variables]

        for v in nc_vars:
            print(v, ds.variables[v].dtype, ds.variables[v].dims)

        sys.exit()
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

        out.createVariable('latitude', 'f8', ('y', 'x'))
        out.variables['latitude'][:,:] = ds.variables['latitude'][:,:]

        out.createVariable('longitude', 'f8', ('y', 'x'))
        out.variables['longitude'][:,:] = ds.variables['longitude'][:,:]

        out.createVariable('SWdown', 'f8', ('time', 'y', 'x'))
        out.variables['SWdown'][:,:,:] = ds.variables['SWdown'][:,:]
        #w_nc_var.setncatts({'long_name': u"Latitude",\
        #            'units': u"degrees_north", 'level_desc': u'Surface',\
        #            'var_desc': u"Air temperature departure",\
        #            'statistic': u'Mean\nM'})
        #w_nc_fid.variables['air_dep'][:] = departure

        #nc_vars.remove("time")
        #nc_vars.remove("x")
        #nc_vars.remove("y")
        #nc_vars.remove("z")



        #print(nc_dims)
        #print(nc_vars)
        #print(nc_attrs)


        #sys.exit()

        ds.close()

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
    M.setup_ini_spin_met_file(start_yr_spin, end_yr_spin, start_met_yr,
                              end_met_yr)
