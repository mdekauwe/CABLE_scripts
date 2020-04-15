#!/usr/bin/env python

"""
Generate a single pixel netcdf timeseries from spatial inputs for CABLE testing

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (21.10.2018)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import netCDF4 as nc
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import datetime
import glob as glob

def main(vars, start_yr, end_yr,  data_type, row, col):

    LWdown_out = np.zeros(0)
    PSurf_out = np.zeros(0)
    Qair_out = np.zeros(0)
    Rainf_out = np.zeros(0)
    Snowf_out = np.zeros(0)
    SWdown_out = np.zeros(0)
    Tair_out = np.zeros(0)
    Wind_out = np.zeros(0)

    for yr in range(start_yr, end_yr +1):
        print(yr)

        for i, var in enumerate(vars):

            fname = "%s_%d_%d_%d.nc" % (var, yr, row, col)
            ds = xr.open_dataset(fname)

            vals = ds[var][:,0,0].values
            lat = ds["lat"][0].values
            lon = ds["lon"][0].values

            if var == "LWdown":
                LWdown_out = np.append(LWdown_out, vals)
            elif var == "PSurf":
                PSurf_out = np.append(PSurf_out, vals)
            elif var == "Qair":
                Qair_out = np.append(Qair_out, vals)
            elif var == "Rainf":
                Rainf_out = np.append(Rainf_out, vals)
            elif var == "Snowf":
                Snowf_out = np.append(Snowf_out, vals)
            elif var == "SWdown":
                SWdown_out = np.append(SWdown_out, vals)
            elif var == "Tair":
                Tair_out = np.append(Tair_out, vals)
            elif var == "Wind":
                Wind_out = np.append(Wind_out, vals)

            ds.close()


    ofname = "AWAP_single_pixel_%.2f_%.2f.nc" % (lat, lon)

    times = []
    secs = 0.0
    n_timesteps = len(Tair_out)
    for i in range(n_timesteps):
        times.append(secs)
        secs += 1800.

    #plt.plot(Tair_out)
    #plt.show()

    ndim = 1

    # create file and write global attributes
    f = nc.Dataset(ofname, 'w', format='NETCDF4')
    f.description = '%s met data, created by Martin De Kauwe' % (data_type)
    f.history = "Created by: %s" % (os.path.basename(__file__))
    f.creation_date = "%s" % (datetime.datetime.now())
    f.contact = "mdekauwe@gmail.com"

    # set dimensions
    f.createDimension('time', n_timesteps)
    f.createDimension('z', ndim)
    f.createDimension('y', ndim)
    f.createDimension('x', ndim)
    #f.Conventions = "CF-1.0"

    # create variables
    time = f.createVariable('time', 'f8', ('time',))
    time.units = "seconds since %s-01-01 00:00:00" % (start_yr)
    time.long_name = "time"
    time.calendar = "standard"

    z = f.createVariable('z', 'f8', ('z',))
    z.long_name = "z"
    z.long_name = "z dimension"

    y = f.createVariable('y', 'f8', ('y',))
    y.long_name = "y"
    y.long_name = "y dimension"

    x = f.createVariable('x', 'f8', ('x',))
    x.long_name = "x"
    x.long_name = "x dimension"

    latitude = f.createVariable('latitude', 'f8', ('y', 'x',))
    latitude.units = "degrees_north"
    latitude.missing_value = -9999.
    latitude.long_name = "Latitude"

    longitude = f.createVariable('longitude', 'f8', ('y', 'x',))
    longitude.units = "degrees_east"
    longitude.missing_value = -9999.
    longitude.long_name = "Longitude"

    SWdown = f.createVariable('SWdown', 'f8', ('time', 'y', 'x',))
    SWdown.units = "W/m^2"
    SWdown.missing_value = -9999.
    SWdown.long_name = "Surface incident shortwave radiation"
    SWdown.CF_name = "surface_downwelling_shortwave_flux_in_air"

    Tair = f.createVariable('Tair', 'f8', ('time', 'z', 'y', 'x',))
    Tair.units = "K"
    Tair.missing_value = -9999.
    Tair.long_name = "Near surface air temperature"
    Tair.CF_name = "surface_temperature"

    Rainf = f.createVariable('Rainf', 'f8', ('time', 'y', 'x',))
    Rainf.units = "mm/s"
    Rainf.missing_value = -9999.
    Rainf.long_name = "Rainfall rate"
    Rainf.CF_name = "precipitation_flux"

    Snowf = f.createVariable('Snowf', 'f8', ('time', 'y', 'x',))
    Snowf.units = "mm/s"
    Snowf.missing_value = -9999.
    Snowf.long_name = "Snowfall rate"
    Snowf.CF_name = "snowfall_flux"

    Qair = f.createVariable('Qair', 'f8', ('time', 'z', 'y', 'x',))
    Qair.units = "kg/kg"
    Qair.missing_value = -9999.
    Qair.long_name = "Near surface specific humidity"
    Qair.CF_name = "surface_specific_humidity"

    Wind = f.createVariable('Wind', 'f8', ('time', 'z', 'y', 'x',))
    Wind.units = "m/s"
    Wind.missing_value = -9999.
    Wind.long_name = "Scalar windspeed" ;
    Wind.CF_name = "wind_speed"

    PSurf = f.createVariable('PSurf', 'f8', ('time', 'y', 'x',))
    PSurf.units = "Pa"
    PSurf.missing_value = -9999.
    PSurf.long_name = "Surface air pressure"
    PSurf.CF_name = "surface_air_pressure"

    LWdown = f.createVariable('LWdown', 'f8', ('time', 'y', 'x',))
    LWdown.units = "W/m^2"
    LWdown.missing_value = -9999.
    LWdown.long_name = "Surface incident longwave radiation"
    LWdown.CF_name = "surface_downwelling_longwave_flux_in_air"

    # write data to file
    x[:] = ndim
    y[:] = ndim
    z[:] = ndim
    time[:] = times
    latitude[:] = lat
    longitude[:] = lon

    SWdown[:,0,0] = SWdown_out.reshape(n_timesteps, ndim, ndim)
    LWdown[:,0,0] = LWdown_out.reshape(n_timesteps, ndim, ndim)
    Tair[:,0,0,0] = Tair_out.reshape(n_timesteps, ndim, ndim, ndim)
    Rainf[:,0,0] = Rainf_out.reshape(n_timesteps, ndim, ndim)
    Snowf[:,0,0] = Snowf_out.reshape(n_timesteps, ndim, ndim)
    Qair[:,0,0,0] = Qair_out.reshape(n_timesteps, ndim, ndim, ndim)
    Wind[:,0,0,0] = Wind_out.reshape(n_timesteps, ndim, ndim, ndim)
    PSurf[:,0,0] = PSurf_out.reshape(n_timesteps, ndim, ndim)

    f.close()

if __name__ == "__main__":


    vars = ["LWdown","PSurf","Qair","Rainf","SWdown","Snowf","Tair","Wind"]
    start_yr = 2016
    end_yr = 2019
    data_type = "AWAP"

    tmp_fn = glob.glob("SWdown_%d_*_*.nc" % (start_yr))[0]
    row = int(tmp_fn.split("_")[2])
    col = int(tmp_fn.split("_")[3].split(".")[0])

    main(vars, start_yr, end_yr,  data_type, row, col)
