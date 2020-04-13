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

def main(fdir, var, start_yr, end_yr, row, col, data_type, ofname):

    out = np.zeros(0)

    for yr in range(start_yr, end_yr +1):
        print(yr)

        if data_type == "GSWP3":
            fname = "%s/%s.BC.%s.3hrMap.%d.nc" % (var, data_type, var, yr)
        else:
            fname = "%s/%s.%s.3hr.%d.nc" % (var, data_type, var, yr)
        fpath = os.path.join(fdir, fname)
        ds = xr.open_dataset(fpath)
        vals = ds[var][:,row,col].values
        lat = ds["lat"][row].values
        lon = ds["lon"][col].values

        out = np.append(out, vals)

        ds.close()

    out.tofile("%s_%d.csv" % (var, start_yr), sep=',')

if __name__ == "__main__":

    # Expecting PFT to be supplied on cmd line, e.g.
    # $ python extract_single_pixel_timseries_from_AWAP.py "rf"
    if len(sys.argv) < 2:
        raise TypeError("Expecting pft name to be supplied on cmd line!")
    var = sys.argv[1]
    start_yr = int(sys.argv[2])
    end_yr = int(sys.argv[2])


    fdir = "/g/data1a/w35/mgk576/research/AWAP_interpolation/interpolated"
    #vars = ["LWdown","PSurf","Qair","Rainf","SWdown","Snowf","Tair","Wind"]

    row = 409 #410
    col = 793

    data_type = "AWAP"
    ofname = "single_pixel_spinup_%s_%d_%d_%d_%d.nc" % \
                (data_type, row, col, start_yr, end_yr)
    main(fdir, vars, start_yr, end_yr, row, col, data_type, ofname)
