#!/usr/bin/env python

"""
Generate a single pixel netcdf timeseries from spatial inputs for CABLE testing

python extract_single_pixel.py "LWdown"
python extract_single_pixel.py "PSurf"
python extract_single_pixel.py "Qair"
python extract_single_pixel.py "Rainf"
python extract_single_pixel.py "SWdown"
python extract_single_pixel.py "Snowf"
python extract_single_pixel.py "Tair"
python extract_single_pixel.py "Wind"


That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (21.10.2018)"
__email__ = "mdekauwe@gmail.com"

import subprocess
import sys
import os

data_type = "AWAP"
fdir = "/g/data1a/w35/mgk576/research/AWAP_interpolation/interpolated"
vars = ["LWdown","PSurf","Qair","Rainf","SWdown","Snowf","Tair","Wind"]
start_yr = 2016
end_yr = 2019
row = 409 #410
col = 793

if len(sys.argv) < 2:
    raise TypeError("Expecting pft name to be supplied on cmd line!")
var = str(sys.argv[1])

for yr in range(start_yr, end_yr +1):

    if data_type == "GSWP3":
        fname = "%s/%s/%s.BC.%s.3hrMap.%d.nc" % (fdir, var, data_type, var, yr)
    else:
        fname = "%s/%s/%s.%s.3hr.%d.nc" % (fdir, var, data_type, var, yr)
    ofname = "%s_%d.nc" % (var, yr)
    qs_cmd = "ncks -dlat,%d,%d -dlon,%d,%d %s %s" % \
                (row, row, col, col, fname, ofname)
    error = subprocess.call(qs_cmd, shell=True)
    if error is 1:
        print("Job failed to submit")
