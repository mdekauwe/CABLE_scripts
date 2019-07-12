#!/usr/bin/env python
"""
Generate a landsea mask for just SE Australia

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (12.07.2019)"
__email__ = "mdekauwe@gmail.com"


import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import sys

fname = "gswp3_landmask_nomissing.nc"
out_fname = "SE_AUS_gswp3_landmask_nomissing.nc"

ds = xr.open_dataset(fname)

# Subset by SE aus.
ds['landsea'] = ds.landsea.where((ds.lon >= 140) & (ds.lon <= 154) & \
                          (ds.lat >= -40) & (ds.lat <= -28), 1.0, 0.0)
#plt.imshow(ds['landsea'])
#plt.colorbar()
#plt.show()

ds.to_netcdf(out_fname)
ds.close()

dsx = xr.open_dataset(out_fname)
plt.imshow(dsx['landsea'])
plt.colorbar()
plt.show()
