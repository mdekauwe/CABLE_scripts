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

fname = "/Users/mdekauwe/Desktop/gswp3_landmask_nomissing.nc"
out_fname = "/Users/mdekauwe/Desktop/SE_AUS_gswp3_landmask_nomissing.nc"
grid_fname = "/Users/mdekauwe/Desktop/gridinfo_mmy_MD_elev_orig_std_avg-sand_mask.nc"

ds = xr.open_dataset(fname)
ds_grid = xr.open_dataset(grid_fname)
iveg = ds_grid.iveg.values


# Subset by SE aus.
ds['landsea'] = ds.landsea.where((ds.lon >= 140) & (ds.lon <= 154) &
                                 (ds.lat >= -40) & (ds.lat <= -28),
                                 1.0)

# The logic above means we could re-introduce a masked sea pixel, so mask
# where iveg is nan
ds['landsea'] = ds.landsea.where(~np.isnan(iveg), 1.0)

ds.to_netcdf(out_fname)
ds.close()

dsx = xr.open_dataset(out_fname)
plt.imshow(dsx['landsea'])
plt.colorbar()
plt.show()
