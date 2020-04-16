#!/usr/bin/env python

# Belinda's pixels.

import numpy as np


def get_row_col(lat0, lon0):

    ncols = 841
    nrows = 681

    aus_latvals = np.zeros((nrows, ncols))
    aus_lonvals = np.zeros((nrows, ncols))

    # lower left corner
    llcrnrlon = 112.
    llcrnrlat = -44.

    for i in range(nrows):
        for j in range(ncols):

            lat = llcrnrlat + (float(i) * 0.05)
            lon = llcrnrlon + (float(j) * 0.05)
            aus_latvals[i,j] = lat
            aus_lonvals[i,j] = lon

    dist_sq_min = 1.0e30
    for i in range(nrows):
        for j in range(ncols):
            latval = aus_latvals[i,j]
            lonval = aus_lonvals[i,j]
            dist_sq = (latval - lat0)**2 + (lonval - lon0)**2
            if dist_sq < dist_sq_min:
                i_idx, j_idx, dist_sq_min = i, j, dist_sq


    print("%f:%f <-> %d:%d" % (lat0, lon0, i_idx, j_idx))

#lats = [-30.4172,-30.416,-30.417,-30.4915,-30.4919,-30.491,-30.4174,-30.4176]
#lons = [151.616,151.6153,151.6155,151.6438,151.643,151.6425,151.6275,151.6234]

#for lat,lon in zip(lats, lons):
#    get_row_col(lat, lon)

lats = [-31.0927,-31.0927,-30.0]
lons = [150.9320,142.5, 141.4]

for lat,lon in zip(lats, lons):
    get_row_col(lat, lon)
