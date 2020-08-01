#!/usr/bin/env python

"""
Check for equilibrium by splitting into three latitudinal groups

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (01.08.2020)"
__email__ = "mdekauwe@gmail.com"

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import optparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

def cmd_line_parser():

    p = optparse.OptionParser()
    p.add_option("--f1", default="", help="CABLE prev file")
    p.add_option("--f2", default="", help="CABLE new file")
    p.add_option("-n", default="0", help="CO2 concentration")
    options, args = p.parse_args()

    return (options.f1, options.f2, int(options.n))

def split_data(var):
    kg_2_g = 1000.
    var_split = np.array_split(var, 3)
    top = np.nansum(var_split[0].values) #* kg_2_g
    middle = np.nansum(var_split[1].values) #* kg_2_g
    bottom = np.nansum(var_split[2].values) #* kg_2_g

    return (np.array([top, middle, bottom]))


tol_npp = 5E-6 # kg, delta < 10^-4 g C m-2, Xia et al. 2013
tol_plant = 0.01 # delta steady-state carbon (%), Xia et al. 2013
tol_soil = 0.01
tol_pass = 0.5

(fn_prev, fn_new, num) = cmd_line_parser()

ds_prev = xr.open_dataset(fn_prev)
lats = ds_prev.latitude[:,0]
npp_prev = ds_prev.NPP.mean(axis=0).mean(axis=1)
leaf_prev = ds_prev.PlantCarbLeaf.mean(axis=0).mean(axis=1)
wood_prev = ds_prev.PlantCarbWood.mean(axis=0).mean(axis=1)
root_prev = ds_prev.PlantCarbFineRoot.mean(axis=0).mean(axis=1)
fast_prev = ds_prev.SoilCarbFast.mean(axis=0).mean(axis=1)
slow_prev = ds_prev.SoilCarbSlow.mean(axis=0).mean(axis=1)
passive_prev = ds_prev.SoilCarbPassive.mean(axis=0).mean(axis=1)

(npp_prev) = split_data(npp_prev)
(leaf_prev) = split_data(leaf_prev)
(wood_prev) = split_data(wood_prev)
(root_prev) = split_data(root_prev)
(fast_prev) = split_data(fast_prev)
(slow_prev) = split_data(slow_prev)
(passive_prev) = split_data(passive_prev)

ds_new = xr.open_dataset(fn_new)
lats = ds_new.latitude[:,0]
npp_new = ds_new.NPP.mean(axis=0).mean(axis=1)
leaf_new = ds_new.PlantCarbLeaf.mean(axis=0).mean(axis=1)
wood_new = ds_new.PlantCarbWood.mean(axis=0).mean(axis=1)
root_new = ds_new.PlantCarbFineRoot.mean(axis=0).mean(axis=1)
fast_new = ds_new.SoilCarbFast.mean(axis=0).mean(axis=1)
slow_new = ds_new.SoilCarbSlow.mean(axis=0).mean(axis=1)
passive_new = ds_new.SoilCarbPassive.mean(axis=0).mean(axis=1)

(npp_new) = split_data(npp_new)
(leaf_new) = split_data(leaf_new)
(wood_new) = split_data(wood_new)
(root_new) = split_data(root_new)
(fast_new) = split_data(fast_new)
(slow_new) = split_data(slow_new)
(passive_new) = split_data(passive_new)

out_fname = "stability_log.txt"
if os.path.exists(out_fname):
    of = open(out_fname, 'a')
else:
    of = open(out_fname, 'w')
    print("N,equilibrium,delta_npp,delta_plant,delta_soil", file=of)

delta_npp = np.zeros(3)
delta_plant = np.zeros(3)
delta_soil = np.zeros(3)
pass_stabile = np.zeros(3).astype(int)
for i in range(3):

    delta_npp[i] = np.fabs(npp_new[i] - npp_prev[i])

    cplant_new = leaf_new + wood_new + root_new
    cplant_prev = leaf_prev + wood_prev + root_prev
    delta_plant[i] = np.fabs((cplant_new[i] - cplant_prev[i]) / cplant_new[i])

    csoil_new = fast_new + slow_new + passive_new
    csoil_prev = fast_prev + slow_prev + passive_prev
    delta_soil[i] = np.fabs((csoil_new[i] - csoil_prev[i]) / csoil_new[i])

    #print(i, cplant_new[i], cplant_prev[i], delta_plant[i])
    #print(i, csoil_new[i], csoil_prev[i], delta_soil[i])
    #print(i, npp_new[i], npp_prev[i], delta_npp[i])

in_equilibrium = False
if ( (delta_npp[0] < tol_npp) and
     (delta_npp[1] < tol_npp) and
     (delta_npp[2] < tol_npp) and
     (delta_plant[0] < tol_plant) and
     (delta_plant[1] < tol_plant) and
     (delta_plant[2] < tol_plant) and
     (delta_soil[0] < tol_soil) and
     (delta_soil[1] < tol_soil) and
     (delta_soil[2] < tol_soil) ):
    in_equilibrium = True
    print("1")
else:
    print("0")
print(num, in_equilibrium, np.mean(delta_npp), np.mean(delta_plant),
      np.mean(delta_soil), file=of)

of.close()
