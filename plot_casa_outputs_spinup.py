#!/usr/bin/env python

"""
CASA spinup diagnostics...C plant pools, C soil pools and CO2, Ndep and Pdep

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (07.02.2020)"
__email__ = "mdekauwe@gmail.com"

import matplotlib.pyplot as plt
import sys
import datetime as dt
import pandas as pd
import numpy as np
from matplotlib.ticker import FixedLocator
import os
import xarray as xr
import glob

def plot_spinup_state(cycle, cf, cw, cr, ca, cs, cp, co2, ndep, pdep):

    fig = plt.figure(figsize=(15,6))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.3)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(3,3,1)
    ax2 = fig.add_subplot(3,3,2)
    ax3 = fig.add_subplot(3,3,3)
    ax4 = fig.add_subplot(3,3,4)
    ax5 = fig.add_subplot(3,3,5)
    ax6 = fig.add_subplot(3,3,6)
    ax7 = fig.add_subplot(3,3,7)
    ax8 = fig.add_subplot(3,3,8)
    ax9 = fig.add_subplot(3,3,9)

    co2, ndep, pdep

    ax1.set_title("C foliage")
    ax1.plot(cf)

    ax2.set_title("C wood")
    ax2.plot(cw)

    ax3.set_title("C root")
    ax3.plot(cr)

    ax4.set_title("C soil (active)")
    ax4.plot(ca)

    ax5.set_title("C soil (slow)")
    ax5.plot(cs)

    ax6.set_title("C soil (passive)")
    ax6.plot(cp)

    ax7.set_title("CO2")
    ax7.plot(co2)

    ax8.set_title("N deposition")
    ax8.plot(ndep)

    ax9.set_title("P deposition")
    ax9.plot(pdep)

    plot_fname = "%s_spinup_carbon_pools.pdf" % (cycle)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)


def get_data(casa_fname, cable_fname):

    ds = xr.open_dataset(casa_fname)
    csoil = ds.csoil[-1,:,0]
    cplant = ds.cplant[-1,:,0]
    ndep = ds.Nmindep[-1,0]
    pdep = ds.Pdep[-1,0]

    ds = xr.open_dataset(cable_fname, decode_times=False)
    co2 = ds.CO2air[-1,0,0]


    return (cplant, csoil, co2, ndep, pdep)


if __name__ == "__main__":


    site = "AU-Tum"
    year = 1850

    cf = []
    cw = []
    cr = []

    ca = []
    cs = []
    cp = []

    co2 = []
    ndep = []
    pdep = []

    #for cycle in ["C", "CN", "CNP"]:
    #for cycle in ["C"]:
    for cycle in ["C", "CN"]:

        files = glob.glob("outputs/%s_%s_out_casa_spin_*.nc" % (site, cycle))
        nums = sorted([int(f.split(".")[0].split("_")[-1]) for f in files])

        for i in range(nums[0], nums[-1]):
            print(cycle, i, nums[-1])
            casa_fn = "outputs/%s_%s_out_casa_spin_%d.nc" % (site, cycle, i)
            cable_fn = "outputs/%s_%s_out_cable_spin_%d.nc" % (site, cycle, i)

            (cplant, csoil, co2x, ndepx, pdepx) = get_data(casa_fn, cable_fn)

            cf.append(cplant[0])
            cw.append(cplant[1])
            cr.append(cplant[2])

            ca.append(csoil[0])
            cs.append(csoil[1])
            cp.append(csoil[2])

            co2.append(co2x)
            ndep.append(ndepx)
            pdep.append(pdepx)

        plot_spinup_state(cycle, cf, cw, cr, ca, cs, cp, co2, ndep, pdep)
