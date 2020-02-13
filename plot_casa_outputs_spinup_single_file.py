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

def plot_spinup_state(ds):

    cf = ds.cplant[:,0,0]
    cw = ds.cplant[:,1,0]
    cr = ds.cplant[:,2,0]
    ca = ds.csoil[:,0,0]
    cs = ds.csoil[:,1,0]
    cp = ds.csoil[:,2,0]
    gpp = ds.Cgpp[:,0]
    npp = ds.Cnpp[:,0]
    nep = ds.Cnep[:,0]

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

    ax1.set_title("$C_{\mathrm{f}}$ (g C m$^{-2}$ d$^{-1}$)")
    ax1.plot(cf)

    ax2.set_title("$C_{\mathrm{w}}$ (g C m$^{-2}$ d$^{-1}$)")
    ax2.plot(cw)

    ax3.set_title("$C_{\mathrm{r}}$ (g C m$^{-2}$ d$^{-1}$)")
    ax3.plot(cr)

    ax4.set_title("$C_{\mathrm{active}}$ (g C m$^{-2}$ d$^{-1}$)")
    ax4.plot(ca)

    ax5.set_title("$C_{\mathrm{slow}}$ (g C m$^{-2}$ d$^{-1}$)")
    ax5.plot(cs)

    ax6.set_title("$C_{\mathrm{passive}}$ (g C m$^{-2}$ d$^{-1}$)")
    ax6.plot(cp)

    ax7.set_title("GPP (g C m$^{-2}$ d$^{-1}$)")
    ax7.plot(gpp)

    ax8.set_title("NPP (g C m$^{-2}$ d$^{-1}$)")
    ax8.plot(npp)

    ax9.set_title("NEP (g C m$^{-2}$ d$^{-1}$)")
    ax9.plot(nep)

    #plot_fname = "%s_spinup_carbon_pools.pdf" % (cycle)
    plot_fname = "spinup_carbon_pools_single_file.png"
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, plot_fname), dpi=150, bbox_inches='tight',
                pad_inches=0.1)

if __name__ == "__main__":

    fname = "outputs/AU-Tum_C_out_casa_spin_1.nc"
    ds = xr.open_dataset(fname, decode_times=False)
    plot_spinup_state(ds)
