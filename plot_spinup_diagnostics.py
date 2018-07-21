#!/usr/bin/env python

"""
Plot some diagnostics to see if the model has hit steady-state

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (17.07.2018)"
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

def plot_plant(tag, cycle, cf, cw, cr, nf, nw, nr, pf, pw, pr):

    fig = plt.figure(figsize=(15,10))
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

    ax1 = fig.add_subplot(3,4,1)
    ax2 = fig.add_subplot(3,4,2)
    ax3 = fig.add_subplot(3,4,3)
    ax4 = fig.add_subplot(3,4,4)

    ax5 = fig.add_subplot(3,4,5)
    ax6 = fig.add_subplot(3,4,6)
    ax7 = fig.add_subplot(3,4,7)
    ax8 = fig.add_subplot(3,4,8)

    ax9 = fig.add_subplot(3,4,9)
    ax10 = fig.add_subplot(3,4,10)
    ax11 = fig.add_subplot(3,4,11)
    ax12 = fig.add_subplot(3,4,12)

    ax1.set_title("Cf")
    ax1.plot(cf, label="Cf")

    ax2.set_title("Cw")
    ax2.plot(cw)

    ax3.set_title("Cr")
    ax3.plot(cr)

    ax4.set_title("Cplant")
    ax4.plot(cf + cw + cr)

    ax5.set_title("Nf")
    ax5.plot(nf, label="nf")

    ax6.set_title("Nr")
    ax6.plot(nw)

    ax7.set_title("Nr")
    ax7.plot(nr)

    ax8.set_title("Nplant")
    ax8.plot(nf+nw+nr)

    ax9.set_title("Pf")
    ax9.plot(pf, label="pf")

    ax10.set_title("Pr")
    ax10.plot(pw)

    ax11.set_title("Pr")
    ax11.plot(pr)

    ax12.set_title("Pplant")
    ax12.plot(pf + pw + pr)

    plot_fname = "%s_spinup_plant.pdf" % (tag)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)

def plot_soil(tag, cycle, cfast, cslow, cpassive, nfast, nslow, npassive,
              pfast, pslow, ppassive):

    fig = plt.figure(figsize=(15,10))
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

    ax1 = fig.add_subplot(3,4,1)
    ax2 = fig.add_subplot(3,4,2)
    ax3 = fig.add_subplot(3,4,3)
    ax4 = fig.add_subplot(3,4,4)

    ax5 = fig.add_subplot(3,4,5)
    ax6 = fig.add_subplot(3,4,6)
    ax7 = fig.add_subplot(3,4,7)
    ax8 = fig.add_subplot(3,4,8)

    ax9 = fig.add_subplot(3,4,9)
    ax10 = fig.add_subplot(3,4,10)
    ax11 = fig.add_subplot(3,4,11)
    ax12 = fig.add_subplot(3,4,12)

    ax1.set_title("C fast")
    ax1.plot(cfast, label="cfast")

    ax2.set_title("C slow")
    ax2.plot(cslow)

    ax3.set_title("C passive")
    ax3.plot(cpassive)

    ax4.set_title("C soil")
    ax4.plot(cfast + cslow + cpassive)

    ax5.set_title("N fast")
    ax5.plot(nfast)

    ax6.set_title("N slow")
    ax6.plot(nslow)

    ax7.set_title("N passive")
    ax7.plot(npassive)

    ax8.set_title("N plant")
    ax8.plot(nfast+nslow+npassive)

    ax9.set_title("P fast")
    ax9.plot(pfast, label="pfast")

    ax10.set_title("P slow")
    ax10.plot(pslow)

    ax11.set_title("P passive")
    ax11.plot(ppassive)

    ax12.set_title("P soil")
    ax12.plot(pfast+pslow+ppassive)

    plot_fname = "%s_spinup_soil.pdf" % (tag)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)

def open_file(fname):
    return xr.open_dataset(fname)



if __name__ == "__main__":

    for cycle in ["C", "CN", "CNP"]:

        print(cycle)

        if cycle == "C":
            experiment_id = "Cumberland_C"
            tag = "C"
        elif cycle == "CN":
            experiment_id = "Cumberland_CN"
            tag = "CN"
        elif cycle == "CNP":
            experiment_id = "Cumberland_CNP"
            tag = "CNP"

        fname = "outputs/%s_out_CASA_zero.nc" % (experiment_id)
        zero = open_file(fname)

        cf = zero.cplant[:,0,0].values
        cw = zero.cplant[:,1,0].values
        cr = zero.cplant[:,2,0].values

        nf = zero.nplant[:,0,0].values
        nw = zero.nplant[:,1,0].values
        nr = zero.nplant[:,2,0].values

        pf = zero.pplant[:,0,0].values
        pw = zero.pplant[:,1,0].values
        pr = zero.pplant[:,2,0].values

        cfast = zero.csoil[:,0,0].values
        cslow = zero.csoil[:,1,0].values
        cpassive = zero.csoil[:,2,0].values

        nfast = zero.nsoil[:,0,0].values
        nslow = zero.nsoil[:,1,0].values
        npassive = zero.nsoil[:,2,0].values

        pfast = zero.psoil[:,0,0].values
        pslow = zero.psoil[:,1,0].values
        ppassive = zero.psoil[:,2,0].values

        files = glob.glob("outputs/%s_out_CASA_ccp*.nc" % (experiment_id))
        for fname in sorted(files):
            ds = open_file(fname)

            cf = np.append(cf, ds.cplant[:,0,0].values)
            cw = np.append(cw, ds.cplant[:,1,0].values)
            cr = np.append(cr, ds.cplant[:,2,0].values)

            nf = np.append(nf, ds.nplant[:,0,0].values)
            nw = np.append(nw, ds.nplant[:,1,0].values)
            nr = np.append(nr, ds.nplant[:,2,0].values)

            pf = np.append(pf, ds.pplant[:,0,0].values)
            pw = np.append(pw, ds.pplant[:,1,0].values)
            pr = np.append(pr, ds.pplant[:,2,0].values)

            cfast = np.append(cfast, ds.csoil[:,0,0].values)
            cslow = np.append(cslow, ds.csoil[:,1,0].values)
            cpassive = np.append(cpassive, ds.csoil[:,2,0].values)

            nfast = np.append(nfast, ds.nsoil[:,0,0].values)
            nslow = np.append(nslow, ds.nsoil[:,1,0].values)
            npassive = np.append(npassive, ds.nsoil[:,2,0].values)

            pfast = np.append(pfast, ds.psoil[:,0,0].values)
            pslow = np.append(pslow, ds.psoil[:,1,0].values)
            ppassive = np.append(ppassive, ds.psoil[:,2,0].values)

        plot_plant(tag, cycle, cf, cw, cr, nf, nw, nr, pf, pw, pr)
        plot_soil(tag, cycle, cfast, cslow, cpassive, nfast, nslow, npassive,
                  pfast, pslow, ppassive)
