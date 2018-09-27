#!/usr/bin/env python

"""
Plot some CASA outputs for the spinup

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

def plot_carbon_fluxes(tag, cycle, ds, type, window=6):

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

    ax1 = fig.add_subplot(2,3,1)
    ax2 = fig.add_subplot(2,3,2)
    ax3 = fig.add_subplot(2,3,3)
    ax4 = fig.add_subplot(2,3,4)
    ax5 = fig.add_subplot(2,3,5)
    ax6 = fig.add_subplot(2,3,6)

    ax1.set_title("GPP")
    gpp = ds.Cgpp[:,0].to_dataframe()
    ax1.plot(gpp.rolling(window=window).mean())
    ax1.set_xticks(gpp.index.to_pydatetime())
    ax1.locator_params(tight=True, nbins=6)

    ax2.set_title("CUE")
    ds['cue'] = ds.Cnpp[:,0]/ds.Cgpp[:,0]
    cue = ds.cue.to_dataframe()
    ax2.plot(cue.rolling(window=window).mean())
    ax2.set_xticks(cue.index.to_pydatetime())
    ax2.locator_params(tight=True, nbins=6)

    ax3.set_title("LAI")
    lai = ds.glai[:,0].to_dataframe()
    ax3.plot(lai.rolling(window=window).mean())
    ax3.set_xticks(lai.index.to_pydatetime())
    ax3.locator_params(tight=True, nbins=6)

    ax4.set_title("Vcmax")
    vcmax = ds.vcmax[:,0].to_dataframe() * 1E6
    ax4.plot(vcmax.rolling(window=window).mean())
    ax4.set_xticks(vcmax.index.to_pydatetime())
    ax4.locator_params(tight=True, nbins=6)

    ax5.set_title("Allocation")
    af = ds.fracCalloc[:,0].squeeze(dim=["land"], drop=True).to_dataframe()
    aw = ds.fracCalloc[:,1].squeeze(dim=["land"], drop=True).to_dataframe()
    ar = ds.fracCalloc[:,2].squeeze(dim=["land"], drop=True).to_dataframe()
    ax5.plot(af.rolling(window=window).mean(), label="Af")
    ax5.plot(aw.rolling(window=window).mean(), label="Aw")
    ax5.plot(ar.rolling(window=window).mean(), label="Ar")
    ax5.legend(numpoints=1, loc="best")
    ax5.set_xticks(af.index.to_pydatetime())
    ax5.locator_params(tight=True, nbins=6)

    ax6.set_title("C labile")
    clabile = ds.clabile[:,0].to_dataframe()
    ax6.plot(clabile.rolling(window=window).mean())
    ax6.set_xticks(clabile.index.to_pydatetime())
    ax6.locator_params(tight=True, nbins=6)

    plot_fname = "%s_%s_carbon_fluxes.pdf" % (tag, type)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)


def plot_nitrogen_fluxes(tag, cycle, ds, type, window=6):

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

    ax1 = fig.add_subplot(3,3,1)
    ax2 = fig.add_subplot(3,3,2)
    ax3 = fig.add_subplot(3,3,3)
    ax4 = fig.add_subplot(3,3,4)
    ax5 = fig.add_subplot(3,3,5)
    ax6 = fig.add_subplot(3,3,6)
    ax7 = fig.add_subplot(3,3,7)
    ax8 = fig.add_subplot(3,3,8)
    ax9 = fig.add_subplot(3,3,9)

    ax1.set_title("N fixation (g N/m^2/year)")
    nfix = ds.Nminfix[:,0].to_dataframe()
    ax1.plot(nfix.rolling(window=window).mean())
    ax1.set_xticks(nfix.index.to_pydatetime())
    ax1.locator_params(tight=True, nbins=6)

    ax2.set_title("N deposition (g N/m^2/year)")
    ndep = ds.Nmindep[:,0].to_dataframe()
    ax2.plot(ndep.rolling(window=window).mean())
    ax2.set_xticks(ndep.index.to_pydatetime())
    ax2.locator_params(tight=True, nbins=6)

    ax3.set_title("N loss")
    nloss = ds.Nminloss[:,0].to_dataframe()
    ax3.plot(nloss.rolling(window=window).mean())
    ax3.set_xticks(nloss.index.to_pydatetime())
    ax3.locator_params(tight=True, nbins=6)

    ax4.set_title("N leach")
    nleach = ds.Nminloss[:,0].to_dataframe()
    ax4.plot(nleach.rolling(window=window).mean())
    ax4.set_xticks(nleach.index.to_pydatetime())
    ax4.locator_params(tight=True, nbins=6)

    ax5.set_title("N uptake")
    nup = ds.Nupland[:,0].to_dataframe()
    ax5.plot(nup.rolling(window=window).mean())
    ax5.set_xticks(nup.index.to_pydatetime())
    ax5.locator_params(tight=True, nbins=6)

    ax6.set_title("N gross mineralisation")
    ngross = ds.Nsmin[:,0].to_dataframe()
    ax6.plot(ngross.rolling(window=window).mean())
    ax6.set_xticks(ngross.index.to_pydatetime())
    ax6.locator_params(tight=True, nbins=6)

    ax7.set_title("N net mineralisation")
    nmin = ds.Nsnet[:,0].to_dataframe()
    ax7.plot(nmin.rolling(window=window).mean())
    ax7.set_xticks(nmin.index.to_pydatetime())
    ax7.locator_params(tight=True, nbins=6)

    ax8.set_title("N immobilisation")
    nimmob = ds.Nsimm[:,0].to_dataframe()
    ax8.plot(nimmob.rolling(window=window).mean())
    ax8.set_xticks(nimmob.index.to_pydatetime())
    ax8.locator_params(tight=True, nbins=6)

    ax9.set_title("Leaf N:C")
    ds['leaf_nc'] = ds.nplant[:,0,0]/ds.cplant[:,0,0]
    leaf_nc = ds.leaf_nc.to_dataframe()
    ax9.plot(leaf_nc.rolling(window=window).mean())
    ax9.set_xticks(leaf_nc.index.to_pydatetime())
    ax9.locator_params(tight=True, nbins=6)

    plot_fname = "%s_%s_nitrogen_fluxes_spinup.pdf" % (tag, type)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)

def plot_phosphorus_fluxes(tag, cycle, ds, type, window=6):

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

    ax1 = fig.add_subplot(3,3,1)
    ax2 = fig.add_subplot(3,3,2)
    ax3 = fig.add_subplot(3,3,3)
    ax4 = fig.add_subplot(3,3,4)
    ax5 = fig.add_subplot(3,3,5)
    ax6 = fig.add_subplot(3,3,6)
    ax7 = fig.add_subplot(3,3,7)
    ax8 = fig.add_subplot(3,3,8)


    ax1.set_title("P deposition (g N/m^2/year)")
    pdep = ds.Pdep[:,0].to_dataframe()
    ax1.plot(pdep.rolling(window=window).mean())
    ax1.set_xticks(pdep.index.to_pydatetime())
    ax1.locator_params(tight=True, nbins=6)

    ax2.set_title("P loss")
    ploss = ds.Ploss[:,0].to_dataframe()
    ax2.plot(ploss.rolling(window=window).mean())
    ax2.set_xticks(ploss.index.to_pydatetime())
    ax2.locator_params(tight=True, nbins=6)

    ax3.set_title("P leach")
    pleach = ds.Pleach[:,0].to_dataframe()
    ax3.plot(pleach.rolling(window=window).mean())
    ax3.set_xticks(pleach.index.to_pydatetime())
    ax3.locator_params(tight=True, nbins=6)

    ax4.set_title("P uptake")
    pup = ds.Pupland[:,0].to_dataframe()
    ax4.plot(pup.rolling(window=window).mean())
    ax4.set_xticks(pup.index.to_pydatetime())
    ax4.locator_params(tight=True, nbins=6)

    ax5.set_title("P gross mineralisation")
    pgross = ds.Psmin[:,0].to_dataframe()
    ax5.plot(pgross.rolling(window=window).mean())
    ax5.set_xticks(pgross.index.to_pydatetime())
    ax5.locator_params(tight=True, nbins=6)

    ax6.set_title("P net mineralisation")
    pnmin = ds.Psnet[:,0].to_dataframe()
    ax6.plot(pnmin.rolling(window=window).mean())
    ax6.set_xticks(pnmin.index.to_pydatetime())
    ax6.locator_params(tight=True, nbins=6)

    ax7.set_title("P immobilisation")
    pimmob = ds.Psimm[:,0].to_dataframe()
    ax7.plot(pimmob.rolling(window=window).mean())
    ax7.set_xticks(pimmob.index.to_pydatetime())
    ax7.locator_params(tight=True, nbins=6)

    ax8.set_title("Leaf N:P")
    ds['leaf_np'] = ds.nplant[:,0,0] / ds.cplant[:,0,0]
    leaf_np = ds.leaf_np.to_dataframe()
    ax8.plot(leaf_np.rolling(window=window).mean())
    ax8.set_xticks(leaf_np.index.to_pydatetime())
    ax8.locator_params(tight=True, nbins=6)

    plot_fname = "%s_%s_phosphorus_fluxes.pdf" % (tag, type)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)


def open_file(fname):

    ds = xr.open_dataset(fname)
    dates = pd.date_range('1/1/1750', periods=len(ds.time), freq='M')
    ds['time'] = dates

    return ds



if __name__ == "__main__":

    type = "transient"

    for cycle in ["C", "CN", "CNP"]:
    #for cycle in ["CN"]:

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


        fname = "outputs/%s_out_CASA_transient.nc" % (experiment_id)
        ds = open_file(fname)

        plot_carbon_fluxes(tag, cycle, ds, type)
        plot_nitrogen_fluxes(tag, cycle, ds, type)
        plot_phosphorus_fluxes(tag, cycle, ds, type)
