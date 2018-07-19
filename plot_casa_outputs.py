#!/usr/bin/env python

"""
Plot some CASA outputs for the simulation

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

def plot_carbon_fluxes(tag, cycle, ds):

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
    ax1.plot(ds.Cgpp[:,0], label="Cf")

    ax2.set_title("CUE")
    ax2.plot(ds.Cnpp[:,0]/ds.Cgpp[:,0])

    ax3.set_title("LAI")
    ax3.plot(ds.glai[:,0])

    ax4.set_title("Vcmax")
    ax4.plot(ds.vcmax[:,0]*1E6)

    ax5.set_title("Allocation")
    ax5.plot(ds.fracCalloc[:,0], label="Af")
    ax5.plot(ds.fracCalloc[:,1], label="Aw")
    ax5.plot(ds.fracCalloc[:,2], label="Ar")
    ax5.legend(numpoints=1, loc="best")

    ax6.set_title("C labile")
    ax6.plot(ds.clabile[:,0])

    plot_fname = "%s_simulation_carbon_fluxes.pdf" % (tag)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)


def plot_nitrogen_fluxes(tag, cycle, ds):

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
    ax1.plot(ds.Nminfix[:,0], label="Cf")

    ax2.set_title("N deposition (g N/m^2/year)")
    ax2.plot(ds.Nmindep[:,0])

    ax3.set_title("N loss")
    ax3.plot(ds.Nminloss[:,0])

    ax4.set_title("N leach")
    ax4.plot(ds.Nminleach[:,0])

    ax5.set_title("N uptake")
    ax5.plot(ds.Nupland[:,0])

    ax6.set_title("N gross mineralisation")
    ax6.plot(ds.Nsmin[:,0])

    ax7.set_title("N net mineralisation")
    ax7.plot(ds.Nsnet[:,0])

    ax8.set_title("N immobilisation")
    ax8.plot(ds.Nsimm[:,0])

    ax9.set_title("Leaf N:C")
    ax9.plot(ds.nplant[:,0,0]/ds.cplant[:,0,0])


    plot_fname = "%s_simulation_nitrogen_fluxes.pdf" % (tag)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)

def plot_phosphorus_fluxes(tag, cycle, ds):

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

    ax1.set_title("P fixation (g N/m^2/year)")
    ax1.plot(ds.Nminfix[:,0], label="Cf")

    ax2.set_title("P deposition (g N/m^2/year)")
    ax2.plot(ds.Pdep[:,0])

    ax3.set_title("P loss")
    ax3.plot(ds.Ploss[:,0])

    ax4.set_title("P leach")
    ax4.plot(ds.Pleach[:,0])

    ax5.set_title("P uptake")
    ax5.plot(ds.Pupland[:,0])

    ax6.set_title("P gross mineralisation")
    ax6.plot(ds.Psmin[:,0])

    ax7.set_title("P net mineralisation")
    ax7.plot(ds.Psnet[:,0])

    ax8.set_title("P immobilisation")
    ax8.plot(ds.Psimm[:,0])

    ax9.set_title("Leaf N:P")
    ax9.plot(ds.nplant[:,0,0]/ds.cplant[:,0,0])


    plot_fname = "%s_simulation_phosphorus_fluxes.pdf" % (tag)
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

        fname = "outputs/%s_out_casa.nc" % (experiment_id)
        simulation = open_file(fname)

        plot_carbon_fluxes(tag, cycle, simulation)
        plot_nitrogen_fluxes(tag, cycle, simulation)
        plot_phosphorus_fluxes(tag, cycle, simulation)
