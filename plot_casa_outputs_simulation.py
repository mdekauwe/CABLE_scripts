#!/usr/bin/env python

"""
Plot some CASA/CABLE outputs for the simulation

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (21.07.2018)"
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

def plot_carbon_fluxes(cycle, ds, ds_cable):

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

    ax1.set_title("GPP/NPP (g C m$^{-2}$ d$^{-1}$)")
    ax1.plot(ds.time, ds.Cgpp[:,0], label="GPP")
    ax1.plot(ds.time, ds.Cnpp[:,0], label="NPP")
    ax1.legend(numpoints=1, loc="best")

    ax2.set_title("CUE (-)")
    ax2.plot(ds.time, ds.Cnpp[:,0]/ds.Cgpp[:,0])

    ax3.set_title("LAI (m$^{2}$ m$^{-2}$)")
    ax3.plot(ds.time, ds.glai[:,0])

    ax4.set_title("Vcmax ($\mathrm{\mu}$mol m$^{-2}$ s$^{-1}$)")
    ax4.plot(ds.time, ds.vcmax[:,0]*1E6)

    ax5.set_title("Allocation (-)")
    ax5.plot(ds.time, ds.fracCalloc[:,0], label="Af")
    ax5.plot(ds.time, ds.fracCalloc[:,1], label="Aw")
    ax5.plot(ds.time, ds.fracCalloc[:,2], label="Ar")
    ax5.legend(numpoints=1, loc="best")

    ax6.set_title("C labile (g C m$^{-2}$ d$^{-1}$)")
    ax6.plot(ds.time, ds.clabile[:,0])

    ax7.set_title("CO$_2$ ($\mathrm{\mu}$mol mol$^{-1}$)")
    ax7.plot(ds_cable.time, ds_cable.CO2air[:,0,0])

    ax8.set_title("Cplant (g C m$^{-2}$ d$^{-1}$)")
    ax8.plot(ds.time, ds.cplant[:,0], label="Cf")
    ax8.plot(ds.time, ds.cplant[:,1], label="Cw")
    ax8.plot(ds.time, ds.cplant[:,2], label="Cr")
    ax8.legend(numpoints=1, loc="best")

    ax9.set_title("Csoil (g C m$^{-2}$ d$^{-1}$)")
    ax9.plot(ds.time, ds.csoil[:,0], label="Active")
    ax9.plot(ds.time, ds.csoil[:,1], label="Slow")
    ax9.plot(ds.time, ds.csoil[:,2], label="Passive")
    ax9.legend(numpoints=1, loc="best")

    plot_fname = "%s_simulation_carbon_fluxes_and_pools.pdf" % (cycle)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)


def plot_nitrogen_fluxes(cycle, ds):

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

    ax1.set_title("N deposition (g N/m^2/year)")
    ax1.plot(ds.time, ds.Nmindep[:,0])

    ax2.set_title("N fixation (g N/m^2/year)")
    ax2.plot(ds.time, ds.Nminfix[:,0], label="Nfix")

    ax3.set_title("N loss (g N/m^2/year)")
    ax3.plot(ds.time, ds.Nminloss[:,0])

    ax4.set_title("N leach")
    ax4.plot(ds.time, ds.Nminleach[:,0])

    ax5.set_title("N uptake")
    ax5.plot(ds.time, ds.Nupland[:,0])

    ax6.set_title("N gross mineralisation")
    ax6.plot(ds.time, ds.Nsmin[:,0])

    ax7.set_title("N net mineralisation")
    ax7.plot(ds.time, ds.Nsnet[:,0])

    ax8.set_title("N immobilisation")
    ax8.plot(ds.time, ds.Nsimm[:,0])

    ax9.set_title("Plant N:C")
    ax9.plot(ds.time, ds.nplant[:,0]/ds.cplant[:,0], label="Leaf")
    ax9.plot(ds.time, ds.nplant[:,1]/ds.cplant[:,1], label="Wood")
    ax9.plot(ds.time, ds.nplant[:,2]/ds.cplant[:,2], label="Root")
    ax9.legend(numpoints=1, loc="best")


    plot_fname = "%s_simulation_nitrogen_fluxes.pdf" % (cycle)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)

def plot_phosphorus_fluxes(cycle, ds):

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
    ax1.plot(ds.time, ds.Pdep[:,0])

    ax2.set_title("P loss")
    ax2.plot(ds.time, ds.Ploss[:,0])

    ax3.set_title("P leach")
    ax3.plot(ds.time, ds.Pleach[:,0])

    ax4.set_title("P uptake")
    ax4.plot(ds.time, ds.Pupland[:,0])

    ax5.set_title("P gross mineralisation")
    ax5.plot(ds.time, ds.Psmin[:,0])

    ax6.set_title("P net mineralisation")
    ax6.plot(ds.time, ds.Psnet[:,0])

    ax7.set_title("P immobilisation")
    ax7.plot(ds.time, ds.Psimm[:,0])

    ax8.set_title("Plant P:N")
    ax8.plot(ds.time, ds.pplant[:,0]/ds.nplant[:,0], label="Leaf")
    ax8.plot(ds.time, ds.pplant[:,1]/ds.nplant[:,1], label="Wood")
    ax8.plot(ds.time, ds.pplant[:,2]/ds.nplant[:,2], label="Root")
    ax8.legend(numpoints=1, loc="best")

    plot_fname = "%s_simulation_phosphorus_fluxes.pdf" % (cycle)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)

def plot_cnp_states(cycle, ds):

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


    ax1.set_title("Carbon")
    ax1.plot(ds.time, ds.cplant[:,0], label="Cf")
    ax1.plot(ds.time, ds.cplant[:,1], label="Cr")
    ax1.plot(ds.time, ds.cplant[:,2], label="Cw")
    ax1.legend(numpoints=1, loc="best")

    ax2.set_title("N plant")
    ax2.plot(ds.time, ds.nplant[:,0], label="Cf")
    ax2.plot(ds.time, ds.nplant[:,1], label="Cr")
    ax2.plot(ds.time, ds.nplant[:,2], label="Cw")
    ax2.legend(numpoints=1, loc="best")

    ax3.set_title("P plant")
    ax3.plot(ds.time, ds.pplant[:,0], label="Nf")
    ax3.plot(ds.time, ds.pplant[:,1], label="Nr")
    ax3.plot(ds.time, ds.pplant[:,2], label="Nw")
    ax3.legend(numpoints=1, loc="best")

    ax4.set_title("C soil")
    ax4.plot(ds.time, ds.csoil[:,0], label="Cf")
    ax4.plot(ds.time, ds.csoil[:,1], label="Cr")
    ax4.plot(ds.time, ds.csoil[:,2], label="Cw")
    ax4.legend(numpoints=1, loc="best")

    ax5.set_title("N soil")
    ax5.plot(ds.time, ds.nsoil[:,0], label="Nf")
    ax5.plot(ds.time, ds.nsoil[:,1], label="Nr")
    ax5.plot(ds.time, ds.nsoil[:,2], label="Nw")
    ax5.legend(numpoints=1, loc="best")

    ax6.set_title("P soil")
    ax6.plot(ds.time, ds.psoil[:,0], label="Pf")
    ax6.plot(ds.time, ds.psoil[:,1], label="Pr")
    ax6.plot(ds.time, ds.psoil[:,2], label="Pw")
    ax6.legend(numpoints=1, loc="best")

    plot_fname = "%s_simulation_state.pdf" % (cycle)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)

def open_casa_and_add_time(fname, start_date):
    ds = xr.open_dataset(fname)
    N = len(ds.Nupland)
    first_date = pd.to_datetime(start_date)
    time = first_date + pd.to_timedelta(np.arange(N), 'D')
    ds['time'] = time

    return (ds)

def open_cable_and_add_time(fname, start_date):
    ds = xr.open_dataset(fname, decode_times=False)
    N = len(ds.CO2air)
    first_date = pd.to_datetime(start_date)
    time = first_date + pd.to_timedelta(np.arange(N), 'D')
    ds['time'] = time

    return (ds)


if __name__ == "__main__":


    #for cycle in ["C", "CN", "CNP"]:
    for cycle in ["C"]:

        fname = "*_%s_out_casa_simulation.nc" % (cycle)
        fname = glob.glob(os.path.join("outputs", fname))[0]
        ds_casa = open_casa_and_add_time(fname, start_date="01/01/2002")

        fname = "*_%s_out_cable_simulation.nc" % (cycle)
        fname = glob.glob(os.path.join("outputs", fname))[0]
        ds_cable = open_cable_and_add_time(fname, start_date="01/01/2002")


        plot_carbon_fluxes(cycle, ds_casa, ds_cable)
        plot_nitrogen_fluxes(cycle, ds_casa)
        plot_phosphorus_fluxes(cycle, ds_casa)
        plot_cnp_states(cycle, ds_casa)
