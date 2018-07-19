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

def plot_plant(tag, cycle, zero, ccp1, ccp2, ccp3, ccp4, transient, simulation):

    tol = 1E-02 #5E-03
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


    zero, ccp1, ccp2, ccp3, ccp4, transient, simulation

    cf = np.hstack((zero.cplant[:,0,0].values,
                    ccp1.cplant[:,0,0].values,
                    ccp2.cplant[:,0,0].values,
                    ccp3.cplant[:,0,0].values,
                    ccp4.cplant[:,0,0].values))
    cw = np.hstack((zero.cplant[:,1,0].values,
                    ccp1.cplant[:,1,0].values,
                    ccp2.cplant[:,1,0].values,
                    ccp3.cplant[:,1,0].values,
                    ccp4.cplant[:,1,0].values))
    cr = np.hstack((zero.cplant[:,2,0].values,
                    ccp1.cplant[:,2,0].values,
                    ccp2.cplant[:,2,0].values,
                    ccp3.cplant[:,2,0].values,
                    ccp4.cplant[:,2,0].values))

    ax1.set_title("Cf")
    ax1.plot(cf, label="Cf")

    ax2.set_title("Cr")
    ax2.plot(cw)

    ax3.set_title("Cr")
    ax3.plot(cr)

    ax4.set_title("Cplant")
    ax4.plot(cf+cw+cr)

    nf = np.hstack((zero.nplant[:,0,0].values,
                    ccp1.nplant[:,0,0].values,
                    ccp2.nplant[:,0,0].values,
                    ccp3.nplant[:,0,0].values,
                    ccp4.nplant[:,0,0].values))
    nw = np.hstack((zero.nplant[:,1,0].values,
                    ccp1.nplant[:,1,0].values,
                    ccp2.nplant[:,1,0].values,
                    ccp3.nplant[:,1,0].values,
                    ccp4.nplant[:,1,0].values))
    nr = np.hstack((zero.nplant[:,2,0].values,
                    ccp1.nplant[:,2,0].values,
                    ccp2.nplant[:,2,0].values,
                    ccp3.nplant[:,2,0].values,
                    ccp4.nplant[:,2,0].values))

    ax5.set_title("Nf")
    ax5.plot(nf, label="nf")

    ax6.set_title("Nr")
    ax6.plot(nw)

    ax7.set_title("Nr")
    ax7.plot(nr)

    ax8.set_title("Nplant")
    ax8.plot(nf+nw+nr)

    pf = np.hstack((zero.pplant[:,0,0].values,
                    ccp1.pplant[:,0,0].values,
                    ccp2.pplant[:,0,0].values,
                    ccp3.pplant[:,0,0].values,
                    ccp4.pplant[:,0,0].values))
    pw = np.hstack((zero.pplant[:,1,0].values,
                    ccp1.pplant[:,1,0].values,
                    ccp2.pplant[:,1,0].values,
                    ccp3.pplant[:,1,0].values,
                    ccp4.pplant[:,1,0].values))
    pr = np.hstack((zero.pplant[:,2,0].values,
                    ccp1.pplant[:,2,0].values,
                    ccp2.pplant[:,2,0].values,
                    ccp3.pplant[:,2,0].values,
                    ccp4.pplant[:,2,0].values))

    ax9.set_title("Pf")
    ax9.plot(pf, label="pf")

    ax10.set_title("Pr")
    ax10.plot(pw)

    ax11.set_title("Pr")
    ax11.plot(pr)

    ax12.set_title("Pplant")
    ax12.plot(pf+pw+pr)

    plot_fname = "%s_spinup_plant.pdf" % (tag)
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)

    cplant = cf + cw + cr
    cplant = cplant[-50:]
    #delta = np.diff(cplant)

    for i,val in enumerate(cplant[1:]):
        if np.fabs(cplant[i] - val) < tol:
            print("C plant (%s): steady-state" % (cycle))


def plot_soil(tag, cycle, zero, ccp1, ccp2, ccp3, ccp4, transient, simulation):

    tol = 1E-02 #5E-03
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


    zero, ccp1, ccp2, ccp3, ccp4, transient, simulation

    cfast = np.hstack((zero.csoil[:,0,0].values,
                      ccp1.csoil[:,0,0].values,
                      ccp2.csoil[:,0,0].values,
                      ccp3.csoil[:,0,0].values,
                      ccp4.csoil[:,0,0].values))
    cslow = np.hstack((zero.csoil[:,1,0].values,
                      ccp1.csoil[:,1,0].values,
                      ccp2.csoil[:,1,0].values,
                      ccp3.csoil[:,1,0].values,
                      ccp4.csoil[:,1,0].values))
    cpassive = np.hstack((zero.csoil[:,2,0].values,
                         ccp1.csoil[:,2,0].values,
                         ccp2.csoil[:,2,0].values,
                         ccp3.csoil[:,2,0].values,
                         ccp4.csoil[:,2,0].values))

    ax1.set_title("C fast")
    ax1.plot(cfast, label="cfast")

    ax2.set_title("C slow")
    ax2.plot(cslow)

    ax3.set_title("C passive")
    ax3.plot(cpassive)

    ax4.set_title("Cplant")
    ax4.plot(cfast+cslow+cpassive)

    nfast = np.hstack((zero.nsoil[:,0,0].values,
                       ccp1.nsoil[:,0,0].values,
                       ccp2.nsoil[:,0,0].values,
                       ccp3.nsoil[:,0,0].values,
                       ccp4.nsoil[:,0,0].values))
    nslow = np.hstack((zero.nsoil[:,1,0].values,
                       ccp1.nsoil[:,1,0].values,
                       ccp2.nsoil[:,1,0].values,
                       ccp3.nsoil[:,1,0].values,
                       ccp4.nsoil[:,1,0].values))
    npassive = np.hstack((zero.nsoil[:,2,0].values,
                          ccp1.nsoil[:,2,0].values,
                          ccp2.nsoil[:,2,0].values,
                          ccp3.nsoil[:,2,0].values,
                          ccp4.nsoil[:,2,0].values))

    ax5.set_title("N fast")
    ax5.plot(nfast)

    ax6.set_title("N slow")
    ax6.plot(nslow)

    ax7.set_title("N passive")
    ax7.plot(npassive)

    ax8.set_title("N plant")
    ax8.plot(nfast+nslow+npassive)

    pfast = np.hstack((zero.psoil[:,0,0].values,
                       ccp1.psoil[:,0,0].values,
                       ccp2.psoil[:,0,0].values,
                       ccp3.psoil[:,0,0].values,
                       ccp4.psoil[:,0,0].values))
    pslow = np.hstack((zero.psoil[:,1,0].values,
                       ccp1.psoil[:,1,0].values,
                       ccp2.psoil[:,1,0].values,
                       ccp3.psoil[:,1,0].values,
                       ccp4.psoil[:,1,0].values))
    ppassive = np.hstack((zero.psoil[:,2,0].values,
                          ccp1.psoil[:,2,0].values,
                          ccp2.psoil[:,2,0].values,
                          ccp3.psoil[:,2,0].values,
                          ccp4.psoil[:,2,0].values))

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

    csoil = cfast+cslow+cpassive
    csoil = csoil[-50:]

    for i,val in enumerate(csoil[1:]):
        if np.fabs(csoil[i] - val) < tol:
            print("C soil (%s): steady-state" % (cycle))

def open_file(fname):
    return xr.open_dataset(fname)



if __name__ == "__main__":

    biogeoC = False
    biogeoCN = False
    biogeoCNP = True

    for i in ["C", "CN", "CNP"]:
        if i == "C":
            biogeoC = False
            biogeoCN = False
            biogeoCNP = True
        elif i == "CN":
            biogeoC = False
            biogeoCN = True
            biogeoCNP = False
        elif i == "CNP":
            biogeoC = False
            biogeoCN = False
            biogeoCNP = True

        if biogeoC:
            experiment_id = "Cumberland_C"
            tag = "C"
        elif biogeoCN:
            experiment_id = "Cumberland_CN"
            tag = "CN"
        elif biogeoCNP:
            experiment_id = "Cumberland_CNP"
            tag = "CNP"

        fname = "outputs/%s_out_CASA_zero.nc" % (experiment_id)
        zero = open_file(fname)
        fname = "outputs/%s_out_CASA_ccp1.nc" % (experiment_id)
        ccp1 = open_file(fname)
        fname = "outputs/%s_out_CASA_ccp2.nc" % (experiment_id)
        ccp2 = open_file(fname)
        fname = "outputs/%s_out_CASA_ccp3.nc" % (experiment_id)
        ccp3 = open_file(fname)
        fname = "outputs/%s_out_CASA_ccp4.nc" % (experiment_id)
        ccp4 = open_file(fname)
        fname = "outputs/%s_out_casa_transient.nc" % (experiment_id)
        transient = open_file(fname)
        fname = "outputs/%s_out_casa.nc" % (experiment_id)
        simulation = open_file(fname)

        plot_plant(tag, zero, ccp1, ccp2, ccp3, ccp4, transient, simulation)
        plot_soil(tag, zero, ccp1, ccp2, ccp3, ccp4, transient, simulation)
