#!/usr/bin/env python

"""
Plot a heap of CASA outputs

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

def main(fname, plot_fname=None):

    ds = xr.open_dataset(fname)
    
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

    ax1 = fig.add_subplot(4,4,1)
    ax2 = fig.add_subplot(4,4,2)
    ax3 = fig.add_subplot(4,4,3)
    ax4 = fig.add_subplot(4,4,4)

    ax5 = fig.add_subplot(4,4,5)
    ax6 = fig.add_subplot(4,4,6)
    ax7 = fig.add_subplot(4,4,7)
    ax8 = fig.add_subplot(4,4,8)

    ax9 = fig.add_subplot(4,4,9)
    ax10 = fig.add_subplot(4,4,10)
    ax11 = fig.add_subplot(4,4,11)
    ax12 = fig.add_subplot(4,4,12)

    ax13 = fig.add_subplot(4,4,13)
    ax14 = fig.add_subplot(4,4,14)
    ax15 = fig.add_subplot(4,4,15)
    ax16 = fig.add_subplot(4,4,16)

    ax1.set_title("Cplant")
    ax1.plot(ds.cplant[:,0,0], label="0")
    ax1.plot(ds.cplant[:,1,0], label="1")
    ax1.plot(ds.cplant[:,2,0], label="2")
    ax1.legend(numpoints=1, loc="best")

    ax2.set_title("Nplant")
    ax2.plot(ds.nplant[:,0,0])
    ax2.plot(ds.nplant[:,1,0])
    ax2.plot(ds.nplant[:,2,0])

    ax3.set_title("Pplant")
    ax3.plot(ds.pplant[:,0,0])
    ax3.plot(ds.pplant[:,1,0])
    ax3.plot(ds.pplant[:,2,0])

    ax4.set_title("Cgpp")
    ax4.plot(ds.Cgpp[:,0])

    ax5.set_title("Csoil")
    ax5.plot(ds.csoil[:,0,0])
    ax5.plot(ds.csoil[:,1,0])
    ax5.plot(ds.csoil[:,2,0])

    ax6.set_title("Nsoil")
    ax6.plot(ds.nsoil[:,0,0])
    ax6.plot(ds.nsoil[:,1,0])
    ax6.plot(ds.nsoil[:,2,0])

    ax7.set_title("Psoil")
    ax7.plot(ds.nsoil[:,0,0])
    ax7.plot(ds.nsoil[:,1,0])
    ax7.plot(ds.nsoil[:,2,0])

    ax8.set_title("CUE")
    ax8.plot(ds.Cnpp[:,0]/ds.Cgpp[:,0])

    ax9.set_title("glai")
    ax9.plot(ds.glai[:,0])

    ax10.set_title("Calloc")
    ax10.plot(ds.fracCalloc[:,0])
    ax10.plot(ds.fracCalloc[:,1])
    ax10.plot(ds.fracCalloc[:,2])

    ax11.set_title("vcmax")
    ax11.plot(ds.vcmax[:,0]*1E6)

    ax12.set_title("Nupland")
    ax12.plot(ds.Nupland[:,0])

    ax13.set_title("Nlittermin")
    ax13.plot(ds.Nlittermin[:,0])

    ax14.set_title("clabile")
    ax14.plot(ds.clabile[:,0])

    ax15.set_title("N mineralisation")
    ax15.plot(ds.Nsmin[:,0])

    ax16.set_title("N immobilisation")
    ax16.plot(ds.Nsimm[:,0])


    plt.show()

    if plot_fname is None:
        plt.show()
    else:
        plot_dir = "plots"
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                    pad_inches=0.1)



if __name__ == "__main__":

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-f", "--fname", dest="fname",
                      action="store", help="filename",
                      type="string")
    parser.add_option("-p", "--plot_fname", dest="plot_fname", action="store",
                      help="Benchmark plot filename", type="string")
    (options, args) = parser.parse_args()

    main(options.fname, options.plot_fname)
