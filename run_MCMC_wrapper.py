#!/usr/bin/env python

"""
Wrapper script to do MCMC param estimation on CABLE. For info on MCMC see lib
documentation -> https://mc3.readthedocs.io/en/latest/mcmc_tutorial.html

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (27.12.2019)"
__email__ = "mdekauwe@gmail.com"


import os
import sys
import glob
import shutil
import subprocess
import multiprocessing as mp
import numpy as np
import mc3
import matplotlib.pyplot as plt
import netCDF4 as nc
import uuid
#import tempfile
import pandas as pd

from cable_utils import adjust_nml_file
from cable_utils import get_svn_info
from cable_utils import change_LAI, change_params
from cable_utils import add_attributes_to_output_file

from run_cable_site import RunCable

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def run_and_unpack_cable(params, param_names):

    #------------- Change stuff ------------- #
    met_dir = "../../met_data/plumber_met"
    log_dir = "logs"
    output_dir = "outputs"
    restart_dir = "restart_files"
    namelist_dir = "namelists"
    aux_dir = "../../src/CABLE-AUX/"
    cable_src = "../../src/trunk/trunk"
    mpi = False
    num_cores = 1 # set to a number, if None it will use all cores...!
    met_subset = ['TumbaFluxnet.1.4_met.nc']

    # MCMC
    adjust_params = True

    print("\n")
    print("===========")
    print(params)
    print("===========")
    print("\n")

    # ------------------------------------------- #
    C = RunCable(met_dir=met_dir, log_dir=log_dir, output_dir=output_dir,
                 restart_dir=restart_dir, aux_dir=aux_dir,
                 namelist_dir=namelist_dir, met_subset=met_subset,
                 cable_src=cable_src, mpi=mpi, num_cores=num_cores,
                 adjust_params=adjust_params)



    #temp = tempfile.NamedTemporaryFile(suffix=".nc")
    #out_fname = temp.name
    ##temp = tempfile.NamedTemporaryFile(suffix=".txt")
    #out_log_fname = temp.name

    id = str(uuid.uuid4())
    out_fname = os.path.join(output_dir, "%s.nc" % (id))
    out_log_fname = os.path.join(log_dir, "%s.txt" % (id))

    C.main(param_names=param_names, param_values=params, out_fname=out_fname,
           out_log_fname=out_log_fname)

    f = nc.Dataset(out_fname)
    time = nc.num2date(f.variables['time'][:],
                       f.variables['time'].units)
    df = pd.DataFrame(f.variables['Qle'][:,0,0], columns=['Qle'])
    df['date'] = time
    df = df.set_index('date')
    df = df.between_time('5:00', '20:00')
    #df = df[df.index.year == 2002]
    df = df.resample("D").agg("mean")
    mod = df.Qle.values
    #plt.plot(mod)
    #plt.show()
    #sys.exit()
    #ds = xr.open_dataset(out_fname, decode_times=False)
    #mod = ds.Qle.values[:,0,0]
    #print("model:", np.mean(mod), np.min(mod), np.max(mod))
    if os.path.exists(out_fname):
        os.remove(out_fname)
    if os.path.exists(out_log_fname):
        os.remove(out_log_fname)

    return mod

obs_dir = "../../flux_files/plumber"
fn = os.path.join(obs_dir, 'TumbaFluxnet.1.4_flux.nc')
f = nc.Dataset(fn)
time = nc.num2date(f.variables['time'][:],
                   f.variables['time'].units)
df = pd.DataFrame(f.variables['Qle'][:,0,0], columns=['Qle'])
df['date'] = time
df = df.set_index('date')
df = df.between_time('5:00', '20:00')
#df = df[df.index.year == 2002]
df = df.resample("D").agg("mean")
obs = df.Qle.values
#plt.plot(obs)
#plt.show()

#ds = xr.open_dataset(fn)
#obs = ds.Qle.values[:,0,0]
#uncert = np.sqrt(np.abs(obs))
uncert = 0.1 * np.abs(obs)
#plt.plot(uncert)
#plt.show()
#sys.exit()
#print("obs:", np.mean(obs), np.min(obs), np.max(obs))

# Define the modeling function as a callable, comparing Qle.
func = run_and_unpack_cable

# Array of initial-guess values of fitting parameters:
#param_names = ["g1", "vcmax"]
#params = np.array([2.0, 50.0])

param_names = ["vcmax"]
params = np.array([50.0])


# Lower and upper boundaries for the MCMC exploration:
#pmin = np.array([0.0, 10.0])   # kPa^0.5, umol/m2/s
#pmax = np.array([8.0, 120.0])  # kPa^0.5, umol/m2/s
#pstep = np.array([1.0, 1.0])

pmin = np.array([10.0])   # kPa^0.5, umol/m2/s
pmax = np.array([120.0])  # kPa^0.5, umol/m2/s
pstep = np.array([1.0])

# Parameter prior probability distributions:
# uniform priors
#prior = np.array([0.0, 0.0])
#priorlow = np.array([0.0, 0.0])
#priorup = np.array([0.0, 0.0])

prior = np.array([0.0])
priorlow = np.array([0.0])
priorup = np.array([0.0])

# Parameter names:
#pnames = ['g1', 'vcmax']
#texnames = [r'$g_{1}$', r'$V_{cmax}$']

pnames = ['vcmax']
texnames = [r'$V_{cmax}$']

# List of additional arguments of func (if necessary):
indparams = [param_names]

# Sampler algorithm, choose from: 'snooker', 'demc' or 'mrw'.
sampler = 'snooker'

# MCMC setup:
nsamples = 5000
burnin = 500
nchains = 6 # set to a multiple of ncups
ncpu = 3
thinning = 1

# MCMC initial draw, choose from: 'normal' or 'uniform'
kickoff = 'normal'

# DEMC snooker pre-MCMC sample size:
hsize = 10

# Optimization before MCMC, choose from: 'lm' or 'trf':
# Levenberg-Marquardt = lm
leastsq = 'lm'
chisqscale = False

# MCMC Convergence:
grtest  = True
grbreak = 1.01
grnmin  = 0.5

# Carter & Winn (2009) Wavelet-likelihood method:
wlike = False

fgamma = 1.0  # Scale factor for DEMC's gamma jump.
fepsilon = 0.0  # Jump scale factor for DEMC's "e" distribution

# Logging:
log = 'MCMC_CABLE.log'
#log = None

# File outputs:
savefile = 'MCMC_CABLE.npz'
plots = True
rms = True

# Run the MCMC:
mc3_output = mc3.sample(data=obs, uncert=uncert, func=func, params=params,
                        indparams=indparams, pmin=pmin, pmax=pmax,
                        pstep=pstep, pnames=pnames, texnames=texnames,
                        prior=prior, priorlow=priorlow, priorup=priorup,
                        sampler=sampler, nsamples=nsamples,  nchains=nchains,
                        ncpu=ncpu, burnin=burnin, thinning=thinning,
                        leastsq=leastsq, chisqscale=chisqscale,
                        grtest=grtest, grbreak=grbreak, grnmin=grnmin,
                        hsize=hsize, kickoff=kickoff,
                        wlike=wlike, log=log,
                        plots=plots, savefile=savefile, rms=rms)
