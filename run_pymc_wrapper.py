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
import pymc3 as pm
import matplotlib.pyplot as plt
import netCDF4 as nc
import uuid
import pandas as pd
import theano
import theano.tensor as t

from run_cable_site import RunCable

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

@theano.compile.ops.as_op(itypes=[t.dscalar, t.dscalar], otypes=[t.dvector])
def run_and_unpack_cable(g1, vcmax):
    params = np.array([g1, vcmax])
    param_names = ["g1", "vcmax"]

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

    C.main(param_names=param_names, param_values=params,
           out_fname=out_fname, out_log_fname=out_log_fname)

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

    return np.array(mod)



# set up Observations...
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
uncert = 0.1 * np.abs(obs) # not using, letting pymc fit this below...


niter = 1000
with pm.Model() as model:
    g1 = pm.Uniform('g1', lower=0, upper=8)
    vcmax = pm.Uniform('vcmax', lower=10., upper=120)
    sigma = pm.Uniform('sigma', lower=0, upper=20) # fit error?

    # define likelihood, i.e. call CABLE...
    #mod = pm.Deterministic('mod', run_and_unpack_cable(g1, vcmax))
    mod = run_and_unpack_cable(g1, vcmax)
    y_obs = pm.Normal('Y_obs', mu=mod, sd=sigma, observed=obs)

    # inference
    #step = pm.NUTS() # Hamiltonian MCMC with No U-Turn Sampler
    #step = pm.Metropolis()
    #trace = pm.sample(niter, step=step, progressbar=True)
    #pm.traceplot(trace)

    #trace = pm.sample(10, chains=1, step=pm.Metropolis())
    trace = pm.sample(10, chains=1)
_ = pm.traceplot(trace)
pm.summary(trace)
