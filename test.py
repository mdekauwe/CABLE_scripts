#!/usr/bin/env python

"""
Wrapper script to do MCMC param estimation on CABLE. This is using PYMC3,
everything needs to be cast to float 64

THEANO_FLAGS='floatX=float64' ./run_pymc_wrapper.py

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
import theano.tensor as tt

from run_cable_site import RunCable


class LogLike(tt.Op):

    """
    Specify what type of object will be passed and returned to the Op when it is
    called. In our case we will be passing it a vector of values (the parameters
    that define our model) and returning a single "scalar" value (the
    log-likelihood)
    """
    itypes = [tt.dvector] # expects a vector of parameter values when called
    #otypes = [tt.dscalar] # outputs a single scalar value (the log likelihood)
    otypes = [tt.dvector] # outputs a vector for the log likelihood)


    def __init__(self, loglike, obs, sigma):
        """
        Initialise the Op with various things that our log-likelihood function
        requires. Below are the things that are needed in this particular
        example.

        Parameters
        ----------
        loglike:
            The log-likelihood (or whatever) function we've defined
        data:
            The "observed" data that our log-likelihood function takes in
        x:
            The dependent variable (aka 'x') that our model requires
        sigma:
            The noise standard deviation that our function requires.
        """

        # add inputs as class attributes
        self.likelihood = loglike
        self.obs = obs
        self.sigma = sigma

    def perform(self, node, inputs, outputs):
        # the method that is used when calling the Op
        theta, = inputs  # this will contain my variables

        # call the log-likelihood function
        logl = self.likelihood(theta, self.obs, self.sigma)

        outputs[0][0] = np.array(logl) # output the log-likelihood

def my_likelihood(theta, obs, sigma):
    """
    A Gaussian log-likelihood function for a model with parameters given
    in theta
    """

    model = run_and_unpack_cable(theta)

    return -(0.5/sigma**2)*np.sum((obs - model)**2)
    #return np.sum(-(0.5/sigma**2)*np.sum((obs - model)**2))

def my_loglikelihood(theta, obs, sigma):

    model = run_and_unpack_cable(theta)

    return -np.sum(0.5 * \
            (np.log(2. * np.pi * sigma ** 2.) + ((obs - model) / sigma) ** 2))

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

#@theano.compile.ops.as_op(itypes=[tt.dscalar, tt.dscalar], otypes=[tt.dvector])
#@theano.as_op(itypes=[tt.dscalar, tt.dscalar], otypes=[tt.dvector])
def run_and_unpack_cable(theta):

    g1, vcmax = theta  # unpack line gradient and y-intercept

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
    mcmc_tag = id
    C.main(param_names=param_names, param_values=params,
           out_fname=out_fname, out_log_fname=out_log_fname,
           mcmc_tag=mcmc_tag)

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
    mod = mod.astype(np.float64)
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
obs = obs.astype(np.float64)
#plt.plot(obs)
#plt.show()

#ds = xr.open_dataset(fn)
#obs = ds.Qle.values[:,0,0]
#uncert = np.sqrt(np.abs(obs))
uncert = 0.1 * np.abs(obs) # not using, letting pymc fit this below...
uncert = uncert.astype(np.float64)

# create our Op
logl = LogLike(my_likelihood, obs, uncert)
#logl = LogLike(my_loglikelihood, obs, uncert)
with pm.Model() as model:
    g1 = pm.Uniform('g1', lower=0.0, upper=8.0)
    vcmax = pm.Uniform('vcmax', lower=10.0, upper=120.0)
    #sigma = pm.Uniform('sigma', lower=0.0, upper=20.0) # fit error?

    theta = tt.as_tensor_variable([g1, vcmax])

    # use a DensityDist (use a lamdba function to "call" the Op)
    pm.DensityDist('likelihood', lambda v: logl(v), observed={'v': theta})

    #trace = pm.sample(10, chains=1, step=pm.Metropolis())
    trace = pm.sample(10, chains=1, step=pm.NUTS())
    #trace = pm.sample(10, chains=1, step=pm.Slice())
    #trace = pm.sample(10, chains=1)
_ = pm.traceplot(trace)
pm.summary(trace)
