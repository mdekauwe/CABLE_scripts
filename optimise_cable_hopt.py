#!/usr/bin/env python

"""
Run CABLE either for a single site, a subset, or all the flux sites pointed to
in the met directory

- Only intended for biophysics
- Set mpi = True if doing a number of flux sites

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (02.08.2018)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import shutil
import subprocess
import multiprocessing as mp
import numpy as np
import pandas as pd
from datetime import datetime
import xarray as xr
import matplotlib.pyplot as plt

from cable_utils import adjust_nml_file
from cable_utils import get_svn_info
from cable_utils import change_LAI
from cable_utils import add_attributes_to_output_file
from cable_utils import change_iveg, change_traits

from hyperopt import hp
from hyperopt import rand, tpe
from hyperopt import Trials
from hyperopt import fmin
from functools import partial

from summary_stats import rmse, bias, nash_sutcliffe, willmott_agreement_indx


def objective(params, b_plant, c_plant, vcmax, iveg, cable_src, zse, lai,
              met_subset, met_dir, log_dir, output_dir, restart_dir,
              namelist_dir, aux_dir, df_obs):

    #bch = params['bch']
    #hyds = params['hyds']
    Kmax = params['Kmax']
    bch = 11.056 #14.1244#9.4148
    hyds = 4.54E-06

    C = RunCable(met_dir=met_dir, log_dir=log_dir, output_dir=output_dir,
                 restart_dir=restart_dir, aux_dir=aux_dir,
                 namelist_dir=namelist_dir, met_subset=met_subset,
                 cable_src=cable_src, mpi=mpi, num_cores=num_cores,
                 co2_conc=410., fwsoil="profitmax", lai_dir="lai_data",
                 zse=zse, bch=bch, b_plant=b_plant, c_plant=c_plant,
                 vcmax=vcmax, Kmax=Kmax, iveg=iveg, hyds=hyds)
    C.main()

    site = met_subset[0].split(".")[0].split("_")[0]
    out_fname = os.path.join(output_dir, "%s_out.nc" % (site))
    shutil.move(out_fname, os.path.join(output_dir,
                "hydraulics_bch_%f_%f_%f.nc" % (bch, hyds, Kmax)))

    df = read_nc_file(os.path.join(output_dir,
                      "hydraulics_bch_%f_%f_%f.nc" % (bch, hyds, Kmax)))
    df = df.between_time('5:00', '21:00')

    # daily avgs
    df.Qle[df.Qle < 0] = np.nan # don't average across negative LE values
    method = {"Qle":"mean"}
    df = df.resample("D").agg(method)

    df_obs.Qle[df_obs.Qle < 0] = np.nan # don't average across negative LE values
    method = {"Qle":"mean"}
    df_obs = df_obs.resample("D").agg(method)

    # Just consider the extended growing season
    start_date = '2022-05-01'
    end_date = '2022-09-01'
    mask = (df.index >= start_date) & (df.index <= end_date)
    df = df.loc[mask]
    mask = (df_obs.index >= start_date) & (df_obs.index <= end_date)
    df_obs = df_obs.loc[mask]

    # Make sure dataframes are the same length, as the obs prob have gaps
    df = df[df.index.isin(df_obs.index)]

    m = df["Qle"].values
    o = df_obs["Qle"].values
    m = m[~np.isnan(m)]
    o = o[~np.isnan(o)]
    o_avg = np.mean(o)# Make sure dataframes are the same length, as the obs prob have gaps
    df = df[df.index.isin(df_obs.index)]

    #rmse_norm = rmse(m, o)
    rmse_norm = rmse(m, o) / o_avg # normlaised RMSE (RMSE_n)

    return rmse_norm

class RunCable(object):

    def __init__(self, met_dir=None, log_dir=None, output_dir=None,
                 restart_dir=None, aux_dir=None, namelist_dir=None,
                 nml_fname="cable.nml",
                 grid_fname="gridinfo_CSIRO_1x1.nc",
                 phen_fname="modis_phenology_csiro.txt",
                 cnpbiome_fname="pftlookup_csiro_v16_17tiles.csv",
                 elev_fname="GSWP3_gwmodel_parameters.nc", fwsoil="standard",
                 lai_dir=None, fixed_lai=None, co2_conc=400.0,
                 met_subset=[], cable_src=None, cable_exe="cable", mpi=True,
                 num_cores=None, verbose=True, zse=None, bch=None,
                 hyds=None, Kmax=None, b_plant=None, c_plant=None, vcmax=None,
                 iveg=None):

        self.met_dir = met_dir
        self.log_dir = log_dir
        self.output_dir = output_dir
        self.restart_dir = restart_dir
        self.aux_dir = aux_dir
        self.namelist_dir = namelist_dir
        self.nml_fname = nml_fname
        self.biogeophys_dir = os.path.join(self.aux_dir, "core/biogeophys")
        self.grid_dir = os.path.join(self.aux_dir, "offline")
        self.biogeochem_dir = os.path.join(self.aux_dir, "core/biogeochem/")
        self.grid_fname = os.path.join(self.grid_dir, grid_fname)
        self.phen_fname = os.path.join(self.biogeochem_dir, phen_fname)
        self.cnpbiome_fname = os.path.join(self.biogeochem_dir, cnpbiome_fname)
        self.elev_fname = elev_fname
        self.co2_conc = co2_conc
        self.met_subset = met_subset
        self.cable_src = cable_src
        self.cable_exe = os.path.join(cable_src, "offline/%s" % (cable_exe))
        self.verbose = verbose
        self.mpi = mpi
        self.num_cores = num_cores
        self.lai_dir = lai_dir
        self.fixed_lai = fixed_lai
        self.fwsoil = fwsoil
        self.zse = zse
        self.bch = bch
        self.hyds = hyds
        self.Kmax = Kmax
        self.b_plant = b_plant
        self.c_plant = c_plant
        self.vcmax = vcmax
        self.iveg = iveg

    def main(self):

        (met_files, url, rev) = self.initialise_stuff()

        # Setup multi-processor jobs
        if self.mpi:
            if self.num_cores is None: # use them all!
                self.num_cores = mp.cpu_count()
            chunk_size = int(np.ceil(len(met_files) / float(self.num_cores)))
            pool = mp.Pool(processes=self.num_cores)
            processes = []

            for i in range(self.num_cores):
                start = chunk_size * i
                end = chunk_size * (i + 1)
                if end > len(met_files):
                    end = len(met_files)

                # setup a list of processes that we want to run
                p = mp.Process(target=self.worker,
                               args=(met_files[start:end], url, rev, ))
                processes.append(p)

            # Run processes
            for p in processes:
                p.start()
        else:
            self.worker(met_files, url, rev)

    def worker(self, met_files, url, rev):

        for fname in met_files:
            site = os.path.basename(fname).split(".")[0].split("_")[0]

            base_nml_fn = os.path.join(self.grid_dir, "%s" % (self.nml_fname))
            nml_fname = "cable_%s.nml" % (site)
            shutil.copy(base_nml_fn, nml_fname)

            (out_fname, out_log_fname) = self.clean_up_old_files(site)

            # Add LAI to met file?
            if self.fixed_lai is not None or self.lai_dir is not None:
                fname = change_LAI(fname, site, fixed=self.fixed_lai,
                                   lai_dir=self.lai_dir)

            if self.fwsoil == "hydraulics":
                #fname = change_iveg(fname, site, 2)
                fname_new = change_iveg(fname, site, 19)
            elif self.fwsoil == "standard":
                fname_new = change_iveg(fname, site, 4)
            elif self.fwsoil == "profitmax":

                fname_new = change_traits(fname, id, site, self.Kmax,
                                          self.b_plant, self.c_plant,
                                          self.vcmax, self.iveg, self.zse,
                                          self.bch, self.hyds)

            replace_dict = {
                            "filename%met": "'%s'" % (fname_new),
                            "filename%out": "'%s'" % (out_fname),
                            "filename%log": "'%s'" % (out_log_fname),
                            "filename%restart_out": "' '",
                            "filename%type": "'%s'" % (self.grid_fname),
                            "output%restart": ".FALSE.",
                            "fixedCO2": "%.2f" % (self.co2_conc),
                            "casafile%phen": "'%s'" % (self.phen_fname),
                            "casafile%cnpbiome": "'%s'" % (self.cnpbiome_fname),
                            "cable_user%GS_SWITCH": "'medlyn'",
                            "cable_user%GW_MODEL": ".FALSE.",
                            "cable_user%or_evap": ".FALSE.",
                            "redistrb": ".FALSE.",
                            "spinup": ".FALSE.",
                            "cable_user%litter": ".TRUE.",
                            "cable_user%FWSOIL_SWITCH": "'%s'" % (self.fwsoil),
            }
            adjust_nml_file(nml_fname, replace_dict)

            self.run_me(nml_fname)

            add_attributes_to_output_file(nml_fname, out_fname, url, rev)
            shutil.move(nml_fname, os.path.join(self.namelist_dir, nml_fname))

            if self.fixed_lai is not None or self.lai_dir is not None:
                os.remove("%s_tmp.nc" % (site))

            os.remove("new_sumbal")

    def initialise_stuff(self):

        if not os.path.exists(self.restart_dir):
            os.makedirs(self.restart_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        if not os.path.exists(self.namelist_dir):
            os.makedirs(self.namelist_dir)

        # Run all the met files in the directory
        if len(self.met_subset) == 0:
            met_files = glob.glob(os.path.join(self.met_dir, "*.nc"))
        else:
            met_files = [os.path.join(self.met_dir, i) for i in self.met_subset]

        cwd = os.getcwd()
        (url, rev) = get_svn_info(cwd, self.cable_src)

        # delete local executable, copy a local copy and use that
        local_exe = "cable"
        if os.path.isfile(local_exe):
            os.remove(local_exe)
        shutil.copy(self.cable_exe, local_exe)
        self.cable_exe = local_exe

        return (met_files, url, rev)

    def clean_up_old_files(self, site):
        out_fname = os.path.join(self.output_dir, "%s_out.nc" % (site))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        out_log_fname = os.path.join(self.log_dir, "%s_log.txt" % (site))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        return (out_fname, out_log_fname)

    def run_me(self, nml_fname):
        # run the model
        if self.verbose:
            cmd = './%s %s' % (self.cable_exe, nml_fname)
            error = subprocess.call(cmd, shell=True)
            if error == 1:
                print("Job failed to submit")
                raise
        else:
            # No outputs to the screen: stout and stderr to dev/null
            cmd = './%s %s > /dev/null 2>&1' % (self.cable_exe, nml_fname)
            error = subprocess.call(cmd, shell=True)
            if error == 1:
                print("Job failed to submit")


def read_nc_file(fname):

    vars_to_keep = ['Qle']

    ds = xr.open_dataset(fname, decode_times=False)
    time_jump = int(ds.time[1].values) - int(ds.time[0].values)

    if time_jump == 3600:
        freq = "H"
    elif time_jump == 1800:
        freq = "30MIN"
    else:
        raise("Time problem")

    units, reference_date = ds.time.attrs['units'].split('since')

    try:
        ds = ds[vars_to_keep].squeeze(dim=["x","y","patch"], drop=True)
    except:
        ds = ds[vars_to_keep].squeeze(dim=["x","y"], drop=True)


    df = pd.DataFrame(ds['Qle'].values[:], columns=['Qle'])
    #df = ds.to_dataframe()

    start = reference_date.strip()[0:19].replace("-","/")
    df['dates'] = pd.date_range(start=start, periods=len(df), freq=freq)
    df = df.set_index('dates')

    return df


if __name__ == "__main__":


    flux_dir = "/Users/xj21307/research/Alice_Holt/data"
    flux_fname = os.path.join(flux_dir, "alice_holt_flux_2022.nc")
    df_obs = read_nc_file(flux_fname)
    df_obs = df_obs.between_time('5:00', '21:00')


    #------------- Change stuff ------------- #
    met_dir = "/Users/xj21307/research/Alice_Holt/data"
    #met_dir = os.getcwd()
    log_dir = "logs"
    output_dir = "outputs_optimise_hopt"
    restart_dir = "restart_files"
    namelist_dir = "namelists"

    #cable_src = "../../src/trunk_DESICA_PFTs/trunk_DESICA_PFTs/"

    aux_dir = "../../src/CABLE-AUX/"

    mpi = False
    num_cores = 4 # set to a number, if None it will use all cores...!
    # if empty...run all the files in the met_dir
    met_subset = ['UK-Ham_2002-2003_Met.nc']

    # ------------------------------------------- #

    site = met_subset[0].split(".")[0].split(".")[0].split("_")[0]

    fname = "lai_data/UK-Ham_500m_composite.csv"
    date_parser=lambda x: datetime.strptime(x, '%d/%m/%Y')
    df = pd.read_csv(fname, index_col='date', date_parser=date_parser)
    df = df[df.index.year == 2022]
    df.loc["2022-12-28 00:00:00"] = 0.0 # pad missing values
    df.loc["2022-12-29 00:00:00"] = 0.0
    df.loc["2022-12-30 00:00:00"] = 0.0
    df.loc["2022-12-31 00:00:00"] = 0.0

    lai = df.LAI
    lai.to_csv("lai_data/%s_lai.csv" % (site))
    # Evaluation of LandscapeDNDC Model Predictions of CO2 and N2O Fluxes from an Oak Forest in SE England, Table 2
    zse = np.array([0.02, 0.08, 0.08, 0.2, 0.38, 0.26]) # Experiment with 6 layers = 1 m

    #zse = np.array([0.00478261,0.0126087, 0.02913043, 0.04108696, 0.28804348, 0.62434783]) # Experiment with 6 layers = 1 m
    #zse = np.array([0.011, 0.029, 0.067, 0.0945, 0.6625, 1.436]) # 6 layers = 2.3m
    #zse = np.array([0.01466667,0.03866667,0.08933333,0.126,0.88333333,1.91466667]) # 6 layers = 3m

    cable_src = "../../src/manon_cable/CableHydraulics//"

    b_plant = 4.969847639573251
    c_plant = 7.740170420442604
    vcmax = 50.

    iveg = 2


    # Create the domain space
    #space = {'bch': hp.uniform('bch', 2, 15),
    #         'hyds': hp.uniform('hyds', 5.0E-7,8.5E-3),
    #         'Kmax': hp.uniform('Kmax', 0.1, 2.5)}
    space = {'Kmax': hp.uniform('Kmax', 0.1, 2.5)}

    # Create two trials objects
    tpe_trials = Trials()

    # Wrapper to pass args
    fmin_objective = partial(objective, b_plant=b_plant,
                             c_plant=c_plant, vcmax=vcmax, iveg=iveg,
                             cable_src=cable_src, zse=zse, lai=lai,
                             met_subset=met_subset, met_dir=met_dir,
                             log_dir=log_dir, output_dir=output_dir,
                             restart_dir=restart_dir, namelist_dir=namelist_dir,
                             aux_dir=aux_dir, df_obs=df_obs)


    # Run 2000 evals with the tpe algorithm
    tpe_best = fmin(fn=fmin_objective, space=space, algo=tpe.suggest,
                    trials=tpe_trials, max_evals=2000)

    """
    results = pd.DataFrame({'loss': [x['loss'] for x in tpe_trials.results],
                             'iteration': tpe_trials.idxs_vals[0]['bch'],
                             'bch': tpe_trials.idxs_vals[1]['bch']})

    # Print out information about number of trials
    print('\nNumber of iteration: {}'.format(tpe_trials.best_trial['misc']['idxs']['bch'][0]))

    # Print out information about value of x
    print('\nMiniumum value of bch: {:.4f}'.format(tpe_best['bch']))

    plt.figure(figsize = (8, 6))
    plt.hist(results['bch'], bins = 50, edgecolor = 'k')
    plt.xlabel('bch')
    plt.ylabel('Count')
    #plt.show()
    plt.savefig("/Users/xj21307/Desktop/bch.png")

    results = pd.DataFrame({'loss': [x['loss'] for x in tpe_trials.results],
                             'iteration': tpe_trials.idxs_vals[0]['hyds'],
                             'hyds': tpe_trials.idxs_vals[1]['hyds']})

    # Print out information about number of trials
    print('\nNumber of iteration: {}'.format(tpe_trials.best_trial['misc']['idxs']['hyds'][0]))

    # Print out information about value of x
    print('\nMiniumum value of hyds: {:.10f}'.format(tpe_best['hyds']))

    plt.figure(figsize = (8, 6))
    plt.hist(results['hyds'], bins = 50, edgecolor = 'k')
    plt.xlabel('hyds')
    plt.ylabel('Count')
    #plt.show()
    plt.savefig("/Users/xj21307/Desktop/hyds.png")
    """
    results = pd.DataFrame({'loss': [x['loss'] for x in tpe_trials.results],
                             'iteration': tpe_trials.idxs_vals[0]['Kmax'],
                             'Kmax': tpe_trials.idxs_vals[1]['Kmax']})

    # Print out information about number of trials
    print('\nNumber of iteration: {}'.format(tpe_trials.best_trial['misc']['idxs']['Kmax'][0]))

    # Print out information about value of x
    print('\nMiniumum value of Kmax: {:.4f}'.format(tpe_best['Kmax']))

    plt.figure(figsize = (8, 6))
    plt.hist(results['Kmax'], bins = 50, edgecolor = 'k')
    plt.xlabel('Kmax')
    plt.ylabel('Count')
    #plt.show()
    plt.savefig("/Users/xj21307/Desktop/Kmax.png")
