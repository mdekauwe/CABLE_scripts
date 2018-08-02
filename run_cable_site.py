#!/usr/bin/env python

"""
Run CABLE either for a single site, a subset, or all the flux sites pointed to
in the met directory

- Only intended for biophysics

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (01.08.2018)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import shutil
import subprocess
import multiprocessing as mp
import numpy as np

from cable_utils import adjust_nml_file
from cable_utils import get_svn_info
from cable_utils import change_LAI
from cable_utils import add_missing_options_to_nml_file
from cable_utils import add_attributes_to_output_file

class RunCable(object):

    def __init__(self, met_dir, log_dir, output_dir, restart_dir, aux_dir,
                 nml_fname, veg_fname, soil_fname, grid_fname, phen_fname,
                 cnpbiome_fname, lai_fname, fixed_lai, co2_conc, met_subset,
                 cable_src, mpi, verbose):

        self.met_dir = met_dir
        self.log_dir = log_dir
        self.output_dir = output_dir
        self.restart_dir = restart_dir
        self.aux_dir = aux_dir
        self.nml_fname = nml_fname
        self.biogeophys_dir = os.path.join(self.aux_dir, "core/biogeophys")
        self.grid_dir = os.path.join(self.aux_dir, "offline")
        self.biogeochem_dir = os.path.join(self.aux_dir, "core/biogeochem/")
        self.veg_fname = os.path.join(self.biogeophys_dir, veg_fname)
        self.soil_fname = os.path.join(self.biogeophys_dir, soil_fname)
        self.grid_fname = os.path.join(self.grid_dir, grid_fname)
        self.phen_fname = os.path.join(self.biogeochem_dir, phen_fname)
        self.cnpbiome_fname = os.path.join(self.biogeochem_dir, cnpbiome_fname)
        self.co2_conc = co2_conc
        self.met_subset = met_subset
        self.cable_src = cable_src
        self.cable_exe = os.path.join(self.cable_src, "offline/cable")
        self.verbose = verbose
        self.mpi = False
        self.lai_fname = lai_fname
        self.fixed_lai = fixed_lai

    def main(self):

        (met_files, url, rev) = self.initialise_stuff()

        if self.mpi:
            # Setup multi-processor jobs
            num_cpus = mp.cpu_count()
            chunk_size = int(np.ceil(len(met_files) / float(num_cpus)))

            pool = mp.Pool(processes=num_cpus)
            processes = []
            for i in range(num_cpus):
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
            site = os.path.basename(fname).split(".")[0]
            print(site)
            (out_fname, out_log_fname) = self.clean_up_old_files(site)

            if self.fixed_lai is not None or self.lai_fname is not None:
                (fname) = change_LAI(fname, fixed=self.fixed_lai,
                                     lai_fname=self.lai_fname)

            replace_dict = {
                            "filename%met": "'%s'" % (fname),
                            "filename%out": "'%s'" % (out_fname),
                            "filename%log": "'%s'" % (out_log_fname),
                            "filename%restart_out": "' '",
                            "filename%type": "'%s'" % (self.grid_fname),
                            "filename%veg": "'%s'" % (self.veg_fname),
                            "filename%soil": "'%s'" % (self.soil_fname),
                            "output%restart": ".FALSE.",
                            "fixedCO2": "%.2f" % (self.co2_conc),
                            "casafile%phen": "'%s'" % (self.phen_fname),
                            "casafile%cnpbiome": "'%s'" % (self.cnpbiome_fname),
                            "cable_user%FWSOIL_SWITCH": "'Haverd2013'",
                            "cable_user%GS_SWITCH": "'medlyn'",
                            "cable_user%or_evap": ".TRUE.",
            }
            adjust_nml_file(self.nml_fname, replace_dict)
            self.run_me()
            add_attributes_to_output_file(self.nml_fname, out_fname, url, rev)

            if self.fixed_lai is not None or self.lai_fname is not None:
                os.remove(fname)

    def initialise_stuff(self):

        if not os.path.exists(self.restart_dir):
            os.makedirs(self.restart_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        base_nml_fn = os.path.join(self.grid_dir, "%s" % (self.nml_fname))
        shutil.copy(base_nml_fn, self.nml_fname)
        add_missing_options_to_nml_file(self.nml_fname)

        # Run all the met files in the directory
        if len(met_subset) == 0:
            met_files = glob.glob(os.path.join(self.met_dir, "*.nc"))
        else:
            met_files = [os.path.join(self.met_dir, i) for i in self.met_subset]

        cwd = os.getcwd()
        (url, rev) = get_svn_info(cwd, self.cable_src)

        return (met_files, url, rev)

    def clean_up_old_files(self, site):
        out_fname = os.path.join(self.output_dir, "%s_out.nc" % (site))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        out_log_fname = os.path.join(self.log_dir, "%s_log.txt" % (site))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        return (out_fname, out_log_fname)

    def run_me(self):
        # run the model
        if self.verbose:
            os.system("%s" % (self.cable_exe))
        else:
            # No outputs to the screen, stout and stderr to dev/null
            os.system("%s > /dev/null 2>&1" % (self.cable_exe))



if __name__ == "__main__":

    #------------- Change stuff ------------- #
    met_dir = "/Users/mdekauwe/research/CABLE_runs/met_data/plumber_met"
    log_dir = "logs"
    output_dir = "outputs"
    restart_dir = "restart_files"
    aux_dir = "../../src/CMIP6-MOSRS/CABLE-AUX/"
    nml_fname = "cable.nml"
    veg_fname = "def_veg_params_zr_clitt_albedo_fix.txt"
    soil_fname = "def_soil_params.txt"
    grid_fname = "gridinfo_CSIRO_1x1.nc"
    phen_fname = "modis_phenology_csiro.txt"
    cnpbiome_fname = "pftlookup_csiro_v16_17tiles.csv"
    co2_conc = 380.0
    cable_src = "../../src/CMIP6-MOSRS/CMIP6-MOSRS"
    verbose = False
    mpi = False
    # if empty...run all the files in the met_dir
    met_subset = []#['TumbaFluxnet.1.4_met.nc']
    lai_fname = None
    fixed_lai = None
    # ------------------------------------------- #

    C = RunCable(met_dir, log_dir, output_dir, restart_dir, aux_dir,
                 nml_fname, veg_fname, soil_fname, grid_fname, phen_fname,
                 cnpbiome_fname, lai_fname, fixed_lai, co2_conc, met_subset,
                 cable_src, mpi, verbose)
    C.main()
