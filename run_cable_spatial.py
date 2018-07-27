#!/usr/bin/env python

"""
Run CABLE spatially.

This script sets various things within a qsub script and then submits the run.

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (27.07.2018)"
__email__ = "mdekauwe@gmail.com"

import subprocess
import sys
import os
import glob
import shutil
import tempfile
import optparse

class RunCable(object):

    def __init__(self, met_fname, log_fname):

        self.met_fname = met_fname
        self.log_fname = log_fname
        self.log_dir = log_dir
        self.nml_fn = nml_fn
        self.veg_param_fn = veg_param_fn
        self.grid_fn = grid_fn
        self.aux_dir = aux_dir
        self.veg_dir = os.path.join(self.aux_dir, "core/biogeophys")
        self.grid_dir = os.path.join(self.aux_dir, "offline")
        self.cable_exe = cable_exe
        self.verbose = verbose
        self.soil_fn = soil_fn

    def create_nml_file(self, log_fname):


        replace_dict = {
                        "filename%met": "'%s'" % (self.met_fname),
                        "filename%type": "'%s'" % (os.path.join(self.aux_dir, "offline/CABLE_UNSW_GSWP3_gridinfo_0.5x0.5.nc")),
                        "output%averaging": "'monthly'",

        }
        self.adjust_param_file(replace_dict)

    def adjust_nml_file(self, log_fname):

        out_log_fname = os.path.join(self.log_dir, log_fname)

        replace_dict = {
                        "filename%log: "'%s'" % (out_log_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_out": "' '",
                        "filename%type": "'%s'" % (os.path.join(self.grid_dir, self.grid_fn)),
                        "filename%veg": "'%s'" % (os.path.join(self.veg_dir, self.veg_param_fn)),
                        "filename%soil": "'%s'" % (os.path.join(self.veg_dir, self.soil_fn)),
                        "output%restart": ".FALSE.",
                        "fixedCO2": "380.0",
                        "output%averaging": "'%s'" % (average),

        }
        self.adjust_param_file(replace_dict)

    def run_me(start_yr, end_yr):

        qs_cmd = 'qsub -v start_yr=%d,end_yr=%d %s' % (start_yr, end_yr, template_fn)

        error = subprocess.call(qs_cmd, shell=True)
        if error is 1:
            raise("Job failed to submit")

    def adjust_param_file(self, replacements):
        """ adjust model parameters in the file and save over the original.

        Parameters:
        ----------
        fname : string
            parameter filename to be changed.
        replacements : dictionary
            dictionary of replacement values.

        """
        fin = open(self.nml_fn, 'r')
        param_str = fin.read()
        fin.close()
        new_str = self.replace_keys(param_str, replacements)
        fd, path = tempfile.mkstemp()
        os.write(fd, str.encode(new_str))
        os.close(fd)
        shutil.copy(path, self.nml_fn)
        os.remove(path)

if __name__ == "__main__":

    met_fname = "/g/data1/wd9/MetForcing/Global/GSWP3_2017/"
    log_dir = "logs"
    soil_fname = "def_soil_params.txt"
    veg_fname = "def_veg_params.txt"
    aux_dir = "/g/data1/w35/mrd561/CABLE/CABLE_AUX-dev/"
    co2_fname = "Annual_CO2_concentration_until_2010.txt"
    cable_aux_path = "/g/data1/w35/mrd561/CABLE/CABLE_AUX-dev/offline"

    if not os.path.exists(restart_dir):
        os.makedirs(restart_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    p = optparse.OptionParser()
    p.add_option('--person', '-p', default="world")
    p.add_option("-s", action="store_true", default=False,
                   help="Setup namelist file")
    p.add_option("-a", action="store_true", default=False,
                   help="Adjust namelist file")
    options, arguments = p.parse_args()

    C = RunCable(met_fname, log_fname, log_dir, aux_dir,
                 nml_fn, veg_fn, soil_fn, grid_fn, cable_exe, verbose)
    if options.s:
        C.create_nml_file()
    elif options.a:
        C.adjust_nml_file()
    elif options.r:
        start_yr = 1950
        end_yr = 1951
        C.run_me(start_yr, end_yr)
