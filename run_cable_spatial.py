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

    def __init__(self, met_path, log_dir, output_dir, aux_dir, soil_fname,
                 veg_fname, co2_fname, grid_fname):

        self.met_path = met_path
        self.log_dir = log_dir
        self.output_dir = output_dir
        self.veg_dir = os.path.join(self.aux_dir, "core/biogeophys")
        self.grid_dir = os.path.join(self.aux_dir, "offline")
        self.soil_fname = soil_fname
        self.veg_fname = os.path.join(self.aux_dir, veg_fname)
        self.soil_fname = os.path.join(self.aux_dir, soil_fname)
        self.grid_fname = os.path.join(self.grid_dir, grid_fname)
        self.co2_fname = co2_fname
        self.grid_fname = grid_fname
        self.cable_exe = exe
        self.verbose = verbose

    def setup_nml_file(self):

        replace_dict = {
                        "filename%met": "'%s'" % (self.met_path),
                        "filename%type": "'%s'" % (self.grid_fname),
                        "filename%veg": "'%s'" % (self.veg_fname),
                        "filename%soil": "'%s'" % (self.soil_fname),
                        "output%averaging": "'monthly'",
                        "spinup": ".FALSE.",
                        "gswpfile%mask": "./surface_data/gswp3_landmask_nomissing.nc",

        }
        self.adjust_param_file(replace_dict)

    def adjust_nml_file(self, log_fname, out_fname, restart_in_fname,
                        restart_out_fname, year, co2_conc):

        out_log_fname = os.path.join(self.log_dir, log_fname)
        out_fname = os.path.join(self.output_dir, out_fname)
        restart_in_fname = os.path.join(self.restart_dir, restart_in_fname)
        restart_out_fname = os.path.join(self.restart_dir, restart_out_fname)

        replace_dict = {
                        "filename%log: "'%s'" % (out_log_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_in_fname": "'%s'" % (restart_in_fname),
                        "filename%restart_out_fname": "'%s'" % (restart_out_fname),
                        "fixedCO2": "%d" % (co2_conc),
                        "gswpfile%rainf": "/gswp/Rainf/GSWP3.BC.Rainf.3hrMap.%s.nc" % (year),
                        "gswpfile%snowf": "/gswp/Snowf/GSWP3.BC.Snowf.3hrMap.%s.nc" % (year),
                        "gswpfile%LWdown": "/gswp/LWdown/GSWP3.BC.LWdown.3hrMap.%s.nc" % (year),
                        "gswpfile%SWdown": "/gswp/SWdown/GSWP3.BC.SWdown.3hrMap.%s.nc" % (year),
                        "gswpfile%PSurf": "/gswp/PSurf/GSWP3.BC.PSurf.3hrMap.%s.nc" % (year),
                        "gswpfile%Qair": "/gswp/Qair/GSWP3.BC.Qair.3hrMap.%s.nc" % (year),
                        "gswpfile%Tair": "/gswp/Tair/GSWP3.BC.Tair.3hrMap.%s.nc" % (year),
                        "gswpfile%wind": "/gswp/Wind/GSWP3.BC.Wind.3hrMap.%s.nc" % (year),

        }
        self.adjust_param_file(replace_dict)

    def run_me(start_yr, end_yr):

        qs_cmd = 'qsub -v start_yr=%d,end_yr=%d %s' % \
                    (start_yr, end_yr, template_fn)

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

def cmd_line_parser():

    p = optparse.OptionParser()
    p.add_option('--person', '-p', default="world")
    p.add_option("-s", action="store_true", default=False,
                   help="Setup namelist file")
    p.add_option("-a", action="store_true", default=False,
                   help="Adjust namelist file")
    p.add_option("-y", default="0", help="year")
    p.add_option("-l", default="", help="log filename")
    p.add_option("-o", default="", help="out filename")
    p.add_option("-i", default="", help="restart in filename")
    p.add_option("-r", default="", help="restart out filename")
    p.add_option("-c", default="400.0", help="CO2 concentration")

    return p.parse_args()

if __name__ == "__main__":

    #------------- Change stuff ------------- #
    met_path = "/g/data1/wd9/MetForcing/Global/GSWP3_2017/"
    log_dir = "logs"
    soil_fname = "def_soil_params.txt"
    veg_fname = "def_veg_params.txt"
    aux_dir = "/g/data1/w35/mrd561/CABLE/CABLE_AUX-dev/"
    co2_fname = "Annual_CO2_concentration_until_2010.txt"
    grid_fname = "CABLE_UNSW_GSWP3_gridinfo_0.5x0.5.nc"
    output_dir = "outputs"
    start_yr = 1950
    end_yr = 1951
    #------------- Change stuff ------------- #

    if not os.path.exists(restart_dir):
        os.makedirs(restart_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    options, args = p.cmd_line_parser()

    C = RunCable(met_path, log_dir, output_dir, aux_dir, soil_fname,
                 veg_fname, co2_fname, grid_fname)

    # qsub script is adjusting namelist file
    if options.a:
        log_fname = options.l
        out_fname = options.o
        restart_in_fname = options.i
        restart_out_fname = options.r
        year = int(options.yr)
        co2_conc = int(options.c)
        C.adjust_nml_file(log_fname, out_fname, restart_in_fname,
                          restart_out_fname, year, co2_conc)

    # Setup initial namelist file and submit qsub job
    else:
        C.setup_nml_file()
        C.run_me(start_yr, end_yr)
