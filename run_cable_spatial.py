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

from cable_utils import adjust_nml_file

class RunCable(object):

    def __init__(self, met_dir=None, log_dir=None, output_dir=None,
                 restart_dir=None, aux_dir=None, cable_src=None,
                 namelist_dir="namelists",
                 soil_fname="def_soil_params.txt",
                 veg_fname="def_veg_params_zr_clitt_albedo_fix.txt",
                 co2_fname="Annual_CO2_concentration_until_2010.txt",
                 grid_fname="CABLE_UNSW_GSWP3_gridinfo_0.5x0.5.nc",
                 mask_fname="gswp3_landmask_nomissing.nc",
                 nml_fname="cable.nml",
                 qsub_template_fname="run_cable_spatial_template.sh",
                 cable_exe="cable"):

        self.met_dir = met_dir
        self.log_dir = log_dir
        self.output_dir = output_dir
        self.aux_dir = aux_dir
        self.restart_dir = restart_dir
        self.veg_dir = os.path.join(self.aux_dir, "core/biogeophys")
        self.grid_dir = os.path.join(self.aux_dir, "offline")
        self.soil_fname = soil_fname
        self.biogeophys_dir = os.path.join(self.aux_dir, "core/biogeophys")
        self.biogeochem_dir = os.path.join(self.aux_dir, "core/biogeochem/")
        self.grid_dir = os.path.join(self.aux_dir, "offline")
        self.veg_fname = os.path.join(self.aux_dir, veg_fname)
        self.soil_fname = os.path.join(self.aux_dir, soil_fname)
        self.grid_fname = os.path.join(self.grid_dir, grid_fname)
        self.mask_fname = os.path.join(self.aux_dir,
                                       "offline/%s" % (mask_fname))
        self.namelist_dir = namelist_dir
        self.co2_fname = co2_fname
        self.qsub_template_fname = qsub_template_fname
        self.cable_src = cable_src
        self.cable_exe = os.path.join(cable_src, "offline/%s" % (cable_exe))
        self.nml_fname  = os.path.join(self.grid_dir, "%s" % (nml_fname))


    def initialise_stuff(self):

        if not os.path.exists(self.restart_dir):
            os.makedirs(self.restart_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        if not os.path.exists(self.namelist_dir):
            os.makedirs(self.namelist_dir)

        # delete local executable, copy a local copy and use that
        local_exe = "cable"
        if os.path.isfile(local_exe):
            os.remove(local_exe)
        shutil.copy(self.cable_exe, local_exe)
        self.cable_exe = local_exe

    def setup_nml_file(self, year):

        replace_dict = {
                        "filename%met": "'%s'" % (self.met_dir),
                        "filename%type": "'%s'" % (self.grid_fname),
                        "filename%veg": "'%s'" % (self.veg_fname),
                        "filename%soil": "'%s'" % (self.soil_fname),
                        "gswpfile%mask": "'%s'" % (self.mask_fname),
                        "gswpfile%rainf": "'%s'" % (os.path.join(self.met_dir, "Rainf/GSWP3.BC.Rainf.3hrMap.%s.nc" % (year))),
                        "gswpfile%snowf": "'%s'" % (os.path.join(self.met_dir, "Snowf/GSWP3.BC.Snowf.3hrMap.%s.nc" % (year))),
                        "gswpfile%LWdown": "'%s'" % (os.path.join(self.met_dir, "LWdown/GSWP3.BC.LWdown.3hrMap.%s.nc" % (year))),
                        "gswpfile%SWdown": "'%s'" % (os.path.join(self.met_dir, "SWdown/GSWP3.BC.SWdown.3hrMap.%s.nc" % (year))),
                        "gswpfile%PSurf": "'%s'" % (os.path.join(self.met_dir, "PSurf/GSWP3.BC.PSurf.3hrMap.%s.nc" % (year))),
                        "gswpfile%Qair": "'%s'" % (os.path.join(self.met_dir, "Qair/GSWP3.BC.Qair.3hrMap.%s.nc" % (year))),
                        "gswpfile%Tair": "'%s'" % (os.path.join(self.met_dir, "Tair/GSWP3.BC.Tair.3hrMap.%s.nc" % (year))),
                        "gswpfile%wind": "'%s'" % (os.path.join(self.met_dir, "Wind/GSWP3.BC.Wind.3hrMap.%s.nc" % (year))),
                        "output%averaging": "'monthly'",
                        "spinup": ".FALSE.",
                        "cable_user%FWSOIL_SWITCH": "'Haverd2013'",
                        "cable_user%GS_SWITCH": "'medlyn'",
                        "cable_user%GW_MODEL": ".FALSE.",
                        "cable_user%or_evap": ".TRUE.",
        }
        adjust_nml_file(self.nml_fname, replace_dict)


    def run_me(self, start_yr, end_yr):

        qs_cmd = 'qsub -v start_yr=%d,end_yr=%d,co2_fname=%s %s' % \
                    (start_yr, end_yr, self.co2_fname, self.qsub_template_fname)
        error = subprocess.call(qs_cmd, shell=True)
        if error is 1:
            print("Job failed to submit")
            sys.exit()

    def create_new_nml_file(self, log_fname, out_fname, restart_in_fname,
                            restart_out_fname, year, co2_conc):

        out_log_fname = os.path.join(self.log_dir, log_fname)
        out_fname = os.path.join(self.output_dir, out_fname)
        restart_in_fname = os.path.join(self.restart_dir, restart_in_fname)
        restart_out_fname = os.path.join(self.restart_dir, restart_out_fname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "filename%restart_in": "'%s'" % (restart_in_fname),
                        "filename%restart_out": "'%s'" % (restart_out_fname),
                        "fixedCO2": "%f" % (co2_conc),
                        "gswpfile%rainf": "'%s'" % (os.path.join(self.met_dir, "Rainf/GSWP3.BC.Rainf.3hrMap.%s.nc" % (year))),
                        "gswpfile%snowf": "'%s'" % (os.path.join(self.met_dir, "Snowf/GSWP3.BC.Snowf.3hrMap.%s.nc" % (year))),
                        "gswpfile%LWdown": "'%s'" % (os.path.join(self.met_dir, "LWdown/GSWP3.BC.LWdown.3hrMap.%s.nc" % (year))),
                        "gswpfile%SWdown": "'%s'" % (os.path.join(self.met_dir, "SWdown/GSWP3.BC.SWdown.3hrMap.%s.nc" % (year))),
                        "gswpfile%PSurf": "'%s'" % (os.path.join(self.met_dir, "PSurf/GSWP3.BC.PSurf.3hrMap.%s.nc" % (year))),
                        "gswpfile%Qair": "'%s'" % (os.path.join(self.met_dir, "Qair/GSWP3.BC.Qair.3hrMap.%s.nc" % (year))),
                        "gswpfile%Tair": "'%s'" % (os.path.join(self.met_dir, "Tair/GSWP3.BC.Tair.3hrMap.%s.nc" % (year))),
                        "gswpfile%wind": "'%s'" % (os.path.join(self.met_dir, "Wind/GSWP3.BC.Wind.3hrMap.%s.nc" % (year))),

        }
        adjust_nml_file(self.nml_fname, replace_dict)

        # save copy as we go for debugging - remove later
        shutil.copyfile(self.nml_fname, os.path.join(self.namelist_dir,
                                                     "cable_%d.nml" % (year)))

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
    met_dir = "/g/data1/wd9/MetForcing/Global/GSWP3_2017/"
    log_dir = "logs"
    output_dir = "outputs"
    restart_dir = "restarts"
    aux_dir = "/g/data1/w35/mgk576/research/CABLE_runs/src/trunk/CABLE-AUX"
    cable_src = "../../src/trunk/trunk/"
    start_yr = 1950
    end_yr = 1951
    # ------------------------------------------- #

    options, args = cmd_line_parser()

    C = RunCable(met_dir=met_dir, log_dir=log_dir, output_dir=output_dir,
                 restart_dir=restart_dir, aux_dir=aux_dir, cable_src=cable_src)
    C.initialise_stuff()

    # qsub script is adjusting namelist file
    if options.a:
        log_fname = options.l
        out_fname = options.o
        restart_in_fname = options.i
        restart_out_fname = options.r
        year = int(options.y)
        co2_conc = float(options.c)
        print(out_fname)
        C.create_new_nml_file(log_fname, out_fname, restart_in_fname,
                              restart_out_fname, year, co2_conc)

    # Setup initial namelist file and submit qsub job
    else:
        C.setup_nml_file(start_yr)
        C.run_me(start_yr, end_yr)
