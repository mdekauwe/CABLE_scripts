#!/usr/bin/env python

"""
Run CABLE spatially, i.e. global GSWP3 run.

Ensure you run it with "-s" if you want to do a spin-up run, i.e.

./run_cable_spatial.py -s

Once this is complete just run without the "-s" flag

The script does a few things internally:
- creates a qsub script.
- after the spin-up step it renames the final restart file to be the first
  simulation year restart file so that we can run a longer simulation
- submits the qsub script

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (14.07.2019)"
__email__ = "mdekauwe@gmail.com"

import subprocess
import sys
import os
import glob
import shutil
import tempfile
import optparse

from cable_utils import adjust_nml_file
from cable_utils import generate_spatial_qsub_script


def cmd_line_parser():

    p = optparse.OptionParser()
    p.add_option("-s", action="store_true", default=False,
                   help="Spinup model")
    p.add_option("-a", action="store_true", default=False,
                   help="Adjust namelist file")
    p.add_option("-t", action="store_true", default=False,
                   help="Sort restart files")
    p.add_option("-y", default="1900", help="year")
    p.add_option("-l", default="", help="log filename")
    p.add_option("-o", default="", help="out filename")
    p.add_option("-i", default="", help="restart in filename")
    p.add_option("-r", default="", help="restart out filename")
    p.add_option("-c", default="400.0", help="CO2 concentration")
    p.add_option("-n", default=None, help="nml_fname")
    options, args = p.parse_args()

    return (options.l, options.o, options.i,  options.r, int(options.y),
            float(options.c), options.n, options.s, options.a, options.t)


class RunCable(object):

    def __init__(self, met_dir=None, log_dir=None, output_dir=None,
                 restart_dir=None, aux_dir=None, cable_src=None, nml_fname=None,
                 spin_up=False, qsub_fname=None,
                 spinup_dir="spinup_restart",
                 namelist_dir="namelists",
                 soil_fname="def_soil_params.txt",
                 veg_fname="def_veg_params_zr_clitt_albedo_fix.txt",
                 co2_fname="Annual_CO2_concentration_until_2010.txt",
                 grid_fname="SE_AU_AWAP_NVIS_iveg_csiro_soil_gimms_lai_grid.nc",
                 mask_fname="SE_AUS_AWAP_csiro_soil_landmask.nc",
                 met_data="GSWP3",
                 cable_exe="cable-mpi", walltime=None, mem="64GB", ncpus="48"):

        self.met_dir = met_dir
        self.log_dir = log_dir
        self.output_dir = output_dir
        self.aux_dir = aux_dir
        self.restart_dir = restart_dir
        self.spinup_dir = spinup_dir
        self.grid_dir = os.path.join(self.aux_dir, "offline")
        self.soil_fname = soil_fname
        self.biogeophys_dir = os.path.join(self.aux_dir, "core/biogeophys")
        self.biogeochem_dir = os.path.join(self.aux_dir, "core/biogeochem/")
        self.veg_fname = os.path.join(self.biogeophys_dir, veg_fname)
        self.soil_fname = os.path.join(self.biogeophys_dir, soil_fname)
        self.grid_fname = os.path.join("SE_AUS_AWAP_grid_mask_files/grid", (grid_fname))
        self.mask_fname = os.path.join("SE_AUS_AWAP_grid_mask_files/mask", (mask_fname))
        self.namelist_dir = namelist_dir
        self.co2_fname = co2_fname
        self.qsub_fname = qsub_fname
        self.cable_src = cable_src
        self.cable_exe = os.path.join(cable_src, "offline/%s" % (cable_exe))
        self.spin_up = spin_up
        self.met_data = met_data

        if nml_fname is None:
            nml_fname = "cable.nml"
            base_nml_file = os.path.join(self.grid_dir, "%s" % (nml_fname))
            shutil.copyfile(base_nml_file, nml_fname)
            self.nml_fname = nml_fname
        else:
            self.nml_fname = nml_fname

        # qsub stuff
        self.walltime = walltime
        self.mem = mem
        self.ncpus = ncpus

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
        local_exe = os.path.basename(self.cable_exe)
        if os.path.isfile(local_exe):
            os.remove(local_exe)
        shutil.copy(self.cable_exe, local_exe)

    def setup_nml_file(self):

        replace_dict = {
                        "filename%met": "''",  # not needed for GSWP3 run
                        "filename%type": "'%s'" % (self.grid_fname),
                        "filename%veg": "'%s'" % (self.veg_fname),
                        "filename%soil": "'%s'" % (self.soil_fname),
                        "gswpfile%mask": "'%s'" % (self.mask_fname),
                        "output%averaging": "'monthly'",
                        "spinup": ".FALSE.",
                        "cable_user%FWSOIL_SWITCH": "'standard'",
                        "cable_user%GS_SWITCH": "'medlyn'",
                        "cable_user%GW_MODEL": ".FALSE.",
                        "cable_user%or_evap": ".FALSE.",
                        "cable_user%GSWP3": ".TRUE.",
                        "cable_user%MetType": "'gswp3'",
                        "verbose": ".FALSE.",
        }
        adjust_nml_file(self.nml_fname, replace_dict)

    def run_qsub_script(self, start_yr, end_yr):

        # Create a qsub script for simulations if missing, there is one of spinup
        # and one for simulations, so two qsub_fnames
        if not os.path.isfile(self.qsub_fname):
            generate_spatial_qsub_script(self.qsub_fname, self.walltime,
                                         self.mem, self.ncpus,
                                         spin_up=self.spin_up)

        # Run qsub script
        qs_cmd = 'qsub -v start_yr=%d,end_yr=%d,co2_fname=%s %s' % \
                    (start_yr, end_yr, self.co2_fname, self.qsub_fname)
        error = subprocess.call(qs_cmd, shell=True)
        if error is 1:
            raise("Job failed to submit\n")

    def create_new_nml_file(self, log_fname, out_fname, restart_in_fname,
                            restart_out_fname, year, co2_conc):

        out_log_fname = os.path.join(self.log_dir, log_fname)
        out_fname = os.path.join(self.output_dir, out_fname)

        # i.e. no restart file for first spinup year
        if restart_in_fname == "missing":
            restart_in_fname = ""
        else:
            restart_in_fname = os.path.join(self.restart_dir, restart_in_fname)
        restart_out_fname = os.path.join(self.restart_dir, restart_out_fname)

        if self.met_data == "GSWP3":
            rainf_fn = os.path.join(self.met_dir, "Rainf/GSWP3.BC.Rainf.3hrMap.%s.nc" % (year))
            snowf_fn = os.path.join(self.met_dir, "Snowf/GSWP3.BC.Snowf.3hrMap.%s.nc" % (year))
            lwdown_fn = os.path.join(self.met_dir, "LWdown/GSWP3.BC.LWdown.3hrMap.%s.nc" % (year))
            swdown_fn = os.path.join(self.met_dir, "SWdown/GSWP3.BC.SWdown.3hrMap.%s.nc" % (year))
            psurf_fn = os.path.join(self.met_dir, "PSurf/GSWP3.BC.PSurf.3hrMap.%s.nc" % (year))
            qair_fn = os.path.join(self.met_dir, "Qair/GSWP3.BC.Qair.3hrMap.%s.nc" % (year))
            wind_fn = os.path.join(self.met_dir, "Wind/GSWP3.BC.Wind.3hrMap.%s.nc" % (year))
            tair_fn = (os.path.join(self.met_dir, "Tair/GSWP3.BC.Tair.3hrMap.%s.nc" % (year)))
        elif self.met_data == "AWAP":
            rainf_fn = os.path.join(self.met_dir, "Rainf/AWAP.Rainf.3hr.%s.nc" % (year))
            snowf_fn = os.path.join(self.met_dir, "Snowf/AWAP.Snowf.3hr.%s.nc" % (year))
            lwdown_fn = os.path.join(self.met_dir, "LWdown/AWAP.LWdown.3hr.%s.nc" % (year))
            swdown_fn = os.path.join(self.met_dir, "SWdown/AWAP.SWdown.3hr.%s.nc" % (year))
            psurf_fn = os.path.join(self.met_dir, "PSurf/AWAP.PSurf.3hr.%s.nc" % (year))
            qair_fn = os.path.join(self.met_dir, "Qair/AWAP.Qair.3hr.%s.nc" % (year))
            wind_fn = os.path.join(self.met_dir, "Wind/AWAP.Wind.3hr.%s.nc" % (year))
            tair_fn = (os.path.join(self.met_dir, "Tair/AWAP.Tair.3hr.%s.nc" % (year)))

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "filename%restart_in": "'%s'" % (restart_in_fname),
                        "filename%restart_out": "'%s'" % (restart_out_fname),
                        "fixedCO2": "%f" % (co2_conc),
                        "ncciy": "%s" % (year), # 0 for not using gswp; 4-digit year input for year of gswp met
                        "CABLE_USER%YearStart": "0", # needs to be 0 so the ncciy is set
                        "CABLE_USER%YearEnd": "0",   # needs to be 0 so the ncciy is set
                        "gswpfile%rainf": "'%s'" % (rainf_fn),
                        "gswpfile%snowf": "'%s'" % (snowf_fn),
                        "gswpfile%LWdown": "'%s'" % (lwdown_fn),
                        "gswpfile%SWdown": "'%s'" % (swdown_fn),
                        "gswpfile%PSurf": "'%s'" % (psurf_fn),
                        "gswpfile%Qair": "'%s'" % (qair_fn),
                        "gswpfile%Tair": "'%s'" % (tair_fn),
                        "gswpfile%wind": "'%s'" % (wind_fn),
        }
        adjust_nml_file(self.nml_fname, replace_dict)

        # save copy as we go for debugging - remove later
        shutil.copyfile(self.nml_fname, os.path.join(self.namelist_dir,
                                                     "cable_%d.nml" % (year)))

    def sort_restart_files(self, start_yr, end_yr):

        if not os.path.exists(self.spinup_dir):
            os.makedirs(self.spinup_dir)

            # Copy the last spinup restart file to the backup dir and rename
            # it as if it was the first year
            fn_in = "restart_%d.nc" % (end_yr)
            fn_out = "restart_%d.nc" % (start_yr)

            restart_in_fname = os.path.join(self.restart_dir, fn_in)
            restart_out_fname = os.path.join(self.spinup_dir, fn_out)

            shutil.copyfile(restart_in_fname, restart_out_fname)

            # remove the restart dir and remake it with the equilibrium file
            shutil.rmtree(self.restart_dir, ignore_errors=True)

            if not os.path.exists(self.restart_dir):
                os.makedirs(self.restart_dir)

            fn_in = "restart_%d.nc" % (start_yr)
            restart_in_fname = os.path.join(self.spinup_dir, fn_in)
            restart_out_fname = os.path.join(self.restart_dir, fn_in)
            shutil.copyfile(restart_in_fname, restart_out_fname)



if __name__ == "__main__":

    #------------- Change stuff ------------- #
    #met_data = "AWAP"
    met_data = "GSWP3"
    if met_data == "GSWP3":
        met_dir = "/g/data/wd9/MetForcing/Global/GSWP3_2017/"
    elif met_data == "AWAP":
        met_dir = "/g/data1a/w35/mgk576/research/AWAP_interpolation/interpolated"

    log_dir = "logs"
    output_dir = "outputs"
    restart_dir = "restarts"
    aux_dir = "/g/data/w35/mgk576/research/CABLE_runs/src/CABLE-AUX"
    #cable_src = "../../src/trunk/trunk/"
    cable_src = "../../src/trunk_DESICA_PFTs/trunk_DESICA_PFTs/"
    spinup_start_yr = 1995
    #spinup_end_yr = 1995
    spinup_end_yr = 2000
    run_start_yr = 2000
    run_end_yr = 2010
    # ------------------------------------------- #

    (log_fname, out_fname, restart_in_fname,
     restart_out_fname, year, co2_conc,
     nml_fname, spin_up, adjust_nml, sort_restarts) = cmd_line_parser()

    if spin_up:
        start_yr = spinup_start_yr
        end_yr = spinup_end_yr
        walltime = "4:00:00"
        qsub_fname = "qsub_wrapper_script_spinup.sh"
    else:
        start_yr = run_start_yr
        end_yr = run_end_yr
        walltime = "7:30:00"
        qsub_fname = "qsub_wrapper_script_simulation.sh"

    C = RunCable(met_dir=met_dir, log_dir=log_dir, output_dir=output_dir,
                 restart_dir=restart_dir, aux_dir=aux_dir, spin_up=spin_up,
                 cable_src=cable_src, qsub_fname=qsub_fname, met_data=met_data,
                 nml_fname=nml_fname, walltime=walltime)

    # Sort the restart files out before we run simulations "-t"
    if sort_restarts:
        C.sort_restart_files(spinup_start_yr, spinup_end_yr)
        sys.exit('Restart files fixed up, run simulation')

    # Setup initial namelist file and submit qsub job
    if adjust_nml == False:
        C.initialise_stuff()
        C.setup_nml_file()
        C.run_qsub_script(start_yr, end_yr)

    # qsub script is adjusting namelist file, i.e. for a different year
    else:
        C.create_new_nml_file(log_fname, out_fname, restart_in_fname,
                              restart_out_fname, year, co2_conc)
