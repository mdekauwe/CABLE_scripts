#!/usr/bin/env python

"""
Run CABLE spatially with CNP turned on

To write...



That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (04.02.2020)"
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
    p.add_option("--s1", action="store_true", default=False,
                   help="Spinup model: first round")
    p.add_option("--s2", action="store_true", default=False,
                   help="Spinup model: second round")
    p.add_option("-a", action="store_true", default=False,
                   help="Adjust namelist file")
    p.add_option("-t", action="store_true", default=False,
                   help="Sort restart files")
    p.add_option("-y", default="1900", help="year")
    p.add_option("--rn", default="1", help="rst_num")
    p.add_option("-l", default="", help="log filename")
    p.add_option("-o", default="", help="out filename")
    p.add_option("-i", default="", help="restart in filename")
    p.add_option("-r", default="", help="restart out filename")
    p.add_option("--ci", default="", help="casa restart in filename")
    p.add_option("--cr", default="", help="casa restart out filename")
    p.add_option("-c", default="400.0", help="CO2 concentration")
    p.add_option("-n", default=None, help="nml_fname")
    options, args = p.parse_args()

    return (options.l, options.o, options.i,  options.r, options.ci,
            options.cr, int(options.y),float(options.c), options.n, options.s1,
            options.s2, options.a, options.t, int(options.rn))


class RunCable(object):

    def __init__(self, met_dir=None, log_dir=None, output_dir=None,
                 restart_dir=None, aux_dir=None, cable_src=None, nml_fname=None,
                 spinup_dir="spinup_files",
                 namelist_dir="namelists",
                 soil_fname="def_soil_params.txt",
                 veg_fname="def_veg_params_zr_clitt_albedo_fix.txt",
                 co2_fname="Annual_CO2_concentration_until_2010.txt",
                 grid_fname="gridinfo_mmy_MD_elev_orig_std_casa_avg-sand_mask.nc",
                 phen_fname="modis_phenology_csiro.txt",
                 #cnpbiome_fname="pftlookup_csiro_v16_17tiles_Ticket2.csv",
                 cnpbiome_fname="pftlookup_csiro_v16_17tiles-cheng-m02.csv",
                 mask_fname="gswp3_landmask_nomissing.nc",
                 biogeochem="C", co2_conc=400.0, co2_fixed=284.7,
                 ndep_fixed=0.79, pdep_fixed=0.144,
                 experiment_name="GSWP3_CNP",
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
        self.grid_fname = grid_fname
        self.mask_fname = os.path.join("land_sea_mask/%s" % (mask_fname))
        self.phen_fname = os.path.join(self.biogeochem_dir, phen_fname)
        #self.cnpbiome_fname = os.path.join(self.biogeochem_dir, cnpbiome_fname)
        self.cnpbiome_fname = os.path.join("CNP_param_file", cnpbiome_fname)

        self.namelist_dir = namelist_dir
        self.co2_fname = co2_fname
        self.cable_src = cable_src
        self.cable_exe = os.path.join(cable_src, "offline/%s" % (cable_exe))
        self.co2_fixed = co2_fixed    # umol mol-1
        self.ndep_fixed = ndep_fixed  # kg N ha-1 yr-1
        self.pdep_fixed = pdep_fixed  # kg N ha-1 yr-1
        self.biogeochem_cyc = biogeochem

        if self.biogeochem_cyc == "C":
            self.biogeochem_id = 1
            self.vcmax_feedback = ".FALSE."
        elif self.biogeochem_cyc == "CN":
            self.biogeochem_id = 2
            self.vcmax_feedback = ".TRUE."
        elif self.biogeochem_cyc == "CNP":
            self.biogeochem_id = 3
            self.vcmax_feedback = ".TRUE."

        self.experiment_id = "%s_%s" % (experiment_name, self.biogeochem_cyc)

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

    def setup_nml_file(self, start_yr):

        # set directory paths...
        out_fname = "%s_out_cable_spin.nc" % (self.experiment_id)
        out_fname = os.path.join(self.output_dir, out_fname)

        out_log_fname = "%s_log_spin.txt" % (self.experiment_id)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)

        cable_rst_ifname = ""
        cable_rst_ofname = os.path.join(self.restart_dir,
                                        "cable_restart_%d.nc" % (start_yr))

        casa_rst_ifname = ""
        casa_rst_ofname = os.path.join(self.restart_dir,
                                       "casa_restart_%d.nc" % (start_yr))

        replace_dict = {
                        "filename%met": "''",  # not needed for GSWP3 run
                        "filename%type": "'%s'" % (self.grid_fname),
                        "filename%veg": "'%s'" % (self.veg_fname),
                        "filename%soil": "'%s'" % (self.soil_fname),
                        "gswpfile%mask": "'%s'" % (self.mask_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_in": "'%s'" % (cable_rst_ifname),
                        "filename%restart_out": "'%s'" % (cable_rst_ofname),
                        "casafile%cnpipool": "'%s'" % (casa_rst_ifname),
                        "casafile%cnpepool": "'%s'" % (casa_rst_ofname),
                        "fixedCO2": "%.2f" % (self.co2_fixed),
                        "casafile%phen": "'%s'" % (self.phen_fname),
                        "casafile%cnpbiome": "'%s'" % (self.cnpbiome_fname),
                        "cable_user%RunIden": "'%s'" % (self.experiment_id),
                        "cable_user%vcmax": "'standard'",
                        "l_vcmaxFeedbk": "%s" % (self.vcmax_feedback),
                        "l_laiFeedbk": ".TRUE.", # prognoistic LAI
                        "icycle": "%d" % (self.biogeochem_id),
                        "output%averaging": "'annually'",
                        "cable_user%CASA_OUT_FREQ": "'annually'",
                        "output%casa": ".TRUE.",
                        "leaps": ".FALSE.",
                        "cable_user%CASA_fromZero": ".TRUE.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "cable_user%CASA_NREP": "0",
                        "spinup": ".FALSE.",
                        "spincasa": ".FALSE.",
                        "output%restart": ".TRUE.",
                        "cable_user%FWSOIL_SWITCH": "'standard'",
                        "cable_user%GS_SWITCH": "'medlyn'",
                        #"cable_user%limit_labile": ".TRUE.",
                        "cable_user%FWSOIL_SWITCH": "'standard'",
                        "cable_user%GS_SWITCH": "'medlyn'",
                        "cable_user%GW_MODEL": ".FALSE.",
                        "cable_user%or_evap": ".FALSE.",
                        "cable_user%GSWP3": ".TRUE.",
                        "cable_user%MetType": "'gswp3'",
                        "verbose": ".TRUE.",
        }
        adjust_nml_file(self.nml_fname, replace_dict)

    def run_qsub_script(self, qsub_fname, start_yr, end_yr, spin_up):

        # Create a qsub script for simulations if missing, there is one of spinup
        # and one for simulations, so two qsub_fnames
        if not os.path.isfile(qsub_fname):
            generate_spatial_qsub_script(qsub_fname, self.walltime,
                                         self.mem, self.ncpus,
                                         spin_up=spin_up, CNP=True)


        # Run qsub script
        qs_cmd = 'qsub -v start_yr=%d,end_yr=%d,co2_fname=%s %s' % \
                    (start_yr, end_yr, self.co2_fname, qsub_fname)

        error = subprocess.call(qs_cmd, shell=True)
        if error is 1:
            raise("Job failed to submit\n")

    def create_new_nml_file(self, log_fname, out_fname, cable_rst_ifname,
                            cable_rst_ofname, casa_rst_ifname, casa_rst_ofname,
                            year, co2_conc):

        out_log_fname = os.path.join(self.log_dir, log_fname)
        out_fname = os.path.join(self.output_dir, out_fname)

        # i.e. no restart file for first spinup year
        if cable_restart_ifname == "missing":
            cable_rst_ifname = ""
            casa_rst_ifname = ""
        else:
            cable_rst_ifname = os.path.join(self.restart_dir, cable_rst_ifname)
            casa_rst_ifname = os.path.join(self.restart_dir, casa_rst_ifname)

        cable_rst_ofname = os.path.join(self.restart_dir, cable_rst_ofname)
        casa_rst_ofname = os.path.join(self.restart_dir, casa_rst_ofname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "filename%restart_in": "'%s'" % (cable_rst_ifname),
                        "filename%restart_out": "'%s'" % (cable_rst_ofname),
                        "casafile%cnpipool": "'%s'" % (casa_rst_ifname),
                        "casafile%cnpepool": "'%s'" % (casa_rst_ofname),

                        "fixedCO2": "%f" % (co2_conc),
                        "ncciy": "%s" % (year), # 0 for not using gswp; 4-digit year input for year of gswp met
                        "CABLE_USER%YearStart": "0", # needs to be 0 so the ncciy is set
                        "CABLE_USER%YearEnd": "0",   # needs to be 0 so the ncciy is set
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

    def sort_restart_files(self, start_yr, end_yr, number):

        if not os.path.exists(self.spinup_dir):
            os.makedirs(self.spinup_dir)

        # Copy the last spinup restart file to the backup dir and rename
        # it as if it was the first year
        fn_in = "cable_restart_%d.nc" % (end_yr)
        fn_out = "cable_restart_%d_%d.nc" % (start_yr, number)
        restart_in_fname = os.path.join(self.restart_dir, fn_in)
        restart_out_fname = os.path.join(self.spinup_dir, fn_out)
        shutil.copyfile(restart_in_fname, restart_out_fname)

        # Copy the output file to check for stability
        fn_in = "cable_out_%d.nc" % (end_yr)
        fn_out = "cable_out_s_%d.nc" % (number)
        old_fname = os.path.join(self.output_dir, fn_in)
        new_fname = os.path.join(self.spinup_dir, fn_out)
        shutil.copyfile(old_fname, new_fname)

        # remove the restart dir and remake it with the equilibrium file
        shutil.rmtree(self.restart_dir, ignore_errors=True)
        if not os.path.exists(self.restart_dir):
            os.makedirs(self.restart_dir)

        fn_in = "cable_restart_%d_%d.nc" % (start_yr, number)
        fn_out = "cable_restart_%d.nc" % (start_yr)
        restart_in_fname = os.path.join(self.spinup_dir, fn_in)
        restart_out_fname = os.path.join(self.restart_dir, fn_out)
        shutil.copyfile(restart_in_fname, restart_out_fname)



if __name__ == "__main__":

    #------------- Change stuff ------------- #
    met_dir = "/g/data/wd9/MetForcing/Global/GSWP3_2017/"
    log_dir = "logs"
    output_dir = "outputs"
    restart_dir = "restarts"
    #aux_dir = "/g/data/w35/mgk576/research/CABLE_runs/src/CABLE-AUX"
    aux_dir = "../../src/CABLE-AUX"
    #cable_src = "../../src/trunk/trunk/"
    cable_src = "../../src/trunk_cnp_spatial/trunk_cnp_spatial/"
    spinup_start_yr = 1901 # GSWP3 starts in 1901
    spinup_end_yr = 1930
    run_start_yr = 1924
    run_end_yr = 1950
    biogeochem = "C"
    experiment_name = "GSWP3_CNP"
    # ------------------------------------------- #

    (log_fname, out_fname, cable_rst_ifname,
     cable_rst_ofname, casa_rst_ifname,
     casa_rst_ofname, year, co2_conc,
     nml_fname, spin_up_step_one, spin_up_step_two,
     adjust_nml, sort_restarts,
     rst_num) = cmd_line_parser()

    C = RunCable(met_dir=met_dir, log_dir=log_dir, output_dir=output_dir,
                 restart_dir=restart_dir, aux_dir=aux_dir,
                 cable_src=cable_src, nml_fname=nml_fname, walltime=walltime,
                 biogeochem=biogeochem, experiment_name=experiment_name)

    # First spin up round....
    if spin_up_step_one:
        start_yr = spinup_start_yr
        end_yr = spinup_end_yr
        walltime = "6:00:00"
        qsub_fname = "qsub_wrapper_script_spinup.sh"
    # Second spin up round, here we're reusing the restart files...
    elif spin_up_step_two:
        start_yr = spinup_start_yr
        end_yr = spinup_end_yr
        walltime = "6:00:00"
        qsub_fname = "qsub_wrapper_script_simulation.sh"
        C.sort_restart_files(spinup_start_yr, spinup_end_yr, rst_num)
    else:
        start_yr = run_start_yr
        end_yr = run_end_yr
        walltime = "6:00:00"
        qsub_fname = "qsub_wrapper_script_simulation.sh"

    # Setup initial namelist file and submit qsub job
    if adjust_nml == False:
        C.initialise_stuff()
        C.setup_nml_file(start_yr)
        C.run_qsub_script(qsub_fname, start_yr, end_yr, spin_up=True)

    # qsub script is adjusting namelist file, i.e. for a different year
    else:
        C.create_new_nml_file(log_fname, out_fname, cable_rst_ifname,
                              cable_rst_ofname, casa_rst_ifname,
                              casa_rst_ofname, year, co2_conc)
