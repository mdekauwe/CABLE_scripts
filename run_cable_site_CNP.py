#!/usr/bin/env python

"""
CABLE site run with full biogeochemistry (CNP) and POP
=======================================================

- Spin-up: using pre-industrial CO2, NDEP, PDEP. Currently this is using
           values from AmazonFACE experiment, i.e. CO2=284.7;
           NDEP-0.79 kg N ha-1 yr-1; PDEP=0.144 kg P ha-1 yr-1
- Transient: 1851-1998, varying CO2, NDEP and PDEP with recycled met data
- Simulation: actual experiment CO2, NDEP/PDEP & met

Options to turn use C/CN/CNP with POP on/off

That's all folks.

"""

__author__ = "Martin De Kauwe, Vanessa Haverd"
__version__ = "1.0 (21.07.2018)"
__email__ = "mdekauwe@gmail.com, Vanessa.Haverd@csiro.au"

import os
import sys
import glob
import shutil
import tempfile
import pandas as pd
import xarray as xr
import numpy as np
import subprocess

from cable_utils import adjust_nml_file
from cable_utils import get_svn_info
from cable_utils import get_years
from cable_utils import add_attributes_to_output_file
from cable_utils import check_steady_state

class RunCable(object):

    def __init__(self, experiment_id=None, met_dir=None, dump_dir=None,
                 driver_dir=None, log_dir=None, output_dir=None,
                 co2_ndep_dir=None, restart_dir=None, aux_dir=None,
                 nml_fname="cable.nml", site_nml_fname="site.nml",
                 veg_fname="def_veg_params_zr_clitt_albedo_fix.txt",
                 soil_fname="def_soil_params.txt",
                 grid_fname="gridinfo_CSIRO_1x1.nc",
                 phen_fname="modis_phenology_csiro.txt",
                 cnpbiome_fname="pftlookup.csv", met_fname=None,
                 co2_ndep_fname="AmaFACE_co2npdepforcing_1850_2100_AMB.csv",
                 cable_src=None, biogeochem="C", use_pop=False, use_sli=False,
                 use_clim=False, verbose=True):

        self.experiment_id = experiment_id
        self.met_dir = met_dir
        self.dump_dir = dump_dir
        self.driver_dir = driver_dir
        self.log_dir = log_dir
        self.output_dir = output_dir
        self.co2_ndep_dir = co2_ndep_dir
        self.restart_dir = restart_dir
        self.aux_dir = aux_dir
        self.biogeophys_dir = os.path.join(self.aux_dir, "core/biogeophys")
        self.grid_dir = os.path.join(self.aux_dir, "offline")
        self.biogeochem_dir = os.path.join(self.aux_dir, "core/biogeochem/")
        self.nml_fname = nml_fname
        # set as we go, but we need access to this to add svn info
        self.out_fname = None
        self.out_fname_CASA = None
        self.site_nml_fname = site_nml_fname
        self.veg_fname = os.path.join(self.biogeophys_dir, veg_fname)
        self.soil_fname = os.path.join(self.biogeophys_dir, soil_fname)

        #self.veg_fname = os.path.join(self.driver_dir, veg_fname)
        #self.soil_fname = os.path.join(self.driver_dir, soil_fname)
        self.grid_fname = os.path.join(self.grid_dir, grid_fname)
        self.phen_fname = os.path.join(self.biogeochem_dir, phen_fname)
        #self.cnpbiome_fname = os.path.join(self.biogeochem_dir, cnpbiome_fname)
        self.cnpbiome_fname = os.path.join(self.driver_dir, cnpbiome_fname)
        self.met_fname = os.path.join(self.met_dir, met_fname)
        self.co2_ndep_fname = os.path.join(co2_ndep_dir, co2_ndep_fname)
        self.cable_restart_fname = os.path.join(self.restart_dir, \
                                    "%s_cable_rst.nc" % (self.experiment_id))
        self.casa_restart_fname = os.path.join(self.restart_dir, \
                                    "%s_casa_rst.nc" % (self.experiment_id))
        self.pop_restart_fname = os.path.join(self.restart_dir, \
                                    "%s_pop_rst.nc" % (self.experiment_id))
        self.cable_src = cable_src
        self.cable_exe = os.path.join(self.cable_src, "offline/cable")
        self.biogeochem = biogeochem
        self.use_pop = use_pop
        self.use_sli = use_sli
        self.use_clim = use_clim
        self.verbose = verbose
        self.nyear_spinup = 30

        if biogeochem == "C":
            self.biogeochem = 1
            #self.vcmax = "standard"
            self.vcmax = "Walker2014"
            self.vcmax_feedback = ".TRUE."
        elif biogeochem == "CN":
            self.biogeochem = 2
            self.vcmax = "Walker2014"
            self.vcmax_feedback = ".TRUE."
        elif biogeochem == "CNP":
            self.biogeochem = 3
            self.vcmax = "Walker2014"
            self.vcmax_feedback = ".TRUE."
        else:
            raise ValueError("Unknown biogeochemistry option: C, CN, CNP")

        if self.use_pop:
            self.pop_flag = ".TRUE."
        else:
            self.pop_flag = ".FALSE."

        if self.use_sli:
            self.soil_flag = "sli"
        else:
            self.soil_flag = "default"

        if self.use_clim:
            self.clim_flag = ".TRUE."
            self.climate_restart_fname = os.path.join(self.restart_dir, \
                                    "%s_climate_rst.nc" % (self.experiment_id))
        else:
            self.clim_flag = ".FALSE."
            self.climate_restart_fname = ""

    def main(self, SPIN_UP=False, TRANSIENT=False, SIMULATION=False):

        num = 1
        not_in_equilibrium = True

        (st_yr, en_yr,
         st_yr_trans, en_yr_trans,
         st_yr_spin, en_yr_spin) = get_years(self.met_fname, self.nyear_spinup)

        (url, rev) = self.initial_setup(st_yr_spin, en_yr_spin, st_yr, en_yr)

        if SPIN_UP == True:

            # initial spin
            print("First Spinup\n")
            self.run_me()
            self.clean_up(url, rev, end=False, tag="zero")

            while not_in_equilibrium:
                print("Spinup stage %d\n" % (num))
                self.logfile="log_ccp%d" % (num)
                self.setup_re_spin(number=num)
                self.run_me()
                self.clean_up(url, rev, end=False, tag="ccp%d" % (num))

                print("Analytical stage %d\n" % (num))
                self.logfile="log_sa%d" % (num)
                self.setup_analytical_spin(number=num, st_yr_spin=st_yr_spin,
                                           en_yr_spin=en_yr_spin )
                self.run_me()
                #self.clean_up(url, rev, end=False, tag="saa%d" % (num))

                not_in_equilibrium = check_steady_state(self.experiment_id,
                                                        self.output_dir, num)

                num += 1

            # one final spin
            print("Final spin\n")
            self.logfile="log_ccp%d" % (num)
            self.setup_re_spin(number=num)
            self.run_me()
            self.clean_up(url, rev, end=False, tag="ccp%d" % (num))

        if TRANSIENT == True:
            print("Transient")
            self.setup_transient(st_yr_trans, en_yr_trans, st_yr, en_yr)
            self.run_me()
            self.clean_up(url, rev, end=False, tag="transient")

        if SIMULATION == True:
            print("Simulation")
            self.setup_simulation(st_yr, en_yr)
            self.run_me()

        self.clean_up(url, rev, end=True)

    def initial_setup(self, st_yr_spin, en_yr_spin, st_yr, en_yr):
        """
        Setup CABLE namelist file for spinup from zero
        """

        if not os.path.exists(self.restart_dir):
            os.makedirs(self.restart_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        if not os.path.exists(self.dump_dir):
            os.makedirs(self.dump_dir)


        base_nml_fn = os.path.join(self.grid_dir, "%s" % (self.nml_fname))

        if os.path.isfile(self.nml_fname):
            os.remove(self.nml_fname)
        shutil.copy(base_nml_fn, self.nml_fname)
        #shutil.copyfile(os.path.join(self.driver_dir, "cable.nml"),
        #                self.nml_fname)
        shutil.copyfile(os.path.join(self.driver_dir, "site.nml"),
                        self.site_nml_fname)

        self.out_fname = os.path.join(self.output_dir,
                                 "%s_out_cable_zero.nc" % (self.experiment_id))
        if os.path.isfile(self.out_fname):
            os.remove(self.out_fname)

        self.out_fname_CASA = os.path.join(self.output_dir,
                                 "%s_out_CASA_zero.nc" % (self.experiment_id))
        if os.path.isfile(self.out_fname_CASA):
            os.remove(self.out_fname_CASA)

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_zero.txt" % (self.experiment_id))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        replace_dict = {
                        "RunType": '"spinup"',
                        "CO2NDepFile": "'%s'" % (self.co2_ndep_fname),
                        "spinstartyear": "%d" % (st_yr),
                        "spinendyear": "%d" % (en_yr),
                        "spinCO2": "284.7",
                        "spinNdep": "0.79",
                        "spinPdep": "0.144",
        }
        adjust_nml_file(self.site_nml_fname, replace_dict)

        replace_dict = {
                        "filename%met": "'%s'" % (self.met_fname),
                        "filename%out": "'%s'" % (self.out_fname),
                        "casafile%out": "'%s'" % (self.out_fname_CASA),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_out": "'%s'" % (self.cable_restart_fname),
                        "cable_user%climate_restart_out": "'%s'" % (self.climate_restart_fname),
                        "cable_user%POP_restart_out": "'%s'" % (self.pop_restart_fname),
                        "casafile%cnpepool": "'%s'" % (self.casa_restart_fname),
                        "filename%restart_in": "''" ,
                        "cable_user%climate_restart_in": "''" ,
                        "cable_user%POP_restart_in": "''",
                        "filename%type": "'%s'" % self.grid_fname,
                        "filename%veg": "'%s'" % (self.veg_fname),
                        "filename%soil": "'%s'" % (self.soil_fname),
                        "output%restart": ".TRUE.",
                        "casafile%phen": "'%s'" % (self.phen_fname),
                        "casafile%cnpbiome": "'%s'" % (self.cnpbiome_fname),
                        "cable_user%RunIden": "'%s'" % (self.experiment_id),
                        "cable_user%POP_out": "'ini'",
                        "cable_user%POP_rst": "'./'",
                        "cable_user%POP_fromZero": ".T.",
                        "cable_user%CASA_fromZero": ".T.",
                        "cable_user%CLIMATE_fromZero": ".T.",
                        "cable_user%CALL_CLIMATE": "%s" % (self.clim_flag),
                        "cable_user%SOIL_STRUC": "'%s'" % (self.soil_flag),
                        "cable_user%vcmax": "'%s'" % (self.vcmax),
                        "cable_user%YearStart": "%d" % (st_yr_spin),
                        "cable_user%YearEnd": "%d" % (en_yr_spin),
                        "cable_user%CASA_SPIN_STARTYEAR": "%d" % (st_yr_spin),
                        "cable_user%CASA_SPIN_ENDYEAR": "%d" % (en_yr_spin),
                        "cable_user%CALL_POP": "%s" % (self.pop_flag),
                        "output%averaging": "'monthly'",
                        "icycle": "%d" % (self.biogeochem),
                        "l_vcmaxFeedbk": "%s" % (self.vcmax_feedback),
                        "l_laiFeedbk": ".TRUE.", # prognoistic LAI
                        "spinup": ".FALSE.",
                        "cable_user%FWSOIL_SWITCH": "'Haverd2013'",
                        "cable_user%GS_SWITCH": "'medlyn'",
                        "cable_user%limit_labile": ".F.",
                        "cable_user%SSNOW_POTEV": "'P-M'",
                        "cable_user%CASA_NREP": "0", # number of times to repeat CASA forcing
                        "delsoilM":"0.01",
                        "delsoilT":"0.1",
                        "output%casa": ".TRUE.",
                        "output%grid": "'land'",
                        "leaps": ".FALSE.",
                        "cable_user%litter": ".true.",
                        "cable_user%CASA_OUT_FREQ": "'monthly'",
                        "cable_user%MetType": "'site'",
                        "cable_user%PHENOLOGY_SWITCH": "'modis'",
        }
        adjust_nml_file(self.nml_fname, replace_dict)
        cwd = os.getcwd()
        (url, rev) = get_svn_info(cwd, self.cable_src)

        # delete local executable, copy a local copy and use that
        local_exe = "cable"
        if os.path.isfile(local_exe):
            os.remove(local_exe)
        shutil.copy(self.cable_exe, local_exe)
        self.cable_exe = local_exe

        return (url, rev)

    def setup_re_spin(self, number=None):
        """
        Adjust the CABLE namelist file with the various flags for another spin
        """
        out_log_fname = "%s_log_ccp%d.txt" % (self.experiment_id, number)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        self.out_fname = "%s_out_cable_ccp%d.nc" % (self.experiment_id, number)
        self.out_fname = os.path.join(self.output_dir, self.out_fname)
        if os.path.isfile(self.out_fname):
            os.remove(self.out_fname)

        self.out_fname_CASA = "%s_out_CASA_ccp%d.nc" % \
                                (self.experiment_id, number)
        self.out_fname_CASA = os.path.join(self.output_dir, self.out_fname_CASA)
        if os.path.isfile(self.out_fname_CASA):
            os.remove(self.out_fname_CASA)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_in": "'%s'" % (self.cable_restart_fname),
                        "cable_user%climate_restart_in": "'%s'" % (self.climate_restart_fname),
                        "cable_user%POP_restart_in": "'%s'" % (self.pop_restart_fname),
                        "casafile%cnpipool": "'%s'" % (self.casa_restart_fname),
                        "cable_user%POP_fromZero": ".F.",
                        "cable_user%CASA_fromZero": ".F.",
                        "cable_user%POP_rst": "'./'",
                        "cable_user%CLIMATE_fromZero": ".F.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".TRUE.",
                        "cable_user%CASA_NREP": "0",
                        "cable_user%SOIL_STRUC": "'%s'" % (self.soil_flag),
                        "icycle": "%d" % (self.biogeochem),
                        "leaps": ".TRUE.",
                        "spincasa": ".FALSE.",
                        #"casafile%c2cdumppath": "'./'",
                        "output%restart": ".TRUE.",
                        "filename%out": "'%s'" % (self.out_fname),
                        "casafile%out": "'%s'" % (self.out_fname_CASA),
        }
        adjust_nml_file(self.nml_fname, replace_dict)

    def setup_analytical_spin(self, number, st_yr_spin, en_yr_spin):
        """
        Adjust the CABLE namelist file with the various flags for the
        analytical spin step
        """
        out_log_fname = "%s_log_analytic_%d.txt" % (self.experiment_id, number)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        self.out_fname_CASA = "%s_out_CASA_analytic_%d.nc" % \
                            (self.experiment_id, number)
        self.out_fname_CASA = os.path.join(self.output_dir, self.out_fname_CASA)
        if os.path.isfile(self.out_fname_CASA):
            os.remove(self.out_fname_CASA)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "icycle": "%d" % (self.biogeochem + 10), # Need to add 10 for spinup
                        "cable_user%CASA_DUMP_READ": ".TRUE.",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "cable_user%CASA_NREP": "1",
                        "cable_user%SOIL_STRUC": "'default'", # THIS needs to be default, we don't turn sli on here, presumably because it is too slow but does this make sense?
                        "leaps": ".FALSE.",
                        "spincasa": ".TRUE.",
                        #"casafile%c2cdumppath": "'./'",
                        "cable_user%CASA_SPIN_STARTYEAR": "%d" % (st_yr_spin),
                        "cable_user%CASA_SPIN_ENDYEAR": "%d" % (en_yr_spin),
                        "casafile%out": "'%s'" % (self.out_fname_CASA),
        }
        adjust_nml_file(self.nml_fname, replace_dict)

    def setup_transient(self, st_yr_trans, en_yr_trans, st_yr, en_yr):
        """
        Adjust the CABLE namelist file for the transient run, i.e. 1850 to XXXX
        """
        replace_dict = {
                        "RunType": '"transient"',
                        "CO2NDepFile": "'%s'" % (self.co2_ndep_fname),
                        "spinstartyear": "%d" % (st_yr),
                        "spinendyear": "%d" % (en_yr),
           }
        adjust_nml_file(self.site_nml_fname, replace_dict)

        out_log_fname = "%s_log_transient.txt" % (self.experiment_id)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        self.out_fname = "%s_out_cable_transient.nc" % (self.experiment_id)
        self.out_fname = os.path.join(self.output_dir, self.out_fname)
        if os.path.isfile(self.out_fname):
            os.remove(self.out_fname)

        self.out_fname_CASA = "%s_out_casa_transient.nc" % (self.experiment_id)
        self.out_fname_CASA = os.path.join(self.output_dir, self.out_fname_CASA)
        if os.path.isfile(self.out_fname_CASA):
            os.remove(self.out_fname_CASA)

        replace_dict = {
                        "filename%out": "'%s'" % (self.out_fname),
                        "casafile%out": "'%s'" % (self.out_fname_CASA),
                        "filename%log": "'%s'" % (out_log_fname),
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "cable_user%CASA_NREP": "0",
                        "cable_user%SOIL_STRUC": "'%s'" % (self.soil_flag),
                        "output%restart": ".TRUE.",
                        "output%averaging": "'monthly'",
                        "spinup": ".FALSE.",
                        "icycle": "%d" % (self.biogeochem),
                        "cable_user%POPLUC": ".F.",
                        "cable_user%YearStart": "%d" % (st_yr_trans),
                        "cable_user%YearEnd": "%d" % (en_yr_trans),
        }
        adjust_nml_file(self.nml_fname, replace_dict)

    def setup_simulation(self, st_yr, en_yr):
        """
        Adjust the CABLE namelist file for the experiment years
        """
        replace_dict = {
                        "RunType": '"historical"',
                        "CO2NDepFile": "'%s'" % (self.co2_ndep_fname),
                        "spinstartyear": "%d" % (st_yr),
                        "spinendyear": "%d" % (en_yr)
        }
        adjust_nml_file(self.site_nml_fname, replace_dict)

        out_log_fname = "%s_log_simulation.txt" % (self.experiment_id)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        self.out_fname = "%s_out_cable.nc" % (self.experiment_id)
        self.out_fname = os.path.join(self.output_dir, self.out_fname)

        self.out_fname_CASA = "%s_out_casa.nc" % (self.experiment_id)
        self.out_fname_CASA = os.path.join(self.output_dir, self.out_fname_CASA)

        if os.path.isfile(self.out_fname):
            os.remove(self.out_fname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "output%averaging": "'daily'",
                        "icycle": "%d" % (self.biogeochem),
                        "cable_user%YearStart": "%d" % (st_yr),
                        "cable_user%YearEnd": "%d" % (en_yr),
                        "filename%out": "'%s'" % (self.out_fname),
                        "cable_user%POPLUC": ".F.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "cable_user%SOIL_STRUC": "'%s'" % (self.soil_flag),
                        "spincasa": ".FALSE.",
                        "spinup": ".FALSE.",
                        "output%averaging": "'all'",
                        "casafile%out": "'%s'" % (self.out_fname_CASA),
        }
        adjust_nml_file(self.nml_fname, replace_dict)

    def run_me(self):

        # run the model
        if self.verbose:
            cmd = './%s' % (self.cable_exe)
            error = subprocess.call(cmd, shell=True)
            if error is 1:
                print("Error running CABLE")
                raise
        else:
            # No outputs to the screen: stout and stderr to dev/null
            cmd = './%s > /dev/null 2>&1' % (self.cable_exe)
            error = subprocess.call(cmd, shell=True)
            if error is 1:
                print("Error running CABLE")


    def clean_up(self, url, rev, end=True, tag=None):
        """
        Move restart files to a directory and delete various files we no longer
        need that CABLE spits out as it spins up.
        """

        add_attributes_to_output_file(self.nml_fname, self.out_fname,
                                      url, rev)
        add_attributes_to_output_file(self.nml_fname, self.out_fname_CASA,
                                      url, rev)

        if end:
            for f in glob.glob("c2c_*_dump.nc"):
                shutil.move(f, os.path.join(self.dump_dir, f))
            f = "cnpfluxOut.csv"
            if os.path.isfile(f):
                os.remove(f)
            f = "new_sumbal"
            if os.path.isfile(f):
                os.remove(f)
            #for f in glob.glob("*.out"):
            #    os.remove(f)
            for f in glob.glob("restart_*.nc"):
                os.remove(f)
        else:
            old = self.cable_restart_fname
            new = "%s_%s.nc" % (old[:-3], tag)
            if os.path.isfile(old):
                shutil.copyfile(old, new)

            old = self.casa_restart_fname
            new = "%s_%s.nc" % (old[:-3], tag)
            if os.path.isfile(old):
                shutil.copyfile(old, new)

            old = self.climate_restart_fname
            new = "%s_%s.nc" % (old[:-3], tag)
            if os.path.isfile(old):
                shutil.copyfile(old, new)

            if self.use_pop:
                old = self.pop_restart_fname
                new = "%s_%s.nc" % (old[:-3], tag)
                if os.path.isfile(old):
                    shutil.copyfile(old, new)

if __name__ == "__main__":

    #------------- Change stuff ------------- #
    site = "Cumberland"
    met_dir = "met"
    dump_dir = "dump"
    driver_dir = "driver_files"
    log_dir = "logs"
    output_dir = "outputs"
    co2_ndep_dir = "met"
    restart_dir = "restart_files"
    aux_dir = "../../src/CMIP6-MOSRS_CNP/CABLE-AUX/"
    #aux_dir = "../../src/NESP2pt9_TRENDYv7//CABLE-AUX/"
    met_fname = "AU_Cum_2014_2017_met.nc"
    cable_src = "../../src/CMIP6-MOSRS_CNP/CMIP6-MOSRS_CNP"
    #cable_src = "../../src/NESP2pt9_TRENDYv7/NESP2pt9_TRENDYv7"
    use_pop = False
    verbose = True
    use_sli = False
    use_clim = False
    # ------------------------------------------- #

    for biogeochem in ["C", "CN", "CNP"]:
    #for biogeochem in ["C"]:

        experiment_id = "%s_%s" % (site, biogeochem)
        C = RunCable(experiment_id=experiment_id, met_dir=met_dir,
                     dump_dir=dump_dir, driver_dir=driver_dir, log_dir=log_dir,
                     output_dir=output_dir, co2_ndep_dir=co2_ndep_dir,
                     restart_dir=restart_dir, aux_dir=aux_dir,
                     met_fname=met_fname, cable_src=cable_src,
                     biogeochem=biogeochem, use_pop=use_pop, use_sli=use_sli,
                     use_clim=use_clim, verbose=verbose)
        C.main(SPIN_UP=True, TRANSIENT=True, SIMULATION=True)
