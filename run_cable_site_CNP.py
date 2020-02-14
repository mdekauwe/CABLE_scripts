#!/usr/bin/env python

"""
Run CABLE casa (CNP) for either a single site, a subset of sites, or all the
flux sites in the met directory.

Steps involve:

1. Spin up carbon plant & soil pools to equilibrium at pre-industrial (1850).
2. There is then a transient simulation from 1850 to the start of the FLUXNET
  file using time varying (annual) CO2, Ndep and Pdep.
3. Run the simulation with time varying CO2, Ndep and Pdep.

The script automatically generates the required met forcing(i.e. with embedded
CO2, Ndep/Pdep.

To spin the model up to equilibrium there are 4 steps:

1. An initial spin to generate the first restart file and setup the nml file.
   This step involves initialising from the PFT-level plant & soil C pools
   given in the "cnpbiome_fname" file.
2. At this point the code checks to see that the plant C pools are stabilised,
   if not, cable is re-run until the plant pools have stabilised.
3. CABLE then attempts to use the analytical solution following Xia et al. 2013
   GCB to solve the steady-state soil and litter pools. While speeding things
   up, this does not lead to an instant solution and so CABLE is then re-run
   until the soil pools have stabilised. Currently I have set the check based on
   the total soil C pools, but this means that the passive pools is not quite
   at steady-state from initial testing. However, it is likely to be good enough
   for most applications. If you wish to ensure the passive pools is stabilised
   simply set "check_passive" in the steady state function call.
4. CABLE then repeats step 4, but allowing the labile P and mineral N pools to
   freely vary.

If you've already run the C version, you can speed things up by using the
restart file from this as the initial condition for the CN run. You can of course
do likewise for a CNP run, using the CN restart file. You simply need to set
"dont_have_restart" to be False and supply the restart file number.

NB. Set mpi = True if doing a number of flux sites (*needs checking...)

ps. my suggestion is that you use this script as an experimental record, so
embed flags you want either directly into the "replace_dict" dict...or better
still, pass them via the "sci_config" dict.


That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (05.02.2020)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import shutil
import subprocess
import multiprocessing as mp
import numpy as np
import xarray as xr
import pandas as pd

from cable_utils import adjust_nml_file
from cable_utils import get_svn_info
from cable_utils import change_LAI
from cable_utils import add_attributes_to_output_file
from cable_utils import check_steady_state
from generate_cable_met_files import GenerateMetFiles


class RunCable(object):

    def __init__(self, met_dir=None, log_dir=None, output_dir=None,
                 restart_dir=None, aux_dir=None, namelist_dir=None,
                 nml_fname="cable.nml", dump_dir=None,
                 veg_fname="def_veg_params_zr_clitt_albedo_fix.txt",
                 soil_fname="def_soil_params.txt",
                 grid_fname="gridinfo_CSIRO_1x1.nc",
                 phen_fname="modis_phenology_csiro.txt",
                 #cnpbiome_fname="pftlookup_csiro_v16_17tiles_Ticket2.csv",
                 cnpbiome_fname="pftlookup_csiro_v16_17tiles-cheng-m02.csv",
                 elev_fname="GSWP3_gwmodel_parameters.nc",
                 biogeochem="C", co2_ndep_dir=None,
                 lai_dir=None, fixed_lai=None, co2_conc=400.0, co2_fixed=284.7,
                 ndep_fixed=0.79, pdep_fixed=0.144,met_subset=[],
                 cable_src=None, cable_exe="cable", mpi=True,
                 num_cores=None, verbose=True):

        self.met_dir = met_dir
        self.dump_dir = dump_dir
        self.log_dir = log_dir
        self.output_dir = output_dir
        self.restart_dir = restart_dir
        self.aux_dir = aux_dir
        self.namelist_dir = namelist_dir
        self.co2_ndep_dir = co2_ndep_dir
        self.nml_fname = nml_fname
        self.biogeophys_dir = os.path.join(self.aux_dir, "core/biogeophys")
        self.grid_dir = os.path.join(self.aux_dir, "offline")
        self.biogeochem_dir = os.path.join(self.aux_dir, "core/biogeochem/")
        self.veg_fname = os.path.join(self.biogeophys_dir, veg_fname)
        self.soil_fname = os.path.join(self.biogeophys_dir, soil_fname)
        self.grid_fname = os.path.join(self.grid_dir, grid_fname)
        self.phen_fname = os.path.join(self.biogeochem_dir, phen_fname)
        #self.cnpbiome_fname = os.path.join(self.biogeochem_dir, cnpbiome_fname)
        self.cnpbiome_fname = os.path.join("CNP_param_file", cnpbiome_fname)
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
        self.biogeochem_cyc = biogeochem
        self.co2_fixed = co2_fixed    # umol mol-1
        self.ndep_fixed = ndep_fixed  # kg N ha-1 yr-1
        self.pdep_fixed = pdep_fixed  # kg N ha-1 yr-1

        if self.biogeochem_cyc == "C":
            self.biogeochem_id = 1
            self.vcmax_feedback = ".FALSE."
        elif self.biogeochem_cyc == "CN":
            self.biogeochem_id = 2
            self.vcmax_feedback = ".TRUE."
        elif self.biogeochem_cyc == "CNP":
            self.biogeochem_id = 3
            self.vcmax_feedback = ".TRUE."

    def main(self, sci_config, dont_have_restart=True, restart_num=None):

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
                               args=(met_files[start:end], url, rev,
                                     sci_config, dont_have_restart,
                                     restart_num, ))
                processes.append(p)

            # Run processes
            for p in processes:
                p.start()
        else:
            self.worker(met_files, url, rev, sci_config, dont_have_restart,
                        restart_num,)

    def worker(self, met_files, url, rev, sci_config, dont_have_restart,
               restart_num):

        num = 0

        for fname in met_files:

            site = os.path.basename(fname).split(".")[0].split("_")[0]
            self.experiment_id = "%s_%s" % (site, self.biogeochem_cyc)
            print("\n%s\n" % (site))

            (st_yr, en_yr) = self.get_met_years(fname)

            (fname_spin,
             fname_trans,
             fname_sim) = self.generate_met_files(site, st_yr, en_yr, fname)

            self.setup_nml_file(site)
            self.inital_spin_setup(site, fname_spin, sci_config, num)

            # First phase: spin up with static CO2, Ndep, Pdep. Here we run the
            # model once to read in the initial plant and soil C pools from the
            # PFT level file and setup the restart file.
            if dont_have_restart:

                # Create initial restart file
                print("\n===============================================\n")
                print("Generate initial restart file\n")
                print("===============================================\n\n")
                self.run_me()
                self.clean_up(num, tag="spin")
                num += 1
            else:
                self.setup_inital_restart_file(restart_num, site)
                # logic is based on that final increment above that the user
                # won't do...
                num = restart_num + 1




            # Second phase: find steady-state NPP
            not_stabilised = True
            while not_stabilised:

                print("\n===============================================\n")
                print("Find steady-state NPP: %d\n" % (num))
                print("===============================================\n\n")
                self.setup_spin(number=num)
                self.run_me()
                self.clean_up(num, tag="spin")

                not_stabilised = check_steady_state(self.experiment_id,
                                                    self.restart_dir,
                                                    self.output_dir,
                                                    num, check_npp=True,
                                                    debug=True)
                num += 1


            # Third phase: bring plant biomass pools into equilibrium
            not_stabilised = True
            while not_stabilised:

                print("\n===============================================\n")
                print("Bring plant C pools into equilibrium: %d\n" % (num))
                print("===============================================\n\n")
                self.setup_spin(number=num)
                self.run_me()
                self.clean_up(num, tag="spin")

                not_stabilised = check_steady_state(self.experiment_id,
                                                    self.restart_dir,
                                                    self.output_dir,
                                                    num, check_plant=True,
                                                    debug=True)
                num += 1


            # Fourht phase: bring soil pools into equilibrium using analytical
            # solution
            not_stabilised = True
            while not_stabilised:

                print("\n===================================================\n")
                print("Bring soil C pools into equilibrium: %d\n" % (num))
                print("===================================================\n\n")
                self.setup_analytical_spin(st_yr, en_yr, number=num)
                self.run_me()
                self.clean_up(num, tag="spin_analytic")

                not_stabilised = check_steady_state(self.experiment_id,
                                                    self.restart_dir,
                                                    self.output_dir,
                                                    num, check_soil=True,
                                                    debug=True)

                # Not found a steady-state solution run further simulations...
                for i in range(10):
                    self.setup_spin(number=num)
                    self.run_me()
                    self.clean_up(num, tag="spin")

                    num += 1

            # Fourth phase: bring soil pools into equilibrium using analytical
            # solution, but without restricting N and P pools
            for i in range(3):

                # run another simulation...
                self.setup_spin(number=num, labile=False)
                self.run_me()
                self.clean_up(num, tag="spin")

                num += 1

            not_stabilised = True
            while not_stabilised:

                print("\n===================================================\n")
                print("Bring soil C pools into equilibrium: \n")
                print("Urestricted labile P/mineral N: %d\n" % (num))
                print("===================================================\n\n")
                self.setup_analytical_spin(st_yr, en_yr, number=num)
                self.run_me()
                self.clean_up(num, tag="spin_analytic")

                not_stabilised = check_steady_state(self.experiment_id,
                                                    self.restart_dir,
                                                    self.output_dir,
                                                    num, check_soil=True,
                                                    debug=True)

                # Not found a steady-state solution run further simulations...
                for i in range(5):
                    self.setup_spin(number=num)
                    self.run_me()
                    self.clean_up(num, tag="spin")

                    num += 1

                num += 1

            # Run transient simulation from 1850
            print("\n===================================================\n")
            print("Historical\n")
            print("===================================================\n\n")
            (out_fname) = self.setup_simulation(fname_trans, historical=True,
                                                number=num-1)
            self.run_me()
            self.clean_up(number=None, tag="historical")

            print("\n===================================================\n")
            print("Simulation\n")
            print("===================================================\n\n")
            (out_fname) = self.setup_simulation(fname_sim, historical=False,
                                                number=num-1)
            self.run_me()
            self.clean_up(number=None, tag="simulation")

            add_attributes_to_output_file(self.nml_fname, out_fname, url, rev)
            ofname = os.path.join(self.namelist_dir, self.nml_fname)
            shutil.move(self.nml_fname, ofname)

    def setup_nml_file(self, site):
        # get a clean namelist file

        base_nml_fn = os.path.join(self.grid_dir, "%s" % (self.nml_fname))
        nml_fname = "cable_%s.nml" % (site)
        shutil.copy(base_nml_fn, nml_fname)
        self.nml_fname = nml_fname

    def generate_met_files(self, site, st_yr, en_yr, met_fname):

        local_met_dir = "met"
        if not os.path.exists(local_met_dir):
            os.makedirs(local_met_dir)

        co2_ndep_fname = "AmaFACE_co2npdepforcing_1850_2100_AMB.csv"
        co2_ndep_fname = os.path.join(self.co2_ndep_dir, co2_ndep_fname)

        G = GenerateMetFiles(site, met_fname, co2_ndep_fname)

        fname_spin = os.path.join(local_met_dir, "%s_met_spin.nc" % (site))
        if not os.path.isfile(fname_spin):
            G.create_spin_file(fname_spin, self.co2_fixed, self.ndep_fixed,
                               self.pdep_fixed)

        fname_trans = os.path.join(local_met_dir, "%s_met_trans.nc" % (site))
        if not os.path.isfile(fname_trans):
            G.create_transient_file(fname_trans)

        fname_sim = "%s_met_simulation.nc" % (site)
        fname_sim = os.path.join(local_met_dir, fname_sim)
        if not os.path.isfile(fname_sim):
            G.create_simulation_file(fname_sim)

        return (fname_spin, fname_trans, fname_sim)

    def setup_inital_restart_file(self, number, site):

        if self.biogeochem_cyc == "CN":
            old_experiment_id = "%s_%s" % (site, "C")
        elif self.biogeochem_cyc == "CNP":
            old_experiment_id = "%s_%s" % (site, "CN")

        cable_rst_ofname = "%s_cable_rst_%d.nc" % (old_experiment_id, number)
        cable_rst_ofname = os.path.join(self.restart_dir, cable_rst_ofname)
        new_cable_rst = "%s_cable_rst_%d.nc" % (self.experiment_id, number)
        new_cable_rst = os.path.join(self.restart_dir, new_cable_rst)
        shutil.copy(cable_rst_ofname, new_cable_rst)

        casa_rst_ofname = "%s_casa_rst_%d.nc" % (old_experiment_id, number)
        casa_rst_ofname = os.path.join(self.restart_dir, casa_rst_ofname)
        new_casa_rst = "%s_casa_rst_%d.nc" % (self.experiment_id, number)
        new_casa_rst = os.path.join(self.restart_dir, new_casa_rst)
        shutil.copy(casa_rst_ofname, new_casa_rst)

    def initialise_stuff(self):

        if not os.path.exists(self.restart_dir):
            os.makedirs(self.restart_dir)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        if not os.path.exists(self.namelist_dir):
            os.makedirs(self.namelist_dir)

        if not os.path.exists(self.dump_dir):
            os.makedirs(self.dump_dir)

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

    def inital_spin_setup(self, site, met_fname, sci_config, number=None):
        """
        Initial setup for CASA spinup from zero to generate restart file
        """

        # set directory paths...
        out_fname = "%s_out_cable_spin_%d.nc" % (self.experiment_id, number)
        out_fname = os.path.join(self.output_dir, out_fname)

        out_log_fname = "%s_log_spin_%d.txt" % (self.experiment_id, number)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)

        cable_rst_ofname = "%s_cable_rst_%d.nc" % (self.experiment_id, number)
        cable_rst_ofname = os.path.join(self.restart_dir, cable_rst_ofname)

        casa_rst_ofname = "%s_casa_rst_%d.nc" % (self.experiment_id, number)
        casa_rst_ofname = os.path.join(self.restart_dir, casa_rst_ofname)

        replace_dict = {
                        "filename%met": "'%s'" % (met_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_in": "' '",
                        "filename%restart_out": "'%s'" % (cable_rst_ofname),
                        "casafile%cnpipool": "' '",
                        "casafile%cnpepool": "'%s'" % (casa_rst_ofname),
                        "filename%type": "'%s'" % (self.grid_fname),
                        "filename%veg": "'%s'" % (self.veg_fname),
                        "filename%soil": "'%s'" % (self.soil_fname),
                        "fixedCO2": "%.2f" % (self.co2_conc),
                        "casafile%phen": "'%s'" % (self.phen_fname),
                        "casafile%cnpbiome": "'%s'" % (self.cnpbiome_fname),
                        "cable_user%RunIden": "'%s'" % (self.experiment_id),
                        "cable_user%vcmax": "'standard'",
                        "l_vcmaxFeedbk": "%s" % (self.vcmax_feedback),
                        "l_laiFeedbk": ".TRUE.", # prognoistic LAI
                        "icycle": "%d" % (self.biogeochem_id),
                        "cable_user%CASA_OUT_FREQ": "'annually'",
                        "output%casa": ".TRUE.",
                        "leaps": ".TRUE.",
                        "cable_user%CASA_fromZero": ".TRUE.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".TRUE.",
                        "cable_user%CASA_NREP": "0",
                        "spinup": ".FALSE.",
                        "spincasa": ".FALSE.",
                        "output%restart": ".TRUE.",
                        "cable_user%FWSOIL_SWITCH": "'standard'",
                        "cable_user%GS_SWITCH": "'medlyn'",
                        "cable_user%limit_labile": ".TRUE.",
        }
        # Make sure the dict isn't empty
        if bool(sci_config):
            replace_dict = merge_two_dicts(replace_dict, sci_config)

        adjust_nml_file(self.nml_fname, replace_dict)

    def setup_spin(self, number=None, labile=True):
        """
        Adjust the CABLE namelist file with the various flags for another spin
        """

        out_fname = "%s_out_cable_spin_%d.nc" % (self.experiment_id, number)
        out_fname = os.path.join(self.output_dir, out_fname)

        out_log_fname = "%s_log_spin_%d.txt" % (self.experiment_id, number)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)

        cable_rst_ifname = "%s_cable_rst_%d.nc" % (self.experiment_id, number-1)
        cable_rst_ifname = os.path.join(self.restart_dir, cable_rst_ifname)

        cable_rst_ofname = "%s_cable_rst_%d.nc" % (self.experiment_id, number)
        cable_rst_ofname = os.path.join(self.restart_dir, cable_rst_ofname)

        casa_rst_ifname = "%s_casa_rst_%d.nc" % (self.experiment_id, number-1)
        casa_rst_ifname = os.path.join(self.restart_dir, casa_rst_ifname)

        casa_rst_ofname = "%s_casa_rst_%d.nc" % (self.experiment_id, number)
        casa_rst_ofname = os.path.join(self.restart_dir, casa_rst_ofname)

        # Restrict build up of labile P and mineral N pools during spinup
        if labile:
            restrict_labile = ".TRUE."
        else:
            restrict_labile = ".FALSE."

        replace_dict = {
                        "filename%out": "'%s'" % (out_fname),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_in": "'%s'" % (cable_rst_ifname),
                        "filename%restart_out": "'%s'" % (cable_rst_ofname),
                        "casafile%cnpipool": "'%s'" % (casa_rst_ifname),
                        "casafile%cnpepool": "'%s'" % (casa_rst_ofname),
                        "cable_user%CASA_fromZero": ".FALSE.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".TRUE.",
                        "cable_user%CASA_NREP": "0",
                        "output%restart": ".TRUE.",
                        "spincasa": ".FALSE.",
                        "icycle": "%d" % (self.biogeochem_id),
                        "leaps": ".TRUE.",
                        "cable_user%limit_labile": "%s" % (restrict_labile),
        }
        adjust_nml_file(self.nml_fname, replace_dict)

    def setup_analytical_spin(self, st_yr, en_yr, labile=True, number=None):
        """
        Adjust the CABLE namelist file with the various flags for the
        analytical spin step (to stabilise the soil C pools)
        """

        out_fname = "%s_out_cable_aspin_%d.nc" % (self.experiment_id, number)
        out_fname = os.path.join(self.output_dir, out_fname)

        out_log_fname = "%s_log_aspin_%d.txt" % (self.experiment_id, number)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)

        # This step doesn't write an output CABLE rst file.
        cable_rst_ifname = "%s_cable_rst_%d.nc" % (self.experiment_id, number-1)
        cable_rst_ifname = os.path.join(self.restart_dir, cable_rst_ifname)

        casa_rst_ifname = "%s_casa_rst_%d.nc" % (self.experiment_id, number-1)
        casa_rst_ifname = os.path.join(self.restart_dir, casa_rst_ifname)

        # This will overwrite the current restart file for CASA but that is
        # fine.
        casa_rst_ofname = "%s_casa_rst_%d.nc" % (self.experiment_id, number)
        casa_rst_ofname = os.path.join(self.restart_dir, casa_rst_ofname)

        replace_dict = {
                        "filename%out": "'%s'" % (out_fname),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_in": "'%s'" % (cable_rst_ifname),
                        "casafile%cnpipool": "'%s'" % (casa_rst_ifname),
                        "casafile%cnpepool": "'%s'" % (casa_rst_ofname),
                        "icycle": "%d" % (self.biogeochem_id + 10), # Need to add 10 for spinup
                        "cable_user%CASA_DUMP_READ": ".TRUE.",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "cable_user%CASA_NREP": "1",
                        "spincasa": ".TRUE.",
                        "cable_user%CASA_SPIN_STARTYEAR": "%d" % (st_yr),
                        "cable_user%CASA_SPIN_ENDYEAR": "%d" % (en_yr),
                        "output%restart": ".TRUE.",
                        "leaps": ".TRUE.",

        }
        adjust_nml_file(self.nml_fname, replace_dict)


    def setup_simulation(self, met_fname, historical=True, number=None):
        """
        Adjust the CABLE namelist file for the experiment years
        """

        out_fname = "%s_out_cable_simulation.nc" % (self.experiment_id)
        out_fname = os.path.join(self.output_dir, out_fname)

        out_log_fname = "%s_log_simulation.txt" % (self.experiment_id)
        out_log_fname = os.path.join(self.log_dir, out_log_fname)

        cable_rst_ifname = "%s_cable_rst_%d.nc" % (self.experiment_id, number)
        cable_rst_ifname = os.path.join(self.restart_dir, cable_rst_ifname)

        casa_rst_ifname = "%s_casa_rst_%d.nc" % (self.experiment_id, number)
        casa_rst_ifname = os.path.join(self.restart_dir, casa_rst_ifname)

        if historical:

            out_fname = "%s_out_cable_historical.nc" % (self.experiment_id)
            out_fname = os.path.join(self.output_dir, out_fname)

            cable_rst_ofname = "%s_cable_rst_%d.nc" % (self.experiment_id, number)
            cable_rst_ofname = os.path.join(self.restart_dir, cable_rst_ofname)

            casa_rst_ofname = "%s_casa_rst_%d.nc" % (self.experiment_id, number)
            casa_rst_ofname = os.path.join(self.restart_dir, casa_rst_ofname)

            replace_dict = {
                            "filename%met": "'%s'" % (met_fname),
                            "filename%out": "'%s'" % (out_fname),
                            "filename%log": "'%s'" % (out_log_fname),
                            "filename%restart_in": "'%s'" % (cable_rst_ifname),
                            "filename%restart_out": "'%s'" % (cable_rst_ofname),
                            "casafile%cnpipool": "'%s'" % (casa_rst_ifname),
                            "casafile%cnpepool": "'%s'" % (casa_rst_ofname),
                            "icycle": "%d" % (self.biogeochem_id),
                            "cable_user%CASA_DUMP_READ": ".FALSE.",
                            "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                            "spincasa": ".FALSE.",
                            "output%restart": ".TRUE.",
                            "output%averaging": "'daily'",
                            "leaps": ".TRUE.",
                            "cable_user%limit_labile": ".FALSE.",
            }
        else:
            replace_dict = {
                        "filename%met": "'%s'" % (met_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_in": "'%s'" % (cable_rst_ifname),
                        "filename%restart_out": "' '",
                        "casafile%cnpipool": "'%s'" % (casa_rst_ifname),
                        "casafile%cnpepool": "' '",
                        "icycle": "%d" % (self.biogeochem_id),
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "spincasa": ".FALSE.",
                        "output%restart": ".FALSE.",
                        "output%averaging": "'daily'",
                        "leaps": ".TRUE.",
                        "cable_user%limit_labile": ".FALSE.",

        }
        adjust_nml_file(self.nml_fname, replace_dict)

        return (out_fname)

    def clean_up(self, number, tag):

        # CASA out file is hardwired!
        fname = glob.glob("*_casa_out.nc")
        if len(fname) > 0:
            if number is None:
                out_fname = "%s_out_casa_%s.nc" % \
                                    (self.experiment_id, tag)
            else:
                out_fname = "%s_out_casa_%s_%d.nc" % \
                                    (self.experiment_id, tag, number)
            shutil.move(fname[0], os.path.join(self.output_dir, out_fname))

        f = "cnpfluxOut.csv"
        if os.path.isfile(f):
            os.remove(f)

        f = "new_sumbal"
        if os.path.isfile(f):
            os.remove(f)

        if tag == "simulation":
            for f in glob.glob("c2c_*_dump.nc"):
                shutil.move(f, os.path.join(self.dump_dir, f))

    def run_me(self):
        # run the model
        if self.verbose:
            cmd = './%s %s' % (self.cable_exe, self.nml_fname)
            error = subprocess.call(cmd, shell=True)
            if error is 1:
                print("Job failed to submit")
                raise
        else:
            # No outputs to the screen: stout and stderr to dev/null
            cmd = './%s %s > /dev/null 2>&1' % (self.cable_exe, self.nml_fname)
            error = subprocess.call(cmd, shell=True)
            if error is 1:
                print("Job failed to submit")

    def get_met_years(self, met_fname):

        ds = xr.open_dataset(met_fname)
        st_yr = pd.to_datetime(ds.time[0].values).year
        en_yr = pd.to_datetime(ds.time[-1].values).year #- 1

        return (st_yr, en_yr)

def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

if __name__ == "__main__":

    #------------- Change stuff ------------- #
    met_dir = "met"
    co2_ndep_dir = "co2_ndep"
    dump_dir = "dump"
    log_dir = "logs"
    output_dir = "outputs"
    restart_dir = "restarts"
    namelist_dir = "namelists"
    aux_dir = "../../src/CABLE-AUX/"
    #cable_src = "../../src/trunk/trunk"
    cable_src = "../../src/trunk_analytical/trunk_analytical/"
    mpi = False
    num_cores = 1 # set to a number, if None it will use all cores...!
    # if empty...run all the files in the met_dir
    met_subset = ['AU-Tum_2002-2016_OzFlux_Met.nc']
    sci_config = {}
    #biogeochem = "CNP"
    # ------------------------------------------- #

    """
    dont_have_restart = True
    for biogeochem in ["C", "CN", "CNP"]:
        C = RunCable(met_dir=met_dir, log_dir=log_dir, output_dir=output_dir,
                     dump_dir=dump_dir, restart_dir=restart_dir,
                     aux_dir=aux_dir, namelist_dir=namelist_dir,
                     met_subset=met_subset, cable_src=cable_src, mpi=mpi,
                     num_cores=num_cores, biogeochem=biogeochem,
                     co2_ndep_dir=co2_ndep_dir)
        C.main(sci_config)
    """

    #"""
    dont_have_restart = True
    for biogeochem in ["C"]:
        C = RunCable(met_dir=met_dir, log_dir=log_dir, output_dir=output_dir,
                     dump_dir=dump_dir, restart_dir=restart_dir,
                     aux_dir=aux_dir, namelist_dir=namelist_dir,
                     met_subset=met_subset, cable_src=cable_src, mpi=mpi,
                     num_cores=num_cores, biogeochem=biogeochem,
                     co2_ndep_dir=co2_ndep_dir)
        C.main(sci_config)
    #"""

    """
    dont_have_restart = False
    num = 511
    for biogeochem in ["CN"]:
        C = RunCable(met_dir=met_dir, log_dir=log_dir, output_dir=output_dir,
                     dump_dir=dump_dir, restart_dir=restart_dir,
                     aux_dir=aux_dir, namelist_dir=namelist_dir,
                     met_subset=met_subset, cable_src=cable_src, mpi=mpi,
                     num_cores=num_cores, biogeochem=biogeochem,
                     co2_ndep_dir=co2_ndep_dir)
        C.main(sci_config, dont_have_restart, num)
    """
