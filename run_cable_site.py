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
import tempfile
from adjust_namelist_files import adjust_nml_file
from add_missing_options_to_nml import add_missing_options_to_nml_file

class RunCable(object):

    def __init__(self, met_files, log_dir, output_dir, restart_dir, aux_dir,
                 nml_fname, veg_fname, soil_fname, grid_fname, phen_fname,
                 cnpbiome_fname, co2_conc, cable_exe, verbose):

        self.met_files = met_files
        self.log_dir = log_dir
        self.output_dir = output_dir
        self.restart_dir = restart_dir
        self.aux_dir = aux_dir
        self.cnp_biome_dir = cnp_biome_dir
        self.nml_fname = nml_fname
        self.veg_fname = veg_fname
        self.soil_fname = soil_fname
        self.grid_fname = grid_fname
        self.phen_fname = phen_fname
        self.cnpbiome_fname = cnpbiome_fname
        self.co2_conc = co2_conc
        self.cable_exe = cable_exe
        self.verbose = verbose
        self.biogeophys_dir = os.path.join(self.aux_dir, "core/biogeophys")
        self.grid_dir = os.path.join(self.aux_dir, "offline")
        self.biogeochem_dir = os.path.join(self.aux_dir, "core/biogeochem/")


    def main(self):

        for fname in self.met_files:
            site = os.path.basename(fname).split(".")[0]
            (out_fname, out_log_fname) = self.clean_up_old_files(site)

            replace_dict = {
                            "filename%met": "'%s'" % (fname),
                            "filename%out": "'%s'" % (out_fname),
                            "filename%log": "'%s'" % (out_log_fname),
                            "filename%restart_out": "' '",
                            "filename%type": "'%s'" % (os.path.join(self.grid_dir, self.grid_fname)),
                            "filename%veg": "'%s'" % (os.path.join(self.biogeophys_dir, self.veg_fname)),
                            "filename%soil": "'%s'" % (os.path.join(self.biogeophys_dir, self.soil_fname)),
                            "output%restart": ".FALSE.",
                            "fixedCO2": "%f" % (self.co2_conc),
                            "casafile%phen": "'%s'" % (os.path.join(self.biogeochem_dir, self.phen_fname)),
                            "casafile%cnpbiome": "'%s'" % (os.path.join(self.biogeochem_dir, self.cnpbiome_fname)),
            }
            adjust_nml_file(self.nml_fname, replace_dict)
            self.run_me()

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
            os.system("%s 1>&2" % (self.cable_exe))



if __name__ == "__main__":

    cwd = os.getcwd()

    #------------- Change stuff ------------- #
    met_dir = "../../met_data/plumber_met/"
    log_dir = "logs"
    output_dir = "outputs"
    restart_dir = "restart_files"
    aux_dir = "../../src/trunk/CABLE-AUX/"
    nml_fname = "cable.nml"
    veg_fname = "def_veg_params_zr_clitt_albedo_fix.txt"
    soil_fname = "def_soil_params.txt"
    grid_fname = "gridinfo_CSIRO_1x1.nc"
    phen_fname = "modis_phenology_csiro.txt"
    cnpbiome_fname = "pftlookup_csiro_v16_17tiles.csv"
    co2_conc = 380.0
    cable_exe = "../../src/trunk/CABLE_trunk/offline/cable"
    verbose = True
    all_met_files = False
    subset = ['TumbaFluxnet.1.4_met.nc']
    # ------------------------------------------- #

    if not os.path.exists(restart_dir):
        os.makedirs(restart_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    base_nml_fn = os.path.join(aux_dir, "offline/%s" % (nml_fname))
    shutil.copy(base_nml_fn, nml_fn)
    add_missing_options_to_nml_file(nml_fn)

    if all_met_files:
        met_files = glob.glob(os.path.join(met_dir, "*.nc"))
    else:
        met_files = [os.path.join(met_dir, i) for i in subset]

    C = RunCable(met_files, log_dir, output_dir, restart_dir, aux_dir,
                 nml_fname, veg_fname, soil_fname, grid_fname, phen_fname,
                 cnpbiome_fname, co2_conc, cable_exe, verbose)
    C.main()
