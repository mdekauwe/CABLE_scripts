#!/usr/bin/env python

"""
Site run - just biophysics
==========================

Run CABLE either for a single site, a subset, or all the flux sites pointed to
in the met directory. The script creates cable.nml file.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (14.09.2017)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import shutil
import tempfile

class RunCable(object):

    def __init__(self, met_files, output_dir, log_dir, aux_dir,
                 nml_fn, veg_param_fn, soil_fn, grid_fn, cable_exe, verbose):

        self.met_files = met_files
        self.output_dir = output_dir
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

    def main(self):

        for fname in self.met_files:
            site = os.path.basename(fname).split(".")[0]
            (out_fname, out_log_fname) = self.clean_up_old_files(site)

            replace_dict = {
                            "filename%met": "'%s'" % (fname),
                            "filename%out": "'%s'" % (out_fname),
                            "filename%log": "'%s'" % (out_log_fname),
                            "filename%restart_out": "' '",
                            "filename%type": "'%s'" % (os.path.join(self.grid_dir, self.grid_fn)),
                            "filename%veg": "'%s'" % (os.path.join(self.veg_dir, self.veg_param_fn)),
                            "filename%soil": "'%s'" % (os.path.join(self.veg_dir, self.soil_fn)),
                            "output%restart": ".FALSE.",
                            "fixedCO2": "380.0",
                            "casafile%phen": "'%s'" % (os.path.join(self.aux_dir, "core/biogeochem/modis_phenology_csiro.txt")),
                            "casafile%cnpbiome": "'%s'" % (os.path.join(self.aux_dir, "core/biogeochem/pftlookup_csiro_v16_17tiles.csv")),
            }
            self.adjust_param_file(replace_dict)
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

    def replace_keys(self, text, replacements_dict):
        """ Function expects to find GDAY input file formatted key = value.

        Parameters:
        ----------
        text : string
            input file data.
        replacements_dict : dictionary
            dictionary of replacement values.

        Returns:
        --------
        new_text : string
            input file with replacement values

        """
        lines = text.splitlines()
        for i, row in enumerate(lines):
            # skip blank lines
            if not row.strip():
                continue
            if "=" not in row:
                lines[i] = row
                continue
            elif not row.startswith("&"):
                key = row.split("=")[0]
                val = row.split("=")[1]
                lines[i] = " ".join((key.rstrip(), "=",
                                     replacements_dict.get(key.strip(),
                                     val.lstrip())))

        return '\n'.join(lines) + '\n'

def add_missing_options_to_nml_file(fname, line_start=60):
    # Some of the flags we may wish to change are missin from the default
    # file so we can't adjust them via this script...add them

    f = open(fname, "r")
    contents = f.readlines()
    f.close()

    arg = "   cable_user%GS_SWITCH = 'medlyn'\n"
    contents.insert(line_start, arg)
    line_start += 1

    arg = "   cable_user%GW_MODEL = .FALSE.\n"
    contents.insert(line_start, arg)
    line_start += 1

    arg = "   cable_user%or_evap = .TRUE.\n"
    contents.insert(line_start, arg)
    line_start += 1

    tmp_fname = "tmp.nml"
    f = open(tmp_fname, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()

    shutil.move(tmp_fname, fname)

if __name__ == "__main__":

    cwd = os.getcwd()
    met_dir = "../../met_data/plumber_met/"
    log_dir = "logs"
    output_dir = "outputs"
    restart_dir = "restart_files"
    aux_dir = "../../src/trunk/CABLE-AUX/"
    cable_exe = "../../src/trunk/CABLE_trunk/offline/cable"

    if not os.path.exists(restart_dir):
        os.makedirs(restart_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    base_nml_fn = os.path.join(aux_dir, "offline/cable.nml")
    nml_fn = "cable.nml"
    shutil.copy(base_nml_fn, nml_fn)
    add_missing_options_to_nml_file(nml_fn)

    veg_fn = "def_veg_params_zr_clitt_albedo_fix.txt"
    soil_fn = "def_soil_params.txt"
    base_veg_param_fn = os.path.join(aux_dir, veg_fn)
    shutil.copy(base_nml_fn, nml_fn)
    grid_fn = "gridinfo_CSIRO_1x1.nc"
    verbose = True
    all_met_files = False

    if all_met_files:
        met_files = glob.glob(os.path.join(met_dir, "*.nc"))
    else:
        subset = ['TumbaFluxnet.1.4_met.nc']
        met_files = [os.path.join(met_dir, i) for i in subset]

    C = RunCable(met_files, output_dir, log_dir, aux_dir,
                 nml_fn, veg_fn, soil_fn, grid_fn, cable_exe, verbose)
    C.main()
