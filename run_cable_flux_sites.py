#!/usr/bin/env python

"""
Run CABLE either for a subset(single) flux set, or all the sites in a directory.
This script also allows you to edit the namelist internal.

Todo
----
 * Pass site dictionary so that the namelist can vary by site.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (14.10.2017)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import shutil
import tempfile

class RunCable(object):

    def __init__(self, driver_dir, met_files, output_dir, nml_fn, veg_param_fn,
                 log_dir, exe, aux_dir, verbose):

        self.driver_dir = driver_dir
        self.met_files = met_files
        self.output_dir = output_dir
        self.nml_fn = nml_fn
        self.veg_param_fn = veg_param_fn
        self.log_dir = log_dir
        self.cable_exe = exe
        self.aux_dir = aux_dir
        self.verbose = verbose

    def main(self):

        for fname in self.met_files:
            site = os.path.basename(fname).split(".")[0]

            out_fname = os.path.join(self.output_dir, "%s_out.nc" % (site))
            if os.path.isfile(out_fname):
                os.remove(out_fname)

            out_log_fname = os.path.join(self.log_dir, "%s_log.nc" % (site))
            if os.path.isfile(out_log_fname):
                os.remove(out_log_fname)

            replace_dict = {
                            "filename%met": "'%s'" % (fname),
                            "filename%out": "'%s'" % (out_fname),
                            "filename%log": "'%s'" % (out_log_fname),
                            "filename%restart_out": "' '",
                            "filename%type": "'%sgridinfo_CSIRO_1x1.nc'" % (self.driver_dir),
                            "filename%veg": "'%s%s'" % (self.driver_dir, veg_param_fn),
                            "filename%soil": "'%sdef_soil_params.txt'" % (self.driver_dir),
                            "output%restart": ".FALSE.",
                            "fixedCO2": "380.0",
                            "casafile%phen": "'%s'" % (os.path.join(self.aux_dir, "core/biogeochem/modis_phenology_csiro.txt")),
                            "casafile%cnpbiome": "'%s'" % (os.path.join(self.aux_dir, "core/biogeochem/pftlookup_csiro_v16_17tiles.csv")),
            }
            self.adjust_param_file(replace_dict)

            # run the model
            if self.verbose:
                os.system("%s" % (self.cable_exe))
            else:
                os.system("%s >& %s" % (self.cable_exe, out_log_fname))

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
                lines[i] = " ".join((key, "=",
                                     replacements_dict.get(key.strip(), val)))

        return '\n'.join(lines) + '\n'

if __name__ == "__main__":

    driver_dir = "../../driver_files/"
    met_dir = "../../met_data/plumber_met/"
    output_dir = "outputs"
    base_nml_fn = os.path.join(driver_dir, "cable.nml")
    cwd = os.getcwd()
    nml_fn = os.path.join(cwd, "cable.nml")
    shutil.copy(base_nml_fn, nml_fn)
    veg_param_fn = "def_veg_params.txt"
    log_dir = "log_files"
    exe = "../../src/trunk/CABLE_trunk/offline/cable"
    #exe = "../../src/tag-2.3.4/CABLE-2.3.4_tag/offline/cable"
    aux_dir = "../../src/tag-2.3.4/CABLE-AUX/"
    verbose = True
    all_met_files = False

    if all_met_files:
        met_files = glob.glob(os.path.join(met_dir, "*.nc"))
    else:
        subset = ['TumbaFluxnet.1.4_met.nc']
        met_files = [os.path.join(met_dir, i) for i in subset]

    C = RunCable(driver_dir, met_files, output_dir, nml_fn, veg_param_fn,
                 log_dir, exe, aux_dir, verbose)
    C.main()
