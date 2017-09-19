#!/usr/bin/env python

"""
Tumbarumba
==========

- Model spin-up: using K34 tower info, CO2=284.7; NDEP-0.79 kg N ha-1 yr-1;
                 PDEP=0.144 kg P ha-1 yr-1
- Transient: 1851-1998, varying CO2 and NDEP, but just recycling the intact
             16-year met forcing.
- CNP + POP switched on.

During the spinup, we are recycling in 30 year chunks.

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (19.09.2017)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import shutil
import tempfile

class RunCable(object):

    def __init__(self, site, driver_dir, output_dir, restart_dir, met_fname,
                 co2_ndep_fname, nml_fn, site_nml_fn, veg_param_fn, log_dir,
                 exe, aux_dir, verbose, SPIN_UP=False):

        self.site = site
        self.driver_dir = driver_dir
        self.output_dir = output_dir
        self.restart_dir = restart_dir
        self.met_fname = met_fname
        self.co2_ndep_fname = co2_ndep_fname
        self.nml_fn = nml_fn
        self.site_nml_fn = site_nml_fn
        self.veg_param_fn = veg_param_fn
        self.log_dir = log_dir
        self.cable_exe = exe
        self.aux_dir = aux_dir
        self.verbose = verbose
        self.SPIN_UP = SPIN_UP

    def main(self):

        if self.SPIN_UP == True:

            #"""
            self.setup_ini_spin()
            self.run_me()
            self.clean_up(ini=True)
            #"""
            """
            restart_fname = "%s_casa_rst.nc" % (self.site)
            self.setup_re_spin(restart_fname, number=1)
            self.run_me()
            self.clean_up(ini=False, number=1)
            """
            #restart_fname = "%s_casa_rst.nc" % (self.site)
            #self.re_spin(restart_fname, analytical_spin=True)
            #self.clean_up(ini=False, number=1)

            #"""

            """
            restart_fname = "%s_casa_rst.nc" % (self.site)
            for i,n in enumerate([1,2]):
                self.re_spin(restart_fname)
                self.clean_up(ini=False, number=n)
            """

    def adjust_nml_file(self, fname, replacements):
        """ adjust CABLE NML file and write over the original.

        Parameters:
        ----------
        fname : string
            parameter filename to be changed.
        replacements : dictionary
            dictionary of replacement values.

        """
        fin = open(fname, 'r')
        param_str = fin.read()
        fin.close()
        new_str = self.replace_keys(param_str, replacements)
        fd, path = tempfile.mkstemp()
        os.write(fd, str.encode(new_str))
        os.close(fd)
        shutil.copy(path, fname)
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

    def run_me(self):
        # run the model
        if self.verbose:
            os.system("%s" % (self.cable_exe))
        else:
            os.system("%s 1>&2" % (self.cable_exe))

    def clean_up(self, ini=True, number=None):

        # is this robust?
        tag = self.site[:-2]

        if ini:
            for f in glob.glob("*.out"):
                os.remove(f)
            os.remove("new_sumbal")
            os.remove("cnpfluxOut.csv")
            os.remove("%s_1822_1851_casa_out.nc" % (tag))

            fromx = "pop_%s_ini.nc" % (tag)
            to = "pop_%s_ini_zero.nc" % (tag)
            to = os.path.join(self.restart_dir, to)
            shutil.copyfile(fromx, to)

            fromx = "%s_climate_rst.nc" % (tag)
            to = "%s_climate_rst_zero.nc" % (tag)
            to = os.path.join(self.restart_dir, to)
            shutil.copyfile(fromx, to)

            fromx = "%s_casa_rst.nc" % (tag)
            to = "%s_casa_rst_zero.nc" % (tag)
            to = os.path.join(self.restart_dir, to)
            shutil.copyfile(fromx, to)
        else:
            fromx = "pop_%s_ini.nc" % (tag)
            to = "pop_%s_ini_ccp%d.nc" % (tag, number)
            to = os.path.join(self.restart_dir, to)
            shutil.copyfile(fromx, to)

            fromx = "%s_climate_rst.nc" % (tag)
            to = "%s_climate_rst_ccp%d.nc" % (tag, number)
            to = os.path.join(self.restart_dir, to)
            shutil.copyfile(fromx, to)

            fromx = "%s_casa_rst.nc" % (tag)
            to = "%s_casa_rst_ccp%d.nc" % (tag, number)
            to = os.path.join(self.restart_dir, to)
            shutil.copyfile(fromx, to)

            for f in glob.glob("c2c_*_dump.nc"):
                os.remove(f)

            for f in glob.glob("*_casa_out"):
                os.remove(f)

            for f in glob.glob("*.out"):
                os.remove(f)

        return

    def setup_ini_spin(self):
        shutil.copyfile(os.path.join(self.driver_dir, "site.nml"),
                        self.site_nml_fn)
        shutil.copyfile(os.path.join(self.driver_dir, "cable.nml"),
                        self.nml_fn)

        # Replace with Vanessa's starting file, this is hardwired until,
        # cable.nml reflects Vanessa's inputs
        vanessa_nml_fn = "ancillary_files/cable.nml.cable_casa_POP_from_zero"
        shutil.copyfile(vanessa_nml_fn, self.nml_fn)

        out_fname = os.path.join(self.output_dir, "%s_spin.nc" % (site))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_zero.nc" % (site))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        replace_dict = {
                        "RunType": '"spinup"',
                        "CO2NDepFile": "'%s'" % (self.co2_ndep_fname),
                        "spinstartyear": "2002",
                        "spinendyear": "2003",
                        "spinCO2": "284.7",
                        "spinNdep": "0.79",
                        "spinPdep": "0.144",
        }
        self.adjust_nml_file(self.site_nml_fn, replace_dict)

        restart_fname = "%s_restart.nc" % (site)
        replace_dict = {
                        "filename%met": "'%s'" % (self.met_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_out": "'%s'" % (restart_fname),
                        "filename%type": "'%s'" % (os.path.join(self.aux_dir, "offline/gridinfo_CSIRO_1x1.nc")),
                        "filename%veg": "'%s%s'" % (self.driver_dir, veg_param_fn),
                        "filename%soil": "'%sdef_soil_params.txt'" % (self.driver_dir),
                        "output%restart": ".TRUE.",
                        "fixedCO2": "380.0",
                        "casafile%phen": "'%s'" % (os.path.join(self.aux_dir, "core/biogeochem/modis_phenology_csiro.txt")),
                        "casafile%cnpbiome": "'%s'" % (os.path.join(self.driver_dir, "pftlookup_csiro_v16_17tiles_Cumberland.csv")),
                        "cable_user%RunIden": "'%s'" % (self.site),
                        "cable_user%POP_out": "'rst'",
                        "cable_user%POP_rst": "'./'",
                        "cable_user%POP_fromZero": ".T.",
                        "cable_user%CASA_fromZero": ".T.",

        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_re_spin(self, restart_fname, number=None):

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_ccp_%d.nc" % (site, number))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_in": "'%s'" % (restart_fname),
                        "filename%restart_out": "'%s'" % (restart_fname),
                        "cable_user%POP_fromZero": ".F.",
                        "cable_user%CASA_fromZero": ".F.",
                        "cable_user%CLIMATE_fromZero": ".F.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".TRUE.",
                        "cable_user%CASA_NREP": "0",
                        "cable_user%SOIL_STRUC": "'sli'",
                        "icycle": "2",
                        "cable_user%POP_out": "'rst'",
                        "leaps": ".TRUE.",
                        "spincasa": ".FALSE.",
                        "casafile%cnpipool": "' '",
                        "casafile%c2cdumppath": "' '",
                        #"cable_user%CASA_SPIN_STARTYEAR": "'%s'" % (str(start_yr)),
                        #"cable_user%CASA_SPIN_ENDYEAR": "'%s'" % (str(end_yr)),
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def analytical_spin(self, restart_fname, number=None, analytical_spin=False):

        if analytical_spin:
            out_log_fname = os.path.join(self.log_dir,
                                         "%s_log_analytic.nc" % (site))
            if os.path.isfile(out_log_fname):
                os.remove(out_log_fname)
            replace_dict = {
                            "filename%log": "'%s'" % (out_log_fname),
                            "icycle": "12",
                            "filename%restart_out": "'%s'" % (restart_fname),
                            "cable_user%POP_fromZero": ".F.",
                            "cable_user%POP_fromZero": ".F.",
                            "cable_user%CLIMATE_fromZero": ".F.",
                            "cable_user%CASA_DUMP_READ": ".TRUE.",
                            "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                            "cable_user%CASA_NREP": "1",
                            "cable_user%POP_out": "'ini'",
                            "cable_user%SOIL_STRUC": "'default'",
                            "leaps": ".FALSE.",
                            "spincasa": ".TRUE.",
                            "casafile%cnpipool": 'poolcnpIn.csv',
                            "casafile%c2cdumppath": "'./'",
            }
            self.adjust_nml_file(self.nml_fn, replace_dict)




if __name__ == "__main__":

    site = "TumbaFluxnet"

    cwd = os.getcwd()
    driver_dir = "../../driver_files/"
    met_dir = "../../met_data/plumber_met/"
    co2_ndep_dir = "../../met_data/co2_ndep"
    output_dir = "outputs"
    restart_dir = "restart_files"
    nml_fn = "cable.nml"
    site_nml_fn = "site.nml"

    met_fname = os.path.join(met_dir, '%s.1.4_met.nc' % (site))
    co2_ndep_fname = os.path.join(co2_ndep_dir,
                                  "AmaFACE_co2npdepforcing_1850_2100_AMB.csv")
    veg_param_fn = "def_veg_params_zr_clitt_fixed.txt"
    log_dir = "logs"

    exe = "../../src/CABLE_SLI_JV_ratio/CABLE-trunk_checks_extract_sli_optimise_JVratio_vanessa/offline/cable"
    aux_dir = "../../src/CABLE-AUX/"
    verbose = True
    SPIN_UP = True

    C = RunCable(site, driver_dir, output_dir, restart_dir, met_fname,
                 co2_ndep_fname, nml_fn, site_nml_fn, veg_param_fn, log_dir,
                 exe, aux_dir, verbose, SPIN_UP=True)
    C.main()
