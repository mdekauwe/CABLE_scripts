#!/usr/bin/env python

"""
Site run with full CNP (and Pop)
================================

- Model spin-up: using K34 tower info, CO2=284.7; NDEP-0.79 kg N ha-1 yr-1;
                 PDEP=0.144 kg P ha-1 yr-1
- Transient: 1851-1998, varying CO2 and NDEP, but just recycling the met data
- Historical: actual met (dates correspond to simulated dates) and actual CO2
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
import netCDF4 as nc

class RunCable(object):

    def __init__(self, site, driver_dir, output_dir, restart_dir, met_fname,
                 co2_ndep_fname, nml_fn, site_nml_fn, veg_param_fn, log_dir,
                 exe, aux_dir, verbose):

        self.site = site
        self.driver_dir = driver_dir
        self.output_dir = output_dir
        self.restart_dir = restart_dir
        self.met_fname = met_fname
        self.co2_ndep_fname = co2_ndep_fname
        self.nml_fn = nml_fn
        self.site_nml_fn = site_nml_fn
        self.veg_param_fn = veg_param_fn
        self.restart_fname = "%s_cable_rst.nc" % (self.site)
        self.log_dir = log_dir
        self.cable_exe = exe
        self.aux_dir = aux_dir
        self.verbose = verbose

    def main(self, SPIN_UP=False, TRANSIENT=False, SIMULATION=False):

        f = nc.Dataset(self.met_fname)
        time = nc.num2date(f.variables['time'][:],
                           f.variables['time'].units)

        start_yr = time[0].year
        end_yr = time[-1].year

        if SPIN_UP == True:

            # Initial spin
            self.setup_ini_spin()
            self.run_me()
            self.clean_up(ini=True, tag="zero")

            # 3 sets of spins & analytical spins
            for num in range(1, 4):
                self.setup_re_spin(number=num)
                self.run_me()
                self.clean_up(re_spin=True, tag="ccp%d" % (num))

                self.setup_analytical_spin(number=num)
                self.run_me()
                self.clean_up(analytical=True, tag="saa%d" % (num))

            # one final spin
            num += 1
            self.setup_re_spin(number=num)
            self.run_me()
            self.clean_up(re_spin=True, tag="ccp%d" % (num))

        if TRANSIENT == True:
            self.setup_transient()
            self.run_me()
            self.clean_up(transient=True, tag="transient")

        if SIMULATION == True:
            self.setup_simulation(start_yr, end_yr)
            self.run_me()
            self.clean_up(tag="simulation")

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
                lines[i] = " ".join((key.rstrip(), "=",
                                     replacements_dict.get(key.strip(),
                                     val.lstrip())))

        return '\n'.join(lines) + '\n'

    def setup_ini_spin(self):
        shutil.copyfile(os.path.join(self.driver_dir, "site.nml"),
                        self.site_nml_fn)
        shutil.copyfile(os.path.join(self.driver_dir, "cable.nml"),
                        self.nml_fn)

        # Replace with Vanessa's starting file, this is hardwired until,
        # cable.nml reflects Vanessa's inputs
        vanessa_nml_fn = "ancillary_files/cable.nml.cable_casa_POP_from_zero"
        shutil.copyfile(vanessa_nml_fn, self.nml_fn)

        out_fname = os.path.join(self.output_dir,
                                 "%s_out_cable_zero.nc" % (site))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_zero.txt" % (site))
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

        replace_dict = {
                        "filename%met": "'%s'" % (self.met_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_out": "'%s'" % (self.restart_fname),
                        "filename%type": "'%s'" % (os.path.join(self.aux_dir, "offline/gridinfo_CSIRO_1x1.nc")),
                        "filename%veg": "'%s%s'" % (self.driver_dir, veg_param_fn),
                        "filename%soil": "'%sdef_soil_params.txt'" % (self.driver_dir),
                        "output%restart": ".TRUE.",
                        "casafile%phen": "'%s'" % (os.path.join(self.aux_dir, "core/biogeochem/modis_phenology_csiro.txt")),
                        "casafile%cnpbiome": "'%s'" % (os.path.join(self.driver_dir, "pftlookup_csiro_v16_17tiles_Cumberland.csv")),
                        "cable_user%RunIden": "'%s'" % (self.site),
                        "cable_user%POP_out": "'rst'",
                        "cable_user%POP_rst": "'./'",
                        "cable_user%POP_fromZero": ".T.",
                        "cable_user%CASA_fromZero": ".T.",

        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_re_spin(self, number=None):

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_ccp%d.txt" % (site, number))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        out_fname = os.path.join(self.output_dir,
                                 "%s_out_cable_ccp%d.nc" % (site, number))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_in": "'%s'" % (self.restart_fname),
                        "filename%restart_out": "'%s'" % (self.restart_fname),
                        "cable_user%POP_fromZero": ".F.",
                        "cable_user%CASA_fromZero": ".F.",
                        "cable_user%POP_out": "'rst'",
                        "cable_user%POP_rst": "'./'",
                        "cable_user%CLIMATE_fromZero": ".F.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".TRUE.",
                        "cable_user%CASA_NREP": "0",
                        "cable_user%SOIL_STRUC": "'sli'",
                        "icycle": "2",
                        "cable_user%POP_out": "'rst'",
                        "leaps": ".TRUE.",
                        "spincasa": ".FALSE.",
                        "casafile%cnpipool": "''",
                        "casafile%c2cdumppath": "' '",
                        "output%restart": ".TRUE.",
                        "filename%out": "'%s'" % (out_fname),
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_analytical_spin(self, number):

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_analytic_%d.txt" % (site, number))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "icycle": "12",
                        "filename%restart_out": "'%s'" % (self.restart_fname),
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
                        "casafile%cnpipool": "'poolcnpIn.csv'",
                        "casafile%c2cdumppath": "'./'",
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_transient(self):
        replace_dict = {
                        "RunType": '"transient"',
        }
        self.adjust_nml_file(self.site_nml_fn, replace_dict)

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_transient.txt" % (site))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        out_fname = os.path.join(self.output_dir,
                                 "%s_out_cable_transient.nc" % (site))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "output%averaging": "'monthly'",
                        "spinup": ".FALSE.",
                        "icycle": "2",
                        "cable_user%YearStart": "1852",
                        "cable_user%YearEnd": "2001",
                        "casafile%cnpipool": "''",
                        "cable_user%POP_out": "'epi'",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "POPLUC": ".T.",
                        "filename%out": "'%s'" % (out_fname),
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_simulation(self, start_yr, end_yr):
        replace_dict = {
                        "RunType": '"historical"',
        }
        self.adjust_nml_file(self.site_nml_fn, replace_dict)

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_simulation.txt" % (site))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        out_fname = os.path.join(self.output_dir,
                                 "%s_out_cable.nc" % (site))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "spinup": ".FALSE.",
                        "output%averaging": "'all'",
                        "icycle": "3",
                        "cable_user%YearStart": "%d" % (start_yr),
                        "cable_user%YearEnd": "%d" % (end_yr),
                        "casafile%cnpipool": "''",
                        "cable_user%POP_out": "'epi'",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "filename%out": "'%s'" % (out_fname),
                        "POPLUC": ".F.",
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def run_me(self):
        # run the model
        if self.verbose:
            os.system("%s" % (self.cable_exe))
        else:
            os.system("%s 1>&2" % (self.cable_exe))

    def clean_up(self, ini=False, re_spin=False, analytical=False,
                  transient=False, tag=None):

        if ini or re_spin:
            for f in glob.glob("*.out"):
                os.remove(f)
            os.remove("new_sumbal")
            os.remove("cnpfluxOut.csv")
            os.remove(glob.glob("%s_*_casa_out.nc" % (site))[0])

        if analytical:
            for f in glob.glob("c2c_*_dump.nc"):
                os.remove(f)
            os.remove("cnpfluxOut.csv")

        if transient:
            for f in glob.glob("*.out"):
                os.remove(f)
            for f in glob.glob("*_out.nc"):
                os.remove(f)
            for f in glob.glob("*_out.nc"):
                os.remove(f)
            for f in glob.glob("%s_*_pop_rst.nc" % (self.site)):
                os.remove(f)
            os.remove("new_sumbal")
            os.remove("cnpfluxOut.csv")
            os.remove("cnpspinlast5.txt")

        fromx = "pop_%s_ini.nc" % (self.site)
        to = os.path.join(self.restart_dir,
                         "pop_%s_ini_%s.nc" % (self.site, tag))
        shutil.copyfile(fromx, to)

        fromx = "%s_climate_rst.nc" % (self.site)
        to = os.path.join(self.restart_dir,
                          "%s_climate_rst_%s.nc" % (self.site, tag))
        shutil.copyfile(fromx, to)

        fromx = "%s_casa_rst.nc" % (self.site)
        to = os.path.join(self.restart_dir,
                          "%s_casa_rst_%s.nc" % (self.site, tag))
        shutil.copyfile(fromx, to)

        fromx = "%s_cable_rst.nc" % (self.site)
        to = os.path.join(self.restart_dir,
                          "%s_cable_rst_%s.nc" % (self.site, tag))
        shutil.copyfile(fromx, to)


if __name__ == "__main__":

    site = "TumbaFluxnet"

    cwd = os.getcwd()
    driver_dir = "../../driver_files/"
    met_dir = "../../met_data/plumber_met/"
    co2_ndep_dir = "../../met_data/co2_ndep"
    log_dir = "logs"
    output_dir = "outputs"
    restart_dir = "restart_files"
    nml_fn = "cable.nml"
    site_nml_fn = "site.nml"

    if not os.path.exists(restart_dir):
        os.makedirs(restart_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    met_fname = os.path.join(met_dir, '%s.1.4_met.nc' % (site))
    co2_ndep_fname = os.path.join(co2_ndep_dir,
                                  "AmaFACE_co2ndepforcing_1850_2015_AMB.csv")
    veg_param_fn = "def_veg_params_zr_clitt_fixed.txt"


    exe = "../../src/CABLE_SLI_JV_ratio/CABLE-trunk_checks_extract_sli_optimise_JVratio_vanessa/offline/cable"
    aux_dir = "../../src/CABLE-AUX/"
    verbose = False

    SPIN_UP = True
    TRANSIENT = True
    SIMULATION = True
    C = RunCable(site, driver_dir, output_dir, restart_dir, met_fname,
                 co2_ndep_fname, nml_fn, site_nml_fn, veg_param_fn, log_dir,
                 exe, aux_dir, verbose)
    C.main(SPIN_UP, TRANSIENT, SIMULATION)
