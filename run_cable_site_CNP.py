#!/usr/bin/env python

"""
CABLE site run with full CNP (and POP)
================================

- Model spin-up: using pre-industrial CO2, NDEP, PDEP. Currently this is using
                 values from AmazonFACE experiment, i.e. CO2=284.7;
                 NDEP-0.79 kg N ha-1 yr-1; PDEP=0.144 kg P ha-1 yr-1
- Transient: 1851-1998, varying CO2 and NDEP and Pdep, but just recycling the
             met data
- Historical: actual met (dates correspond to simulated dates) and actual CO2
- CNP + POP switched on.

During the spinup, we are recycling in the input meteorological forcing file in
chunks.

That's all folks.

"""

__author__ = "Martin De Kauwe, Vanessa Haverd"
__version__ = "1.0 (16.07.2018)"
__email__ = "mdekauwe@gmail.com, Vanessa.Haverd@csiro.au"

import os
import sys
import glob
import shutil
import tempfile
import netCDF4 as nc
import math


class RunCable(object):

    def __init__(self, site, driver_dir, param_dir, output_dir, restart_dir,
                 met_fname, co2_ndep_fname, nml_fn, site_nml_fn, veg_param_fn,
                 log_dir, exe, aux_dir, biogeochem, pop_on, verbose):

        self.site = site
        self.driver_dir = driver_dir
        self.param_dir = param_dir
        self.output_dir = output_dir
        self.restart_dir = restart_dir
        self.met_fname = met_fname
        self.co2_ndep_fname = co2_ndep_fname
        self.nml_fn = nml_fn
        self.site_nml_fn = site_nml_fn
        self.veg_param_fn = veg_param_fn
        self.restart_fname = "%s_cable_rst.nc" % (self.site)
        self.casa_restart_fname = "%s_casa_rst.nc" % (self.site)
        self.pop_restart_fname = "%s_pop_rst.nc" % (self.site)
        self.climate_restart_fname = "%s_climate_rst.nc" % (self.site)
        self.log_dir = log_dir
        self.cable_exe = exe
        self.aux_dir = aux_dir
        self.verbose = verbose
        self.nyear_spinup = nyear_spinup
        if biogeochem == "C":
            self.biogeochem = 1
        elif biogeochem == "CN":
            self.biogeochem = 2
        elif biogeochem == "CNP":
            self.biogeochem = 3
        else:
            raise ValueError("Unknown biogeochemistry option: C, CN, CNP")
        if pop_on:
            self.call_pop = ".TRUE."
        else:
            self.call_pop = ".FALSE."

    def main(self, SPIN_UP=False, TRANSIENT=False, SIMULATION=False):

        (st_yr, en_yr,
         st_yr_trans, en_yr_trans,
         st_yr_spin, en_yr_spin) = self.get_years()

        if SPIN_UP == True:

            # Initial spin
            self.setup_ini_spin(st_yr_spin, en_yr_spin, st_yr, en_yr)
            self.run_me()

            # 3 sets of spins & analytical spins
            for num in range(1, 4):
                self.logfile="log_ccp%d" % (num)
                self.setup_re_spin(number=num)
                self.run_me()
                self.clean_up(re_spin=True, tag="ccp%d" % (num))

                self.logfile="log_sa%d" % (num)
                self.setup_analytical_spin(number=num, st_yr_spin=st_yr_spin,
                                           en_yr_spin=en_yr_spin )
                self.run_me()
                self.clean_up(analytical=True, tag="saa%d" % (num))

            # one final spin
            num += 1
            self.logfile="log_ccp%d" % (num)
            self.setup_re_spin(number=num)
            self.run_me()
            self.clean_up(re_spin=True, tag="ccp%d" % (num))

        if TRANSIENT == True:
            self.setup_transient(st_yr_trans, en_yr_trans, st_yr, en_yr)
            self.run_me()
            self.clean_up(transient=True, tag="transient")

        if SIMULATION == True:
            self.setup_simulation(st_yr, en_yr)
            self.run_me()
            self.clean_up(tag="simulation")

    def get_years(self):
        f = nc.Dataset(self.met_fname)
        time = nc.num2date(f.variables['time'][:],
                           f.variables['time'].units)

        st_yr = time[0].year

        # PALS met files final year tag only has a single 30 min, so need to
        # end at the previous year, which is the real file end
        en_yr = time[-1].year - 1

        # length of met record
        nrec = en_yr - st_yr + 1

        # number of times met dat is recycled in transient simulation from ~1850
        # to yearstart-1
        nloop_transient = math.ceil((st_yr - 1 - 1850) / nrec) - 1

        # number of times met data is recycled with a spinup run of nyear_spinup
        nloop_spin = math.ceil( self.nyear_spinup/ nrec)

        st_yr_transient = st_yr - 1 - nloop_transient * nrec + 1
        en_yr_transient = st_yr_transient + nloop_transient * nrec - 1

        en_yr_spin = st_yr_transient - 1
        st_yr_spin = en_yr_spin - nloop_spin * nrec + 1

        return (st_yr, en_yr, st_yr_transient, en_yr_transient,
                st_yr_spin, en_yr_spin)

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

    def setup_ini_spin(self, st_yr_spin, en_yr_spin, st_yr, en_yr):
        shutil.copyfile(os.path.join(self.driver_dir, "site.nml"),
                        self.site_nml_fn)
        shutil.copyfile(os.path.join(self.driver_dir, "cable.nml"),
                        self.nml_fn)

        out_fname = os.path.join(self.output_dir,
                                 "%s_out_cable_zero.nc" % (site))
        out_fname_CASA = os.path.join(self.output_dir,
                                 "%s_out_CASA_zero.nc" % (site))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_zero.txt" % (site))
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
        self.adjust_nml_file(self.site_nml_fn, replace_dict)

        replace_dict = {
                        "filename%met": "'%s'" % (self.met_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "casafile%out": "'%s'" % (out_fname_CASA),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_out": "'%s%s'" % (self.restart_dir,self.restart_fname),
                        "cable_user%climate_restart_out": "'%s%s'" % (self.restart_dir,self.climate_restart_fname),
                        "cable_user%POP_restart_out": "'%s%s'" % (self.restart_dir,self.pop_restart_fname),
                        "casafile%cnpepool": "'%s%s'" % (self.restart_dir,self.casa_restart_fname),
                        "filename%restart_in": "''" ,
                        "cable_user%climate_restart_in": "''" ,
                        "cable_user%POP_restart_in": "''",
                        "filename%type": "'%s'" % (os.path.join(self.aux_dir, "offline/gridinfo_CSIRO_1x1.nc")),
                        "filename%veg": "'%s%s'" % (self.param_dir, veg_param_fn),
                        "filename%soil": "'%s%s'" % (self.driver_dir, soil_param_fn),
                        "output%restart": ".TRUE.",
                        "casafile%phen": "'%s'" % (os.path.join(self.aux_dir, "core/biogeochem/modis_phenology_csiro.txt")),
                        "casafile%cnpbiome": "'%s'" % (os.path.join(self.param_dir, bgc_param_fn)),
                        "cable_user%RunIden": "'%s'" % (self.site),
                        "cable_user%POP_out": "'ini'",
                        "cable_user%POP_rst": "'./'",
                        "cable_user%POP_fromZero": ".T.",
                        "cable_user%CASA_fromZero": ".T.",
                        "cable_user%CLIMATE_fromZero": ".T.",
                        "cable_user%YearStart": "%d" % (st_yr_spin),
                        "cable_user%YearEnd": "%d" % (en_yr_spin),
                        "cable_user%CASA_SPIN_STARTYEAR": "%d" % (st_yr_spin),
                        "cable_user%CASA_SPIN_ENDYEAR": "%d" % (en_yr_spin),
                        "cable_user%CALL_POP": "%s" % (self.call_pop),
                        "output%averaging": "'monthly'",
                        "icycle": "%d" % (self.biogeochem),
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_re_spin(self, number=None):

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_ccp%d.txt" % (site, number))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        out_fname = os.path.join(self.output_dir,
                                 "%s_out_cable_ccp%d.nc" % (site, number))
        out_fname_CASA = os.path.join(self.output_dir,
                                 "%s_out_CASA_ccp%d.nc" % (site, number))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_out": "'%s%s'" % (self.restart_dir,self.restart_fname),
                        "cable_user%climate_restart_out": "'%s%s'" % (self.restart_dir,self.climate_restart_fname),
                        "cable_user%POP_restart_out": "'%s%s'" % (self.restart_dir,self.pop_restart_fname),
                        "cable_user%cnpepool": "'%s%s'" % (self.restart_dir,self.casa_restart_fname),
                        "filename%restart_in": "'%s%s'" % (self.restart_dir,self.restart_fname),
                        "cable_user%climate_restart_in": "'%s%s'" % (self.restart_dir,self.climate_restart_fname),
                        "cable_user%POP_restart_in": "'%s%s'" % (self.restart_dir,self.pop_restart_fname),
                        "casafile%cnpipool": "'%s%s'" % (self.restart_dir,self.casa_restart_fname),
                        "cable_user%POP_fromZero": ".F.",
                        "cable_user%CASA_fromZero": ".F.",
                        "cable_user%POP_out": "'ini'",
                        "cable_user%POP_rst": "'./'",
                        "cable_user%CLIMATE_fromZero": ".F.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".TRUE.",
                        "cable_user%CASA_NREP": "0",
                        "cable_user%SOIL_STRUC": "'sli'",
                        "icycle": "%d" % (self.biogeochem),
                        "leaps": ".TRUE.",
                        "spincasa": ".FALSE.",
                        "casafile%c2cdumppath": "' '",
                        "output%restart": ".TRUE.",
                        "filename%out": "'%s'" % (out_fname),
                        "casafile%out": "'%s'" % (out_fname_CASA),
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_analytical_spin(self, number, st_yr_spin, en_yr_spin):

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_analytic_%d.txt" % (site, number))
        out_fname_CASA = os.path.join(self.output_dir,
                                 "%s_out_CASA_analytic_%d.nc" % (site, number))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "icycle": "%d" % (self.biogeochem + 10), # Need to add 10 for spinup
                        "filename%restart_out": "'%s%s'" % (self.restart_dir,self.restart_fname),
                        "cable_user%climate_restart_out": "'%s%s'" % (self.restart_dir,self.climate_restart_fname),
                        "cable_user%POP_restart_out": "'%s%s'" % (self.restart_dir,self.pop_restart_fname),
                        "casafile%cnpepool": "'%s%s'" % (self.restart_dir,self.casa_restart_fname),
                        "filename%restart_in": "'%s%s'" % (self.restart_dir,self.restart_fname),
                        "cable_user%climate_restart_in": "'%s%s'" % (self.restart_dir,self.climate_restart_fname),
                        "cable_user%POP_restart_in": "'%s%s'" % (self.restart_dir,self.pop_restart_fname),
                        "casafile%cnpipool": "'%s%s'" % (self.restart_dir,self.casa_restart_fname),
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
                        "casafile%c2cdumppath": "'./'",
                        "cable_user%CASA_SPIN_STARTYEAR": "%d" % (st_yr_spin),
                        "cable_user%CASA_SPIN_ENDYEAR": "%d" % (en_yr_spin),
                        "casafile%out": "'%s'" % (out_fname_CASA),
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_transient(self, st_yr_trans, en_yr_trans, st_yr, en_yr):
        shutil.copyfile(os.path.join(self.driver_dir, "site.nml"),
                        self.site_nml_fn)
        shutil.copyfile(os.path.join(self.driver_dir, "cable.nml"),
                        self.nml_fn)

        replace_dict = {
                        "RunType": '"transient"',
                        "CO2NDepFile": "'%s'" % (self.co2_ndep_fname),
                        "spinstartyear": "%d" % (st_yr),
                        "spinendyear": "%d" % (en_yr),
           }
        self.adjust_nml_file(self.site_nml_fn, replace_dict)

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_transient.txt" % (site))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        out_fname = os.path.join(self.output_dir,
                                 "%s_out_cable_transient.nc" % (site))
        out_fname_CASA = os.path.join(self.output_dir,
                                 "%s_out_casa_transient.nc" % (site))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        replace_dict = {
                        "filename%met": "'%s'" % (self.met_fname),
                        "filename%out": "'%s'" % (out_fname),
                        "filename%log": "'%s'" % (out_log_fname),
                        "filename%restart_out": "'%s%s'" % (self.restart_dir,self.restart_fname),
                        "cable_user%climate_restart_out": "'%s%s'" % (self.restart_dir,self.climate_restart_fname),
                        "cable_user%POP_restart_out": "'%s%s'" % (self.restart_dir,self.pop_restart_fname),
                        "casafile%cnpepool": "'%s%s'" % (self.restart_dir,self.casa_restart_fname),
                        "filename%restart_in": "'%s%s'" % (self.restart_dir,self.restart_fname),
                        "cable_user%climate_restart_in": "'%s%s'" % (self.restart_dir,self.climate_restart_fname),
                        "cable_user%POP_restart_in": "'%s%s'" % (self.restart_dir,self.pop_restart_fname),
                        "casafile%cnpipool": "'%s%s'" % (self.restart_dir,self.casa_restart_fname),
                        "cable_user%POP_fromZero": ".F.",
                        "cable_user%CASA_fromZero": ".F.",
                        "cable_user%POP_out": "'ini'",
                        "cable_user%POP_rst": "'./'",
                        "cable_user%CLIMATE_fromZero": ".F.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "cable_user%CASA_NREP": "0",
                        "cable_user%SOIL_STRUC": "'sli'",
                        "filename%type": "'%s'" % (os.path.join(self.aux_dir, "offline/gridinfo_CSIRO_1x1.nc")),
                        "filename%veg": "'%s%s'" % (self.param_dir, veg_param_fn),
                        "filename%soil": "'%s%s'" % (self.driver_dir, soil_param_fn),
                        "output%restart": ".TRUE.",
                        "casafile%phen": "'%s'" % (os.path.join(self.aux_dir, "core/biogeochem/modis_phenology_csiro.txt")),
                        "casafile%cnpbiome": "'%s'" % (os.path.join(self.param_dir, bgc_param_fn)),
                        "cable_user%RunIden": "'%s'" % (self.site),
                        "output%averaging": "'monthly'",
                        "spinup": ".FALSE.",
                        "icycle": "%d" % (self.biogeochem),
                        "cable_user%POP_out": "'ini'",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "POPLUC": ".T.",
                        "filename%out": "'%s'" % (out_fname),
                        "cable_user%YearStart": "%d" % (st_yr_trans),
                        "cable_user%YearEnd": "%d" % (en_yr_trans),
                        "casafile%out": "'%s'" % (out_fname_CASA),
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def setup_simulation(self, st_yr, en_yr):
        shutil.copyfile(os.path.join(self.driver_dir, "site.nml"),
                        self.site_nml_fn)
        shutil.copyfile(os.path.join(self.driver_dir, "cable.nml"),
                        self.nml_fn)

        replace_dict = {
                        "RunType": '"historical"',
                        "CO2NDepFile": "'%s'" % (self.co2_ndep_fname),
                        "spinstartyear": "%d" % (st_yr),
                        "spinendyear": "%d" % (en_yr)
        }
        self.adjust_nml_file(self.site_nml_fn, replace_dict)

        out_log_fname = os.path.join(self.log_dir,
                                     "%s_log_simulation.txt" % (site))
        if os.path.isfile(out_log_fname):
            os.remove(out_log_fname)

        out_fname = os.path.join(self.output_dir,
                                 "%s_out_cable.nc" % (site))
        out_fname_CASA = os.path.join(self.output_dir,
                                 "%s_out_casa.nc" % (site))
        if os.path.isfile(out_fname):
            os.remove(out_fname)

        replace_dict = {
                        "filename%log": "'%s'" % (out_log_fname),
                        "spinup": ".FALSE.",
                        "filename%restart_out": "'%s%s'" % (self.restart_dir,self.restart_fname),
                        "cable_user%climate_restart_out": "'%s%s'" % (self.restart_dir,self.climate_restart_fname),
                        "cable_user%POP_restart_out": "'%s%s'" % (self.restart_dir,self.pop_restart_fname),
                        "casafile%cnpepool": "'%s%s'" % (self.restart_dir,self.casa_restart_fname),
                        "filename%restart_in": "'%s%s'" % (self.restart_dir,self.restart_fname),
                        "cable_user%climate_restart_in": "'%s%s'" % (self.restart_dir,self.climate_restart_fname),
                        "cable_user%POP_restart_in": "'%s%s'" % (self.restart_dir,self.pop_restart_fname),
                        "casafile%cnpipool": "'%s%s'" % (self.restart_dir,self.casa_restart_fname),
                        "output%averaging": "'daily'",
                        "icycle": "%d" % (self.biogeochem),
                        "cable_user%YearStart": "%d" % (st_yr),
                        "cable_user%YearEnd": "%d" % (en_yr),
                        "cable_user%POP_out": "'ini'",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "filename%out": "'%s'" % (out_fname),
                        "POPLUC": ".F.",
                        "filename%met": "'%s'" % (self.met_fname),
                        "filename%type": "'%s'" % (os.path.join(self.aux_dir, "offline/gridinfo_CSIRO_1x1.nc")),
                        "filename%veg": "'%s%s'" % (self.param_dir, veg_param_fn),
                        "filename%soil": "'%s%s'" % (self.driver_dir, soil_param_fn),
                        "casafile%phen": "'%s'" % (os.path.join(self.aux_dir, "core/biogeochem/modis_phenology_csiro.txt")),
                        "casafile%cnpbiome": "'%s'" % (os.path.join(self.param_dir, bgc_param_fn)),
                        "cable_user%RunIden": "'%s'" % (self.site),
                        "cable_user%POP_out": "'ini'",
                        "cable_user%POP_rst": "'./'",
                        "cable_user%POP_fromZero": ".F.",
                        "cable_user%CASA_fromZero": ".F.",
                        "cable_user%CLIMATE_fromZero": ".F.",
                        "cable_user%CASA_DUMP_READ": ".FALSE.",
                        "cable_user%CASA_DUMP_WRITE": ".FALSE.",
                        "cable_user%SOIL_STRUC": "'sli'",
                        "spincasa": ".FASLE.",
                        "l_laiFeedbk": ".TRUE.",
                        "output%averaging": "'all'",
                        "casafile%out": "'%s'" % (out_fname_CASA),
        }
        self.adjust_nml_file(self.nml_fn, replace_dict)

    def run_me(self):
        # run the model
        if self.verbose:
            os.system("%s > log" % (self.cable_exe ))
            os.system("mv log $s" % (self.logfile))
        else:
            os.system("%s 1>&2" % (self.cable_exe))

    def clean_up(self, ini=False, re_spin=False, analytical=False,
                  transient=False, tag=None):

        if ini or re_spin:
            for f in glob.glob("*.out"):
                if (os.path.isfile(f)):
                    os.remove(f)
            if (os.path.isfile("new_sumbal")):
                os.remove("new_sumbal")
            if (os.path.isfile("cnpfluxOut.csv")):
                os.remove("cnpfluxOut.csv")

        if analytical:
            if (os.path.isfile("cnpfluxOut.csv")):
                os.remove("cnpfluxOut.csv")

        if transient:
            for f in glob.glob("*.out"):
                if (os.path.isfile(f)):
                    os.remove(f)
            for f in glob.glob("*_out.nc"):
                os.remove(f)
            for f in glob.glob("*_out.nc"):
                os.remove(f)
            for f in glob.glob("%s_*_pop_rst.nc" % (self.site)):
                os.remove(f)
            os.remove("new_sumbal")
            os.remove("cnpfluxOut.csv")

        fromx = self.restart_dir + self.restart_fname
        to = fromx[:-3] + "_" + tag + ".nc"
        shutil.copyfile(fromx, to)

        fromx = self.restart_dir + self.casa_restart_fname
        to = fromx[:-3] + "_" + tag + ".nc"
        shutil.copyfile(fromx, to)

        fromx = self.restart_dir + self.climate_restart_fname
        to = fromx[:-3] + "_" + tag + ".nc"
        shutil.copyfile(fromx, to)

        fromx = self.restart_dir + self.pop_restart_fname
        to = fromx[:-3] + "_" + tag + ".nc"
        shutil.copyfile(fromx, to)



if __name__ == "__main__":

    site = "Cumberland"

    cwd = os.getcwd()
    driver_dir = "driver_files/"
    param_dir = "driver_files/"
    met_dir = "met"
    co2_ndep_dir = "met"
    aux_dir = "../../src/NESP2pt9_TRENDYv7/CABLE-AUX/"
    log_dir = "logs/"
    output_dir = "outputs/"
    restart_dir = "restart_files/"
    nml_fn = "cable.nml"
    site_nml_fn = "site.nml"
    #met_fname = os.path.join(met_dir, '%s.1.4_met.nc' % (site)) # PLUMBER sites
    met_fname = os.path.join(met_dir, 'AU_Cum_2014_2017_met.nc')
    co2_ndep_fname = os.path.join(co2_ndep_dir,
                                  "AmaFACE_co2npdepforcing_1850_2100_AMB.csv")
    veg_param_fn = "def_veg_params.txt"
    bgc_param_fn = "pftlookup.csv"
    soil_param_fn = "def_soil_params.txt"   # only used when soilparmnew = .FALSE. in cable.nml
    exe = "../../src/NESP2pt9_TRENDYv7/NESP2pt9_TRENDYv7/offline/cable"
    verbose = False
    nyear_spinup = 5
    biogeochem = "C" # C, CN, CNP
    pop_on = False

    if not os.path.exists(restart_dir):
        os.makedirs(restart_dir)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    C = RunCable(site, driver_dir, param_dir, output_dir, restart_dir,
                 met_fname, co2_ndep_fname, nml_fn, site_nml_fn, veg_param_fn,
                 log_dir, exe, aux_dir, biogeochem, pop_on, verbose)
    C.main(SPIN_UP=True, TRANSIENT=True, SIMULATION=True)
