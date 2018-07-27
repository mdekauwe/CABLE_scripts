#!/usr/bin/env python

"""
Run CABLE spatially.

This script sets various things within a qsub script and then submits the run.

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (27.07.2018)"
__email__ = "mdekauwe@gmail.com"

import subprocess
import sys
import os
from shutil import copyfile

def main(debug, template_fn, cable_aux_path, met_path, co2_file, gw, average,
         start_yr, end_yr):

    qs_cmd = ( 'qsub -v cable_aux_path=%s,met_path=%s,co2_file=%s,gw=%s,'
               'average=%s,start_yr=%d,end_yr=%d %s' % \
                (cable_aux_path, met_path, co2_file, gw, average,
                 start_yr, end_yr, template_fn) )
    if debug:
        print(qs_cmd)
    else:
        error = subprocess.call(qs_cmd, shell=True)
        if error is 1:
            raise("Job failed to submit")


if __name__ == "__main__":

    debug = True
    template_dir = "qsub_scripts"
    template_fn = os.path.join(template_dir, "run_cable_spatial_template.sh")
    cable_aux_path = "/g/data1/w35/mrd561/CABLE/CABLE_AUX-dev/offline"
    met_path = "/g/data1/wd9/MetForcing/Global/GSWP3_2017/"
    co2_file = "Annual_CO2_concentration_until_2010.txt"
    gw = "TRUE"
    average = "monthly"
    start_yr = 1950
    end_yr = 1951
    main(debug, template_fn, cable_aux_path, met_path, co2_file, gw, average,
         start_yr, end_yr)
