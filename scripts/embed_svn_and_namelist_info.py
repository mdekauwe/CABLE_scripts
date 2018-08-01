#!/usr/bin/env python

"""
Add SVN info and cable namelist file to the output file

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (01.08.2018)"
__email__ = "mdekauwe@gmail.com"


import os
import netCDF4

def add_attributes_to_output_file(fname, url, rev):

    # Add SVN info to output file
    nc = netCDF4.Dataset(fname, 'r+')
    nc.setncattr('cable_branch', url)
    nc.setncattr('svn_revision_number', rev)

    # Add namelist to output file
    fp = open(self.nml_fname, "r")
    namelist = fp.readlines()
    fp.close()

    for i, row in enumerate(namelist):
        # skip blank lines
        if not row.strip():
            continue
        if "=" not in row:
            continue
        elif not row.startswith("&"):
            key = str(row.split("=")[0])
            val = str(row.split("=")[1])
            nc.setncattr(key, val)
    nc.close()
