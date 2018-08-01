#!/usr/bin/env python

"""
Add SVN info and cable namelist file to the output file

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (01.08.2018)"
__email__ = "mdekauwe@gmail.com"


import os
import sys
import netCDF4

def add_attributes_to_output_file(nml_fname, fname, url, rev):

    # Add SVN info to output file
    nc = netCDF4.Dataset(fname, 'r+')
    nc.setncattr('cable_branch', url)
    nc.setncattr('svn_revision_number', rev)

    # Add namelist to output file
    fp = open(nml_fname, "r")
    namelist = fp.readlines()
    fp.close()

    for i, row in enumerate(namelist):
        # skip blank lines
        if not row.strip():
            continue
        # skip comment lines
        if "=" not in row:
            continue
        elif not row.startswith("&"):
            key = str(row.strip().split("=")[0]).rstrip()
            val = str(row.strip().split("=")[1]).rstrip()
            nc.setncattr(key, val)

    nc.close()
