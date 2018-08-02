#!/usr/bin/env python

"""
Various CABLE utilities

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (02.08.2018)"
__email__ = "mdekauwe@gmail.com"


import os
import sys
import netCDF4
import shutil
import tempfile
import pandas as pd

def adjust_nml_file(fname, replacements):
    """
    Adjust the params/flags in the CABLE namelise file. Note this writes
    over whatever file it is given!

    Parameters:
    ----------
    replacements : dictionary
        dictionary of replacement values.
    """
    f = open(fname, 'r')
    param_str = f.read()
    f.close()
    new_str = replace_keys(param_str, replacements)
    fd, path = tempfile.mkstemp()
    os.write(fd, str.encode(new_str))
    os.close(fd)
    shutil.copy(path, fname)
    os.remove(path)

def replace_keys(text, replacements_dict):
    """ Function expects to find CABLE namelist file formatted key = value.

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


def add_missing_options_to_nml_file(fname, line_start=None):
    """
    Some of the flags we may wish to change are missin from the default
    file so we can't adjust them via this script...add them
    """

    if line_start is None:
        line_start = sum(1 for line in open(fname)) - 1

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

def get_svn_info(here, there):
    """
    Add SVN info and cable namelist file to the output file
    """

    os.chdir(there)
    os.system("svn info > tmp_svn")
    fname = 'tmp_svn'
    fp = open(fname, "r")
    svn = fp.readlines()
    fp.close()
    os.remove(fname)

    url = [i.split(":", 1)[1].strip() \
            for i in svn if i.startswith('URL')]
    rev = [i.split(":", 1)[1].strip() \
            for i in svn if i.startswith('Revision')]
    os.chdir(here)

    return url, rev


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
        print(row)
        # skip blank lines
        if not row.strip():
            continue
        # skip comment lines
        if "=" not in row:
            continue
        if row.startswith("!"):
            continue
        elif not row.startswith("&"):
            key = str(row.strip().split("=")[0]).rstrip()
            val = str(row.strip().split("=")[1]).rstrip()
            nc.setncattr(key, val)

    nc.close()

def ncdump(nc_fid):
    '''
    ncdump outputs dimensions, variables and their attribute information.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions

    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables

    return nc_attrs, nc_dims, nc_vars

def change_LAI(met_fname=None, fixed=None, lai_fname=None):

    new_met_fname = "tmp.nc"
    if fixed is not None:
        lai = fixed
    else:
        lai = pd.read_csv(lai_fname)

    shutil.copyfile(met_fname, new_met_fname)

    nc = netCDF4.Dataset(new_met_fname, 'r+')
    (nc_attrs, nc_dims, nc_vars) = ncdump(nc)

    nc_var = nc.createVariable('LAI', 'f4', ('time', 'y', 'x'))
    nc.setncatts({'long_name': u"Leaf Area Index",})
    nc.variables['LAI'][:,0,0] = lai
    nc.close()  # close the new file

    return new_met_fname
