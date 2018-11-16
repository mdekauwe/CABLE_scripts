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
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

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
    keys = [] # save our keys
    lines = text.splitlines()
    for i, row in enumerate(lines):
        # skip blank lines
        if not row.strip():
            continue
        elif "=" not in row:
            lines[i] = row
            continue
        elif not row.startswith("&"):
            key = row.split("=")[0]
            val = row.split("=")[1]
            lines[i] = " ".join((key.rstrip(), "=",
                                 replacements_dict.get(key.strip(),
                                 val.lstrip())))
            keys.append(key.strip())

    # Make sure our replacements were in the namelist to begin with, it is
    # possible they were missing completely so we will need to add these
    # manually
    fix_end_statement = False
    for key, val in replacements_dict.items():
        if key not in keys:
            key_to_add = " ".join((key.rstrip(), "=", val.strip()))
            # add 3 extra spaces at the front to line things up
            string_length = len(key_to_add) + 3
            key_to_add = key_to_add.rjust(string_length)
            lines[i] = key_to_add
            fix_end_statement = True

            # Need to add an extra element and sort out counter
            lines.append("")
            i += 1

    if fix_end_statement:
        lines[i] = "&end"

    return '\n'.join(lines) + '\n'

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

        # skip blank lines
        if not row.strip():
            continue
        # Lines without key = value statement
        elif "=" not in row:
            continue
        # Skip lines that are just comments as these can have "=" too
        elif row.strip().startswith("!"):
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

def change_LAI(met_fname, site, fixed=None, lai_dir=None):

    new_met_fname = "%s_tmp.nc" % (site)

    if fixed is not None:
        lai = fixed
    else:
        lai_fname = os.path.join(lai_dir, "%s_lai.csv" % (site))
        df_lai = pd.read_csv(lai_fname)

        if lai_dir is not None:
            ds = xr.open_dataset(met_fname)

            vars_to_keep = ['Tair']
            dfx = ds[vars_to_keep].squeeze(dim=["x","y","z"],
                                          drop=True).to_dataframe()

            out_length = len(dfx.Tair)
            time_idx = dfx.index
            dfx = dfx.reindex(time_idx)
            dfx['year'] = dfx.index.year
            dfx['doy'] = dfx.index.dayofyear

            df_lai_out = pd.DataFrame(columns=['LAI'], index=range(out_length))

            df_non_leap = df_lai.copy()
            df_leap = df_lai.copy()
            extra = pd.DataFrame(columns=["LAI"], index=range(1))
            extra.LAI[0] = df_lai.LAI[364]
            df_leap = df_leap.append(extra)
            df_leap = df_leap.reset_index(drop=True)


            idx = 0
            for yr in np.unique(dfx.index.year):
                #print(yr)
                ndays = int(len(dfx[dfx.index.year == yr]) / 48)

                if ndays == 366:
                    st = idx
                    en = idx + (366 * 48)
                    df_lai_out.LAI[st:en] = np.repeat(df_leap.values, 48)
                else:
                    st = idx
                    en = idx + (365 * 48)
                    df_lai_out.LAI[st:en] = np.repeat(df_non_leap.values, 48)

                idx += (ndays * 48)

    shutil.copyfile(met_fname, new_met_fname)

    nc = netCDF4.Dataset(new_met_fname, 'r+')
    (nc_attrs, nc_dims, nc_vars) = ncdump(nc)

    nc_var = nc.createVariable('LAI', 'f4', ('time', 'y', 'x'))
    nc.setncatts({'long_name': u"Leaf Area Index",})
    if lai_dir is not None:
        nc.variables['LAI'][:,0,0] = df_lai_out.LAI.values.reshape(out_length, 1, 1)
    else:
        nc.variables['LAI'][:,0,0] = lai
    nc.close()  # close the new file

    return new_met_fname

def get_years(met_fname, nyear_spinup):
    """
    Figure out the start and end of the met file, the number of times we
    need to recycle the met data to cover the transient period and the
    start and end of the transient period.
    """
    pre_indust = 1850

    ds = xr.open_dataset(met_fname)

    st_yr = pd.to_datetime(ds.time[0].values).year

    # PALS met files final year tag only has a single 30 min, so need to
    # end at the previous year, which is the real file end
    en_yr = pd.to_datetime(ds.time[-1].values).year - 1

    # length of met record
    nrec = en_yr - st_yr + 1

    # number of times met data is recycled during transient simulation
    nloop_transient = np.ceil((st_yr - 1 - pre_indust) / nrec) - 1

    # number of times met data is recycled with a spinup run of nyear_spinup
    nloop_spin = np.ceil(nyear_spinup / nrec)

    st_yr_transient = st_yr - 1 - nloop_transient * nrec + 1
    en_yr_transient = st_yr_transient + nloop_transient * nrec - 1

    en_yr_spin = st_yr_transient - 1
    st_yr_spin = en_yr_spin - nloop_spin * nrec + 1

    return (st_yr, en_yr, st_yr_transient, en_yr_transient,
            st_yr_spin, en_yr_spin)

def check_steady_state(experiment_id, output_dir, num, debug=True):
    """
    Check whether the plant (leaves, wood and roots) and soil
    (fast, slow and active) carbon pools have reached equilibrium. To do
    this we are checking the state of the last year in the previous spin
    cycle to the state in the final year of the current spin cycle.
    """
    tol = 0.05 # This is quite high, I use 0.005 in GDAY
    g_2_kg = 0.001

    if num == 1:
        prev_cplant = 99999.9
        prev_csoil = 99999.9
    else:
        fname = "%s_out_CASA_ccp%d.nc" % (experiment_id, num-1)
        fname = os.path.join(output_dir, fname)
        ds = xr.open_dataset(fname)
        prev_cplant = ds.cplant[:,:,0].values[-1].sum() * g_2_kg
        prev_csoil = ds.csoil[:,:,0].values[-1].sum() * g_2_kg

    fname = "%s_out_CASA_ccp%d.nc" % (experiment_id, num)
    fname = os.path.join(output_dir, fname)
    ds = xr.open_dataset(fname)
    new_cplant = ds.cplant[:,:,0].values[-1].sum() * g_2_kg
    new_csoil = ds.csoil[:,:,0].values[-1].sum() * g_2_kg

    if ( np.fabs(prev_cplant - new_cplant) < tol and
         np.fabs(prev_csoil - new_csoil) < tol ):
        not_in_equilibrium = False
    else:
        not_in_equilibrium = True

    if debug:
        print("*", num, not_in_equilibrium,
              "*cplant", np.fabs(prev_cplant - new_cplant),
              "*csoil", np.fabs(prev_csoil - new_csoil))

    return not_in_equilibrium
