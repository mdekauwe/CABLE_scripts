#!/usr/bin/env python

"""
Create netcdf files to force CABLE-CNP for FACE MIP type experiment.

Typically this process involves creating:
1. a file for the standard 30 year spin to equillibrium at 1850 (fixed CO2,
   Ndep, Pdep)
2. A transient met file, i.e. 1850-20XX, with time varying CO2, Ndep & Pdep.
3. An expermental run file, with whatever treatment added (i.e. eCO2,
   PPT mainpulation, etc).

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (12.11.2017)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import netCDF4 as nc
import datetime as dt
import numpy as np
import math
import pandas as pd
import calendar

class GenerateMetFiles(object):

    def __init__(self, site,  met_fname, co2_ndep_fname):

        self.site = site
        self.met_fname = met_fname
        self.co2_ndep_fname = co2_ndep_fname
        self.KG_2_G = 1000.0
        self.HA_2_M2 = 10000.0
        self.YR_2_DAY = 365.0

    def create_spin_file(self, ofname, co2_fixed, ndep_fixed, pdep_fixed):
        """
        Read in the met file and dump it back out unmodified except for the
        addition of fixed CO2, fixed Ndep & fixed Pdep

        Todo:
        - Add flag to not add Ndep / Pdep?
        """

        ds = nc.Dataset(self.met_fname)
        out = nc.Dataset(ofname, 'w', format='NETCDF4')

        nc_attrs = ds.ncattrs()
        nc_dims = [dim for dim in ds.dimensions]
        #nc_vars = [var for var in ds.variables]
        nc_vars = ds.variables.keys()

        # Write vars
        for v in nc_vars:
            if len(ds.variables[v].dimensions) == 1:
                out.createDimension(v, ds.variables[v].size)
                out.createVariable(v, ds.variables[v].dtype, (v,))
                ncvar = ds.variables[v]
                if hasattr(ncvar, 'units'):
                    mval = ncvar.units
                    out.variables[v].setncatts({'units': mval})
                out.variables[v][:] = ds.variables[v][:]
            else:
                out.createVariable(v, ds.variables[v].dtype,
                                   ds.variables[v].dimensions)
                if len(ds.variables[v].dimensions) == 2:
                    out.variables[v][:,:] = ds.variables[v][:,:]
                elif len(ds.variables[v].dimensions) == 3:
                    out.variables[v][:,:,:] = ds.variables[v][:,:,:]
                elif len(ds.variables[v].dimensions) == 4:
                    out.variables[v][:,:,:,:] = ds.variables[v][:,:,:,:]
                ncvar = ds.variables[v]
                out = self.write_attributes(v, ncvar, out)

        # Add CO2, NDEP & PDEP
        try:
            out_length = len(ds.variables["Rainf"][:,:])
        except:
            out_length = len(ds.variables["Precip"][:,:])

        #ndim = 1
        #out.createDimension('z', ndim)

        for v in ["CO2air", "Ndep", "Pdep"]:
            out.createVariable(v, 'float32', ('time', 'z', 'y', 'x'))

        out.variables["CO2air"][:,:,:,:] = np.ones((out_length,1,1,1)) * \
                                                   co2_fixed
        out.variables["CO2air"].setncatts({'units': "ppm"})
        out.variables["CO2air"].setncatts({'missing_value': "-9999"})
        out.variables["CO2air"].setncatts({'long_name':
                                           "Atmosphereic CO2 concentration"})

        # kg ha-1 y-1 -> gN m-2 d-1
        conv = self.KG_2_G / self.HA_2_M2 / self.YR_2_DAY
        out.variables["Ndep"][:,:,:,:] = np.ones((out_length,1,1,1)) * \
                                                 (ndep_fixed * conv)
        out.variables["Ndep"].setncatts({'units': "gN/m^2/d^1"})
        out.variables["Ndep"].setncatts({'missing_value': "-9999"})
        out.variables["Ndep"].setncatts({'long_name': "N deposition"})

        out.variables["Pdep"][:,:,:,:] = np.ones((out_length,1,1,1)) *  \
                                                 (pdep_fixed * conv)
        out.variables["Pdep"].setncatts({'units': "gP/m^2/d^1"})
        out.variables["Pdep"].setncatts({'missing_value': "-9999"})
        out.variables["Pdep"].setncatts({'long_name': "P deposition"})

        # write global attributes
        for ncattr in ds.ncattrs():
            if ncattr != "_NCProperties":
                out.setncattr(ncattr, getattr(ds, ncattr))

        ds.close()
        out.close()

    def create_transient_file(self, ofname):
        """
        Read in the met file, generate new met outputs covering a transient
        period from 1850 to the start of the met file. We are also adding
        time varying CO2, Ndep & Pdep.

        Todo:
        - Add flag to not add Ndep / Pdep?
        """
        ds = nc.Dataset(self.met_fname)
        out = nc.Dataset(ofname, 'w', format='NETCDF4')

        #offset = ds.variables['time'][1] - ds.variables['time'][0]
        #time = nc.num2date(ds.variables['time'][:] -offset,
        #                   ds.variables['time'].units)

        time = nc.num2date(ds.variables['time'][:],
                           ds.variables['time'].units)

        pre_industrial = 1850
        start_yr = time[0].year
        end_yr = time[-1].year

        (actual_yrs, yr_seq) = self.get_leap_sequence_of_years(start_yr, end_yr)

        # Figure out the number of timesteps in the new file
        orig_years = np.asarray([i.year for i in time])

        n = 0
        for yr in yr_seq:
            idx = np.argwhere(orig_years==yr)
            n += len(idx)

        nc_attrs = ds.ncattrs()
        nc_dims = [dim for dim in ds.dimensions]
        #nc_vars = [var for var in ds.variables]
        nc_vars = ds.variables.keys()

        # Write vars
        for v in nc_vars:
            if len(ds.variables[v].dimensions) == 1:
                if v == "time":
                    out.createDimension(v, n)
                else:
                    out.createDimension(v, ds.variables[v].size)
                out.createVariable(v, ds.variables[v].dtype, (v,))
                ncvar = ds.variables[v]
                if hasattr(ncvar, 'units'):
                    mval = "seconds since %s-01-01 00:00:00" % (pre_industrial)
                    if v == "time":
                        out.variables[v].setncatts({'units': mval})
                    else:
                        mval = ncvar.units
                        out.variables[v].setncatts({'units': mval})

                if v == "time":
                    jump = ds.variables[v][1]-ds.variables[v][0]
                    tx = []
                    val = 0.0
                    for i in range(n):
                        tx.append(val)
                        val += jump
                    out.variables[v][:] = np.asarray(tx)
                else:
                    out.variables[v][:] = ds.variables[v][:]

            else:
                out.createVariable(v, ds.variables[v].dtype,
                                   ds.variables[v].dimensions)
                if len(ds.variables[v].dimensions) == 2:
                    out.variables[v][:,:] = ds.variables[v][:,:]
                elif len(ds.variables[v].dimensions) == 3:
                    tx = np.zeros((n,1,1))
                    cnt1 = 0
                    cnt2 = 0
                    for yr in yr_seq:
                        idx = np.argwhere(orig_years==yr)
                        st = idx[0][0]
                        en = idx[-1][0] + 1
                        slice = ds.variables[v][st:en,:,:]
                        cnt2 += len(slice)
                        tx[cnt1:cnt2,:,:] = ds.variables[v][st:en,:,:]
                        cnt1 += len(slice)
                    out.variables[v][:,:,:] = tx
                elif len(ds.variables[v].dimensions) == 4:
                    tx = np.zeros((n,1,1,1))
                    cnt1 = 0
                    cnt2 = 0
                    for yr in yr_seq:
                        idx = np.argwhere(orig_years==yr)
                        st = idx[0][0]
                        en = idx[-1][0] + 1
                        slice = ds.variables[v][st:en,:,:,:]
                        cnt2 += len(slice)
                        tx[cnt1:cnt2,:,:,:] = ds.variables[v][st:en,:,:,:]
                        cnt1 += len(slice)
                    out.variables[v][:,:,:,:] = tx
                ncvar = ds.variables[v]
                out = self.write_attributes(v, ncvar, out)


        #ndim = 1
        #out.createDimension('z', ndim)

        for v in ["CO2air", "Ndep", "Pdep"]:
            out.createVariable(v, 'float32', ('time', 'z', 'y', 'x'))

        # Add CO2, NDEP & PDEP
        df = pd.read_csv(self.co2_ndep_fname, sep=';')

        cnt1 = 0
        cnt2 = 0
        cx = np.zeros((n,1,1))
        nx = np.zeros((n,1,1))
        px = np.zeros((n,1,1))
        y = pre_industrial
        for yr in yr_seq:
            idx = np.argwhere(orig_years==yr)
            st = idx[0][0]
            en = idx[-1][0] + 1
            slice = ds.variables["SWdown"][st:en,:,:]
            cnt2 += len(slice)
            cx[cnt1:cnt2,:,:] = df[df.Year == y]["CO2 [ppm]"].values[0]
            nx[cnt1:cnt2,:,:] = df[df.Year == y]["ndepo [kgN/ha/yr]"].values[0]
            px[cnt1:cnt2,:,:] = df[df.Year == y]["pdepo [kgP/ha/yr]"].values[0]
            cnt1 += len(slice)
            y += 1

        out.variables["CO2air"][:,:,:,:] = cx
        out.variables["CO2air"].setncatts({'units': "ppm"})
        out.variables["CO2air"].setncatts({'missing_value': "-9999"})
        out.variables["CO2air"].setncatts({'long_name':
                                           "Atmosphereic CO2 concentration"})

        # kg ha-1 y-1 -> gN m-2 d-1
        conv = self.KG_2_G / self.HA_2_M2 / self.YR_2_DAY
        out.variables["Ndep"][:,:,:,:] = nx * conv
        out.variables["Ndep"].setncatts({'units': "gN/m^2/d^1"})
        out.variables["Ndep"].setncatts({'missing_value': "-9999"})
        out.variables["Ndep"].setncatts({'long_name': "N deposition"})

        out.variables["Pdep"][:,:,:,:] = px * conv
        out.variables["Pdep"].setncatts({'units': "gP/m^2/d^1"})
        out.variables["Pdep"].setncatts({'missing_value': "-9999"})
        out.variables["Pdep"].setncatts({'long_name': "P deposition"})

        # write global attributes
        for ncattr in ds.ncattrs():
            if ncattr != "_NCProperties":
                out.setncattr(ncattr, getattr(ds, ncattr))

        ds.close()
        out.close()

    def create_simulation_file(self, ofname):

        ds = nc.Dataset(self.met_fname)
        out = nc.Dataset(ofname, 'w', format='NETCDF4')

        time = nc.num2date(ds.variables['time'][:],
                           ds.variables['time'].units)

        start_yr = time[0].year
        end_yr = time[-1].year
        num_yrs = start_yr - end_yr
        yr_sequence = np.arange(start_yr, end_yr+1)

        # Figure out the number of timesteps in the new file
        years = np.asarray([i.year for i in time])

        n = 0
        for yr in yr_sequence:
            idx = np.argwhere(years==yr)
            n += len(idx)

        nc_attrs = ds.ncattrs()
        nc_dims = [dim for dim in ds.dimensions]
        #nc_vars = [var for var in ds.variables]
        nc_vars = ds.variables.keys()

        # Write vars
        for v in nc_vars:
            if len(ds.variables[v].dimensions) == 1:
                out.createDimension(v, ds.variables[v].size)
                out.createVariable(v, ds.variables[v].dtype, (v,))
                ncvar = ds.variables[v]
                if hasattr(ncvar, 'units'):
                    mval = ncvar.units
                    out.variables[v].setncatts({'units': mval})
                out.variables[v][:] = ds.variables[v][:]
            else:
                out.createVariable(v, ds.variables[v].dtype,
                                   ds.variables[v].dimensions)
                if len(ds.variables[v].dimensions) == 2:
                    out.variables[v][:,:] = ds.variables[v][:,:]
                elif len(ds.variables[v].dimensions) == 3:
                    out.variables[v][:,:,:] = ds.variables[v][:,:,:]
                elif len(ds.variables[v].dimensions) == 4:
                    out.variables[v][:,:,:,:] = ds.variables[v][:,:,:,:]
                ncvar = ds.variables[v]
                out = self.write_attributes(v, ncvar, out)

        #ndim = 1
        #out.createDimension('z', ndim)

        for v in ["CO2air", "Ndep", "Pdep"]:
            out.createVariable(v, 'float32', ('time', 'z', 'y', 'x'))

        # Add CO2, NDEP & PDEP
        df = pd.read_csv(self.co2_ndep_fname, sep=';')

        cnt1 = 0
        cnt2 = 0
        cx = np.zeros((n,1,1))
        nx = np.zeros((n,1,1))
        px = np.zeros((n,1,1))
        y = start_yr
        for yr in yr_sequence:
            idx = np.argwhere(years==yr)
            st = idx[0][0]
            en = idx[-1][0] + 1
            slice = ds.variables["SWdown"][st:en,:,:]
            cnt2 += len(slice)
            #print(y, st, en, cnt1, cnt2, n)
            cx[cnt1:cnt2,:,:] = df[df.Year == y]["CO2 [ppm]"].values[0]
            nx[cnt1:cnt2,:,:] = df[df.Year == y]["ndepo [kgN/ha/yr]"].values[0]
            px[cnt1:cnt2,:,:] = df[df.Year == y]["pdepo [kgP/ha/yr]"].values[0]
            cnt1 += len(slice)
            y += 1

        out.variables["CO2air"][:,:,:,:] = cx
        out.variables["CO2air"].setncatts({'units': "ppm"})
        out.variables["CO2air"].setncatts({'missing_value': "-9999"})
        out.variables["CO2air"].setncatts({'long_name':
                                           "Atmosphereic CO2 concentration"})

        # kg ha-1 y-1 -> gN m-2 d-1
        conv = self.KG_2_G / self.HA_2_M2 / self.YR_2_DAY
        out.variables["Ndep"][:,:,:,:] = nx * conv
        out.variables["Ndep"].setncatts({'units': "gN/m^2/d^1"})
        out.variables["Ndep"].setncatts({'missing_value': "-9999"})
        out.variables["Ndep"].setncatts({'long_name': "N deposition"})

        out.variables["Pdep"][:,:,:,:] = px * conv
        out.variables["Pdep"].setncatts({'units': "gP/m^2/d^1"})
        out.variables["Pdep"].setncatts({'missing_value': "-9999"})
        out.variables["Pdep"].setncatts({'long_name': "P deposition"})

        # write global attributes
        for ncattr in ds.ncattrs():
            if ncattr != "_NCProperties":
                out.setncattr(ncattr, getattr(ds, ncattr))

        ds.close()
        out.close()



    def check_differences(self, ofname):
        ds1 = nc.Dataset(self.met_fname)
        ds2 = nc.Dataset(ofname)

        nc_vars = [var for var in ds1.variables]
        for v in nc_vars:
            if not np.array_equal(ds1.variables[v], ds2.variables[v]):
                print(v, "Differ")

        #for v in nc_vars:
        #    print(v,  ds1.variables[v].dtype, ds2.variables[v].dtype)

        import matplotlib.pyplot as plt
        plt.plot(ds2.variables["Tair"][0:250,0,0]-273.15)
        plt.plot(ds1.variables["Tair"][0:250,0,0]-273.15, "r.")
        plt.show()

    def write_attributes(self, v, ncvar, out):

        if hasattr(ncvar, 'units'):
            val = ncvar.units
            out.variables[v].setncatts({'units': val})
        if hasattr(ncvar, 'missing_value'):
            val = ncvar.missing_value
            out.variables[v].setncatts({'missing_value': val})
        if hasattr(ncvar, 'long_name'):
            val = ncvar.long_name
            out.variables[v].setncatts({'long_name': val})
        if hasattr(ncvar, 'values'):
            val = ncvar.values
            out.variables[v].setncatts({'values': val})
        if hasattr(ncvar, 'CF_name'):
            val = ncvar.CF_name
            out.variables[v].setncatts({'CF_name': val})

        return out

    def get_leap_sequence_of_years(self, start_yr, end_yr):

        pre_industrial = 1850

        # Set the seed so we can repeat this if required
        np.random.seed(42)

        # length of met record
        nrec = end_yr - start_yr

        req_yrs = start_yr - pre_industrial - 1

        # number of times met data is recycled during transient simulation
        nloop_transient = int(np.ceil(req_yrs / nrec))

        # Create a big sequence of years.
        seq = np.tile(np.tile(np.arange(start_yr, end_yr),
                        nloop_transient), 20)

        i = 0
        yr_seq = np.zeros(0)
        actual_yrs = np.zeros(0)
        for yr in np.arange(pre_industrial, start_yr):
            year = seq[i]

            if calendar.isleap(yr):

                while calendar.isleap(year) == False:
                    i += 1
                    year = seq[i]
            i += 1

            yr_seq = np.append(yr_seq, year)
            actual_yrs = np.append(actual_yrs, yr)
        return actual_yrs, yr_seq


if __name__ == "__main__":

    site = "AU-Tum"
    met_dir = "met"
    co2_ndep_dir = "co2_ndep"
    met_fname = os.path.join(met_dir, "AU-Tum_2002-2016_OzFlux_Met.nc")
    co2_ndep_fn = "AmaFACE_co2npdepforcing_1850_2100_AMB.csv"
    co2_ndep_fname = os.path.join(co2_ndep_dir, co2_ndep_fn)

    G = GenerateMetFiles(site, met_fname, co2_ndep_fname)

    # CREATE SPIN UP FILE
    co2_fixed = 284.7  # umol mol-1
    ndep_fixed = 0.79  # kg N ha-1 yr-1
    pdep_fixed = 0.144 # kg N ha-1 yr-1
    ofname = os.path.join(met_dir, "%s_met_spin.nc" % (site))
    G.create_spin_file(ofname, co2_fixed, ndep_fixed, pdep_fixed)

    # CREATE TRANSIENT FILE
    ofname = os.path.join(met_dir, "%s_met_trans.nc" % (site))
    G.create_transient_file(ofname)

    # CREATE SIMULATION FILE
    ofname = os.path.join(met_dir, "%s_met_simulation.nc" % (site))
    G.create_simulation_file(ofname)
    G.check_differences(ofname)
