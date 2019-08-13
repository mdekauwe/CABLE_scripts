#!/usr/bin/Rscript

# Generate the GSWP3 forcing to run a single CABLE pixel
#
# $ module load gdal
# $ module load proj
# $ module load R/3.3.2
# $ R
# $ install.packages("raster")
#

library(raster)
library(ncdf4)

# clear R environment
rm(list=ls(all=TRUE))

year <- 1995
xx <- 662 + 1 # for R
yy <- 121 + 1 # for R

indir <- "/g/data1/wd9/MetForcing/Global/GSWP3_2017/"
met_vars <- c("LWdown",	"PSurf", "SWdown", "Tair",
              "Qair", "Rainf", "Snowf", "Wind")
units <- c("W/m2", "Pa", "W/m2", "K", "kg/kg", "kg/m2/s", "kg/m2/s", "m/s")

# Read met data
met_data <- lapply(met_vars, function(x) brick(list.files(paste0(indir, "/", x),
                                               pattern=as.character(year),
                                               full.names=TRUE)))

# Extract one pixel
pixel <- lapply(met_data, function(x) x[xx,yy])


# Create lon and lat vectors
lat <- rev(seq(-89.75, 89.75, by=0.5))[yy]
lon <- (seq(0.25, 359.75, by=0.5))[xx]#(seq(0.25, 359.75, by=0.5)-180)[100]

# Create new nc_file
time <- seq(0, by=60*60*3, length.out=length(pixel[[1]]))
time_unit <- "seconds since 1995-01-01 00:00:00"

xd = ncdim_def('x',vals=c(1),units='')
yd = ncdim_def('y',vals=c(1),units='')
zd = ncdim_def('z',vals=c(1),units='')

# Define time dimension:
td = ncdim_def('time', unlim=TRUE, units=time_unit, vals=time)

# First set correct dimensions
# (Tair, Qair, CO2air and Wind need an extra z-dimension)
ind_dim <- which(met_vars %in% c("Tair", "Qair", "PSurf", "Wind"))

dims <- lapply(1:length(met_vars), function(x) if(x %in% ind_dim)
                                                  list(xd,yd,zd,td) else
                                                  list(xd,yd,td))

# Create variable definitions for time series variables
var_defs <- mapply(function(i, dim) ncvar_def(name=met_vars[i],
                                              units=units[i],
                                              dim=dim, missval=-9999,
                                              longname=""),
                   i=1:length(met_vars), dim=dims, SIMPLIFY=FALSE)

# First necessary non-time variables:
# Define latitude:
latdim <- ncvar_def('latitude','degrees_north',dim=list(xd,yd),
                    missval=-9999, longname='Latitude')

# Define longitude:
londim <- ncvar_def('longitude','degrees_east',dim=list(xd,yd),
                    missval=-9999,longname='Longitude')

# Collate variables
all_vars <- append(var_defs, c(list(latdim), list(londim)))

# Create
ncid <- nc_create("single_pixel_gswp3_forcing.nc", vars=all_vars)

# Add variable data to file:
ncvar_put(ncid, latdim, vals=lat)
ncvar_put(ncid, londim, vals=lon)

# Time dependent variables:
lapply(1:length(var_defs), function(x) ncvar_put(nc=ncid,
                                                 varid=var_defs[[x]],
                                                 vals=pixel[[x]]))


nc_close(ncid)
