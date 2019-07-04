#!/bin/bash

# This is the template file, we shouldn't need to touch this, the
# python script should set this off and submit the qsub job
#
# Martin De Kauwe, 4 July 2019
#

#PBS -m ae
#PBS -P w35
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -l mem=32GB
#PBS -l ncpus=16
#PBS -j oe
#PBS -l wd
#PBS -l other=gdata1

module load dot
module add intel-cc
module add intel-fc
module load netcdf/4.3.3.1
module load openmpi
module load subversion/1.9.0

#module load intel-mpi/4.1.1.036

source activate nci

cpus=16
exe="./cable"

# variables passed to this script by the python script
start_yr=$start_yr
prev_yr="$(($start_yr-1))"
end_yr=$end_yr
co2_fname=$co2_fname

year=$start_yr
while [ $year -le $end_yr ]
do

    co2_conc=$(gawk -v yr=$year 'NR==yr' $co2_fname)
    echo $co2_conc
    # adjust and make a new namelist file
    restart_in="restart_$prev_year.nc"
    restart_out="restart_$year.nc"
    outfile="cable_out_$year.nc"
    logfile="cable_log_$year.txt"

    python ./run_cable_spatial.py -a -y $year -l $logfile -o $outfile \
                                  -i $restart_in -r $restart_out -c $co2_conc

    mpirun -n $cpus $exe

    year=$[$year+1]
    prev_yr=$[$prev_yr+1]

done
