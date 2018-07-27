#!/bin/bash

# This is the template file, we shouldn't need to touch this, the
# python script should set this off and submit the qsub job
#
# Martin De Kauwe, 27 July 2018
#

#PBS -m ae
#PBS -M mdekauwe\@gmail.com
#PBS -P dt6
#PBS -q normalbw
#PBS -l walltime=2:00:00
#PBS -l mem=10GB
#PBS -l ncpus=28
#PBS -j oe
#PBS -l wd
#PBS -l other=gdata1

module load intel-mpi/4.1.1.036
module load netcdf/4.2.1.1

cable_aux_path="/g/data1/w35/mrd561//CABLE/CABLE_AUX-dev/offline"
met_path="/g/data1/wd9/MetForcing/Global/GSWP3_2017/"
co2_file="Annual_CO2_concentration_until_2010.txt"
gw=$gw           #"TRUE"
average=$avg     #"monthly"
start_yr=1950
prev_yr="$(($start_yr-1))"
end_yr=1951

# Set output stuff
outdir=$outdir
logfile=$outdir"/cable_log_${year}.txt"
outfile=$outdir"/cable_out_${year}.nc"
restart_in=$outdir"/restart_${prev_year}.nc"
restart_out=$outdir"/restart_${year}.nc"
namelist=$outdir"/cable_${year}.nml"
cpus=28
exe="./cable-mpi"

# Create output directory
if [ ! -d $outdir ]
then
    mkdir -p $outdir
fi

# set data dirs
ln -s $cable_aux_path surface_data
ln -s $met_path gswp

year=$start_yr
while [ $year <= $end_yr ]
do

	#Read CO2 concentration
	#set co2_file=`cat Annual_CO2_concentration_until_2010.txt`
	#set co2=$co2_file[$year]

    python ./create_spatial_nml.py -y $year -e "false" -g $gw -a $average \
                                   -l $logfile -o $outfile -i $restart_in \
                                   -r $restart_out -c $co2
	mpirun -n $cpus $exe
    cp ./cable.nml $namelist

    year=$[$year+1]
	prev_yr=$[$prev_yr+1]

done
