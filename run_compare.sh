#!/bin/bash
# This is an edited version of Ned Haughton's Jenkins script to do a quick
# multi-site benchmarking comparison on model trunk vs a branch of CABLE.
#
# It produces some statistics and seasonal plots.
#
#
# TODO: We need to compile the src as part of this...
#

module load pbs

# To be changed, add cmd line parser?
trunk_src="/short/x45/vxh599/CABLEbranches/trunk"
branch_src="/short/x45/vxh599/CABLEbranches/CMIP6-MOSRS"
project=w35

DATA_DIR=/g/data1/wd9/cable_jenkins/PALS_datasets

while [ $# -gt 0 ] # Until you run out of parameters...
do
    case "$1" in
        -t|--TRUNK)
		      trunk_src="$2"
              ;;
        -b|--BRANCH)
		      branch_src="$2"
              ;;
        -p|--PROJECT)
		      project="$2"
              ;;
	    -u|-U|--usage)
		      echo "Usage: [-t/--TRUNK   link to trunk src]"
              echo "       [-b/--BRANCH  link to testing (branch) src]"
              echo "       [-p/--PROJECT project code]"
              echo "       [-h/--help    print this message]"
              exit
              ;;
        -h|-H|--help)
		      echo "Usage: [-t/--TRUNK   link to trunk src]"
              echo "       [-b/--BRANCH  link to testing (branch) src]"
              echo "       [-p/--PROJECT project code]"
              echo "       [-h/--help    print this message]"
              exit
              ;;
    esac
    shift       # Check next set of parameters.
done

if [ -d cable_trunk ]
then
    rm -rf cable_trunk
fi

if [ -d "cable_branch" ]
then
    rm -rf cable_branch
fi

ln -sf $trunk_src cable_trunk
ln -sf $branch_src cable_branch

SITES=(Amplero Blodgett Bugac ElSaler ElSaler2 Espirra FortPeck Harvard Hesse
       Howard Howlandm Hyytiala Kruger Loobos Merbleue Mopane Palang Sylvania
       Tumba UniMich)
#SITES=(Amplero)

# Set up and submit all jobs
i=0

# Perhaps we don't want to run the control case again?
RUN_TRUNK=1
if [ $RUN_TRUNK -eq 1 ]
then
    vers=(trunk branch)
else
    vers=(branch)
fi

for version in "${vers[@]}"
do
    ./cable_$version/scripts/get_aux_files.sh offline run_${version}/

    pushd run_$version

    for site in "${SITES[@]}"
    do
        # Set the input and ouput file names
        met_file=${DATA_DIR}/met/${site}Fluxnet.1.4_met.nc
        flux_file=${DATA_DIR}/flux/${site}Fluxnet.1.4_flux.nc
        out_file=${site}_cable_${version}.nc

        mkdir -p $site

        ../cable_testing/scripts/edit_cable_nml.sh -o ${site}/cable.nml\
            cable.nml "filename%met='${met_file}'" "filename%out='${out_file}'"

        pushd $site
        ln -sf ../cable ../def_soil_params.txt ../def_veg_params.txt \
            ../gridinfo_CSIRO_1x1.nc ./

        if [ $version == "trunk" ]
        then
            cable_exe=$trunk_src/offline/cable
        else
            cable_exe=$branch_src/offline/cable
        fi

        # Run offline
        qsub -N cj$i -P $project -q express -l walltime=0:08:00,mem=100MB,wd \
            -W block=true -- $cable_exe &

        popd

        i=$(echo $i + 1 | bc)
    done
    popd
done

# wait for all cable scripts to run
wait

# Check that all jobs have run
if [[ 40 -ne $(ls run_*/*/*_cable_*.nc | wc -l) ]]
then
    echo "Some runs appear to have failed!"
    ls run_*/*/*_cable_*.nc
    exit 1
fi

# Analyse all outputs
#for version in trunk branch ; do
#    pushd run_$version
#
#    for site in "${SITES[@]}" ; do
#        flux_file=${DATA_DIR}/flux/${site}Fluxnet.1.4_flux.nc
#        out_file=${site}_cable_${version}.nc
#
#        pushd $site
#
#         Calculate metrics
#        ../../cable_testing/scripts/offline_metrics.py calculate $out_file \
#        $flux_file \
#            --out=${version}_${site}_metrics.csv --name=$site
#
#        popd
#
#    done
#
#    # Concatenate CSVs, skipping headers
#    head -n 1 ${SITES[0]}/${version}_${SITES[0]}_metrics.csv >\
#        ${version}_all_metrics.csv
#    for site in "${SITES[@]}" ; do
#        tail -n +2 ${site}/${version}_${site}_metrics.csv >>\
#            ${version}_all_metrics.csv
#    done
#
#    popd
#
#done
#
#cable_testing/scripts/offline_metrics.py compare \
#    Branch run_branch/branch_all_metrics.csv \
#    Trunk run_trunk/trunk_all_metrics.csv

mkdir -p plots

# Metrics plot
#cable_testing/scripts/offline_metrics.py plot normalised_metrics.csv\
#    normalised_metrics.png

for site in "${SITES[@]}"; do
    # Seasonal cycle plots
    ./cable_testing/scripts/benchmark_seasonal_plot.py \
        -o run_trunk/${site}/${site}_cable_trunk.nc \
        -n run_branch/${site}/${site}_cable_branch.nc \
        -p plots/${site}_seasonal_cycle.png
done

#cable_testing/scripts/offline_metrics.py xml -u $BUILD_URL \
#    -o normalised_metrics.xml normalised_metrics.png plots/*



#chmod g+rw -R ./*
