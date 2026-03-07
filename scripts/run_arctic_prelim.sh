#! /usr/bin/env bash

# run the arctic regions.  The first argument should be the release file, and should specify the release number
# as Release=[NNN]
# the second argument should be the location file and should specify
# ATL14_root=[path]
# If no other arguments are specified, will run RA IS CN CS and SV, otherwise will run the specified region

release_file=$1
shift

loc_file=$1
shift

period_file=$1
shift



if [ ! -f $release_file ]; then
    echo "release file not found"; exit
fi

$(grep -q monthly $period_file) && hemi_suffix="_monthly" || hemi_suffix=""
echo $hemi_suffix

release=`grep Release $release_file | sed s/\=/\ / | awk '{print $NF}'`
root=`grep ATL14_root $loc_file | sed s/\=/\ / | awk '{print $NF}'`
cycles=`grep cycles $release_file | sed s/\=/\ / | awk '{print $NF}'`
version=`grep version $release_file | sed s/\=/\ / | awk '{print $NF}'`

if [ $# -eq 0 ]; then
    regions="RA IS CN CS SV"
else
    regions=$1
fi

for reg in $regions; do
    base=${root}/rel${release}/north${hemi_suffix}/${reg}/

    echo $base

    if $(grep -q monthly $period_file); then
    	ATL14_ref_str=--ATL14_reference_file=$root/rel$release/north/$reg/ATL14_${reg}_${cycles}_100m_${release}_${version}.nc
    	suffix="_monthly"
    fi

    [ -d $base ] && rm -r $base

    [ -d $reg"_prelim" ] && rm -r $reg"_prelim"
    setup_ATL1415_region.py $loc_file $release_file $period_file default_args/north.txt default_args/$reg.txt $ATL14_ref_str

    make_ATL1415_queue.py prelim $base/input_args_$reg.txt
    run_name=${reg}${hemi_suffix}_prelim
    setup_slurm_run.py --run_name $run_name -q 1415_queue_$reg"_prelim.txt" --time 04:00:00 -j 7 -e ATL14

    pushd $run_name
    sbatch slurm_run.sh
    popd

done
