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

if [ ! -f $release_file ]; then
    echo "release file not found"; exit
fi

release=`grep Release $release_file | sed s/\=/\ / | awk '{print $NF}'`
root=`grep ATL14_root $loc_file | sed s/\=/\ / | awk '{print $NF}'`


if [ $# -eq 0 ]; then
    regions="RA IS CN CS SV"
fi

for reg in $regions; do 
    [ -d $root/rel$release/north/$reg/ ] && rm -r /discover/nobackup/projects/icesat2/ATL14_processing/rel$release/north/$reg/

    [ -d $reg"_prelim" ] && rm -r $reg"_prelim"
    setup_ATL1415_region.py $loc_file $release_file default_args/north.txt default_args/$reg.txt

    make_ATL1415_queue.py prelim $root/rel$release/north/$reg/input_args_$reg.txt

    setup_slurm_run.py --run_name $reg"_prelim" -q 1415_queue_$reg"_prelim.txt" --time 04:00:00 -j 7 -e ATL14

    pushd $reg"_prelim"
    sbatch slurm_run.sh
    popd

done
