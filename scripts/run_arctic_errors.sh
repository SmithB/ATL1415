
# add SV here:


release_file=$1

if [ ! -f $release_file ]; then
    echo "release file not found"; exit
fi

release=`grep Release $release_file | sed s/\=/\ / | awk '{print $NF}'`


for reg in SV; do
    echo $reg
#for reg in RA IS CN CS; do
#    [ -d /discover/nobackup/projects/icesat2/ATL14_processing/rel$release/north/$reg/ ] && rm -r /discover/nobackup/projects/icesat2/ATL14_processing/rel$release/north/$reg/

    [ -d $reg"_prelim" ] && rm -r $reg"_prelim"
    #setup_ATL1415_region.py default_args/discover.txt $release_file default_args/north.txt default_args/$reg.txt 
    #setup_ATL1415_region.py default_args/discover.txt default_args/rel_003.txt default_args/north.txt default_args/$reg.txt

    make_ATL1415_queue.py prelim /discover/nobackup/projects/icesat2/ATL14_processing/rel$release/north/$reg/input_args_$reg.txt --errors_only

    setup_slurm_run.py --run_name $reg"_prelim" -q 1415_queue_$reg"_prelim.txt" --time 04:00:00 -j 4 -e ATL14

    pushd $reg"_prelim"
    #sbatch slurm_run.sh
    popd

done
