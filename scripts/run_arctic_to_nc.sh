#! /usr/bin/env bash

# make the netcdf files for the arctic regions.  The first argument should be the release file, and should specify the release number
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

> arctic_2nc_queue.txt

for reg in $regions; do 
    echo "ATL14_write2nc.py @${root}/rel${release}/north/${reg}/input_args_${reg}.txt" >> arctic_2nc_queue.txt
    echo "ATL15_write2nc.py @${root}/rel${release}/north/${reg}/input_args_${reg}.txt" >> arctic_2nc_queue.txt
done

rm -r arctic_2nc
setup_slurm_run.py --run_name arctic_2nc -q arctic_2nc_queue.txt --time 02:00:00 -j 1 -e ATL14
mv arctic_2nc_queue.txt arctic_2nc
pushd arctic_2nc; sbatch slurm_run.sh; popd
