#! /usr/bin/env bash

# make the netcdf files for the antarctic regions.  The first argument should be the release file, and should specify the release number
# as Release=[NNN]
# the second argument should be the location file and should specify
# ATL14_root=[path]
# If no other arguments are specified, will run A1 A2 A3 A4, otherwise will run the specified region(s)

release_file=$1
shift

loc_file=$1
shift

if [ ! -f $release_file ]; then
    echo "release file not found"; exit
fi

release=`grep '^--Release=' $release_file | sed s/\=/\ / | awk '{print $NF}'`
root=`grep '^--ATL14_root=' $loc_file | sed s/\=/\ / | awk '{print $NF}'`

if [ -z "$release" ]; then
    echo "could not find --Release= in $release_file"; exit
fi
if [ -z "$root" ]; then
    echo "could not find --ATL14_root= in $loc_file"; exit
fi

if [ $# -eq 0 ]; then
    regions="A1 A2 A3 A4"
else
    regions=$1
fi

> AA_2nc_queue.txt

for reg in $regions; do
    echo "ATL14_write2nc.py @${root}/rel${release}/south/${reg}/input_args_${reg}.txt" >> AA_2nc_queue.txt
    echo "ATL15_write2nc.py @${root}/rel${release}/south/${reg}/input_args_${reg}.txt" >> AA_2nc_queue.txt
done

setup_slurm_run.py --run_name AA_2nc -q AA_2nc_queue.txt --time 03:00:00 -j 8 -e IS2
mv AA_2nc_queue.txt AA_2nc
pushd AA_2nc
sbatch slurm_run.sh
