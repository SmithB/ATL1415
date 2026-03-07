


> AA_2nc_queue.txt

for reg in A1 A2 A3 A4; do
    echo "ATL14_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel005/south/${reg}/input_args_${reg}.txt" >> AA_2nc_queue.txt
    echo "ATL15_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel005/south/${reg}/input_args_${reg}.txt" >> AA_2nc_queue.txt
done

setup_slurm_run.py --run_name AA_2nc -q AA_2nc_queue.txt --time 03:00:00 -j 8 -e IS2
mv AA_2nc_queue.txt AA_2nc
pushd AA_2nc
sbatch slurm_run.sh
