
release_file=$1

if [ ! -f $release_file ]; then
    echo "release file not found"; exit
fi

release=`grep Release $release_file | sed s/\=/\ / | awk '{print $NF}'`

> arctic_2nc_queue.txt

#for reg in RA IS CS SV CN; do
for reg in SV; do
    echo "ATL14_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel$release/north/${reg}/input_args_${reg}.txt" >> arctic_2nc_queue.txt
    echo "ATL15_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel$release/north/${reg}/input_args_${reg}.txt" >> arctic_2nc_queue.txt
done

rm -r arctic_2nc
setup_slurm_run.py --run_name arctic_2nc -q arctic_2nc_queue.txt --time 02:00:00 -j 1 -e ATL14
mv arctic_2nc_queue.txt arctic_2nc
#pushd arctic_2nc; sbatch slurm_run.sh; popd


#> GL_2nc_queue.txt
#reg='GL'
#echo "ATL14_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel$release/north/${reg}/input_args_${reg}.txt" >> GL_2nc_queue.txt
#echo "ATL15_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel$release/north/${reg}/input_args_${reg}.txt" >> GL_2nc_queue.txt


#rm -r GL_2nc
#setup_slurm_run.py --run_name GL_2nc -q GL_2nc_queue.txt --time 02:00:00 -j 5 -e ATL14
#mv GL_2nc_queue.txt GL_2nc
#pushd GL_2nc ; sbatch slurm_run.sh; 




