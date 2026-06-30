# make the index for ATL11 indexing for North and South
make_ATL11_index.py @default_args/discover.txt @default_args/latest_release.txt

setup_slurm_run.py -r index_ATL11 -q index_queue.txt -e IS2 --lines_per_task 50 -j 2 -t 01:00:00
bash post_queue.txt
