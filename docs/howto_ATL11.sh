# make the index for ATL11 indexing for North and South
# EDIT THIS FOR RELEASE AND VERSION (3 occurrences)
scripts/make_ATL11_index /discover/nobackup/projects/icesat2/ATL14_processing/ 007_0330_v03 /discover/nobackup/bjelley/ATL11_processing/Arctic_007_cycle_03_30_v03/007 /discover/nobackup/bjelley/ATL11_processing/Antarctic_007_cycle_03_30_v03/007

setup_slurm_run.py -r index_ATL11 -q index_queue.txt -e IS2 --lines_per_task 50 -j 2 -t 01:00:00
bash post_queue.txt
