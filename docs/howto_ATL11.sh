# DISCOVER/SLURM ONLY.  There is deliberately NO MAAP counterpart: ATL11 indexing
# stays on discover (Q1/Q13).  The index is then staged to the bucket -- see
# docs/howto_MAAP_staging.sh step S3, which is the only per-release handoff.

# make the index for ATL11 indexing for North and South
make_ATL11_index.py @default_args/discover.txt @default_args/latest_release.txt

setup_slurm_run.py -r index_ATL11 -q index_queue.txt -e IS2 --lines_per_task 50 -j 2 -t 01:00:00
bash post_queue.txt
