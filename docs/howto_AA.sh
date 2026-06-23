# see howto_ATL11.sh for indexing

# setup Antarctica:
setup_ATL1415_region.py default_args/discover.txt default_args/rel_005_0329.txt default_args/south.txt default_args/AA_0329.txt

# make the queue for AA_north
make_ATL1415_queue.py prelim /discover/nobackup/projects/icesat2/ATL14_processing/rel005/south/AA/input_args_AA.txt --min_xy 360000 

# make the slurm run for AA_north
setup_slurm_run.py --run_name AA_prelim_north -q 1415_queue_AA_prelim.txt --time 04:00:00  -e ATL14

# make the field size reports:
scripts/make_field_size_report.py /discover/nobackup/projects/icesat2/ATL14_processing/rel005/south/AA/prelim AA prelim

# make the AA_matched queue for north:
scripts/make_ATL1415_queue.py matched /discover/nobackup/projects/icesat2/ATL14_processing/rel005/south/AA/input_args_AA.txt --min_xy 360000
setup_slurm_run.py --run_name AA"_matched_north" -q 1415_queue_AA"_matched.txt" --time 04:00:00  -e ATL14 --lines_per_task 4

# setup Antarctica south:
make_ATL1415_queue.py prelim /discover/nobackup/projects/icesat2/ATL14_processing/rel005/south/AA_44km/input_args_AA_44km.txt --max_xy 440000
setup_slurm_run.py --run_name AA"_prelim_south" -q 1415_queue_AA"_prelim.txt" --time 06:00:00  -e ATL14

# setup Antarctica matched south
make_ATL1415_queue.py matched /discover/nobackup/projects/icesat2/ATL14_processing/rel005/south/AA_44km/input_args_AA_44km.txt --max_xy 440000
setup_slurm_run.py --run_name AA"_matched_south" -q 1415_queue_AA"_matched.txt" --time 04:00:00  -e ATL14

# mosaic the northern tiles:
make_200km_tiles.py  ~/shared/ATL14_processing/rel005/south/AA AA --last_lag 28
# then run the slurm_run in tile_run_AA
scripts/make_200km_tiles.py  ~/shared/ATL14_processing/rel005/south/AA_44km AA --name AA_south --W 44000 --spacing 40000 --last_lag 28

# now run git_repos/notebooks/make_and_check_links.ipynb

# then mosaic the tiles:
for sector in A1 A2 A3 A4; do bash scripts/make_200km_to_mosaic_jobs.sh /discover/nobackup/projects/icesat2/ATL14_processing/rel005/south/$sector; done

# make the netCDFs
scripts/run_antarctic_tonc.sh

