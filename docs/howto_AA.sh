# see howto_ATL11.sh for indexing

# setup Antarctica:
setup_ATL1415_region.py default_args/discover.txt default_args/latest_release.txt default_args/AA_latest.txt default_args/quarterly.txt --Hemisphere=-1

# make the queue for AA_north
make_ATL1415_queue.py prelim /discover/nobackup/projects/icesat2/ATL14_processing/rel006/south/AA/input_args_AA.txt --min_xy 360000 

# make the slurm run for AA_north
setup_ATL1415_run.py --run_name AA_prelim_north -q 1415_queue_AA_prelim.txt --time 04:00:00  -e ATL14

# make the field size reports:
#scripts/make_field_size_report.py /discover/nobackup/projects/icesat2/ATL14_processing/rel006/south/AA/prelim

# make the AA_matched queue for north:
scripts/make_ATL1415_queue.py matched /discover/nobackup/projects/icesat2/ATL14_processing/rel006/south/AA/input_args_AA.txt --min_xy 360000
setup_ATL1415_run.py --run_name AA"_matched_north" -q 1415_queue_AA"_matched.txt" --time 04:00:00  -e ATL14 --lines_per_task 4

# setup Antarctica south:
make_ATL1415_queue.py prelim /discover/nobackup/projects/icesat2/ATL14_processing/rel006/south/AA_44km/input_args_AA_44km.txt --max_xy 440000
setup_ATL1415_run.py --run_name AA"_prelim_south" -q 1415_queue_AA"_prelim.txt" --time 06:00:00  -e ATL14

# setup Antarctica matched south
make_ATL1415_queue.py matched /discover/nobackup/projects/icesat2/ATL14_processing/rel006/south/AA_44km/input_args_AA_44km.txt --max_xy 440000
setup_ATL1415_run.py --run_name AA"_matched_south" -q 1415_queue_AA"_matched.txt" --time 04:00:00  -e ATL14

# mosaic the northern tiles:
make_200km_tiles.py  ~/shared/ATL14_processing/rel006/south/AA AA -t 2018.75,2026.5
# then run the slurm_run in tile_run_AA
scripts/make_200km_tiles.py  ~/shared/ATL14_processing/rel006/south/AA_44km AA --name AA_south --W 44000 --spacing 40000 -t 2018.75,2026.5

# setup the four antarctic sectors:
setup_AA_sectors.py ~/shared/ATL14_processing/rel006/south

# then mosaic the tiles:
for sector in A1 A2 A3 A4; do make_200km_to_mosaic_jobs.py -b /discover/nobackup/projects/icesat2/ATL14_processing/rel006/south/$sector -rr $sector -t 2018.75,2026.5; done

# make the netCDFs
scripts/run_antarctic_tonc.sh default_args/latest_release.txt default_args/discover.txt

###############################
# monthly
###############################
# setup
setup_ATL1415_region.py default_args/discover.txt default_args/latest_release.txt default_args/AA_latest.txt default_args/monthly.txt --Hemisphere=-1 --ATL14_reference_file "/discover/nobackup/projects/icesat2/ATL14_processing/rel005_0329/south/A*/ATL14_*_0329_100m_005_02.nc"
ln -s /discover/nobackup/projects/icesat2/ATL14_processing/rel006/south_monthly/AA ./AA_monthly

# prelim queue
make_ATL1415_queue.py prelim AA_monthly/input_args_AA.txt
setup_ATL1415_run.py --run_name AAm_prelim -q 1415_queue_AA_prelim.txt --time 04:00:00  -e ATL14

make_ATL1415_queue.py matched AA_monthly/input_args_AA.txt
setup_ATL1415_run.py --run_name AAm_matched -q 1415_queue_AA_matched.txt --time 04:00:00  -e ATL14

# tiles:

make_200km_tiles.py  ~/shared/ATL14_processing/rel006/south_monthly/AA AA -t 2018.75,2026.5 @/discover/nobackup/projects/icesat2/ATL14_processing/rel006/south_monthly/AA/input_args_AA.txt

setup_AA_sectors.py ~/shared/ATL14_processing/rel006/south_monthly --near_pole_radius 0

for sector in A1 A2 A3 A4; do make_200km_to_mosaic_jobs.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel006/south_monthly/$sector/input_args_$sector.txt  -rr $sector -t 2018.75,2026.5; done

run_antarctic_tonc.sh default_args/latest_release.txt default_args/discover.txt
