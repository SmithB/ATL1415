# Run the indexing commands from howto_AA.sh
# make the index for ATL11 indexing for North:
#scripts/make_ATL11_index /discover/nobackup/projects/icesat2/ATL14_processing/ 006_0326_v11 /discover/nobackup/bjelley/ATL11_processing/Arctic_006_cycle_03_26_v11 None

# note- need to make new versions of these three files:
# default_args/rel_004_0325.txt default_args/north.txt default_args/GL_0325.txt  
# check each for the cycle numbers
# setup Greenland:
setup_ATL1415_region.py default_args/discover.txt default_args/rel_005_0329.txt default_args/north.txt default_args/GL_0329.txt

# make the queue
make_ATL1415_queue.py prelim /discover/nobackup/projects/icesat2/ATL14_processing/rel005/north/GL/input_args_GL.txt

# setup the slurm run:
setup_slurm_run.py --run_name GL"_prelim" -q 1415_queue_GL"_prelim.txt" --time 04:00:00 -j 7 -e ATL14

# generate the tile-size report
make_field_size_report.py /discover/nobackup/projects/icesat2/ATL14_processing/rel005/north/GL/prelim GL prelim
# now look at the tile sizes using the check_tiles.ipynb notebook in ATL1415. 

# setup matched
make_ATL1415_queue.py matched /discover/nobackup/projects/icesat2/ATL14_processing/rel005/north/GL/input_args_GL.txt
setup_slurm_run.py --run_name GL"_matched" -q 1415_queue_GL"_matched.txt" --time 04:00:00 -j 8 -e ATL14

# setup the mosaic run:
scripts/make_mosaic_jobs /discover/nobackup/projects/icesat2/ATL14_processing/rel005/north/GL

# run the tonc
ATL14_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel005/north/GL/input_args_GL.txt
ATL15_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel005/north/GL/input_args_GL.txt
