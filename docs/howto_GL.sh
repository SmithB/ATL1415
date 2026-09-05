# DISCOVER/SLURM variant.  The MAAP/DPS counterpart is docs/howto_MAAP_GL.sh.
# This file remains the production path; the MAAP one is still tentative.

# Run the indexing commands from howto_ATL11.txt

# Update default_args/latest_release.txt and default_args/GL_latest.txt.  Each should be a symbolic link to a file containing release-specific details

# setup Greenland:
setup_ATL1415_region.py default_args/discover.txt default_args/quarterly.txt default_args/latest_release.txt default_args/GL_latest.txt --Hemisphere=1
ln -s /discover/nobackup/projects/icesat2/ATL14_processing/rel006/north/GL .


# make the queue
make_ATL1415_queue.py prelim GL/input_args_GL.txt

# setup the slurm run:
setup_ATL1415_run.py --run_name GL"_prelim" -q 1415_queue_GL"_prelim.txt" --time 04:00:00 -j 7 -e ATL14

# generate the tile-size report  # this is now taken care of in ATL11_to_ATL15.py
#make_field_size_report.py /discover/nobackup/projects/icesat2/ATL14_processing/rel006/north/GL/prelim
# now look at the tile sizes using the check_tiles.ipynb notebook in ATL1415. 

# setup matched
make_ATL1415_queue.py matched /discover/nobackup/projects/icesat2/ATL14_processing/rel006/north/GL/input_args_GL.txt
setup_ATL1415_run.py --run_name GL_matched -q 1415_queue_GL_matched.txt --time 04:00:00 -j 8 -e ATL14

# setup the mosaic run:
make_mosaic_jobs.py -b /discover/nobackup/projects/icesat2/ATL14_processing/rel006/north/GL -rr GL

# run the tonc

echo "ATL14_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel006/north/GL/input_args_GL.txt" > GL_2nc_queue.sh
echo "ATL15_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel006/north/GL/input_args_GL.txt" > GL_2nc_queue.sh
setup_slurm_run.py --run_name GL"_2nc" -q GL_2nc_queue.sh --time 04:00:00 -j 8 -e ATL14 ; pushd GL_2nc ; sbatch slurm_run.sh; popd

#-----------------------------------
# monthly:
# need to choose an ATL14 reference file, then execute:
setup_ATL1415_region.py default_args/discover.txt default_args/monthly.txt default_args/latest_release.txt default_args/GL_latest.txt --Hemisphere=1 --ATL14_reference_file /discover/nobackup/projects/icesat2/ATL14_processing/rel005_0329/north/GL/ATL14_GL_0329_100m_005_02.nc

ln -s /discover/nobackup/projects/icesat2/ATL14_processing/rel006/north_monthly/GL ./GL_monthly

make_ATL1415_queue.py prelim GL_monthly/input_args_GL.txt
setup_ATL1415_run.py --run_name GLm_prelim -q 1415_queue_GL_prelim.txt --time 04:00:00 -j 7 -e ATL14

make_ATL1415_queue.py matched GL_monthly/input_args_GL.txt
setup_ATL1415_run.py --run_name GLm_matched -q 1415_queue_GL_matched.txt --time 04:00:00 -j 7 -e ATL14


make_mosaic_jobs.py @GL_monthly/input_args_GL.txt


echo "ATL14_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel006/north_monthly/GL/input_args_GL.txt" > GLm_2nc_queue.sh
echo "ATL15_write2nc.py @/discover/nobackup/projects/icesat2/ATL14_processing/rel006/north_monthly/GL/input_args_GL.txt" >> GLm_2nc_queue.sh
setup_slurm_run.py --run_name GLm_2nc -q GLm_2nc_queue.sh --time 04:00:00 -j 8 -e ATL14 ; pushd GLm_2nc ; sbatch slurm_run.sh; popd
