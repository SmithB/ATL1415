# howto_MAAP_GL.sh -- Greenland on MAAP (per-tile solves on DPS)
#
# ############################################################################
# ##  REWRITTEN 2026-09-19 from what the Iceland run learned.               ##
# ##  GL HAS NEVER RUN ON MAAP.  EVERY STEP IS TENTATIVE until it has.      ##
# ##  The PROCEDURE is the one IS ran end to end, quarterly and monthly     ##
# ##  (docs/howto_MAAP_arctic.sh, which explains each step and cites the   ##
# ##  IS record); this file carries it for one region and marks what GL    ##
# ##  does that IS never did.  Revise the tags as GL goes through.         ##
# ############################################################################
#
# The discover/SLURM variant is docs/howto_GL.sh, which is still the
# production path and is NOT replaced by this file.
#
# Tags:  [ADE] / [DPS];  [OK on IS] the same command ran for Iceland;
#   [UNTESTED] never run in this form;  [NEW FOR GL] exercises something no
#   MAAP job has.  Steps are numbered 0-20, parallel to the arctic howto's
#   ("GL step 7" is "arctic step 7" for Greenland).
#
# WHAT GL DOES THAT IS NEVER DID:
#   a. A GEOTIFF MASK (GreenlandIceMask_2018.1_2026.0_100m_v4.1.tif), read
#      through /vsis3 as a gridded mask.  IS used a vector .db.  The /vsis3
#      read itself was verified on the GL geotiffs 2026-09-06 (the pyTMD
#      AWS_NO_SIGN_REQUEST fix, Transition_to_maap.md); a GL SOLVE never has.
#   b. TIDES: --tide_mask_file and --tide_model=Gr1km-v2, read anonymously
#      from s3://pytmd.  AA's transect proved tides on DPS with CATS2008;
#      Gr1km-v2 has never been read by a job.
#   c. ERROR-SCALING MAPS: --E_d3zdx2dt_scale_file and --E_d2z0dx2_file.
#   d. SCALE: 1483 tile centers against IS's 29.  The submitter, collector
#      and fetcher were built and tested at 29; they make API or S3 calls per
#      job, so expect minutes-to-hours of ADE time per pass, not seconds.
#   e. MEMORY: unmeasured.  IS peaked at 9.5 GiB (16 GiB queue); AA's densest
#      tiles near the pole reached 21.5 GiB (32 GiB queue).  Memory tracks
#      N_fit.  The smoke tiles (step 3) decide the queue.
# All GL masks named in GL_0331.txt are staged (checked 2026-09-19).

conda activate ATL14
cd ~/git_repos/ATL1415
repo=$PWD
reg=GL
rel_file=default_args/latest_release.txt
rel=$(grep '^--Release=' $rel_file | cut -d= -f2)     # 006
cyc=$(grep '^--cycles=' $rel_file | cut -d= -f2)      # 0332
ver=$(grep '^--version=' $rel_file | cut -d= -f2)     # 02
tspan=$(grep '^-t=' $rel_file | cut -d= -f2)          # 2018.75,2026.5
ATL14_root=/home/jovyan/ATL14_processing
s3_root=s3://maap-ops-workspace/ben_smith
ledgers=$ATL14_root/maap_ledgers
runs=$ATL14_root/runs
tile_list=$repo/ATL1415/resources/GL/40km_tile_list.txt     # 1483 centers
# paths [_monthly] -- as the arctic howto's, for GL
paths () {
    region_dir=$ATL14_root/rel$rel/north$1/GL
    s3_run=$s3_root/ATL1415/run_args/rel$rel/north$1/GL
    s3_out=$s3_root/ATL14_processing/rel$rel/north$1/GL
    tag=GL_rel${rel}_${cyc}$1
    L=$ledgers/GL_${cyc}$1
}
paths


# ===========================================================================
# 0. [ADE] [OK on IS]  The DPS build is the commit you mean to run.  (arctic 0)
# ===========================================================================
/srv/conda/envs/notebook/bin/python register_algorithm.py --dry-run   # "push check: OK"
scripts/maap/check_build_id.py $s3_root/ATL1415/run_args/rel006/north/IS/input_args_IS.txt \
    maap-dps-worker-16gb                  # VERDICT: MATCH, maap_pgt=set


# ===========================================================================
# 1. [ADE] [OK]  Point the release symlinks at this release.
# ===========================================================================
ln -sf rel_006_0332.txt default_args/latest_release.txt
ln -sf GL_0331.txt      default_args/GL_latest.txt
# GL_latest.txt does not exist until this makes it.  GL_0331.txt names only
# masks, the tide model and the error-scaling maps -- nothing tied to the
# ATL11 generation -- so it serves cycles 03-32 unchanged, its name
# notwithstanding.  The mask covers 2018.1-2026.0.


# ===========================================================================
# 2. [ADE] [UNTESTED for GL]  Compose and publish the args.  (arctic 2)
# ===========================================================================
setup_ATL1415_region.py default_args/MAAP_dps.txt $rel_file \
    default_args/GL_latest.txt default_args/quarterly.txt --Hemisphere=1
aws s3 cp $region_dir/input_args_GL.txt $s3_run/
# CHECK the GL-specific lines came through as bucket URIs, and the previous
# product as a CMR search (--previous_product=005_0329, no
# --previous_product_top):
grep -E '^(--mask_file|--tide_mask_file|--tide_model|--E_d3zdx2dt_scale_file|--E_d2z0dx2_file|--previous_product)' \
    $region_dir/input_args_GL.txt


# ===========================================================================
# 3. [DPS] [NEW FOR GL]  Two smoke tiles: the densest, and a tide tile.
# ===========================================================================
# Chosen 2026-09-19 against the list and the staged tide mask:
#   200000 -1880000  E200_N-1880, 19 km from Summit: dense interior, no
#                    floating ice (tide mask 0.00).  Sizes the queue.
#   480000 -1040000  E480_N-1040, the 79N glacier tongue: 35% floating by the
#                    tide mask -- a grounding-line tile, the kind that tests
#                    the mask boundary (AA transect).  Exercises Gr1km-v2.
# (Petermann, E-280_N-960, is 32% floating: the alternative.)
# The xy file goes under $ledgers, NOT in the checkout.
printf '200000 -1880000\n480000 -1040000\n' > $ledgers/GL_smoke_xy.txt
scripts/maap/submit_MAAP_jobs.py --xy_file $ledgers/GL_smoke_xy.txt \
    --step prelim --args_url $s3_run/input_args_GL.txt \
    --tile_prefix $s3_out --queue maap-dps-worker-32gb \
    --tag ${tag}_smoke --ledger ${L}_smoke_jobs.csv --dry-run
# (then without --dry-run)
scripts/maap/collect_jobs.py ${L}_smoke_jobs.csv
# GATES: both successful; N_XO > 0; commit = the step 0 build; on E480 the
# log shows the tide read (Gr1km-v2) without error.  THE QUEUE: if the
# interior tile peaks under ~10 GiB, fan out on -16gb; otherwise -32gb.  The
# smoke itself runs on -32gb so it cannot OOM before it has measured.


# ===========================================================================
# 4. [DPS] [UNTESTED for GL]  Fan out.  (arctic 4)
# ===========================================================================
nohup scripts/maap/submit_MAAP_jobs.py --tile_list $tile_list \
    --step prelim --args_url $s3_run/input_args_GL.txt \
    --tile_prefix $s3_out --queue maap-dps-worker-32gb \
    --tag ${tag}_prelim --ledger ${L}_prelim_jobs.csv \
    --max_in_flight 200 > ${L}_prelim_submit.log 2>&1 &
# --queue from step 3.  --max_in_flight holds the unfinished jobs at N,
# polling until one finishes, so the submitter runs for hours: hence nohup,
# with its output in a log beside the ledger.  N=200 is a RECOMMENDATION, not
# a measured limit -- IS's 29 all ran at once; nothing is known at 1483.
# The ledger is written and flushed row by row, so an interrupted submitter
# still leaves a readable ledger (and a re-run needs a NEW --ledger and only
# the centers not yet submitted).
# NEVER register while this runs (arctic f).


# ===========================================================================
# 5-6. [ADE] [OK on IS; UNTESTED at GL scale]  Watch, fetch, check.
# ===========================================================================
scripts/maap/collect_jobs.py ${L}_prelim_jobs.csv > ${L}_prelim_collect.txt
scripts/maap/fetch_tiles.py  ${L}_prelim_jobs.csv $region_dir --step prelim
scripts/check_field_sizes.py $region_dir/prelim @$region_dir/input_args_GL.txt
# As arctic 5-6.  Over 1483 rows both scripts are slow (several calls per
# job) -- save their output rather than scroll it.  fetch_tiles saves the
# no-data centers to $region_dir/prelim/no_data_tiles.txt (step 20); failed
# jobs are listed under NOT FETCHED and retried as arctic 5 says.


# ===========================================================================
# 7-8. [DPS] [OK on IS; UNTESTED for GL]  Matched.  (arctic 7-8)
# ===========================================================================
nohup scripts/maap/submit_MAAP_jobs.py --tile_list $tile_list \
    --step matched --args_url $s3_run/input_args_GL.txt \
    --tile_prefix $s3_out --queue maap-dps-worker-32gb \
    --tag ${tag}_matched --ledger ${L}_matched_jobs.csv \
    --max_in_flight 200 > ${L}_matched_submit.log 2>&1 &
# Skips, by name, centers with no prelim tile.  Then:
scripts/maap/collect_jobs.py ${L}_matched_jobs.csv > ${L}_matched_collect.txt
scripts/maap/fetch_tiles.py  ${L}_matched_jobs.csv $region_dir --step matched
scripts/check_field_sizes.py $region_dir/matched @$region_dir/input_args_GL.txt


# ===========================================================================
# 9. [ADE] [OK on IS; UNTESTED at GL size]  Mosaic.  (arctic 9)
# ===========================================================================
cd $runs
make_mosaic_jobs.py -b $region_dir -rr GL -t $tspan -e ATL14 \
    --run_name GL_${cyc}_mosaic @$repo/default_args/quarterly.txt
cd GL_${cyc}_mosaic
seq 1 $(ls queue | wc -l) | xargs -P 4 -I{} env SLURM_ARRAY_TASK_ID={} bash slurm_run.sh
check_mosaic_outputs.py $runs/GL_${cyc}_mosaic --values
cd $repo
# -P 4, not IS's 12: the z0 task alone is a 100 m grid over Greenland, orders
# of magnitude larger than Iceland's 0.25 GiB tasks, on an ADE with 124 GiB.
# UNTIMED.  Watch the first tasks' memory (scripts/run_with_rusage.py wraps a
# command and reports its peak) and raise -P if there is room.


# ===========================================================================
# 10. [ADE] [OK on IS; UNTESTED for GL]  netCDF, and the comparison.
# ===========================================================================
mkdir -p $runs/GL_${cyc}_nc && cd $runs/GL_${cyc}_nc
ATL14_write2nc.py @$region_dir/input_args_GL.txt > ATL14.log 2>&1
ATL15_write2nc.py @$region_dir/input_args_GL.txt > ATL15.log 2>&1
cd $repo
# No INVALID line in either log; XO rows NOT_SET in four attributes only.
# COMPARE WITH rel005 (Ben's bar: no >10 m errors, no major gaps) before
# monthly builds on it -- method in plan_cycles_03_32.sh T8 "I9g6".
# BROWSE PLOTS, as discover does -- [UNTESTED on MAAP; not planned for IS]:
#   ATL14_browse_plots.py @$region_dir/input_args_GL.txt
#   ATL15_browse_plots.py @$region_dir/input_args_GL.txt


# ===========================================================================
# 11-18. [ADE+DPS] [OK on IS; UNTESTED for GL]  Monthly.  (arctic 11-18)
# ===========================================================================
# Publish the quarterly products (the ATL14 is the reference, one file, read
# by URI):
for f in $region_dir/ATL1[45]_GL_${cyc}_*_${rel}_${ver}.nc; do aws s3 cp $f $s3_out/; done
ref=$s3_out/ATL14_GL_${cyc}_100m_${rel}_${ver}.nc
setup_ATL1415_region.py default_args/MAAP_dps.txt $rel_file \
    default_args/GL_latest.txt default_args/monthly.txt --Hemisphere=1 \
    --ATL14_reference_file=$ref
paths _monthly
aws s3 cp $region_dir/input_args_GL.txt $s3_run/
# Then steps 3-10 again with these paths (tag and ledger carry _monthly):
#   - smoke E200_N-1880: its monthly N_fit must be about its quarterly N_fit
#     -- the only sign of an unreadable reference;
#   - prelim and matched from the same $tile_list;
#   - mosaic with @$repo/default_args/monthly.txt and --run_name
#     GL_${cyc}_monthly_mosaic; INSPECT the inferred lags (1,3,...,84) first;
#   - ATL15_write2nc.py ONLY, into $runs/GL_${cyc}_monthly_nc;
#   - compare with quarterly on RAW delta_h (arctic 17).
# Then publish the monthly ATL15 files to $s3_out.


# ===========================================================================
# 19. [ADE] [SUGGESTION, NO SOFTWARE]  Annotate the build history.  (ogc O12b)
# ===========================================================================


# ===========================================================================
# 20. [ADE] [OK on IS]  Take the no-data centers out of the list.  (arctic 20)
# ===========================================================================
cat $ATL14_root/rel$rel/north{,_monthly}/GL/prelim/no_data_tiles.txt 2>/dev/null \
    | sort -u > $ledgers/GL_no_data.txt
grep -vxFf $ledgers/GL_no_data.txt $tile_list > $ledgers/t && mv $ledgers/t $tile_list
git commit -m "Drop no-data centers from the GL tile list" $tile_list && git push
