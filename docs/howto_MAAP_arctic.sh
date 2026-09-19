# howto_MAAP_arctic.sh -- the arctic regions (RA IS CN CS SV) on MAAP
#
# ############################################################################
# ##  STATUS 2026-09-18: IS HAS RUN STEPS 0-10 at cycles 03-32 with       ##
# ##  COMPLETE lineage (plan_cycles_03_32.sh T5-T8), AND THE MONTHLY        ##
# ##  PRODUCT, step 10b (plan_monthly_on_maap.sh).  Steps 5 and 11 are      ##
# ##  still NEEDS CODE.  The banner below is the 2026-09-05 original.       ##
# ##                                                                        ##
# ##  TENTATIVE.  Written 2026-09-05 BEFORE any of it has been run end to   ##
# ##  end -- no ATL1415 tile has been solved on DPS yet.  This is the plan,  ##
# ##  not a record of a successful run.  Expect steps to move, split and    ##
# ##  change as testing advances; revise this file as that happens.         ##
# ############################################################################
#
# The discover/SLURM variant is docs/howto_arctic.sh, which is still the
# production path and is NOT replaced by this file.
#
# READ docs/howto_MAAP_GL.sh FIRST.  GL is the reference workflow; this file
# marks only what the arctic regions do differently.  Same tags:
#   [ADE] / [DPS]   [OK] [UNTESTED] [NEEDS CODE: x]
# Steps are numbered 1..12.
#
# Prerequisite: docs/howto_MAAP_staging.sh S1-S6; step 5 is gated on S7.
#
# WHY THIS IS THE ONE TO TEST FIRST (Q21):
#   ICELAND (IS) IS THE TEST REGION.  It is small enough to fan out under the
#   ~10 jobs/hr public-queue throttle, it is a real region rather than a toy,
#   and the crossover work was verified against an Iceland box (13 tiles).
#   Do the whole of this file for IS alone before running the five-region loop.
#
# WHAT MAKES THE ARCTIC DIFFERENT:
#   a. FIVE REGIONS, not one.  The discover workflow drives them through
#      scripts/run_arctic_{prelim,matched,mosaic,to_nc}.sh, which loop over
#      "RA IS CN CS SV".  Those greps already read --Release/--ATL14_root/
#      --cycles/--version out of MAAP_dps.txt unchanged; ONLY the
#      `setup_slurm_run.py ...; sbatch` tail has to change.
#      PUT THE MAAP VARIANTS IN scripts/maap/ so the SLURM originals stay put.
#   b. THE MASKS ARE VECTOR .db FILES, not geotiffs (Q12).  masks/RGI_reduced/
#      is on the bucket as real files, and GDAL's SQLite driver reads a .db
#      over /vsis3.  ATL1415/make_mask_from_vector.py called ogr.Open() on the
#      raw URI, which failed; as of 2026-09-05 it routes through
#      pc.io_utils.as_gdal_path().  EXERCISED AGAINST THE BUCKET 2026-09-06:
#      it works -- but only after a SECOND bug had to be fixed, pyTMD's
#      AWS_NO_SIGN_REQUEST (step 2).  That one was not arctic-specific: it
#      broke every mask read in every region.
#   c. the region args files have no cycles suffix: default_args/{RA,IS,CN,CS,SV}.txt.
#   d. all five are northern hemisphere, so rel006/north/<REG>/.

conda activate ATL14
cd ~/git_repos/ATL1415
reg=IS                       # <-- do IS alone first (Q21)
region_dir=/home/jovyan/ATL14_processing/rel006/north/$reg
s3_run=s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/north/$reg
s3_out=s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north/$reg


# ===========================================================================
# 0. [ADE] [OK 2026-09-11]  Rebuild the DPS image if any code has changed.
# ===========================================================================
# DPS DOES NOT RUN THIS WORKING COPY: a build clones on_s3 from GitHub.  Push,
# register, and check the image -- staging S5 and S5b carry the rule, and
# howto_MAAP_ogc O3-O6 the detail:
/srv/conda/envs/notebook/bin/python register_algorithm.py           # refuses unpushed work
/srv/conda/envs/notebook/bin/python scripts/maap/check_build_id.py  # must say MATCH


# ===========================================================================
# 1. [ADE] [OK -- 0332 since 2026-09-17]  Point the release symlink at this release.
# ===========================================================================
ln -sf rel_006_0332.txt default_args/latest_release.txt
# rel_006_0331.txt is kept as the record of the first IS run.  CMR no longer
# lists that ATL11 generation, so 0331 cannot be solved (plan_cycles_03_32.sh).


# ===========================================================================
# 2. [ADE] [OK 2026-09-06]  CHECK THE .db MASK READS FROM THE BUCKET. <-- do this first
# ===========================================================================
# Q12's one-line fix.  If this raises "No such file or directory", nothing else
# in this file will work, and the problem is in ATL1415/make_mask_from_vector.py,
# not in the staging.  Run it again on any fresh account or env: it is the
# cheapest check that the bucket, the credentials and the GDAL path all agree.
#
# THE TEST TILE, given by Ben 2026-09-05: x0,y0 = 1260, -2620 km.  43K points --
# a medium-sized dataset, likely on the edge of the ice sheet, so it exercises
# the mask edge rather than a saturated interior tile.  On the grid: centers sit
# at odd multiples of tile_spacing/2, and 1260 = 63 x 20 km, -2620 = -131 x 20 km.
# The same tile is the DPS smoke test -- howto_MAAP_staging.sh S7.
python - <<'EOF'
import pointCollection as pc
from ATL1415.make_mask_from_vector import make_mask_from_vector
db = 's3://maap-ops-workspace/ben_smith/ATL1415/masks/RGI_reduced/06_rgi60_Iceland_reduced.db'
print(pc.io_utils.as_gdal_path(db))          # expect /vsis3/maap-ops-workspace/...
m = make_mask_from_vector(db, W={'x':6e4, 'y':6e4},
                          ctr={'x':1260000., 'y':-2620000.}, spacing=100,
                          srs_proj4=('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 '
                                     '+k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))
print(m.z.shape, m.z.sum())    # expect a nonzero count, and NOT the full 601x601 --
                               # an edge tile should be part ice, part not
EOF
# RESULT 2026-09-06: PASSES.  601x601, 72638 ice cells = 20.1% -- a real edge
# tile, as predicted.  It did NOT pass on the first attempt, and what stopped
# it was not the .db path at all:
#
#   pyTMD v3.0.9 sets AWS_NO_SIGN_REQUEST=YES in os.environ when pyTMD.io is
#   imported, so that its own reads of the public s3://pytmd stores are
#   anonymous.  The variable is PROCESS-WIDE and GDAL honours it for every
#   /vsis3 read, so every mask, geoid and tide mask read from
#   s3://maap-ops-workspace came back HTTP 403 -- surfacing as ogr.Open()
#   returning None here, and as from_geotif() returning an object with no .z on
#   the GL/AA geotiffs, which were verified to fail the same way.
#   ATL11_to_ATL15.py:47 imports pyTMD at module level, so EVERY ATL1415
#   process had it set: this would have broken every DPS tile job, not just the
#   arctic ones.  Fixed 2026-09-06 in ATL1415/__init__.py, which pops the
#   variable after the pyTMD import that its own first line triggers.  Our tide
#   reads do not need it -- ATL1415/tides.py passes anon= to s3fs explicitly,
#   and s3fs ignores the GDAL variable.  STILL TO DO: report upstream to pyTMD.
#
#   If a process has already made an anonymous read before the variable is
#   cleared, GDAL caches that decision and gdal.VSICurlClearCache() is needed as
#   well.  Clearing at import time, as the fix does, avoids that entirely.


# ===========================================================================
# 3. [ADE] [OK 2026-09-07]  Compose the args file.   (as GL step 2)
# ===========================================================================
setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
    default_args/$reg.txt default_args/quarterly.txt --Hemisphere=1
# RUN 2026-09-06, and RECOMPOSED 2026-09-07 after the Q27 W1/W3/W4 fixes
# (b293807) changed what setup emits: writes
# /home/jovyan/ATL14_processing/rel006/north/IS/input_args_IS.txt, now 970
# bytes rather than the 1093 of the first composition -- shorter because the
# two /discover/... previous-product paths are gone.  It needs --ATL14_root to
# exist first -- see staging S1.  RECOMPOSE AFTER ANY CHANGE TO MAAP_dps.txt OR
# TO setup_ATL1415_region.py, and republish (step 4): the bucket copy is what
# DPS reads, and nothing reconciles the two.
#
# EVERY CLOUD INPUT IN THE COMPOSED FILE IS AN s3:// URI OR A CMR SEARCH, as
# intended: --ATL11_index, --tide_directory, --geoid_file, --mask_file (the
# Iceland .db), and -- since 2026-09-07 -- the previous product, which is now
#   --previous_product_earthaccess
#   --previous_product=005_0329
# in place of the two /discover/... paths this file used to warn about (Q27
# W1/W4 are fixed; --previous_product_top is dropped in cloud mode).  The
# previous product WILL therefore be read on the smoke test, not silently
# skipped: expect the log to name ATL14_IS_0329_100m_005_02.nc and
# ATL15_IS_0329_01km_005_02.nc, found by a bounding-box CMR search.
#
# THAT EXPECTATION HOLDS ONLY IF THE IMAGE CARRIES b293807 -- step 0.  Against
# the 2026-09-04 build the run does not get that far: without 1023306 the
# Iceland .db read 403s first.  If a smoke-test log shows the old silent
# "no previous product" skip, the image is stale, not the code.
# -b is the one local path left, and run.sh overrides it after the args file on
# purpose, so it is not a problem.


# ===========================================================================
# 4. [ADE] [OK 2026-09-07]  Publish the args file.   (as GL step 3)
# ===========================================================================
aws s3 cp $region_dir/input_args_$reg.txt $s3_run/
# RUN 2026-09-06 for IS: the prefix did not exist beforehand and `aws s3 cp`
# made it.  REPUBLISHED 2026-09-07 03:58 UTC with the recomposed file, and
# CONFIRMED ON THE BUCKET 2026-09-08 -- 970 bytes at
# s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/north/IS/input_args_IS.txt
# -- which is the args_file the smoke test submits (staging S7).  Check the
# size, not just the presence: a 1093-byte object there is the superseded
# composition that still names the discover tree.


# ===========================================================================
# 5. [ADE] [SUPERSEDED 2026-09-18 -- the centers are ATL1415/resources/<region>/40km_tile_list.txt, step 6]  Tile centers.
# ===========================================================================
# Same four blockers as GL step 4, but NOT the 1 km mask work (Q6/Q16):
# CONFIRMED 2026-09-06 that none of it applies here.  make_ATL1415_queue.py's
# .db branch (lines 211-219) does not look for a 1 km sibling at all -- it
# requires <mask_base>_40km.tif and takes the tile centers straight from it.
# Those files are already on the bucket beside the .db masks, built there by
# make_lowres_masks.sh with `gdal_rasterize -tr 40000 40000 -at -tap`, and -at
# (all-touched) is the same "any ice in the cell" rule the Q16 (b) answer
# chose -- so the arctic already does it, one rasterization earlier.
#
# WHAT DOES BITE HERE: line 214 tests os.path.isfile(mask_base+'_40km.tif'),
# which is False for an s3:// URI, so this raises
#   OSError: gridded mask file s3://.../06_rgi60_Iceland_reduced_40km.tif not found
# on a file that is sitting right there.  Line 109 has the same bug for
# --tide_mask_file.  Both want pc.io_utils.path_exists(), not os.path.isfile.
make_ATL1415_queue.py prelim $region_dir/input_args_$reg.txt --xy_out ${reg}_prelim_xy.txt
# IS DODGED THIS (plan_IS_run.sh QI1/I1): its 29 centers were read off the 40 km
# tif by hand and frozen in region_files/IS_prelim_xy.txt.  The bugs are
# unfixed, and every region without a frozen list needs them fixed first.


# ===========================================================================
# 6. [DPS] [OK on IS 2026-09-15, plan_IS_run.sh I2]  Fan out.  (as GL step 5)
# ===========================================================================
# ONE NAMED SMOKE TILE FIRST, then the whole list (the smoke tile re-runs as
# part of it; dedup=False).  --tile_prefix makes run.sh write each tile to
# $s3_out/prelim/ as well as to its dps_output prefix -- the matched step
# (9) cannot run without it.  -16gb was enough: IS prelim peaked at 9.09 GiB.
# Keep ledgers OUTSIDE the checkout: an untracked file makes
# register_algorithm.py refuse.
# THE CENTERS COME FROM THE REGION'S TILE LIST, ATL1415/resources/$reg/
# 40km_tile_list.txt (Ben, 2026-09-18; docs/plan_tile_lists.sh TL2) -- no
# longer region_files/${reg}_prelim_xy.txt, which IS used and which is
# retired for fan-outs.  No-data centers come out of it after each run (step 8).
# [UNTESTED on DPS -- --tile_list is new; dry-run verified on IS]
ledgers=~/ATL14_processing/maap_ledgers
tile_list=ATL1415/resources/$reg/40km_tile_list.txt
scripts/maap/submit_MAAP_jobs.py --tile_list $tile_list \
    --step prelim --args_url $s3_run/input_args_$reg.txt \
    --tile_prefix $s3_out --queue maap-dps-worker-16gb \
    --tag ${reg}_rel006_prelim --ledger $ledgers/${reg}_prelim_jobs.csv
# NEVER commit and re-register while a fan-out is queued: IS's 29 silently
# split across two builds that way (I3).  Only collect_jobs.py's per-tile
# commit column shows it.


# ===========================================================================
# 7. [ADE] [OK on IS 2026-09-16, plan_IS_run.sh I3]  Watch.
# ===========================================================================
# The MAAP analogue of `slurm_run_status.py`.  check_MAAP_jobs.py was never
# written and was not needed: collect_jobs.py reports status, wall clock, peak
# RSS, N_ATL11/N_AT/N_XO, iterations and the build each tile ran.
# Gate: N_XO > 0 on every tile, or crossovers are not being read.
scripts/maap/collect_jobs.py $ledgers/${reg}_prelim_jobs.csv
# A failed job's tile is unrecoverable -- it gets no dps_output prefix.


# ===========================================================================
# 8. [ADE] [OK on IS 2026-09-16, plan_IS_run.sh I4]  Collect.  (as GL step 7)
# ===========================================================================
# The deterministic prefix (Q9/QI4) landed as --tile_prefix, so the tiles are
# at $s3_out/prelim/ too.  IS came down with the ledger-driven fetcher, which
# also reaches tiles solved without --tile_prefix:
scripts/maap/fetch_tiles.py $ledgers/${reg}_prelim_jobs.csv $region_dir --step prelim --dry-run
scripts/maap/fetch_tiles.py $ledgers/${reg}_prelim_jobs.csv $region_dir --step prelim
# Check the tiles' field-size reports: dz/dz of the shape -W, -g and -t give,
# prelim sigma_dz == dz/dz, and a report for every tile.  Exit 0 OK, 1
# problems, 2 the check did not happen.  [OK on IS 2026-09-17, plan I5]
scripts/check_field_sizes.py $region_dir/prelim @$region_dir/input_args_$reg.txt
# NO-DATA CENTERS ARE SAVED, NOT PRUNED (Ben, 2026-09-19; plan_tile_lists.sh
# TL8).  The fetch above lists every successful job that left no tile --
# the prelim fit or the uncertainty step found no data -- under "NO DATA", and
# adds their names to $region_dir/prelim/no_data_tiles.txt (merged, so retry
# ledgers add to it).  Failed jobs are NOT in it: they show as FAILED under
# NOT FETCHED, and are real faults until their logs are read.
# Nothing changes the tile list during the run -- matched (step 9) simply
# skips centers without a prelim tile.  CLEANUP, AFTER THE RUN: take the
# saved names out of the list and commit, so the next run does not submit
# them.  [OK on IS 2026-09-19: E1020_N-2580, both ledgers]
grep -vxFf $region_dir/prelim/no_data_tiles.txt $tile_list > t && mv t $tile_list
git -C ~/git_repos/ATL1415 commit -m "Drop no-data centers from the $reg tile list" $tile_list


# ===========================================================================
# 9. [DPS] [OK on IS 2026-09-16, plan_IS_run.sh I6-I8]  Matched.  (as GL step 9)
# ===========================================================================
# FROM THE SAME TILE LIST (AM8: the lists drive prelim AND matched), which
# still holds this run's no-data centers -- the cleanup is after the run.
# submit_MAAP_jobs.py lists $s3_out/prelim/ once and SKIPS, by name, every
# center with no prelim tile (QT3), so matched runs exactly where a prelim
# tile exists -- the rule IS followed by hand (region_files/IS_matched_xy.txt).
# Each matched job fetches its own and its neighbours' prelim tiles from
# --tile_prefix; missing neighbours are logged, not fatal.
scripts/maap/submit_MAAP_jobs.py --tile_list $tile_list \
    --step matched --args_url $s3_run/input_args_$reg.txt \
    --tile_prefix $s3_out --queue maap-dps-worker-16gb \
    --tag ${reg}_rel006_matched --ledger $ledgers/${reg}_matched_jobs.csv
scripts/maap/collect_jobs.py $ledgers/${reg}_matched_jobs.csv
scripts/maap/fetch_tiles.py $ledgers/${reg}_matched_jobs.csv $region_dir --step matched
scripts/check_field_sizes.py $region_dir/matched @$region_dir/input_args_$reg.txt
# Matched tiles have NO sigma_dz, by design: the uncertainties come from the
# prelim tiles.  IS matched peaked at 8.99 GiB; memory tracks N_fit.


# ===========================================================================
# 10. [ADE] [OK on IS 2026-09-16/17, plan_IS_run.sh I9]  Mosaic and netCDF.
# ===========================================================================
# NO run_queue_local.sh: make_mosaic_jobs.py's slurm_run.sh is plain bash
# (the #SBATCH lines are comments), so running it with SLURM_ARRAY_TASK_ID
# set does the same queue -> running -> done bookkeeping as on discover.
# Build the run directory OUTSIDE the checkout (untracked files block
# registration).  -e ATL14: the default, IS2, is discover's env.
# IS: 41 tasks, 0.25 GiB each, under a minute at -P 12.
cd ~/ATL14_processing/runs
make_mosaic_jobs.py -b $region_dir -rr $reg -t 2018.75,2026.5 -e ATL14 \
    --run_name ${reg}_mosaic @$HOME/git_repos/ATL1415/default_args/quarterly.txt
cd ${reg}_mosaic
n_tasks=$(ls queue | wc -l)
seq 1 $n_tasks | xargs -P 12 -I{} env SLURM_ARRAY_TASK_ID={} bash slurm_run.sh
# Check the outputs: exit codes are not enough.  Metadata only by default;
# --values also flags all-NaN fields.
check_mosaic_outputs.py ~/ATL14_processing/runs/${reg}_mosaic --values
#
# netCDF: run the two writers directly, in the ADE (IS: 12 s + 18 s; five
# files).  They need no lineage flags -- as of 28b4f72 they never open
# ATL11; the prelim tiles carry each granule's attributes (/meta/lineage,
# since build 61a19af, plan_lineage_at_solve_time.sh).  EXPECT NO INVALID
# warning.  One is a real fault: a tile solved on an older build.  The
# only NOT_SET values are on the ATL11XO rows, in start/end_orbit and
# start/end_region, which the XO granules do not carry -- that is correct.
mkdir -p ~/ATL14_processing/runs/${reg}_nc && cd ~/ATL14_processing/runs/${reg}_nc
ATL14_write2nc.py @$region_dir/input_args_$reg.txt > ATL14.log 2>&1
ATL15_write2nc.py @$region_dir/input_args_$reg.txt > ATL15.log 2>&1
# ATL15 writes all four resolutions (1, 10, 20, 40 km) in one call.
# To open the products with GDAL, use the notebook env: ATL14's GDAL has no
# netCDF or HDF5 plugin.


# ===========================================================================
# 10b. [ADE+DPS] [OK on IS 2026-09-18, plan_monthly_on_maap.sh M0-M10]  Monthly.
# ===========================================================================
# The monthly product (dt 1/12 yr): the SAME four steps as 6-10, with
# default_args/monthly.txt added and the QUARTERLY ATL14 of the same release,
# cycles and version as a reference DEM.  Only ATL15 is written.
# PREREQUISITE: steps 0-10 done for this region -- monthly subtracts the
# quarterly ATL14 from every point.  No rebuild: the solver already takes a
# reference file by URI.  Check it anyway before the smoke tile (step 0).
# IS measured: prelim 496-1071 s at 1.3-4.1 GiB, matched 85-307 s at
# 1.1-3.9 GiB -- ~3x faster and ~2.3x lighter than quarterly (8x fewer
# unknowns, the same data), so maap-dps-worker-16gb is ample.
q_dir=$HOME/ATL14_processing/rel006/north/$reg
s3_q=s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north/$reg
ref=ATL14_${reg}_0332_100m_006_02.nc

# a. Publish the quarterly ATL14 beside the quarterly tiles, and read it back
#    from the URI the way the solver will (plan M1-M2): pc.grid.mosaic()
#    .from_list(['$s3_q/$ref'], group='', bounds=..., fields=['h','h_sigma'])
#    against the local file.  Compare with equal_nan=True -- plain
#    array_equal says False on NaNs, which is not a difference.
aws s3 cp $q_dir/$ref $s3_q/

# b. Compose and publish the monthly args (plan M4-M5).  They differ from the
#    quarterly args in exactly four lines: -g, --dzdt_lags,
#    --ATL14_reference_file and -b.  No --hemi_suffix line is written; it
#    only sets the directory, rel006/north_monthly/<reg>.
setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
    default_args/$reg.txt default_args/monthly.txt --Hemisphere=1 \
    --ATL14_reference_file=$s3_q/$ref
m_dir=$HOME/ATL14_processing/rel006/north_monthly/$reg
s3_run_m=s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/north_monthly/$reg
s3_out_m=s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north_monthly/$reg
aws s3 cp $m_dir/input_args_$reg.txt $s3_run_m/

# c. Smoke ONE prelim tile, then fan out (plan M6-M7), exactly as step 6 --
#    the SAME --tile_list: one list per region serves both periods -- with
#    --args_url $s3_run_m/input_args_$reg.txt --tile_prefix $s3_out_m and a
#    _monthly_ tag.  Smoke on the quarterly memory high-water tile.
#    THE GATE THAT MATTERS: N_fit the same order as that tile's quarterly
#    N_fit.  An unreadable reference does not raise -- it leaves nothing
#    valid, and _expand_reference_files guards local paths, not URIs.
#    check_field_sizes.py reads -g=...,1/12 (plan M3): expect
#    dz/dz [25, 25, 94] for -W 60000.
#    A center the quarterly ATL14 does not cover has an all-NaN reference and
#    NO DATA for the fit ("smooth_fit: no valid data").  Normally it never
#    gets here: it wrote no quarterly tile either, so it came out of
#    the list at cleanup.  On a build before plan_tile_lists.sh TL1 such a fit
#    FAILED the job (IS E1020_N-2580, 9266c3d7, deterministic -- do not
#    retry); from TL1 on it is a successful job with no tile, which step 8's
#    fetch saves to prelim/no_data_tiles.txt for the cleanup.

# d. Matched, from the prelim tiles that EXIST (plan M8), as step 9.

# e. Mosaic, as step 10, with monthly.txt in place of quarterly.txt (plan M9).
#    No z0 task (skip_z0: the z0 spacing is over 1000 m).  IS: 44 tasks,
#    74 s at -P 12.  INSPECT THE QUEUE before running: make_mosaic_jobs.py
#    INFERS the dzdt lags from -t and -g rather than reading --dzdt_lags --
#    they must come out 1,3,6,12,24,36,48,60,72,84.
cd ~/ATL14_processing/runs
make_mosaic_jobs.py -b $m_dir -rr $reg -t 2018.75,2026.5 -e ATL14 \
    --run_name ${reg}_0332_monthly_mosaic @$HOME/git_repos/ATL1415/default_args/monthly.txt
cat ${reg}_0332_monthly_mosaic/queue/* | grep -oE "lag[0-9]+" | sort -t g -k2 -n -u

# f. netCDF: ATL15 ONLY (plan M10) -- ATL15_<reg>_0332_1mo_{2.5,10,20,40}km_006_02.nc,
#    91 epochs.
ATL15_write2nc.py @$m_dir/input_args_$reg.txt > ATL15.log 2>&1
#    COMPARING WITH QUARTERLY: raw delta_h IS comparable.  Both are 0 at the
#    2020 reference with sigma 0, and IS's raw median difference is -0.05 m.
#    Do NOT difference against the first epoch: it spreads that epoch's
#    errors across every other.  Expect a systematic SEASONAL difference of
#    a few tenths of a metre (IS: Apr -0.37 m, Oct +0.31 m) -- monthly
#    resolves a cycle quarterly smooths -- so many values will exceed 3x
#    sigma without anything being wrong.


# ===========================================================================
# 11. [ADE+DPS] [NEEDS CODE: scripts/maap/run_arctic_*.sh]  All five regions.
# ===========================================================================
# ONLY AFTER IS HAS GONE THROUGH END TO END.  These are the MAAP counterparts
# of scripts/run_arctic_{prelim,matched,mosaic,to_nc}.sh -- same argument
# order (release file, location file, period file, optional region), same
# per-region loop, same directory bookkeeping; the tail submits to DPS or runs
# the queue locally instead of calling sbatch.
bash scripts/maap/run_arctic_prelim.sh  default_args/latest_release.txt \
     default_args/MAAP_dps.txt default_args/quarterly.txt
for j in RA IS CN CS SV; do echo $j; check_MAAP_jobs.py ${j}_prelim_jobs.csv; done

bash scripts/maap/run_arctic_matched.sh default_args/latest_release.txt \
     default_args/MAAP_dps.txt default_args/quarterly.txt
for j in RA IS CN CS SV; do echo $j; check_MAAP_jobs.py ${j}_matched_jobs.csv; done

bash scripts/maap/run_arctic_mosaic.sh  default_args/latest_release.txt \
     default_args/MAAP_dps.txt default_args/quarterly.txt
bash scripts/maap/run_arctic_to_nc.sh   default_args/latest_release.txt \
     default_args/MAAP_dps.txt default_args/quarterly.txt

# monthly: the same four with default_args/monthly.txt.  The run_arctic_*.sh
# scripts already derive hemi_suffix="_monthly" by grepping the period file,
# so that part needs no change (Q17).


# ===========================================================================
# 12. [ADE] [SUGGESTION, NO SOFTWARE]  Annotate the run's build history.
# ===========================================================================
# The run's last step, quarterly or monthly: the git history across the
# builds its tiles ran, annotated with what changed and which tiles were
# rerun on which build.  The procedure is howto_MAAP_ogc O12b; it is not
# repeated here.
