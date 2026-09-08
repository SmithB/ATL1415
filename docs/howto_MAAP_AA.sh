# howto_MAAP_AA.sh -- Antarctica on MAAP (per-tile solves on DPS)
#
# ############################################################################
# ##  TENTATIVE.  Written 2026-09-05 BEFORE any of it has been run end to   ##
# ##  end -- no ATL1415 tile has been solved on DPS yet.  This is the plan,  ##
# ##  not a record of a successful run.  Expect steps to move, split and    ##
# ##  change as testing advances; revise this file as that happens.         ##
# ############################################################################
#
# The discover/SLURM variant is docs/howto_AA.sh, which is still the
# production path and is NOT replaced by this file.
#
# READ docs/howto_MAAP_GL.sh FIRST.  GL is the reference workflow; this file
# marks only what Antarctica does differently, and refers to "GL step N" for
# the parts that are identical.  Tags and numbering follow the same scheme:
#   [ADE] / [DPS]   where it runs
#   [OK] [UNTESTED] [NEEDS CODE: x]
# Steps are numbered 1..14.
#
# Prerequisite: docs/howto_MAAP_staging.sh S1-S6; step 6 is gated on S7.
#
# READ THIS BEFORE DEBUGGING ANY 403 ON A MASK: pyTMD v3.0.9 set
# AWS_NO_SIGN_REQUEST=YES process-wide at import, which made GDAL read every
# /vsis3 object anonymously, so every mask read from s3://maap-ops-workspace
# came back HTTP 403.  Fixed 2026-09-06 in ATL1415/__init__.py.  Verified on
# the GL geotiff masks in both directions; full account in howto_MAAP_arctic.sh
# step 2 and in Transition_to_maap.md, "The pyTMD AWS_NO_SIGN_REQUEST bug".
#
# WHAT MAKES AA DIFFERENT, in one place:
#   a. it is submitted as TWO HALVES on a 400 km line, and stays that way on
#      DPS (Q7).  The halves differ in tile geometry, not just in extent.
#   b. the two halves should target DIFFERENT QUEUES (Q22): the near-pole
#      south tiles are the expensive ones.  This is exactly what the
#      `queue_name` input on the registered algorithm is for -- a per-job
#      override, no re-registration.
#   c. between the tile solves and the mosaic there are three extra ADE
#      stages: 200 km tiles, the four sectors, and per-sector mosaic jobs.
#   d. it is ~20 GiB of previous product across A1-A4, which is why the
#      previous-product read path had to stop downloading whole granules.
#   e. WHETHER THE ADE CAN MOSAIC AA AT ALL IS STILL AN OPEN MEASUREMENT
#      (Q4/Q18).  See step 12.

conda activate ATL14
cd ~/git_repos/ATL1415
region_dir=/home/jovyan/ATL14_processing/rel006/south/AA
region_dir_44=/home/jovyan/ATL14_processing/rel006/south/AA_44km
s3_run=s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/south/AA
s3_out=s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/south/AA
s3_out_44=s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/south/AA_44km

# TWO OUTPUT PREFIXES, ONE PER HALF, and they must stay separate.  The two
# halves solve DIFFERENT TILE SIZES -- 60 km / 40 km spacing for the north
# half, 44 km / 40 km for the south -- but a tile is named for its CENTER
# alone, 'E%d_N%d.h5' (make_ATL1415_queue.py:257), with nothing in the name to
# say which width produced it.  The discover workflow keeps them apart by
# giving each half its own region directory ($region_dir vs $region_dir_44);
# on the bucket that separation has to be made explicitly, or the halves
# overwrite each other.
#
# THEY DELIBERATELY OVERLAP -- confirmed by Ben 2026-09-08.  Step 4 filters
# with --min_xy 360000 (keep if max|xy| >= 360 km) and step 5 with
# --max_xy 440000 (keep if every |xy| <= 440 km), so any tile whose max|xy|
# falls in 360000..440000 is queued by BOTH halves: same center, same
# 'E%d_N%d.h5' filename, two different widths.  DO NOT "FIX" THE LIMITS.
#
# That makes the two prefixes above load-bearing rather than tidy: in the
# overlap band a tile center legitimately has TWO valid solutions, and any
# namespace that holds only one of them silently keeps whichever was written
# last.  Two consequences for code that is not written yet:
#
#   Q9, the deterministic output prefix.  The key CANNOT be the tile center
#   alone.  It has to carry the half (or the width) as well, or the 44 km and
#   60 km solutions of an overlap tile collide wherever they meet -- on the
#   bucket, in a job ledger, or in a requeue check that asks "does the output
#   for this center already exist?"  The answer to that question is only
#   well-posed per half.
#
#   Q8, the matched neighbourhood.  A matched job localizes its tile's prelim
#   output plus its 8 neighbours BY NAME.  For a tile in the overlap band,
#   "the prelim tile at this center" is ambiguous, and a neighbourhood
#   assembled across the two halves would mix 44 km and 60 km fits.  Each
#   half's matched pass must draw only on its own prelim tree -- which is what
#   the discover workflow gets for free from $region_dir vs $region_dir_44,
#   and what step 9 below has to reproduce explicitly.


# ===========================================================================
# 0. [ADE] [UNTESTED]  Rebuild the DPS image if any code has changed.
# ===========================================================================
# DPS DOES NOT RUN THIS WORKING COPY.  It clones repository_url at
# algorithm_version (on_s3) FROM GITHUB at build time and bakes the result into
# a container, so anything uncommitted, unpushed, or committed since the last
# build is simply not on the worker -- and nothing in a job log says so.  A
# stale image fails as a wrong-looking runtime error, not as a version error.
#
#   git -C ~/git_repos/ATL1415 status --short                     # nothing uncommitted
#   git -C ~/git_repos/ATL1415 log --oneline origin/on_s3..on_s3  # empty
#
# If either is non-empty, push, then re-register and wait for the build to go
# green before submitting anything -- staging S5, which carries the rule and
# the record of the one rebuild this has already forced.


# ===========================================================================
# 1. [ADE] [OK]  Point the release symlinks at this release.
# ===========================================================================
ln -sf rel_006_0331.txt default_args/latest_release.txt
ln -sf AA_0331.txt      default_args/AA_latest.txt


# ===========================================================================
# 2. [ADE] [UNTESTED]  Compose the args file.        (as GL step 2)
# ===========================================================================
setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
    default_args/AA_latest.txt default_args/quarterly.txt --Hemisphere=-1

# AA is the region that actually exercises --tide_adjustment, which was being
# silently dropped by the greedy defaults regex until that was fixed in
# setup_ATL1415_region.py (Q15).  Confirm it survived into the composed file:
grep -E '^(--tide_adjustment|--tide_model|--mask_dir)' $region_dir/input_args_AA.txt
#
# THE PREVIOUS PRODUCT IS NOW A CMR SEARCH, not a discover path (Q27 W3/W4,
# b293807).  MAAP_dps.txt carries --previous_product_earthaccess, so setup
# rewrites --previous_product_top into --previous_product=<release>_<cycles>
# and drops the /discover/... tree.  The composed file should therefore contain
#   --previous_product_earthaccess
#   --previous_product=005_0329
# and NO --previous_product_top.  EXPECTATION, from MAAP_dps.txt and the IS run
# of 2026-09-06 -- not yet observed for this region.
grep -E '^--previous_product' $region_dir/input_args_AA.txt
#
# AA IS THE REGION WHERE W3's SIMPLIFICATION SHOWS.  The discover workflow
# globs A1-A4 and reads up to ~20 GiB of previous product; the bounding-box CMR
# search returns only the sectors a tile actually touches, which on the tested
# Antarctic tile was exactly one ATL14_A2 granule.


# ===========================================================================
# 3. [ADE] [UNTESTED]  Publish the args file.        (as GL step 3)
# ===========================================================================
aws s3 cp $region_dir/input_args_AA.txt $s3_run/
# ...and the south half's, composed with the overrides file (see step 5):
setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
    default_args/AA_latest.txt default_args/quarterly.txt default_args/AA_44km.txt \
    --Hemisphere=-1
aws s3 cp $region_dir_44/input_args_AA_44km.txt $s3_run/
# RUN 2026-09-08: 1375 and 1385 bytes, both under $s3_run.


# ===========================================================================
# 3b. [ADE+DPS] [READY, BLOCKED ON STEP 3]  The cost-characterisation transect.
# ===========================================================================
# WHY, and it is not part of the production workflow: the Iceland smoke tile
# (staging S7) took 26.6 minutes for 239613 ATL11 points, and it is close to
# the FLOOR of the cost range -- a small, low-latitude, grounded, tide-free
# tile.  Nothing yet says what the expensive end looks like, and the S6 queue
# request cannot be written without it.  This runs a spread of Antarctic tiles
# and records time, peak memory and input size for each.
#
# 16 tiles, in scripts/maap/AA_queue_xy.txt, chosen against the real masks and
# documented one row each in scripts/maap/AA_queue_manifest.csv:
#
#   10 along azimuth 90 deg through East Antarctica, from the pole outward.
#      One sits INSIDE the ICESat-2 pole hole (ice_frac 0.00, r=100 km): the
#      mask encodes the hole exactly, no ATL11 data exists there, and the tile
#      is a free production test of the empty-tile fix -- it should skip
#      cleanly rather than raise.  One sits on the hole's EDGE at 88 S, where
#      track convergence peaks.  The rest run out to the coastal margin at
#      2420 km, where track density is lowest.
#
#    6 that exercise the TIDE CORRECTION, which no DPS job has touched: three
#      fully floating (Ross, Ronne, Amery: tide_frac 1.00) and three straddling
#      a grounding line (tide_frac 0.32/0.44/0.67).  The grounding-line tiles
#      are the interesting ones -- a tile that is partly floating is where the
#      tide mask boundary actually has to be right.  These also exercise the
#      ANONYMOUS s3://pytmd read, the one credential path of the three that no
#      job has used yet: AA sets --tide_model=CATS2008-v2023 and
#      --tide_adjustment, and IS does not.
#
# Run it (needs the args file from step 3):
scripts/maap/submit_AA_queue.py scripts/maap/AA_queue_xy.txt \
    $s3_run/input_args_AA.txt maap-dps-worker-32gb AA_queue_jobs.csv
scripts/maap/collect_AA_queue.py AA_queue_jobs.csv
#
# The collector joins each job's status to the peak RSS and elapsed time the
# job reports about ITSELF (scripts/run_with_rusage.py, one line per fit /
# error / matched step), plus N_ATL11 and N_fit parsed from its log.  It does
# not rely on getJobMetrics, which returned an empty dict for the one job that
# has succeeded so far.
#
# EXPECT SOME TO FAIL, and that is a result too: a tile that OOMs on
# maap-dps-worker-32gb has told us that its class needs a bigger queue, which
# is exactly what S6 has to ask for.  Re-run those on maap-dps-worker-64gb.
#
# Regenerate the queue with scripts/maap/make_AA_queue.py (which recomputes
# ice_frac and tide_frac from the staged masks) if the tile geometry changes.


# ===========================================================================
# 4. [ADE] [NEEDS CODE: make_ATL1415_queue.py --xy_out]  North-half centers.
# ===========================================================================
# Same four blockers as GL step 4, plus the 1 km grid mask -- Q6/Q16 are
# ANSWERED, so what is missing there is the code, not a decision.  For AA
# the mask problem is worse: AntarcticIceMask_..._240m_v4.1.tif has neither
# '100m' nor '125m' in its name, so make_ATL1415_queue.py raises ValueError
# outright rather than merely failing to find a sibling.
make_ATL1415_queue.py prelim $region_dir/input_args_AA.txt --min_xy 360000 \
    --xy_out AA_north_prelim_xy.txt


# ===========================================================================
# 5. [ADE] [NEEDS CODE, as step 4]  South-half centers, different geometry.
# ===========================================================================
# The south half is a SEPARATE REGION DIRECTORY with a different tile size:
# W=44000, spacing 40000, against the 60 km / 40 km of the north half.  That
# is why it cannot simply be a --max_xy filter on the same queue.
#
# THE 44 km ARGS FILE, which nothing used to compose.  It is used here, at step
# 6 and at step 9, but step 2 composed only input_args_AA.txt -- and the
# discover howto (docs/howto_AA.sh:23) has the identical gap.  Resolved
# 2026-09-08 with default_args/AA_44km.txt, an overrides file carrying just
#   --region=AA_44km
#   -W=44000
# layered AFTER AA_latest.txt and the release file, which set --region=AA and
# -W=60000.  It has to be a FILE: setup_ATL1415_region.py takes only
# defaults_files, --ATL14_reference_file and --Hemisphere on the command line,
# so neither --region nor -W can be overridden there.  --region drives both the
# directory and the args-file name, so this lands exactly where steps 5/6/9
# expect it.  Both files were composed and published 2026-09-08; they differ in
# three lines -- --region, -W and -b -- and in nothing else.
make_ATL1415_queue.py prelim $region_dir_44/input_args_AA_44km.txt --max_xy 440000 \
    --xy_out AA_south_prelim_xy.txt


# ===========================================================================
# 6. [DPS] [NEEDS CODE: scripts/submit_MAAP_jobs.py]  Fan out, two submissions.
# ===========================================================================
# GATED ON THE SMOKE TEST (staging S7).  Two ledgers, and per Q22 two queues:
# the south half gets the larger instance.  Sizing is a guess until a real
# tile is timed -- that is smoke-test question 5.
submit_MAAP_jobs.py --xy_file AA_north_prelim_xy.txt --step prelim \
    --args_url $s3_run/input_args_AA.txt --out_prefix $s3_out/prelim \
    --queue maap-dps-worker-32gb \
    --tag AA_rel006_prelim_north --ledger AA_north_prelim_jobs.csv

submit_MAAP_jobs.py --xy_file AA_south_prelim_xy.txt --step prelim \
    --args_url $s3_run/input_args_AA_44km.txt --out_prefix $s3_out_44/prelim \
    --queue maap-dps-worker-32vcpu-64gb \
    --tag AA_rel006_prelim_south --ledger AA_south_prelim_jobs.csv


# ===========================================================================
# 7. [DPS] [NEEDS CODE: scripts/check_MAAP_jobs.py]  Watch both halves.
# ===========================================================================
for L in AA_north_prelim_jobs.csv AA_south_prelim_jobs.csv; do
    echo $L; check_MAAP_jobs.py $L
done


# ===========================================================================
# 8. [ADE] [NEEDS CODE: deterministic output prefix]  Collect.  (as GL step 7)
# ===========================================================================
aws s3 sync $s3_out/prelim/    $region_dir/prelim/
aws s3 sync $s3_out_44/prelim/ $region_dir_44/prelim/


# ===========================================================================
# 9. [DPS] [NEEDS CODE: run.sh prelim_prefix input]  Matched, both halves.
# ===========================================================================
# As GL step 9, twice.  The discover workflow uses --lines_per_task 4 for AA
# matched; on DPS that has no analogue -- one job is one tile.
make_ATL1415_queue.py matched $region_dir/input_args_AA.txt --min_xy 360000 \
    --xy_out AA_north_matched_xy.txt
make_ATL1415_queue.py matched $region_dir_44/input_args_AA_44km.txt --max_xy 440000 \
    --xy_out AA_south_matched_xy.txt
# ... then two submit_MAAP_jobs.py calls with --step matched, EACH POINTED AT
# ITS OWN HALF'S PRELIM TREE -- north at --prelim_prefix $s3_out/prelim, south
# at --prelim_prefix $s3_out_44/prelim.  NOT one shared prefix, and this is the
# step where getting it wrong is least visible: a matched job localizes its
# tile's own prelim output plus its 8 neighbours by name, so in the deliberate
# 360-440 km overlap band a shared prefix would hand it a neighbourhood mixing
# 44 km and 60 km fits at identical filenames.  The result would be a solved
# tile, not an error.  See the note at the top of this file.
aws s3 sync $s3_out/matched/    $region_dir/matched/
aws s3 sync $s3_out_44/matched/ $region_dir_44/matched/


# ===========================================================================
# 10. [ADE] [NEEDS CODE: run_queue_local.sh]  200 km tiles, both halves.
# ===========================================================================
# make_200km_tiles.py emits a queue directory plus a slurm_run.sh, so it needs
# the local runner for the same reason the mosaic does.
make_200km_tiles.py $region_dir AA -t 2018.75,2026.5
run_queue_local.sh tile_run_AA -P 8

make_200km_tiles.py $region_dir_44 AA --name AA_south --W 44000 --spacing 40000 \
    -t 2018.75,2026.5
run_queue_local.sh tile_run_AA_south -P 8


# ===========================================================================
# 11. [ADE] [UNTESTED]  Set up the four Antarctic sectors.
# ===========================================================================
# Pure bookkeeping over the 200 km tiles; no solve, no cloud read.  Expected to
# work unchanged, but it has never been run in the ADE.
setup_AA_sectors.py /home/jovyan/ATL14_processing/rel006/south


# ===========================================================================
# 12. [ADE] [NEEDS CODE: run_queue_local.sh]  Mosaic, per sector.
# ===========================================================================
# THIS IS THE STEP Q4/Q18 IS ABOUT.  On discover an AA mosaic is a 4-hour,
# 4-task SLURM job per field group.  The ADE has 16 cores, 124.4 GiB RAM and
# no walltime limit -- four times the tasks and no clock -- so the shape of
# the answer is encouraging, but NOTHING HAS BEEN TIMED, because
# pointCollection was not importable in the ADE until staging S1.
#
# TO SETTLE IT: after S1, time one z0 field group for a small region (IS or
# GL) and extrapolate by tile count to AA.  If the answer is no, the mosaic
# stops being an ADE step and this file changes shape.
for sector in A1 A2 A3 A4; do
    make_200km_to_mosaic_jobs.py -b /home/jovyan/ATL14_processing/rel006/south/$sector \
        -rr $sector -t 2018.75,2026.5
    run_queue_local.sh ${sector}_mosaic -P 8
done


# ===========================================================================
# 13. [ADE] [NEEDS CODE: a MAAP variant of run_antarctic_tonc.sh]  netCDF.
# ===========================================================================
# scripts/run_antarctic_tonc.sh ends in `setup_slurm_run.py ...; sbatch`.  Its
# release/root greps already work against MAAP_dps.txt, so only that tail
# changes.  PUT THE VARIANT IN scripts/maap/ so the SLURM original stays
# untouched.
bash scripts/maap/run_antarctic_tonc.sh default_args/latest_release.txt \
    default_args/MAAP_dps.txt


# ===========================================================================
# 14. [ADE+DPS] [NEEDS CODE, as above]  The monthly variant.
# ===========================================================================
# Same 13 steps with monthly.txt, an --ATL14_reference_file, south_monthly in
# every path (Q17), and setup_AA_sectors.py --near_pole_radius 0.
#
# --ATL14_reference_file: W5 (b293807) FIXED THE SILENT FAILURE HERE BUT NOT
# THE USE CASE.  The discover workflow passes a GLOB PATTERN across sectors
# ("rel005_0329/south/A*/ATL14_*_0329_100m_005_02.nc"), and a URI containing a
# wildcard now RAISES ValueError instead of quietly yielding an empty reference
# DEM -- correct, and better than before, but AA monthly still has no cloud
# path.  GL's single-file case works; AA's does not.  Not on the quarterly
# critical path, so this is where it stops for now.
#
# OPEN, and it is a small one: the ValueError says "Name the granules
# explicitly, one --ATL14_reference_file each", but the argument is a plain
# type=path_or_uri with no action='append' (ATL11_to_ATL15.py:987), so a second
# occurrence overwrites the first.  Either the message or the argument is
# wrong.  Decide which when AA monthly is picked up.
