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
# ...and the south half's, DERIVED from it rather than composed (see step 5):
scripts/maap/make_AA_44km_args.py $region_dir/input_args_AA.txt \
    $region_dir_44/input_args_AA_44km.txt
aws s3 cp $region_dir_44/input_args_AA_44km.txt $s3_run/
# RUN 2026-09-08.  The two files differ in exactly two lines, -W and -b; the
# script refuses to run if the source does not carry --region=AA.


# ===========================================================================
# 3b. [ADE+DPS] [SUBMITTED 2026-09-11 on ab84687, running]  The cost-characterisation transect.
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
# CROSSOVERS CHANGED UNDER THIS STEP, 2026-09-09 (e798344).  Every DPS job so
# far has been along-track only -- E220_N20 reported N_AT=935506, N_XO=0 -- and
# a cloud run now reads crossovers too, keyed by --ATL11xo_version out of the
# release args file.  TWO CONSEQUENCES HERE:
#   - REBUILD FIRST (step 0).  DPS bakes the repo in at build time, so a
#     transect run on the current image would measure the OLD behaviour and
#     silently look like a valid cost characterisation.
#   - THE COST NUMBERS MOVE.  Crossovers add points to the fit, so time, peak
#     RSS and N all rise relative to anything measured before the rebuild.  Do
#     not mix pre- and post-rebuild rows in the same table.
# It is also the cheapest confirmation that the fix works on a worker rather
# than only against CMR: N_XO must now be non-zero.  collect_AA_queue.py
# reports it as a column since 2026-09-10 (from the solve's
# "Decimate_data: N_AT=..., N_XO=..." line), so no hand-run grep is needed.
# THE BASELINE IS ZERO EVERYWHERE: the ported collector, run over all twelve
# pre-fix transect jobs, reads N_XO=0 on every one of them -- and reproduces
# AA_cost_results.csv exactly, field for field, while doing it.
# Under the OGC runner the solve's lines are in the job's _stderr.txt, not
# _stdout.txt (howto_MAAP_ogc QD); the collector reads both.
#
# ---------------------------------------------------------------------------
# 3b-i. [OK, 2026-09-11]  VERIFY THE CROSSOVER READ FIRST, on two tiles, not sixteen.
# ---------------------------------------------------------------------------
# >>> BLOCKED 2026-09-10 on docs/howto_MAAP_ogc.sh: submit_AA_queue.py uses
# submitJob, gone in maap-py 5.x.  This step is OGC step O8. <<<
#
# DECIDED 2026-09-10 (Ben): before re-measuring cost, just establish that
# ATL11XO is read on a worker at all.  Two tiles near the pole hole answer it,
# and the full transect can wait until the answer is yes.
#
# Step 0 first, and then S5b -- there is no point submitting anything until the
# image is known to carry the crossover commit:
#     /srv/conda/envs/notebook/bin/python register_algorithm.py
#     /srv/conda/envs/notebook/bin/python scripts/maap/check_build_id.py
#
scripts/maap/submit_AA_queue.py scripts/maap/AA_xo_check_xy.txt \
    $s3_run/input_args_AA.txt $s3_run/input_args_AA_44km.txt \
    maap-dps-worker-32gb AA_xo_check_jobs.csv
#
# FIXED 2026-09-10: this command and step 3b's passed FOUR arguments to a
# five-positional script (xy, args_60km, args_44km, queue, ledger), which put
# the queue name in the 44 km args slot and the ledger name in the queue
# slot -- confirmed with --dry-run: queue=AA_xo_check_jobs.csv,
# args=maap-dps-worker-32gb.  The script now refuses that arrangement.
#
# WHY THESE TWO (scripts/maap/AA_xo_check_xy.txt):
#   220000 20000  the pole-hole EDGE tile, 88 S, where track convergence peaks
#                 -- so crossover density is at its highest anywhere on the
#                 continent, which makes N_XO=0 unambiguous if the fix failed.
#                 It is ALSO the exact tile that produced the pre-fix evidence:
#                 E220_N20 reported N_AT=935506, N_XO=0.  Rerunning it is a
#                 direct before/after on one tile rather than an argument.
#   300000 20000  just outside the hole, fully grounded (ice_frac 1.0), as the
#                 control: whatever E220_N20 does, an ordinary interior tile
#                 should do too.
# Both have max|xy| below 360 km, so halves_for() routes each to the 44 km half
# ONLY -- two centers, two jobs, not four.
#
# WHAT SAYS IT WORKED:
scripts/maap/collect_AA_queue.py AA_xo_check_jobs.csv
# N_XO > 0 on both.  For E220_N20 compare against N_AT=935506, N_XO=0; for
# E300_N20, N_AT=1256488, N_XO=0 -- both read back 2026-09-10 from the
# pre-fix jobs by the collector itself.
# RUN 2026-09-10/11 on build 8935494 -- PASSED:
#   E220_N20  N_XO=171488  N_AT=935506   3.1 h (was 2.2)  21.40 GiB (was 20.93)
#   E300_N20  N_XO=52882   N_AT=1256488  1.3 h (was 1.1)  16.32 GiB (was 15.70)
# Along-track counts unchanged, and N_ATL11 rose by exactly N_XO on both.
# So 3b's cost numbers (AA_cost_results.csv) are pre-crossover and low --
# most of all near the pole; rerun the transect before sizing S6.
#
# EXPECT IT TO COST MORE THAN THE PRE-FIX RUN.  44km_E220_N20 was the most
# expensive tile in the whole transect at 132.5 min and 20.9 GiB on a 32 GiB
# worker -- about 35% headroom -- and crossovers only add points.  If it OOMs,
# that is a RESULT and not a failure: rerun on maap-dps-worker-64gb, and note
# that the S6 queue request has to assume the post-crossover numbers.  N_XO is
# printed by decimate_data early in the fit, so even a job that later dies has
# already answered the question.
#
# Run it (needs the args file from step 3):
scripts/maap/submit_AA_queue.py scripts/maap/AA_queue_xy.txt \
    $s3_run/input_args_AA.txt $s3_run/input_args_AA_44km.txt \
    maap-dps-worker-32gb AA_queue_jobs.csv
scripts/maap/collect_AA_queue.py AA_queue_jobs.csv
#
# SUBMITTED 2026-09-11 ~15:50 UTC, the POST-CROSSOVER rerun, on build ab84687
# (MATCH, maap_pgt=set: howto_MAAP_ogc O6 run 3); queue -32gb.  All 17
# submit_job calls accepted -- 16 centers, E420_N20 in both halves.  Ledger,
# outside the checkout so it cannot block register_algorithm.py:
#   ~/ATL14_processing/maap_ledgers/AA_transect_ab84687_jobs.csv
# Read it with
scripts/maap/collect_AA_queue.py ~/ATL14_processing/maap_ledgers/AA_transect_ab84687_jobs.csv
# First run on a build whose tiles print their BUILD_ID and write /meta
# build_* (howto_MAAP_ogc O11 option 1, O12a), so its commit column should
# read ab84687 on every row -- anything else is a stale worker.  Its rows
# REPLACE AA_cost_results.csv; do not mix the two (see above).
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
# 2026-09-08 with scripts/maap/make_AA_44km_args.py, which DERIVES it from the
# north-half file, changing exactly two lines: -W to 44000 and -b to the
# AA_44km region directory.
#
# IT IS NOT A default_args OVERRIDES FILE, which was tried first and broke all
# four 44 km jobs.  setup_ATL1415_region.py derives both the region directory
# and the args-file name from --region, so producing input_args_AA_44km.txt
# through it means --region=AA_44km -- and --region is NOT A LABEL.
# ATL11_to_ATL15.py:595 loads the gridded mask only for region in ['AA','GL'];
# any other value falls through both branches with mask_data left None, and the
# solve dies at line 619 with "'NoneType' object has no attribute 'z'".  The
# 60 km jobs succeeded at the same tile centers, including E420_N20 in the
# overlap band, which is what isolated it to the args file.  The south half
# therefore KEEPS --region=AA and differs only in geometry and output location.
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
