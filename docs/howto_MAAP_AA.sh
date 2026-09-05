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


# ===========================================================================
# 3. [ADE] [UNTESTED]  Publish the args file.        (as GL step 3)
# ===========================================================================
aws s3 cp $region_dir/input_args_AA.txt $s3_run/


# ===========================================================================
# 4. [ADE] [NEEDS CODE: make_ATL1415_queue.py --xy_out]  North-half centers.
# ===========================================================================
# Same four blockers as GL step 4, plus the 1 km mask (Q6/Q16) -- and for AA
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
    --args_url $s3_run/input_args_AA_44km.txt --out_prefix $s3_out/prelim \
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
aws s3 sync $s3_out/prelim/ $region_dir/prelim/


# ===========================================================================
# 9. [DPS] [NEEDS CODE: run.sh prelim_prefix input]  Matched, both halves.
# ===========================================================================
# As GL step 9, twice.  The discover workflow uses --lines_per_task 4 for AA
# matched; on DPS that has no analogue -- one job is one tile.
make_ATL1415_queue.py matched $region_dir/input_args_AA.txt --min_xy 360000 \
    --xy_out AA_north_matched_xy.txt
make_ATL1415_queue.py matched $region_dir_44/input_args_AA_44km.txt --max_xy 440000 \
    --xy_out AA_south_matched_xy.txt
# ... then two submit_MAAP_jobs.py calls with --step matched and
# --prelim_prefix $s3_out/prelim, exactly as in GL step 9.
aws s3 sync $s3_out/matched/ $region_dir/matched/


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
# --ATL14_reference_file has the glob-over-a-URI problem of Q27 W5; on AA the
# discover workflow passes it a GLOB PATTERN across sectors
# ("rel005_0329/south/A*/ATL14_*_0329_100m_005_02.nc"), which makes it a
# harder case than GL's single file.  Not on the quarterly critical path.
