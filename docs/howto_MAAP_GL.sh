# howto_MAAP_GL.sh -- Greenland on MAAP (per-tile solves on DPS)
#
# ############################################################################
# ##  TENTATIVE.  Written 2026-09-05 BEFORE any of it has been run end to   ##
# ##  end -- no ATL1415 tile has been solved on DPS yet.  This is the plan,  ##
# ##  not a record of a successful run.  Expect steps to move, split and    ##
# ##  change as testing advances; revise this file as that happens.         ##
# ############################################################################
#
# The discover/SLURM variant of this workflow is docs/howto_GL.sh, which is
# still the production path and is NOT replaced by this file.
#
# Every step is tagged:
#   [ADE] / [DPS]   where it runs
#   [OK]            expected to work as written; the pieces exist
#   [UNTESTED]      the pieces exist but this command has not been run on MAAP
#   [NEEDS CODE: x] blocked -- x does not exist yet
# and numbered 1..12 so steps can be referred to ("GL step 5") from
# docs/Transition_to_maap.md and from the other MAAP howtos.
#
# Prerequisite: docs/howto_MAAP_staging.sh S1-S6.  Nothing below works until
# the ADE has an ATL1415 env (S1) and the masks and ATL11 index are on the
# bucket (S2, S3).  Step 5 is gated on the smoke test, S7.
#
# THE TWO STRUCTURAL DIFFERENCES from howto_GL.sh:
#   - the location layer is default_args/MAAP_dps.txt, not discover.txt.  Use
#     it INSTEAD OF MAAP.txt: MAAP.txt names the ~/my-private-bucket mount,
#     which a DPS worker does not have.
#   - there is no hemisphere layer.  north.txt/south.txt exist only to supply
#     --ATL11_index, --ATL11_xover_dir and --Hemisphere; in cloud mode
#     MAAP_dps.txt carries the first, the third is a command-line option, and
#     the second has nothing to point at.  Pass --Hemisphere=1 directly.

conda activate ATL14
cd ~/git_repos/ATL1415
region_dir=/home/jovyan/ATL14_processing/rel006/north/GL
s3_run=s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/north/GL
s3_out=s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north/GL


# ===========================================================================
# 1. [ADE] [OK]  Point the release symlinks at this release.
# ===========================================================================
# Same as the discover workflow: each is a symlink to the release-specific file.
ln -sf rel_006_0331.txt default_args/latest_release.txt
ln -sf GL_0331.txt      default_args/GL_latest.txt


# ===========================================================================
# 2. [ADE] [UNTESTED]  Compose the args file.
# ===========================================================================
# Writes $region_dir/input_args_GL.txt.  Every path in it must resolve ON THE
# WORKER -- that is what MAAP_dps.txt is for.
setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
    default_args/GL_latest.txt default_args/quarterly.txt --Hemisphere=1

# CHECK BEFORE GOING ON: no path in the composed file should be local except
# --ATL14_root, which only setup itself reads.
grep -E '^(--mask_dir|--ATL11_index|--tide_directory|--ATL14_root)' $region_dir/input_args_GL.txt


# ===========================================================================
# 3. [ADE] [UNTESTED]  Publish the args file where DPS can localize it.
# ===========================================================================
# DPS localizes a `file` input by URL.  Whether it accepts an
# s3://maap-ops-workspace/... URL is open until the smoke test (staging S7.3).
aws s3 cp $region_dir/input_args_GL.txt $s3_run/


# ===========================================================================
# 4. [ADE] [NEEDS CODE: make_ATL1415_queue.py --xy_out]  Tile centers.
# ===========================================================================
# DPS wants x0/y0 pairs, NOT the shell command lines make_ATL1415_queue.py
# emits today.  Four separate problems block this step; all four are listed in
# Transition_to_maap.md ("Code that has to exist"):
#   - os.path.isfile() on the --ATL11_index URI, which is always False
#   - the --mask_dir join, which needs join_path_or_uri()/paths.exists()
#   - the greedy defaults regex that silently drops bare store_true flags
#     (--ATL11_earthaccess, --tide_adjustment), already fixed in
#     setup_ATL1415_region.py but still live here
#   - --xy_out itself
# ALSO GATED ON THE 1 km MASK (Q6/Q16): the v4.1 masks have no 1 km version,
# and the tile grid is derived from one.  Settled shape: an explicit
# --grid_mask_file, a tile list frozen per release in region_files/, and a 1 km
# decimation built on the fly with gdal_translate using MAX resampling, cached
# back to the bucket.
make_ATL1415_queue.py prelim $region_dir/input_args_GL.txt --xy_out GL_prelim_xy.txt


# ===========================================================================
# 5. [DPS] [NEEDS CODE: scripts/submit_MAAP_jobs.py]  Fan out one job per tile.
# ===========================================================================
# GATED ON THE SMOKE TEST (staging S7) -- do not fan out thousands of jobs
# before one has been shown to work.
#
# Policy the submitter needs (Q11): --max_in_flight (poll and top up), --rate,
# and record-and-continue on a submitJob error with the failure in the ledger.
# No account limits are known; for testing they are not yet a problem.  The
# ~10 jobs/hr public-queue throttle IS a problem for production -- see S6.
submit_MAAP_jobs.py --xy_file GL_prelim_xy.txt --step prelim \
    --args_url $s3_run/input_args_GL.txt \
    --out_prefix $s3_out/prelim \
    --queue maap-dps-worker-32gb \
    --tag GL_rel006_prelim --ledger GL_prelim_jobs.csv

# The ledger is a CSV the submitter writes: x0,y0,step,job_id,submitted_at.
# Job state is recovered from it rather than from listJobs() (Q10).


# ===========================================================================
# 6. [DPS] [NEEDS CODE: scripts/check_MAAP_jobs.py]  Watch, and requeue.
# ===========================================================================
# The slurm_run_status.py analogue: counts by state, failures grouped by the
# exception parsed from the job's stderr.
check_MAAP_jobs.py GL_prelim_jobs.csv
check_MAAP_jobs.py GL_prelim_jobs.csv --requeue both --dry_run
check_MAAP_jobs.py GL_prelim_jobs.csv --requeue both


# ===========================================================================
# 7. [ADE] [NEEDS CODE: deterministic output prefix]  Collect the tiles.
# ===========================================================================
# DPS scatters output to dps_output/<algo>/<tag>/YYYY/MM/DD/HH/MM/SS/<us>/,
# which the mosaic step cannot address.  The job itself must instead write to
#   $s3_out/prelim/E%d_N%d.h5
# i.e. the tree the ADE already builds, mirrored on the bucket, with output/
# kept for logs only.  The hemi suffix carries the quarterly/monthly split
# (Q9 + Q17): monthly goes to rel006/north_monthly/GL/.
aws s3 sync $s3_out/prelim/ $region_dir/prelim/


# ===========================================================================
# 8. [ADE] [OK]  Look at the tile sizes.
# ===========================================================================
# ATL11_to_ATL15.py now writes the field-size report itself, so
# make_field_size_report.py is no longer needed.  Inspect with the
# check_tiles.ipynb notebook in ATL1415.


# ===========================================================================
# 9. [DPS] [NEEDS CODE: run.sh prelim_prefix input]  The matched solve.
# ===========================================================================
# Steps 4-7 again with --step matched.  THE ONE EXTRA THING: a matched job
# reads the tile PLUS ITS 8 NEIGHBOURS (prior_edge_include), and DPS localizes
# one URL per input.  Per Q8(a), run.sh gets a `prelim_prefix` input and does
# `aws s3 cp` of the 9 keys itself -- one small change to run.sh, one config
# input, nothing staged.  That depends on step 7's deterministic prefix, so
# these two are really one problem and matched cannot be fanned out until a
# job can name its neighbours' prelim tiles by key.
make_ATL1415_queue.py matched $region_dir/input_args_GL.txt --xy_out GL_matched_xy.txt
submit_MAAP_jobs.py --xy_file GL_matched_xy.txt --step matched \
    --args_url $s3_run/input_args_GL.txt \
    --prelim_prefix $s3_out/prelim \
    --out_prefix $s3_out/matched \
    --queue maap-dps-worker-32gb \
    --tag GL_rel006_matched --ledger GL_matched_jobs.csv
check_MAAP_jobs.py GL_matched_jobs.csv
aws s3 sync $s3_out/matched/ $region_dir/matched/


# ===========================================================================
# 10. [ADE] [NEEDS CODE: scripts/run_queue_local.sh]  Mosaic.
# ===========================================================================
# make_mosaic_jobs.py still emits <run_dir>/queue/task_N plus a slurm_run.sh,
# and THERE IS NO sbatch IN THE ADE.  run_queue_local.sh reproduces
# packable_job.txt's bookkeeping -- mv queue->running->done, log to
# active_logs/ then logs/, copy to error_logs/ on nonzero exit or a
# ##TASK_LINE_FAILED## marker -- with `xargs -P` for concurrency, so
# slurm_run_status.py keeps working unchanged on those directories.
#
# ADE hardware (measured 2026-09-04): 16 cores, 124.4 GiB RAM, no walltime
# limit -- four times what a discover mosaic job asks for.  GL is well within
# it; whether AA is has not been timed (Q4/Q18).
make_mosaic_jobs.py -b $region_dir -rr GL -t 2018.75,2026.5 \
    --run_name GL_mosaic @default_args/quarterly.txt
run_queue_local.sh GL_mosaic -P 8


# ===========================================================================
# 11. [ADE] [NEEDS CODE: run_queue_local.sh]  netCDF and browse.
# ===========================================================================
ATL14_write2nc.py @$region_dir/input_args_GL.txt
ATL15_write2nc.py @$region_dir/input_args_GL.txt
# browse plots
ATL14_browse_plots.py @$region_dir/input_args_GL.txt
ATL15_browse_plots.py @$region_dir/input_args_GL.txt


# ===========================================================================
# 12. [ADE+DPS] [NEEDS CODE, as above]  The monthly variant.
# ===========================================================================
# Same 11 steps with default_args/monthly.txt in place of quarterly.txt, an
# --ATL14_reference_file, and north_monthly in every path (Q17: the monthly
# specifier reuses the hemi suffix).
#
# NOTE --ATL14_reference_file has the SAME glob-over-a-URI problem as
# --previous_product (Q27 W5): ATL11_to_ATL15.py line 682 calls
# glob.glob(ATL14_reference_file), which returns [] for an s3:// URI.  It is
# monthly-only, so it is not on the quarterly critical path -- but this step
# does not work until it is fixed.
region_dir_m=/home/jovyan/ATL14_processing/rel006/north_monthly/GL
setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
    default_args/GL_latest.txt default_args/monthly.txt --Hemisphere=1 \
    --ATL14_reference_file <the rel006 quarterly ATL14_GL netCDF>
