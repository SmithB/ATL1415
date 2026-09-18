#! /usr/bin/env bash
# ===========================================================================
# PLAN: run the MONTHLY product (dt = 1/12 yr) on MAAP.  IS first.
# Written 2026-09-18.  TENTATIVE -- NOTHING HERE HAS RUN YET.  Revise as steps
# land; every step carries its own status tag.
# ===========================================================================
# WHY NOW (Ben, 2026-09-18): "Assuming that the differences from the previous
# product do not indicate > 10m errors or major gaps, move on to the monthly
# steps.  Develop a plan for running monthly on maap."
# READ AS: the monthly PRODUCT (docs/howto_arctic.sh lines 36-41, the discover
# "monthly:" block), not a monthly schedule for the quarterly run.
#
# Provenance per claim: STATEMENT = verified 2026-09-18, with how;
# DECIDED = Ben said so; RECOMMENDATION = mine, overridable;
# QUESTION = open, for Ben, not guessed.
# Tags: [ADE] / [DPS] for where a step runs; [NOT STARTED] [NEEDS CODE: x]
# [BLOCKED: x] [DONE].
#
#
# ===========================================================================
# OPEN QUESTIONS FOR BEN -- first, because M1-M10 wait on them
# ===========================================================================
# QM1. Does the 0332 IS quarterly pass your rel005 bar, so that
#      ATL14_IS_0332_100m_006_02.nc can be the monthly reference DEM?
#      (plan_cycles_03_32.sh T8 has the full numbers.)
#        - gaps: none major (0.11% of rel005 ATL14 cells, 0.09% of ATL15);
#        - ATL15 1 km delta_h: median -0.03 m, |d|>10 m on 0.05% of cells;
#        - ATL14 h: |d|>10 m on 3.8% of cells (max 330 m) -- but 99% of those
#          have data_count 0; where there ARE data it is 0.36% (427 cells).
#      RECOMMENDATION: pass -- the large differences are interpolation between
#      tracks, where both products' h_sigma are ~11 m.  The reference is read
#      only AT DATA POINTS (z_ref = ref_dem.interp(data.x, data.y)), which is
#      where the two agree best.
#      A. Pass; M1 onward.            B. Hold; look at the 427 cells first.
# AM1: Pass
#
# QM2. time_coverage_duration is wrong in every product the ADE writes
#      (ATL1415_attrs_meta.py:314; plan_cycles_03_32.sh T8).  Fix it before
#      M10, and re-write the five 0332 quarterly files, or leave it for later?
#      RECOMMENDATION: fix now -- one line plus a test, ADE-only (the writers
#      run in the ADE; the DPS image does not use attrs_meta), so NO rebuild.
#      Re-writing the quarterly files takes ~35 s.  Do it BEFORE M1, so the
#      reference DEM published to the bucket is the final object.
#      A. Fix now, re-write, then M1.   B. Later; monthly carries the bug too.
# AM2: A. Fix now, re-write, then M1
#
# QM3. Scope: IS alone first, as the quarterly run did?
#      RECOMMENDATION: yes.  GL and AA have never run quarterly on MAAP, and a
#      region's monthly run needs its quarterly ATL14 first (M12).
#      A. IS only.   B. Name the regions.
# AM3: A. IS only. 
#
# QM4. [CLOSED 2026-09-18 BY EVIDENCE, not by Ben]  Integer seconds or an ISO
#      8601 duration string for time_coverage_duration?
#      I recommended ISO 8601 to Ben on the grounds that CF/NSIDC expect it.
#      THAT RECOMMENDATION WAS WRONG and is withdrawn.  STATEMENT, from the
#      code and the released product:
#        - rel005 (ATL14_IS_0329_100m_005_02.nc), which NSIDC accepted,
#          carries 2.17e8 -- a NUMBER of seconds (T8, plan_cycles_03_32.sh:433);
#        - main computes end-start over delta_time, also seconds
#          (ATL14_attrs_meta.py:217);
#        - root_info's own default is 0., a float (ATL1415_attrs_meta.py:28);
#        - the two metadata templates carry the placeholder "SET_BY_PGE",
#          which fixes no type.
#      DECIDED: integer seconds, matching the released series.  Ben said "go
#      ahead with the plan" without picking; the evidence picked for him, and
#      the change is one line and a ~30 s re-write if he disagrees.
#
#
# ===========================================================================
# BACKGROUND.  All STATEMENT, 2026-09-18, by reading the code named.
# ===========================================================================
# WHAT MONTHLY IS ON DISCOVER (scripts/run_arctic_{prelim,matched,mosaic,to_nc}.sh):
#   the SAME pipeline -- prelim, matched, mosaic, netCDF -- with two changes:
#   1. default_args/monthly.txt in place of quarterly.txt:
#        --hemi_suffix=_monthly        region dir rel006/north_monthly/IS
#        -g=1250,2500,1/12             z0 1250 m, dz 2500 m, dt one month
#        --dzdt_lags=1,3,6,12,24,36,48,60,72,84
#   2. --ATL14_reference_file = the QUARTERLY ATL14 of the same release,
#      cycles and version:  rel006/north/IS/ATL14_IS_0332_100m_006_02.nc.
#      The solver subtracts it from every point (z -= z_ref), edits
#      |z - z_ref| >= DEM_tol (50), and skips the DEM three-sigma edit on the
#      error pass (ATL11_to_ATL15.py:753, 793-803).
#   And only ATL15 is written (run_arctic_to_nc.sh: no ATL14_write2nc.py
#   for monthly).
#
# WHAT ALREADY WORKS FOR MAAP, by reading -- none of it RUN for monthly:
#   - --ATL14_reference_file takes a URI: _expand_reference_files()
#     (ATL11_to_ATL15.py:265-295) passes a single s3:// name through
#     unglobbed; pc.grid.mosaic().from_list -> from_nc reads a URI since
#     pointCollection PR #53 (Transition_to_maap.md Q27 W2).  A URI WITH A
#     WILDCARD RAISES -- matters for AA (M12), not IS.
#   - setup_ATL1415_region.py copies --ATL14_reference_file through verbatim
#     (line 60), builds rel<R>/north_monthly/<region> from --hemi_suffix
#     (lines 92-93), and does not write --hemi_suffix out (line 195).
#   - run.sh takes the tile prefix per job and the tile spacing from
#     --tile_spacing (40000, from IS.txt, unchanged); a grep finds nothing
#     quarterly-specific in it.  s3_tiles.py's docstring already names
#     .../rel006/north_monthly/IS as a tile prefix.
#   - make_mosaic_jobs.py parses '1/12' (line 283) and sets skip_z0 because
#     the z0 spacing exceeds 1000 m (line 291).
#   - ATL15_write2nc.py parses '1/12' (line 319), picks '1mo' from delta_t
#     (line 182), and names the native grid from the dz spacing (lines
#     170-181): ATL15_IS_0332_1mo_{2.5,10,20,40}km_006_02.nc.
#   THEREFORE NO DPS REBUILD: build 61a19af carries everything the solve
#   needs.  RECOMMENDATION: confirm with check_build_id before M6 anyway.
#
# WHAT DOES NOT WORK YET:
#   - check_field_sizes.py reads -g with float() (lines 77-80), so '1/12'
#     raises CannotCheck.  [NEEDS CODE] -- M3.
#
# SIZES, computed from the args (not measured):
#   tile dz: 60000/2500 + 1 = 25 per side; (2026.5 - 2018.75)*12 + 1 = 94
#     epochs -> dz/dz [25, 25, 94].  z0 at 1250 m: 49 x 49.
#   netCDF: --t_crop=2019,2026.5 keeps 91 of the 94 epochs.
#   Unknowns per tile ~ 2401 + 58750 against quarterly's 361201 + 119072, so
#   about 8x fewer -- but the DATA are the same points (N_fit up to 273382
#   on IS).  RECOMMENDATION: expect less memory than quarterly's 9.52 GiB
#   peak; MEASURE it (M6) before trusting that.
#
#
# ===========================================================================
# M0. [ADE] [DONE 2026-09-18 -- AM1 "Pass"]  Accept the reference DEM.
# ===========================================================================
# ATL14_IS_0332_100m_006_02.nc is the monthly reference DEM.  AM2 was "fix
# now", so M0b ran first and the file published in M1 is the re-written one.
#
# M0b. [ADE] [DONE 2026-09-18]  Fix time_coverage_duration.
#      FIXED in ATL1415_attrs_meta.py:313-317: int((datetime_end -
#      datetime_start).total_seconds()), replacing
#      int((datetime_start-datetime_end).seconds).  Numeric, per QM4.
#      The duration now agrees with the start and end written beside it --
#      start + duration == end -- and keeps the region offset that
#      set_time_range adds to the start, so the three attributes are
#      self-consistent.
#      TEST: tests/test_time_coverage.py, 8 cases -- the invariant over
#      IS/GL/AA, the known 0332 span (236681615 s), the shape of the old bug,
#      a one-year span, the int type, and the METADATA/Extent mirror.
#      Suite 108 passed 2 skipped (was 100 passed before these 8).
#      RE-WRITTEN, into ~/ATL14_processing/runs/IS_0332_nc (logs overwritten):
#      ATL14_write2nc.py 12 s, ATL15_write2nc.py 18 s, both exit 0, NO INVALID
#      warning -- lineage still complete.
#      VERIFIED, ncdump -h diffed against a snapshot of all five headers taken
#      before the re-write: the ONLY attribute changed is
#      time_coverage_duration, 54385 -> 236681615, in all five files.  The
#      rest of each diff is the expected regeneration stamps (date_created,
#      history, identifier_file_uuid, the per-file uuid).  All five file sizes
#      are byte-identical to before (9915264, 18117749, 895146, 701363,
#      638098).
#      DATA UNCHANGED, re-running T8's own checks against the mosaic:
#      ATL14 h (3001,4201), finite(h) == finite(z0) & ice_area>0 exactly,
#      max |h - z0| = 1.219e-4 m (T8 recorded 1.2e-4); ATL15 1 km delta_h
#      (31,301,421), 31 epochs.
#
#
# ===========================================================================
# M1. [ADE] [DONE 2026-09-18]  Publish the reference DEM.
# ===========================================================================
region_dir=/home/jovyan/ATL14_processing/rel006/north/IS
s3_out=s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north/IS
aws s3 cp $region_dir/ATL14_IS_0332_100m_006_02.nc $s3_out/
# Beside the quarterly tiles, the same place as on discover.
# DONE: all five 0332 files uploaded (the recommendation below taken), to a
# prefix that held only matched/ and prelim/ -- nothing was overwritten.
# VERIFIED: every bucket size equals its local size --
#   ATL14 100m 9915264; ATL15 1km 18117749, 10km 895146, 20km 701363,
#   40km 638098.  The ATL14 is the M0b re-write, so the published reference
#   DEM is final.
# RECOMMENDATION (taken): publish the four ATL15 files too, so the bucket holds
# the whole quarterly product; only ATL14 is needed for monthly.
#
#
# ===========================================================================
# M2. [ADE] [DONE 2026-09-18 -- bit-identical]  Read it back from the bucket the way the solver will.
# ===========================================================================
# pc.grid.mosaic().from_list(['<s3_out>/ATL14_IS_0332_100m_006_02.nc'],
#     group='', bounds=<E1340_N-2460 +/- 32 km>, fields=['h','h_sigma'])
# and compare with the same read of the local file: identical h and
# h_sigma.  Seconds in the ADE; it catches a URI-read failure before a DPS job
# spends its time finding it.  A worker's credentials are not the ADE's, so
# M6 is still the real test.
#
#
# ===========================================================================
# M3. [ADE] [DONE 2026-09-18]  check_field_sizes.py fractional -g.
# ===========================================================================
# Parse each -g entry as a/b the way make_mosaic_jobs.py and ATL15_write2nc.py
# already do, and add a test: -W=60000 -g=1250,2500,1/12 -t=2018.75,2026.5
# gives [25, 25, 94].  The quarterly expectation must stay [61, 61, 32].
# DONE: scripts/check_field_sizes.py gains _spacing(), which reads 'a/b' the
# way those two scripts do, and the derivation line now prints dt as the args
# file wrote it ('1/12', not '0.0833333').
#   monthly   -> [25, 25, 94]   (-W=60000 / 2500 + 1 = 25;  7.75 / 1/12 + 1 = 94)
#   quarterly -> [61, 61, 32]   unchanged
# Four tests added (tests/test_check_field_sizes.py, now 31): the monthly
# shape, the quarterly shape unchanged, a fraction that does not divide the
# span (1/7) raising CannotCheck, and a zero denominator raising CannotCheck
# rather than a traceback.  Suite 112 passed 2 skipped.
# REGRESSION CHECKED on the real quarterly tiles: 28 reports, 28 tiles,
# 28 of 28 passed, 0 problems.
#
#
# ===========================================================================
# M4. [ADE] [DONE 2026-09-18]  Compose the monthly args.
# ===========================================================================
setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
    default_args/IS.txt default_args/monthly.txt --Hemisphere=1 \
    --ATL14_reference_file=$s3_out/ATL14_IS_0332_100m_006_02.nc
monthly_dir=/home/jovyan/ATL14_processing/rel006/north_monthly/IS
# Writes $monthly_dir/input_args_IS.txt.  CHECK, against the quarterly file:
# the ONLY differences are -g=1250,2500,1/12, the monthly --dzdt_lags,
# --ATL14_reference_file (an s3:// URI) and -b.  No --hemi_suffix line.
# VERIFIED, sorted diff against rel006/north/IS/input_args_IS.txt -- exactly
# those four and nothing else:
#   + --ATL14_reference_file=s3://.../rel006/north/IS/ATL14_IS_0332_100m_006_02.nc
#   - --dzdt_lags=1,2,4,8,12,16,20,24,28   + --dzdt_lags=1,3,6,12,24,36,48,60,72,84
#   - -g=100,1000,0.25                     + -g=1250,2500,1/12
#   - -b=.../rel006/north/IS               + -b=.../rel006/north_monthly/IS
# --hemi_suffix is absent, as predicted; cycles 0332, Release 006, version 02,
# --t_crop=2019,2026.5 and the ATL11 release all carry through unchanged.
#
#
# ===========================================================================
# M5. [ADE] [DONE 2026-09-18]  Publish the args.
# ===========================================================================
s3_run_m=s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/north_monthly/IS
aws s3 cp $monthly_dir/input_args_IS.txt $s3_run_m/
# Then diff the bucket copy against the local one, as in T7.
# DONE: the prefix was empty; the bucket copy was pulled back and diffed
# against the local file -- IDENTICAL.
#
#
# ===========================================================================
# M6. [DPS] [DONE 2026-09-18 -- ALL FOUR GATES PASS]  Smoke one prelim tile.
# ===========================================================================
# NO REGISTRATION NEEDED, and Ben was told so.  STATEMENT, from git:
# everything committed since the registered build 61a19af is docs,
# region_files job lists, tests/test_time_coverage.py, and
# ATL1415_attrs_meta.py -- and ATL11_to_ATL15.py does not import attrs_meta,
# while set_time_range is called only from inside attrs_meta itself, by the
# two netCDF writers, which run in the ADE.  scripts/check_field_sizes.py is
# ADE-only too.  So the solve code in the image is unchanged and 61a19af is
# the right build to run monthly on.
# RECOMMENDATION unchanged: run check_build_id.py --expect 61a19af first, to
# confirm the image is still that commit before spending a job on it.
# ===========================================================================
s3_out_m=s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north_monthly/IS
echo "1340000 -2460000" > region_files/IS_0332_monthly_smoke_xy.txt
scripts/maap/submit_MAAP_jobs.py --xy_file region_files/IS_0332_monthly_smoke_xy.txt \
    --step prelim --args_url $s3_run_m/input_args_IS.txt \
    --tile_prefix $s3_out_m --queue maap-dps-worker-16gb \
    --tag IS_rel006_0332_monthly_prelim \
    --ledger ~/ATL14_processing/maap_ledgers/IS_0332_monthly_smoke_jobs.csv
# RECOMMENDATION: E1340_N-2460, the quarterly memory high-water tile
# (N_fit 273382, 9.52 GiB), so the smoke also sizes the queue.
# GATES:
#   a. the log names the s3:// reference file and N_fit is the same order as
#      quarterly's (a missing or empty reference would edit away every point
#      -- W5, which _expand_reference_files guards locally, not for a URI);
#   b. field-size report dz/dz [25, 25, 94], sigma_dz the same (M3);
#   c. /meta/lineage present, as in T6;
#   d. wall time and peak memory, which set M7's queue.
#
# BUILD GATE FIRST: check_build_id.py against the monthly args with
#   --expect 61a19af -> VERDICT MATCH.  Build stamp == live git in the image
#   == cwl commit == 61a19af, tree_state=clean, maap_py=5.1.0, and
#   maap_pgt=SET (the one that matters: unset means the worker falls back to
#   earthaccess with no credentials and cannot read ATL11 at all).  origin is
#   past the build at cb2796a -- expected, and NOT a reason to re-register:
#   everything pushed since 61a19af is docs, tests and ADE-only code.
#   TRAP, fallen into and then fixed: `check_build_id.py --help` SUBMITS A JOB
#   (job 2cfaac9e, 32gb queue, 4 min, wrote nothing).  There is no --help, and
#   the positional args_file swallowed it.  Guarded in 7d64824; see
#   howto_MAAP_ogc.sh O5.
#
# RESULT: job f0c16110-4874-4c0e-9947-bcaa41a03840, submitted 16:32:51Z,
#   SUCCESSFUL in 1060 s on maap-dps-worker-16gb, peak 4.05 GiB, 3 iterations,
#   commit 61a19af.  Steps: fit 853 s at 4.05 GiB, error 189 s at 1.88 GiB.
#   Ledger ~/ATL14_processing/maap_ledgers/IS_0332_monthly_smoke_jobs.csv,
#   list region_files/IS_0332_monthly_smoke_xy.txt.  The output prefix was
#   verified EMPTY before submitting, so the tile that appeared is this job's.
#   Submitted --dry-run first: exactly 1 job, right args URL, right prefix.
# a. PASSES.  The log carries
#      arg: --ATL14_reference_file=s3://.../rel006/north/IS/ATL14_IS_0332_100m_006_02.nc
#    and N_fit = 280955 against the QUARTERLY 273382 for this same tile
#    (plan_cycles_03_32.sh:331) -- 2.8% HIGHER, not the collapse a missing or
#    unreadable reference would cause.  dz/dz is 100% finite and sigma_dz
#    98.9%, which a NaN reference could not produce.  The only log warnings
#    are cwltool/GDAL boilerplate; the "cannot kill container" line is runner
#    cleanup after success.
# b. PASSES.  check_field_sizes.py $monthly_dir/prelim @$monthly_dir/input_args_IS.txt
#    -> 1 reports, 1 tiles, 1 of 1 passed, 0 problems, at the derived
#    dz/dz [25, 25, 94] with sigma_dz the same.  M3's code path, on real data.
# c. PASSES.  /meta/lineage present: 16 granules (14 along-track, 2
#    crossover), 0 NOT_SET attributes, 16 distinct uuids.  Fewer than the
#    quarterly product's 79 because that is the union over 28 tiles; this is
#    one tile's own inputs.
# d. MEASURED, and it settles M7's queue: 1060 s and 4.05 GiB, against the
#    quarterly's 3322 s and 9.52 GiB for the SAME tile -- 3.1x faster and
#    2.4x less memory, as the ~8x-fewer-unknowns estimate predicted.
#    maap-dps-worker-16gb has ample headroom; -32gb is NOT needed for M7.
#
#
# ===========================================================================
# M7. [DPS] [DONE 2026-09-18 -- 28/29 successful, 28 tiles; E1020 failed, see QM5]  Prelim fan-out, 29 centers.
# ===========================================================================
# region_files/IS_prelim_xy.txt, unchanged -- the centers come from the mask,
# not the period.  Tag IS_rel006_0332_monthly_prelim, ledger
# IS_0332_monthly_prelim_jobs.csv.  Then collect_jobs.py, fetch_tiles.py
# <ledger> $monthly_dir --step prelim, check_field_sizes.py $monthly_dir/prelim
# @$monthly_dir/input_args_IS.txt.  EXPECT E1020_N-2580 to write no tile again.
# DO NOT RE-REGISTER while jobs are queued (the split-build trap, I3).
# SUBMITTED 2026-09-18: 29/29 on 61a19af, maap-dps-worker-16gb (M6 gate d),
#   ledger ~/ATL14_processing/maap_ledgers/IS_0332_monthly_prelim_jobs.csv.
#   Dry-run first; committed and pushed (98beeb1) BEFORE submitting.
#
# FINDING -- E1020_N-2580 FAILED, and NOT the way the quarterly runs did.
#   STATEMENT (job 9266c3d7, triaged log
#   triaged_job-...-20260918T170628.949933Z_task-7133567b..., _exit_code 1):
#   read N=267958 (AT 267894, XO 64), then `smooth_fit: no valid data` and
#   exit 1 IN THE FIT STEP, 212 s at 0.73 GiB.  Quarterly (0331 and 0332):
#   the FIT succeeded with N_fit 327/349, and only the ERROR step found no
#   data, which the I7a branch turns into delete-the-tile, exit 0.
#   CAUSE, STATEMENT, measured: the reference DEM has ZERO finite h over this
#   tile's 60 km box (0 of 361201 cells), against 87.6% for the smoke tile.
#   The quarterly run never wrote a tile for this center, and at 40 km
#   spacing it has no neighbour to fill from (plan_IS_run.sh, 2026-09-16), so
#   the quarterly ATL14 is empty here.  The monthly solve subtracts z_ref from
#   every point; with z_ref all NaN, nothing is valid.
#   DETERMINISTIC -- same args, same build, same empty reference.  NOT
#   RETRIED, and a retry would fail identically.
#   OUTCOME IS THE SAME AS QUARTERLY: no tile for E1020_N-2580.  So the
#   matched set is 28 either way, and M8 (built from the tiles that EXIST)
#   is not blocked.  Only the job accounting differs: failed, not
#   successful-with-no-tile.
#   E1180_N-2380 (reference coverage only 0.34%, 1230 cells) SUCCEEDED at
#   N_fit 222, against 158 quarterly -- sparse, not broken.
#
# RESULT, all 29 terminal: 28 successful, 1 failed (E1020_N-2580, above).
#   ALL 29 ON ONE BUILD, 61a19af -- no split-build (collect_jobs "Builds that
#   ran these tiles": 61a19af, 29 tiles).
#   Wall 496..1071 s; peak 1.33..4.13 GiB; N_fit 222..280955; 3 iterations
#   on every tile.  Heaviest three are the high-N_fit centers again --
#   E1340_N-2460 (4.13 GiB, 280955), E1340_N-2500 (3.83, 259550),
#   E1300_N-2500 (3.75, 254509) -- so memory still tracks N_fit.  Against
#   quarterly 0332 prelim (1332..3322 s, 4.42..9.52 GiB): ~3x faster, ~2.3x
#   less memory.
# FETCHED: 27 fetched + 1 already local (E1340_N-2460, from M6), 1 FAILED.
#   M7 re-ran the smoke center, so the local copy (M6, job f0c16110) and the
#   bucket copy (M7) came from different jobs, 42 min apart.  STATEMENT,
#   compared: dz/dz, dz/sigma_dz, z0/z0 and z0/sigma_z0 are BIT-IDENTICAL
#   (np.array_equal, equal_nan=True; max |d| 0.0).  The solve is
#   deterministic on this build, so local == bucket in content and the copy
#   was kept.
# CHECKS, all pass:
#   - check_field_sizes.py: 28 reports, 28 tiles, 28 of 28 passed,
#     0 problems, at [25, 25, 94] with sigma_dz the same.
#   - the SAME 28 centers as the quarterly 0332 prelim (set-equal).
#   - lineage on all 28: none missing or empty, 0 NOT_SET, 79 distinct
#     granules = 69 along-track + 10 crossover -- EXACTLY the quarterly
#     product's 79, so the monthly run read the same inputs.
#   - bucket prelim/ and local prelim/ list the same 28 tiles.
#
# QM5. [OPEN -- for Ben]  E1020_N-2580 fails the monthly FIT for want of any
#      reference-DEM coverage.  What should happen to such a center?
#      WHY IT IS A QUESTION, not settled: Ben's 2026-09-16 decision ("tiles
#      that fail on the uncertainty step are not critical") covers the
#      ERROR step.  This is the same no-data exit (smooth_fit.py:485) on the
#      FIT step, which I7a does not handle.  The decision's own boundary says
#      not to widen it without asking.
#        A. Accept it: a failed job, no tile, nothing changed.  M8 proceeds
#           from the 28 tiles that exist.  No code, no rebuild.
#        B. Extend the I7a no-data branch to the fit step, so it exits 0
#           cleanly.  Solver change in ATL11_to_ATL15.py -> rebuild and
#           REGISTRATION, and only AFTER M7 finishes (the split-build trap).
#        C. Drop centers with no reference coverage from the monthly list
#           before submitting -- a pre-flight coverage check per center,
#           ADE-side, no rebuild.
#      RECOMMENDATION: A now, for IS -- the outcome is already correct and
#      IS has exactly one such center.  But C before M12: this is
#      STRUCTURAL, not an IS quirk.  Every center the quarterly product does
#      not cover will fail monthly the same way, and GL and AA have many
#      sparse edge tiles.  A pre-flight check turns a predictable failed job
#      into a center that is simply never submitted.
# AM5:
#
#
# ===========================================================================
# M8. [DPS] [DONE 2026-09-18 -- 28/28 successful, all checks pass]  Matched.
# ===========================================================================
# The list from the monthly prelim tiles that EXIST, bucket and local
# agreeing -> region_files/IS_0332_monthly_matched_xy.txt; submit --step
# matched with the same prefix; collect, fetch, check (--step matched; no
# sigma, by design).
# LIST BUILT 2026-09-18: region_files/IS_0332_monthly_matched_xy.txt, 28
#   centers, from the 28 prelim tiles that exist -- bucket and local listed
#   IDENTICAL first.  E1020_N-2580 absent.  The list is IDENTICAL to the
#   quarterly region_files/IS_0332_matched_xy.txt.
#   NOT BLOCKED BY QM5: E1020 writes no tile under all three of QM5's
#   options, so the matched list is the same whatever Ben answers.
#   Queue maap-dps-worker-16gb: quarterly matched peaked ~9 GiB against a
#   9.52 GiB prelim; monthly prelim peaked 4.13 GiB.
# RESULT: 28/28 SUCCESSFUL, 0 failed, all on 61a19af (no split build).
#   Dry-run first; matched/ prefix verified EMPTY before submitting.
#   Ledger ~/ATL14_processing/maap_ledgers/IS_0332_monthly_matched_jobs.csv.
#   Wall 85..307 s; peak 1.12..3.90 GiB; N_fit 222..280955; 1 iteration each
#   (as quarterly matched).  Heaviest: E1340_N-2460 3.90 GiB (N_fit 280955),
#   E1300_N-2460 3.57, E1300_N-2500 3.52 -- under the prelim peak, as in the
#   quarterly run.
# FETCHED 28/28, 0.40 GiB, to $monthly_dir/matched.
# CHECKS, all pass:
#   - check_field_sizes.py --step matched: 28 reports, 28 tiles, 28 of 28
#     passed, 0 problems -- dz/dz [25, 25, 94], dz/sigma_dz null (by design).
#   - bucket matched/ and local matched/: the same 28 names AND the same
#     sizes, tile for tile.
#
#
# ===========================================================================
# M9. [ADE] [DONE 2026-09-18 -- 44/44, PROBLEMS 0]  Mosaic.
# ===========================================================================
cd ~/ATL14_processing/runs
make_mosaic_jobs.py -b $monthly_dir -rr IS -t 2018.75,2026.5 -e ATL14 \
    --run_name IS_0332_monthly_mosaic @$HOME/git_repos/ATL1415/default_args/monthly.txt
cd IS_0332_monthly_mosaic
seq 1 $(ls queue | wc -l) | xargs -P 12 -I{} env SLURM_ARRAY_TASK_ID={} bash slurm_run.sh
check_mosaic_outputs.py ~/ATL14_processing/runs/IS_0332_monthly_mosaic --values
# No z0 task (skip_z0).  Check sigma_dzdt coverage against the values, as T8
# did: the quarterly tiles cover only 44% of the dzdt cells.
# QUEUE, inspected BEFORE running: 44 tasks = 3 x (avg_dz + 10 avg_dzdt lags)
#   + dz + 10 dzdt lags, exactly as predicted; NO z0 task (skip_z0 fired on the
#   1250 m z0 spacing).  make_mosaic_jobs.py INFERS the lags from -t and -g
#   rather than reading --dzdt_lags (I9c), so they were checked, not assumed:
#   the queue holds lag1,3,6,12,24,36,48,60,72,84 -- identical to the monthly
#   --dzdt_lags.  Values read from matched/*.h5, sigmas from prelim/*.h5
#   (matched carries no sigma, by design).
# RAN: all 44 at -P 12 in 74 s; xargs exit 0; queue 0, running 0, done 44,
#   error_logs 0.  Run dir ~/ATL14_processing/runs/IS_0332_monthly_mosaic.
# check_mosaic_outputs.py --values: 44 files, 136 fields, PROBLEMS 0.
#   Epoch counts follow nt - lag: lag84 has 10 of 94, lag72 22, lag60 34.
# SIGMA COVERAGE -- the quarterly gap is ABSENT here.  STATEMENT, from the
#   --values report:
#   - all 43 dzdt/avg files: sigma finite fraction == value finite fraction
#     EXACTLY (native grid 40.1% == 40.1%);
#   - dz.h5: dz 40.1%, sigma_dz 39.6% (98.8% of value cells have a sigma).
#     count/misfit_rms are 5.2% -- data-bearing cells only, as expected.
#   Quarterly 0332 at 1 km: sigma_dzdt 16.9% of the box vs 41.4% for the
#   values (T8).  HYPOTHESIS, NOT TESTED: the coarser monthly dz grid (2500 m
#   vs 1000 m) puts more data in each cell, so the error propagation reaches
#   all of them.  It would bear on the quarterly gap Ben has flagged but not
#   pursued -- recorded, not chased.
#
#
# ===========================================================================
# M10. [ADE] [NOT STARTED]  netCDF -- ATL15 only.
# ===========================================================================
mkdir -p ~/ATL14_processing/runs/IS_0332_monthly_nc && cd ~/ATL14_processing/runs/IS_0332_monthly_nc
ATL15_write2nc.py @$monthly_dir/input_args_IS.txt > ATL15.log 2>&1
# EXPECT ATL15_IS_0332_1mo_{2.5,10,20,40}km_006_02.nc, 91 epochs
# 2019.00..2026.50, no INVALID warning.  Checks as T8: lineage complete (the
# monthly prelim tiles carry their own), finite product == finite mosaic &
# ice_area > 0.
# RECOMMENDATION, the science check: monthly against the 0332 QUARTERLY at
# 10 km, at the quarterly epochs.  Monthly delta_h is relative to a
# DIFFERENT surface (the quarterly DEM), so compare delta_h differences
# between epochs, or add the reference back, rather than raw delta_h.
# Bar as Ben's: no >10 m errors, no major gaps.
#
#
# ===========================================================================
# M11. [ADE] [NOT STARTED]  Docs.
# ===========================================================================
#   - howto_MAAP_arctic.sh: a monthly section after step 10, citing M1-M10,
#     like the discover howto's monthly block.
#   - plan_cycles_03_32.sh T9 and this file: record what ran.
#
#
# ===========================================================================
# M12. [NOT STARTED, NOT PLANNED]  Other regions.
# ===========================================================================
#   - ORDER: a region's monthly run needs its quarterly ATL14 first.
#   - AA: on discover the reference is a GLOB over quadrants; a URI cannot be
#     globbed (_expand_reference_files raises), so AA needs every granule
#     named, one --ATL14_reference_file each -- and setup_ATL1415_region.py
#     takes a single value today.  [NEEDS CODE] when AA gets there.
