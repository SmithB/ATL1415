#! /usr/bin/env bash
# ===========================================================================
# PLAN: archive the old IS and AA outputs, rerun IS and the AA transect, run a
# south-to-north GL transect -- all on the new code -- and rebuild the timing
# budget (~/ATL14_processing/maap_resource_estimate.txt) from them.
# Written 2026-09-25, before anything was moved or submitted.  TENTATIVE.
# QT1-QT5 ANSWERED by Ben 2026-09-25 (below); A under way.
# Ben: "Archive the old outputs for Antarctica and Iceland in a separate
# directory for reference.  Rerun Iceland the Antarctic transect and run a
# south-to-north Greenland transect, and use these to repopulate the timing
# budget."
# ===========================================================================
# Provenance per claim: STATEMENT = verified 2026-09-25, with how;
# DECIDED = Ben said so; RECOMMENDATION = mine; QUESTION = open.
# Tags: [ADE] / [DPS] / [BEN]; [NOT STARTED] [NEEDS CODE: x] [DONE].
#
# QUESTIONS FOR BEN (answer on these lines):
#   QT1  Archive location and name.  RECOMMENDATION:
#          s3://maap-ops-workspace/ben_smith/ATL14_processing/archive/2026-09-25_spqr/
#          ~/ATL14_processing/archive/2026-09-25_spqr/
#        with the same rel006/north[_monthly]/IS layout below it, and
#        AA_transect_ab84687/ for the AA tiles.  IS is MOVED (the canonical
#        prefix must be empty before the rerun: matched reads neighbours from
#        it, and must not pick up an old tile); the AA tiles are COPIED out of
#        dps_output (the job records stay where the ledgers point).
#        QT1 answer:  Agree with recommendation
#   QT2  Iceland: rerun all four job types (quarterly + monthly, prelim +
#        matched) AND the ADE mosaic + netCDF, then compare the products with
#        the archive?  RECOMMENDATION: yes -- the comparison is the end-to-end
#        test of cholmod + the reach kernel, and the ADE steps are minutes.
#        QT2 answer: Agree with recommendation
#   QT3  Transects: which job types?  The budget needs prelim AND matched for
#        AA and GL, quarterly AND monthly.  A matched job on an isolated
#        transect tile has no neighbours, so it would run fast and mislead.
#        RECOMMENDATION: quarterly + monthly PRELIM on every transect tile,
#        plus ONE 3x3 block per region (the 8 neighbours of AA 60km_E420_N20
#        and of GL E200_N-1880) so one matched job per region per product
#        runs with a full neighbourhood.  Transect tiles go to a TRANSECT
#        prefix (.../ATL14_processing/transects/2026-09-25/<region>), not the
#        production tree.
#        QT3 answer:  Just run quarterly prelim for the transects
#          DECIDED: no monthly, no 3x3 blocks, no matched on AA/GL.  So AA/GL
#          matched and monthly budget lines are ESTIMATES (E2), not measured.
#   QT4  AA's args are still at 0331 (below), and 0331 is gone from CMR.
#        RECOMMENDATION: recompose both AA args files at 0332 the way IS was
#        (setup_ATL1415_region.py with rel_006_0332.txt), adding
#        --solver=cholmod.  The AA mask file covers 2018.00-2026.25 and the
#        0332 t_crop runs to 2026.5 (GL's mask ends 2026.0): is that
#        acceptable for a TIMING run?  (It is a science question for
#        production.)
#        QT4 answer:  The limited masks are OK for timing.
#          (Taken with the recommendation: recompose at 0332 + --solver=cholmod.)
#   QT5  The GL transect (D below, 19 tiles): OK?
#        QT5 answer: OK.
#   QT6  (ANSWERED 2026-09-25: -w; DONE, see B5.)  Should the z0.h5
#        mosaic tasks pass -w (weighted, like the 10 km dz/dzdt tasks)?
#        STATEMENT: make_mosaic_jobs.py writes the z0 tasks (all six matched
#        fields, and sigma_z0 from prelim) WITHOUT -w, but with -p 5000
#        -f 10000.  pointCollection's make_mosaic.py sets pad and feather to
#        None when -w is absent, so they are ignored; where tiles overlap, the
#        value comes from whichever file glob.glob lists last, i.e.
#        directory order = fetch order.  True since the python port (1861c6d).
#        STATEMENT (scratchpad rebuilds, today's code):
#          - archived tiles -> archived z0.h5 exactly (code did not change);
#          - new tiles -> new z0.h5 exactly;
#          - every one of the 28 tiles' z0/z0 matches the archive to <=4e-8 m;
#          - yet the two z0.h5 differ by up to 2053 m, and 585 m inside both
#            ice masks (509k of 1.13M masked cells differ by >1 cm).
#          - WITH -w, new vs archived tiles agree to 3.4e-8 m.
#        STATEMENT: worst ice cell x=1360100 y=-2510000: the two tiles whose
#        cores cover it (E1340_N-2500, E1380_N-2500, 10 km from centre) say
#        837-838 m; E1340_N-2540 and E1380_N-2540, 30 km from centre (their
#        outer edge), say 1423-1426 m.  The archive took 838, the new run 1423.
#        The other 40 mosaics match the archive to <=3.1e-7 (model and sigma).
#        So the old (archived) ATL14, and presumably every z0 made with this
#        script (discover too), carries edge values at random tile seams.
#        RECOMMENDATION: add -w to the z0 tasks (matched fields and
#        sigma_z0), ADE-only (no registration), then rebuild IS z0.h5 and
#        carry on with B5.  The old product cannot then be compared on z0;
#        compare new-weighted vs archive-tiles-weighted instead (3.4e-8 m above).
#        QT6 answer: The -w flag is needed.  Note this as something that has been fixed that will lead to minor differences between new products and archived products.
#          DONE: make_mosaic_jobs.py z0 tasks (6 matched fields + sigma_z0) pass
#          -w; tests/test_make_mosaic_jobs.py (fails without it); suite 181
#          passed 2 skipped.  Noted as fixed in howto_MAAP_arctic step 9.
#          STATEMENT: the differences vs archived products are small in the
#          median but NOT small at seams -- IS z0 on ice: median 1.7 cm,
#          18.7k of 1.13M cells > 10 m, max 412 m (sigma_z0 max 193 m).
#
#
# ===========================================================================
# WHAT EXISTS NOW.  STATEMENT, 2026-09-25 (ls, aws s3 ls).
# ===========================================================================
# IS quarterly  local ~/ATL14_processing/rel006/north/IS: 28 prelim, 28
#   matched, 41 mosaic .h5, 5 .nc, input_args_IS.txt (1.4 GB).  S3
#   .../ATL14_processing/rel006/north/IS: 117 objects, 1.23 GB (prelim/,
#   matched/, the 5 .nc).  Build 61a19af, SPQR, old error kernel.
# IS monthly    local .../north_monthly/IS (1.0 GB) and S3 (116 objects,
#   0.98 GB), same build.
# BOTH IS args files were changed 2026-09-25 (--solver=cholmod added); the
#   archive gets the args as they were when the outputs were made (the
#   pre-change copies are in the session scratchpad; the only difference is
#   that one line).
# AA            NO outputs in the canonical tree, local or S3 -- only
#   rel006/south/AA{,_44km}/input_args_AA*.txt.  The 17-tile cost transect
#   (scripts/maap/AA_cost_results.csv, 2026-09-11, build ab84687, 0331) left
#   its tiles in dps_output; ledger maap_ledgers/AA_transect_ab84687_jobs.csv
#   (and AA_xo_check_jobs.csv for the crossover check jobs: archived too).
# AA args (local == run_args on S3): --cycles=0331, --version=01,
#   --ATL11_release=007_cycle_03_31_v04, --t_crop=2019,2026.25.  The ATL11
#   index on S3 is 0332 only (ATL11_index_0332_007_05, both hemispheres).
# GL            never run on MAAP; no args composed.  docs/howto_MAAP_GL.sh
#   steps 1-2 compose them; its masks are staged.
# Deployed build d59b140 (MATCH 2026-09-25): cholmod opt-in, reach kernel,
#   threads = physical cores, WORKER/cpu lines.  Nothing in flight.
# Local disk: 269 GB free.
#
#
# ===========================================================================
# A. ARCHIVE.  [ADE]
# ===========================================================================
# A DONE 2026-09-25 (VERIFY OK): s3 .../ATL14_processing/archive/2026-09-25_spqr/
#   (370 objects, 5.79 GB) and ~/ATL14_processing/archive/2026-09-25_spqr/
#   (2.3 GB): rel006/north{,_monthly}/IS (MOVED; args as used, i.e. without
#   --solver), AA_transect_ab84687/<identifier>/ and AA_xo_check/<identifier>/
#   (whole job prefixes COPIED from dps_output: tile, report, logs),
#   AA_transect_ab84687/AA_cost_results.csv, AA_args_0331/.  MANIFEST.json in
#   both.  Canonical IS prefixes now hold only the live args files.
# A1. [DONE]  Manifest first: list every object/file to be moved or
#     copied, with sizes, into the archive's MANIFEST.txt.
# A2. [DONE]  S3 IS: aws s3 mv --recursive (never the FUSE mount) each
#     canonical prefix into the archive; the args as used go alongside.
# A3. [DONE]  Local IS: mv the region dirs' contents (not the new args
#     files, which stay for the rerun) into the local archive.
# A4. [DONE]  AA transect: aws s3 cp each of the 17 tiles (and each
#     job's _stderr.txt, the timing record) from its dps_output prefix into
#     AA_transect_ab84687/<half>/.  Plus scripts/maap/AA_cost_results.csv.
# A5. [DONE]  Verify against the manifest: counts and sizes equal,
#     the canonical IS prefixes empty.  Only then does B start.
#
#
# ===========================================================================
# B. ICELAND RERUN.  [DPS + ADE]  docs/howto_MAAP_arctic.sh, IS only.
# ===========================================================================
# B STATUS 2026-09-25 ~04:45Z: B1 and B3 SUBMITTED together (28 + 28, -16gb,
#   canonical prefixes, build d59b140, --solver=cholmod in both args files):
#   ledgers ~/ATL14_processing/maap_ledgers/IS_0332_cholmod_prelim_jobs.csv
#   and IS_0332_monthly_cholmod_prelim_jobs.csv.  DO NOT REGISTER while they
#   (or B2/B4) run.  When done: collect_jobs.py <ledger>; fetch_tiles.py
#   <ledger> ~/ATL14_processing/rel006/north[_monthly]/IS --step prelim;
#   check_field_sizes.py; then B2/B4 matched (same tile list, same prefix).
#   Helper scripts: ~/ATL14_processing/session_tools_2026-09-25/.
# B RESULTS 2026-09-25 ~05:30Z:
#   B1 quarterly prelim: 27/28 successful, fetched, 27/27 field sizes OK; vs
#     the archive: 0 edit flips on every tile, worst model 1.6e-7 m
#     (E1340_N-2500), worst sigma 8.5e-6 rel.  E1300_N-2620 FAILED in MAAP's
#     own get_maap_pgt_token.py (connect timeout to api.maap-project.org; no
#     NSIDC credentials) -> retried, ledger IS_0332_cholmod_prelim_retry_jobs.csv.
#   B3 monthly prelim: 28/28 FAILED -- MY MISTAKE: A moved the quarterly
#     ATL14_IS_0332_100m_006_02.nc, which the monthly args name as
#     --ATL14_reference_file, and monthly was launched before the quarterly
#     ATL14 exists (howto order: quarterly prelim -> matched -> mosaic+nc ->
#     monthly).  FIX: rerun B3 AFTER B5 writes the NEW quarterly ATL14 (not a
#     copy of the old one back).  pointCollection's mosaic.from_list reports
#     the missing file as UnboundLocalError 'temp' (upstream, unfixed).
# B2 RESULT 2026-09-25 ~15:10Z: quarterly matched 28/28.  27 first time;
#   E1340_N-2420 FAILED in MAAP's stage_in (connect timeout to
#   api.maap-project.org, the same infra fault as the prelim retry), retried
#   alone (ledger IS_0332_cholmod_matched_retry_jobs.csv), successful.  All 28
#   fetched, 28/28 field sizes OK; vs the archive: 0 edit flips on every tile,
#   worst model 2.0e-7 m (E1380_N-2460).  Nothing in flight after B2.
# B5 STATUS ~15:30Z: quarterly mosaic run (runs/IS_0332_cholmod_mosaic -- a
#   NEW run name; the old runs/IS_0332_mosaic is untouched): 41/41 tasks, 56 s,
#   no error logs, check_mosaic_outputs --values 0 problems.  40/41 files
#   match the archive; z0.h5 does NOT -> QT6.  PAUSED: no netCDF written,
#   nothing published, monthly (B3) still waits.
# B5 MOSAIC RERUN ~16:00Z with -w (runs/IS_0332_cholmod_w_mosaic): 41/41,
#   42 s, no error logs, 0 problems.  z0.h5 vs the ARCHIVED TILES mosaicked
#   the same way (-w, scratchpad): model 3.4e-8 m, sigma 4.9e-7 rel, no NaN
#   flips.  The other 40 files vs the archive: <=1.5e-7.  NEXT: netCDF
#   (howto step 10), compare, publish (step 11), then B3 monthly.
# B5 netCDF ~17:40Z (runs/IS_0332_cholmod_nc): ATL14 + 4 ATL15, rc 0, no
#   INVALID.  vs archive: ATL15 all <=6.1e-5 (float32 rounding).  ATL14 h /
#   h_sigma differ exactly as QT6 predicts (18,721 cells > 10 m, max 412 m;
#   23 NaN flips); data_count / misfit_* gain or lose 3,632 NaN cells (now
#   weighted too).  PUBLISHED (step 11), S3 sizes == local.  Monthly args
#   (S3 == local) name this ATL14 and --solver=cholmod.
# B3 RERUN SUBMITTED ~17:45Z: 28 monthly prelim, -16gb, on d59b140 (see D),
#   ledger IS_0332_monthly_cholmod_prelim_rerun_jobs.csv; both monthly
#   prefixes were empty.
# B3 RESULT ~18:30Z: 28/28 successful (all d59b140), fetched, 28/28 field
#   sizes OK.  vs the archive, REPORTED CELLS ONLY (cell_area > 0, Ben
#   2026-09-25; zero-area cells never reach a product): up to 22 m and 31%
#   sigma, except E1020_N-2420 (isolated, no seams) at 1.8e-12 m.
#   STATEMENT: E1260_N-2540 rerun locally against the ARCHIVED ATL14 matches
#   the archived tile to 2.1e-12 m, sigma 1e-14, 0/25693 edit flips.  So the
#   monthly differences are ALL from the new reference ATL14 (QT6 z0 -w fix),
#   none from the code.  Checker: session_tools_2026-09-25/compare_reported.py.
# B4 SUBMITTED ~18:35Z: 28 monthly matched, -16gb, ledger
#   IS_0332_monthly_cholmod_matched_jobs.csv.
# B4 RESULT ~20:15Z: 28/28 successful (all d00568c), fetched, 28/28 field
#   sizes OK.  vs archive on reported cells: E1020_N-2420 1.3e-12 m, worst
#   14.7 m (E1380_N-2500) -- the reference-ATL14 change, as for B3.
# B5 MONTHLY DONE ~20:25Z: mosaic runs/IS_0332_cholmod_monthly_mosaic (44/44,
#   41 s, lags 1 3 6 12 24 36 48 60 72 84, no z0 task, 0 problems); ATL15
#   runs/IS_0332_cholmod_monthly_nc (4 files, rc 0, no INVALID).  vs archived
#   monthly ATL15 delta_h: max 9.3 m (2.5 km) / 3.0 (10 km) / 1.2 (20 km) /
#   1.1 (40 km), none > 10 m, no NaN flips; delta_h_sigma max 0.026 m.
#   Monthly - new quarterly (10 km): median -0.067 m; >10 m in 20 values in
#   the same 6 weak cells at the record ends (archive: 22 values, 6 cells).
#   PUBLISHED (step 18), S3 sizes == local.  B (IS rerun) COMPLETE.
# B1. quarterly prelim, 29 centers from ATL1415/resources/IS/40km_tile_list.txt
#     (E1020_N-2580 is out already), -16gb, --tile_prefix = canonical.
# B2. quarterly matched.  B3. monthly prelim.  B4. monthly matched.
# B5. (QT2) mosaic + netCDF, quarterly and monthly; compare every product
#     with the archive: model fields <= 1e-3 m (QC1), sigma <= 5%.
#
#
# ===========================================================================
# C. ANTARCTIC TRANSECT.  [ADE + DPS]
# ===========================================================================
# C STATUS 2026-09-25 ~05:00Z: C1 DONE -- both AA args recomposed at 0332
#   (diff vs 0331: exactly --cycles, --version, --ATL11_release, --t_crop),
#   + --solver=cholmod, published to run_args (S3 == local).  C2 SUBMITTED:
#   13 (60km) + 4 (44km) prelim, -32gb, tile prefixes
#   .../ATL14_processing/transects/2026-09-25/AA and .../AA_44km (separate:
#   E420_N20 is in both halves); ledgers AA_transect_cholmod_{60km,44km}_jobs.csv.
# C1. (QT4) recompose input_args_AA.txt and input_args_AA_44km.txt at 0332
#     + --solver=cholmod; publish to run_args (the 0331 files go to the archive).
# C2. The same 17 tiles as 2026-09-11, same halves, -32gb (their peaks were
#     up to 21.5 GiB), quarterly prelim; C3 monthly prelim; C4 (QT3) the 3x3
#     block around 60km_E420_N20 and its matched job(s).
#     QT3: C3 and C4 DROPPED -- quarterly prelim only.
#
#
# ===========================================================================
# D. GREENLAND SOUTH-TO-NORTH TRANSECT.  [ADE + DPS]
# ===========================================================================
# D STATUS 2026-09-25 ~05:00Z: D1 DONE -- default_args/GL_latest.txt ->
#   GL_0331.txt committed; input_args_GL.txt composed (all masks on S3,
#   checked), + --solver=cholmod, published.  D2 SUBMITTED: 19 prelim, -32gb,
#   prefix .../transects/2026-09-25/GL; ledger GL_transect_cholmod_jobs.csv.
#   FIRST GL SOLVES EVER ON MAAP (Gr1km-v2 tides, geotiff mask, scaling maps).
# D2 RESULT (~05:40Z): 18/19 successful.  E480_N-1040 (79N tongue, the tide
#   tile) FAILED in pyTMD extrapolation: Gr1km-v2's zarr store has
#   inconsistent chunks along x -> "Object has inconsistent chunks" from
#   Dataset.chunks.  Deterministic; reproduced locally.  FIXED 2c64b87
#   (tides.py: ds.unify_chunks(); values unchanged).  NEEDS A REGISTRATION,
#   then retry E480_N-1040 alone -- only when NO jobs are in flight.
# D STATUS ~17:40Z: Ben reports the registration done, but list_algorithms
#   still shows atl1415_tile_solve:on_s3 modified 03:01Z and check_build_id
#   says the image is d59b140 (MATCH, pre-fix).  E480_N-1040 NOT retried:
#   it would fail the same way.  Retry once check_build_id shows >= 2c64b87.
#   pyTMD bug write-up for Ben: ~/ATL14_processing/pyTMD_inconsistent_chunks_bug.txt
# D DONE ~18:10Z: registration deployed (check_build_id MATCH d00568c, built
#   17:46Z; process lastModifiedTime still read 03:01Z -- not a deploy signal).
#   E480_N-1040 retried alone (GL_transect_cholmod_retry_jobs.csv): successful,
#   778 s, peak 9.42 GiB, m5.2xlarge.  GL transect 19/19.
# D1. howto_MAAP_GL steps 1-2: compose + publish input_args_GL.txt (and the
#     monthly one), + --solver=cholmod.
# D2. 19 tiles, -32gb (the howto's choice until GL memory is measured):
#     the middle tile of every 4th row (160 km) from the south tip to the
#     north coast, plus the howto's two smoke tiles (E200_N-1880 near Summit,
#     E480_N-1040 on the 79N tongue):
#       E80_N-3320 E0_N-3160 E-40_N-3000 E-40_N-2840 E40_N-2680 E120_N-2520
#       E160_N-2360 E240_N-2200 E280_N-2040 E240_N-1880 E200_N-1880
#       E200_N-1720 E200_N-1560 E40_N-1400 E0_N-1240 E80_N-1080 E480_N-1040
#       E80_N-920 E200_N-760
#     (all 19 checked present in ATL1415/resources/GL/40km_tile_list.txt).
#     Quarterly prelim; D3 monthly prelim; D4 (QT3) the 3x3 block around
#     E200_N-1880 and its matched job(s).
#     QT3: D3 and D4 DROPPED -- quarterly prelim only.
#
#
# ===========================================================================
# E. THE TIMING BUDGET.  [ADE]
# ===========================================================================
# E NOTE 2026-09-25 (TABLED by Ben): where AA 44km_E220_N20's fit time goes.
#   STATEMENT (DPS log): fit 3488 s, cpu 1317 s (0.38 cores); first solve at
#   ~52 min; 3 cholmod solves ~270 s wall.  So ~3200 s before the solve,
#   mostly idle -> predicted mostly data reads (ATL11/xover/previous product
#   from S3), with ~800-900 cpu-s of pre-solve compute.  NOT profiled: a local
#   cProfile to the first solve (scratchpad prof_e220/run_prof.py) was stopped
#   unfinished.  Error step: 1691 of 1703 s is uncertainty propagation, peak
#   20.0 GiB (was 11.6 on SPQR/old kernel) -- a memory note for the budget.
# E STATUS 2026-09-25 ~20:50Z: E1-E3 DONE, DRAFT FOR BEN'S REVIEW, NOT SENT.
#   ~/ATL14_processing/maap_resource_estimate.txt rewritten (previous kept as
#   maap_resource_estimate_2026-09-19.txt).  Discover logs are GONE (/tmp was
#   wiped), so E2 used old/new MAAP ratios on the same tiles instead of
#   new-MAAP/discover: IS q prelim 4.0x, q matched 4.3x, m prelim 2.0x,
#   m matched 2.4x; AA q prelim 2.7x; GL q prelim measured directly (median
#   12 min).  AA/GL matched + monthly scaled by the IS ratio, marked est.
#   TOTAL ~6,000-7,500 job-hours (was ~18,000); AA quarterly ~2 days at 100
#   concurrent (was 5-6).  Ratio script: ~/ATL14_processing/session_tools_2026-09-25/budget_ratios.py.
# E NOTE 2026-09-25 ~23:30Z: ATL11 READ VOLUME -- the draft's section 1 is WRONG.
#   Ben doubted the "<= ~10 TB" input estimate (tiles are 68-91% copied data).
#   STATEMENT (local, eth0 + per-object s3fs counts, fit to first solve):
#     IS E1340_N-2460  1.71 GB (ATL11 read 3.3x the granules' 0.50 GB)
#     GL E200_N-1880   4.37 GB (2.4x of 1.78 GB)
#     AA 60km_E900_N20 9.90 GB (1.4x of 6.74 GB)
#   CAUSE: pointCollection DOES merge index ranges (query_xy cleanup); the
#   waste is inside a range read -- fields' compressed chunks interleave and
#   open_remote's 256 KiB readahead cache keeps one block, so blocks are
#   refetched.  FIX MEASURED (fs.open cache_type='blockcache'), identical
#   solve inputs (sha of A and b):  IS 141->51 s, 1.71->0.26 GB;  GL 252->109 s,
#   4.37->0.58 GB;  AA 525->279 s, 9.90->1.54 GB.  Whole-file download is
#   faster on IS/GL but pulls the previous ATL14/15 whole (AA 20 GB): NOT rec.
#   RECOMMENDATION: blockcache in pointCollection io_utils.open_remote; then
#   rewrite estimate section 1 (AA ~180 TB as deployed, ~30 TB with the fix).
#   AWAITS BEN.  Scripts: session_tools_2026-09-25/read_fix.py, read_attrib.py.
#   DONE 2026-09-26 (Ben: "add blockcache to pointCollection open_remote"):
#   pointCollection branch blockcache_reads 2d438e1 (pushed; BEN TO MERGE --
#   DPS installs pointCollection from main), + ATL1415 read_ATL11.py lineage
#   reopen now passes DEFAULT_REMOTE_BLOCK_SIZE (s3fs's default block is
#   50 MiB; that reopen alone was ~0.74 GB of the IS tile's 1.0 GB).  With
#   BOTH, committed code paths: IS 141->53 s, 1.71->0.26 GB; AA E900
#   525->262 s, 9.90->1.54 GB; identical solve inputs, all lineage read.
#   CORRECTION to the E NOTE above: the refetching is INSIDE one range read
#   (45 range reads for the IS tile, one per granule and pair; ~30 MB each
#   on readahead, ~7 MB with blockcache), not many ranges.  Estimate section 1
#   rewritten for the fixed code.  Deploy = merge + register + check_build_id.
# E1. collect_jobs on every ledger: time per step, peak RSS, instance type,
#     lifecycle, cores got, steal.  Workers are SPOT and of mixed type, so
#     every number is reported with the instance type it ran on.
# E2. Method as before (maap_resource_estimate.txt): MAAP time / discover
#     time per tile, median per region and job type, times the discover
#     totals -- now with the new code's MAAP times.  Old-vs-new on the SAME
#     tiles (IS, AA transect) is reported alongside.
#     QT3: AA/GL matched and monthly are not run.  They are ESTIMATED as
#     discover time x (IS new MAAP / discover) for the same job type, and
#     labelled as estimates in the budget.
# E3. Rewrite sections 2 (CPU) and 4 (time) and the memory section; keep
#     the 2026-09-19 numbers as history.  Ben reviews before it is sent.
# ===========================================================================
