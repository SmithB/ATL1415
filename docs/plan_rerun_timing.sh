#! /usr/bin/env bash
# ===========================================================================
# PLAN: archive the old IS and AA outputs, rerun IS and the AA transect, run a
# south-to-north GL transect -- all on the new code -- and rebuild the timing
# budget (~/ATL14_processing/maap_resource_estimate.txt) from them.
# Written 2026-09-25, before anything was moved or submitted.  TENTATIVE.
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
#        QT1 answer:
#   QT2  Iceland: rerun all four job types (quarterly + monthly, prelim +
#        matched) AND the ADE mosaic + netCDF, then compare the products with
#        the archive?  RECOMMENDATION: yes -- the comparison is the end-to-end
#        test of cholmod + the reach kernel, and the ADE steps are minutes.
#        QT2 answer:
#   QT3  Transects: which job types?  The budget needs prelim AND matched for
#        AA and GL, quarterly AND monthly.  A matched job on an isolated
#        transect tile has no neighbours, so it would run fast and mislead.
#        RECOMMENDATION: quarterly + monthly PRELIM on every transect tile,
#        plus ONE 3x3 block per region (the 8 neighbours of AA 60km_E420_N20
#        and of GL E200_N-1880) so one matched job per region per product
#        runs with a full neighbourhood.  Transect tiles go to a TRANSECT
#        prefix (.../ATL14_processing/transects/2026-09-25/<region>), not the
#        production tree.
#        QT3 answer:
#   QT4  AA's args are still at 0331 (below), and 0331 is gone from CMR.
#        RECOMMENDATION: recompose both AA args files at 0332 the way IS was
#        (setup_ATL1415_region.py with rel_006_0332.txt), adding
#        --solver=cholmod.  The AA mask file covers 2018.00-2026.25 and the
#        0332 t_crop runs to 2026.5 (GL's mask ends 2026.0): is that
#        acceptable for a TIMING run?  (It is a science question for
#        production.)
#        QT4 answer:
#   QT5  The GL transect (D below, 19 tiles): OK?
#        QT5 answer:
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
#   (and AA_xo_check_jobs.csv for the crossover check jobs -- archive too?).
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
# A1. [NOT STARTED]  Manifest first: list every object/file to be moved or
#     copied, with sizes, into the archive's MANIFEST.txt.
# A2. [NOT STARTED]  S3 IS: aws s3 mv --recursive (never the FUSE mount) each
#     canonical prefix into the archive; the args as used go alongside.
# A3. [NOT STARTED]  Local IS: mv the region dirs' contents (not the new args
#     files, which stay for the rerun) into the local archive.
# A4. [NOT STARTED]  AA transect: aws s3 cp each of the 17 tiles (and each
#     job's _stderr.txt, the timing record) from its dps_output prefix into
#     AA_transect_ab84687/<half>/.  Plus scripts/maap/AA_cost_results.csv.
# A5. [NOT STARTED]  Verify against the manifest: counts and sizes equal,
#     the canonical IS prefixes empty.  Only then does B start.
#
#
# ===========================================================================
# B. ICELAND RERUN.  [DPS + ADE]  docs/howto_MAAP_arctic.sh, IS only.
# ===========================================================================
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
# C1. (QT4) recompose input_args_AA.txt and input_args_AA_44km.txt at 0332
#     + --solver=cholmod; publish to run_args (the 0331 files go to the archive).
# C2. The same 17 tiles as 2026-09-11, same halves, -32gb (their peaks were
#     up to 21.5 GiB), quarterly prelim; C3 monthly prelim; C4 (QT3) the 3x3
#     block around 60km_E420_N20 and its matched job(s).
#
#
# ===========================================================================
# D. GREENLAND SOUTH-TO-NORTH TRANSECT.  [ADE + DPS]
# ===========================================================================
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
#
#
# ===========================================================================
# E. THE TIMING BUDGET.  [ADE]
# ===========================================================================
# E1. collect_jobs on every ledger: time per step, peak RSS, instance type,
#     lifecycle, cores got, steal.  Workers are SPOT and of mixed type, so
#     every number is reported with the instance type it ran on.
# E2. Method as before (maap_resource_estimate.txt): MAAP time / discover
#     time per tile, median per region and job type, times the discover
#     totals -- now with the new code's MAAP times.  Old-vs-new on the SAME
#     tiles (IS, AA transect) is reported alongside.
# E3. Rewrite sections 2 (CPU) and 4 (time) and the memory section; keep
#     the 2026-09-19 numbers as history.  Ben reviews before it is sent.
# ===========================================================================
