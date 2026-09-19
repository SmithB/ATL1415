#! /usr/bin/env bash
# ===========================================================================
# PLAN: move the rel006 runs from ATL11 cycles 03-31 to cycles 03-32.
# Written 2026-09-17.  T1 is DONE.  QT1-QT5 ANSWERED by Ben the same day
# (AT1-AT5, at the foot of this file) and folded into the steps below.
# ===========================================================================
# WHY NOW (Ben, 2026-09-17, AL1 in docs/plan_lineage_at_solve_time.sh): "After
# the current round of code changes are complete we will transision to cycles
# 03-32."  And it is not optional --
# STATEMENT, CMR, 2026-09-17: the 0331_007_04 generation the IS run used is NO
# LONGER LISTED.  A bounding-box search over Iceland returns 55 granules, all
# 0332_007_05; an exact granule_name search for ATL11_023003_0331_007_04.h5
# returns nothing.  read_ATL11_at filters to --ATL11_release, so on today's
# args EVERY tile would find zero along-track data.  Nothing can be solved --
# IS, GL or AA -- until T2 lands.
#
# Provenance per claim: STATEMENT = verified, with how; DECIDED = Ben said so;
# RECOMMENDATION = mine, overridable; QUESTION = open, not guessed.
#
# WHAT DOES NOT CHANGE, all STATEMENT 2026-09-17:
#   - CROSSOVERS.  Cycles 1 and 2 are still ATL11XO 007_03 (14 of them over
#     the Iceland box); the 007_04 crossovers that exist are cycles 30 and 31
#     only, which this run does not read ([[xover-cycles-1-2-by-design]]).
#     --ATL11xo_version=007_cycle_03_30_v03 stays as it is.
#     DECIDED (Ben, AT5): "Leave these as they are" -- the crossover generation
#     and the previous product (005_0329) both stay put.
#   - masks, geoid, tides, the previous product (005_0329, AT5), the tile centers
#     (region_files/IS_prelim_xy.txt, from the 40 km mask), the queues.
#   - the code: this is an ARGS and DATA change.  The rebuild in T5 is needed
#     for the lineage work (plan_lineage_at_solve_time.sh L1-L5), not for this.
#
#
# ===========================================================================
# T1. [ADE] [DONE 2026-09-17]  Stage the 0332_007_05 index.
# ===========================================================================
# Ben delivered /tmp/ATL11_007_cycle_03_32_v05_index.tgz (203518176 bytes).
# ARCHIVE LAYOUT, read before extracting: 8102 members, no top-level directory
# of its own --
#     north/index/ATL11_<rgt><region>_0332_007_05.h5   3944 granule indexes
#     south/index/ATL11_<rgt><region>_0332_007_05.h5   4156 granule indexes
#     north/index/GeoIndex.h5, south/index/GeoIndex.h5
# THE DELIVERED LAYOUT IS NOT THE ONE THE CLOUD READS.  index_path_for_granule
# (pointCollection/scripts/query_ATL11_cloud.py) builds
#   <index_root>/ATL11_index_<cycles>_<release>_<version>/<granule basename>
# which for this generation is ATL11_index_0332_007_05/ -- flat, both
# hemispheres together, exactly as ATL11_index_0331_007_04/ already is
# (8100 objects, 2136377505 bytes, listed).
#
# WHAT WAS DONE, and it is all verified:
#   a. extracted to /home/jovyan/ATL11_index_staging/ (2.5 GB, 67 s);
#   b. checked the two hemispheres share NO basename (comm over the two
#      listings: nothing in common), so flattening loses nothing -- the same
#      property the 0331 staging relied on;
#   c. flattened with mv into ATL11_index_0332_007_05/: 8100 files, 2.1 GB;
#   d. compared one index against its 0331 counterpart: same structure, one
#      'index' group, and file_0 = '../ATL11_000103_0332_007_05.h5:pair1' --
#      the same relative form with the :pair1 suffix that
#      read_ATL11_granule_cloud_items strips before matching;
#   e. uploaded to
#        s3://maap-ops-workspace/ben_smith/ATL11_index/ATL11_index_0332_007_05/
#      8100 objects, 2136402126 bytes.  EVERY file's size on the bucket was
#      diffed against the local copy: no differences;
#   f. END-TO-END CHECK, the one that matters: of the 55 ATL11 granules CMR
#      returns for the Iceland box at --ATL11_release=007_cycle_03_32_v05,
#      ALL 55 resolve to an index that exists on the bucket (path_exists on
#      index_path_for_granule).  0 missing.  A missing index is a hard error
#      in read_ATL11_granule_cloud_items, so this is the check that says the
#      staging is usable;
#   g. hemisphere manifest written, mirroring the 0331 one but NOT overwriting
#      it -- the 0331 files stay where they are:
#        .../ATL11_index/hemisphere_manifest/0332_007_05/{README,north,south}.txt
#   h. GeoIndex.h5 (one per hemisphere) was NOT uploaded: it is the
#      whole-archive index for LOCAL reads, the cloud path never opens it, and
#      the 0331 staging left it out too.  Both are kept in the staging dir.
# STILL ON THE ADE: /home/jovyan/ATL11_index_staging/ (2.5 GB, the extracted
# archive plus the flattened tree).
# DECIDED (Ben, AT4): "Delete once T6 has solved" -- so after T6's smoke tile
# succeeds, delete BOTH that staging directory and the superseded 0331 index
# on the bucket:
#   rm -r /home/jovyan/ATL11_index_staging
#   aws s3 rm --recursive s3://maap-ops-workspace/ben_smith/ATL11_index/ATL11_index_0331_007_04/
# NOT BEFORE T6: until a tile has actually been solved on 0332, the 0331 index
# is the only staged one, and re-staging it would mean a new archive from Ben.
# The hemisphere_manifest/ files are tiny; leave both generations' manifests.
# DONE 2026-09-18 by Ben, after T6 passed.  Verified: the staging dir is gone,
# ATL11_index_0331_007_04/ lists 0 objects, and the bucket's ATL11_index/
# holds only ATL11_index_0332_007_05/ (8100 objects, 2136402126 bytes -- as
# uploaded in e) and hemisphere_manifest/.
# [[maap-bucket-is-mountpoint-s3]]: this was `aws s3 cp`, not the mount.
#
#
# ===========================================================================
# T2. [ADE] [DONE 2026-09-17]  The release args file.
# ===========================================================================
# STATEMENT: default_args/latest_release.txt is a SYMLINK to rel_006_0331.txt,
# and that file carries every value this transition touches:
#     --cycles=0331            names the products (ATL14_IS_0331_..._006_01.nc)
#     --ATL11_release=007_cycle_03_31_v04      the CMR/index generation filter
#     -t=2018.75,2026.5        the solved time span; sets the epoch count
#     --t_crop=2019,2026.25    the span written to the netCDFs
# DECIDED (Ben 2026-09-17, AT1 and AT2).  Copy it to
# default_args/rel_006_0332.txt, change these, and repoint the symlink --
# leaving rel_006_0331.txt intact as the record of what the first IS run used:
#     --cycles=0332
#     --ATL11_release=007_cycle_03_32_v05
#     --version=02                              (AT2 "Use release 02")
#     -t=2018.75,2026.5        UNCHANGED        } AT1
#     --t_crop=2019,2026.5     was 2019,2026.25 }
# WHAT THAT MEANS, so nobody reads a surprise as a fault:
#   - -t is UNCHANGED, so every tile still has 32 epochs
#     (nt = (2026.5-2018.75)/0.25 + 1) and every field size and field-size
#     report stays as it is.  scripts/check_field_sizes.py keeps deriving
#     [61, 61, 32] for IS.
#   - --t_crop NOW REACHES THE TOP OF THE SOLVED SPAN.  The netCDFs get 31
#     epochs rather than the 30 the 0331 products had, and the last of them is
#     the final solved epoch, at the edge of the fit rather than one dt inside
#     it.  Edge epochs are the least constrained ones, so expect the last
#     delta_h slice to be noisier than its neighbours; that is the chosen
#     trade, not a defect.
#   - the products become ATL14_IS_0332_100m_006_02.nc and
#     ATL15_IS_0332_3mo_<res>_006_02.nc.
#   - STATEMENT, for the record: the data reach 2026.489
#     (ATL11_023003_0332_007_05.h5, ancillary_data/end_delta_time = 267933889 s
#     by the solver's conversion delta_time/24/3600/365.25 + 2018), so the
#     solved span ends essentially at the data.
#
# DONE 2026-09-17: default_args/rel_006_0332.txt is rel_006_0331.txt with
# EXACTLY four lines changed (--cycles, --version, --ATL11_release, --t_crop;
# `diff` shows those and nothing else -- -t, --ATL11xo_version, the E_* terms,
# the masks and --previous_product_top are untouched), and
# default_args/latest_release.txt now points at it.  rel_006_0331.txt is kept
# as the record of the first IS run.
# ALSO REPOINTED, because they would have put the symlink back: the `ln -sf`
# line in howto_MAAP_arctic.sh step 1, howto_MAAP_GL.sh step 1 and
# howto_MAAP_AA.sh step 1.
# NOT DONE HERE: nothing is recomposed or republished -- that is T3, and until
# it runs the composed input_args_IS.txt on the bucket is still the 0331 one.
#
#
# ===========================================================================
# T3. [ADE] [DONE 2026-09-17 -- RE-RUN IT AFTER T4, see below]
#     Recompose and republish the args.
# ===========================================================================
# As arctic howto steps 3 and 4.  The composed file is what DPS reads, and
# nothing reconciles it with default_args:
setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
    default_args/$reg.txt default_args/quarterly.txt --Hemisphere=1
aws s3 cp $region_dir/input_args_$reg.txt $s3_run/
# CHECK, not just run: the new input_args_IS.txt must show --cycles=0332,
# --ATL11_release=007_cycle_03_32_v05, the new --version and --t_crop, and an
# UNCHANGED --ATL11_index (the root; the generation subdirectory comes from
# --ATL11_release) and --ATL11xo_version.
#
# DONE 2026-09-17 for IS.  The recomposed file is 969 bytes (was 970) and
# differs from the 0331 composition in EXACTLY the four lines T2 changed:
# --cycles=0332, --version=02, --ATL11_release=007_cycle_03_32_v05,
# --t_crop=2019,2026.5.  Everything else is byte-identical, including
# --ATL11_earthaccess, --ATL11_index (the root), --ATL11xo_version
# (007_cycle_03_30_v03), --previous_product_earthaccess/=005_0329, the Iceland
# .db mask, -W, -g and -t.
# CONSEQUENCES CHECKED, not assumed: the tile shape the args imply is still
# [61, 61, 32] (scripts/check_field_sizes.py derives it from -W, -g, -t), and
# --t_crop=2019,2026.5 gives 31 netCDF epochs against the 0331 run's 30.
# PUBLISHED to
#   s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/north/IS/input_args_IS.txt
# and the bucket copy was read back and diffed against the local file:
# identical.
#
# THE BUCKET NOW SAYS 0332 WHILE THE TILES ON IT ARE STILL 0331.  Nothing
# reads the args by itself, so this is inert until a job is submitted -- but
# it is one more reason T4 has to happen before any submission.
#
# RE-RUN THIS STEP AFTER T4.  T4 deletes the region directory, and
# input_args_IS.txt lives in it, so the local copy goes with it; the bucket
# copy survives (different prefix).  Re-running setup_ATL1415_region.py
# recreates both the directory and the file, and republishing is then a
# no-op-but-harmless upload of the same bytes.
#
#
# ===========================================================================
# T4. [ADE+bucket] [DONE 2026-09-17 -- deleted, and T3 re-run]
#     Clear the 0331 outputs first.
# ===========================================================================
# STATEMENT, and it is the trap in this transition: the region directory
# carries the RELEASE, not the cycle range --
# setup_ATL1415_region.py:96-98 builds rel<Release>/<hemi>/<region>, so a 0332
# run writes to the same rel006/north/IS as the 0331 run did, and a tile is
# named for its center alone (E1300_N-2500.h5).  So the 0332 prelim and
# matched tiles COLLIDE with the ones on disk and on the bucket
# (s3://.../ATL14_processing/rel006/north/IS/{prelim,matched}/), and so do the
# mosaics (dz.h5, z0.h5, ...).  Only the netCDFs differ by name (0331 vs 0332).
# DECIDED (Ben, AT3): DELETE the 0331 outputs -- do not keep them aside.  They
# are not releasable (invalid lineage) and 0332 supersedes them.
#   local:  rm -r ~/ATL14_processing/rel006/north/IS
#   bucket: aws s3 rm --recursive s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north/IS/
#           (delete with the CLI, not through the mountpoint-s3 mount)
# then RE-RUN T3, which recreates the region directory and the composed args
# file that the local delete takes with it.  The published args under
# .../ATL1415/run_args/... are a DIFFERENT prefix and are not touched by the
# delete above -- do not widen the delete to reach them.
# DO IT BEFORE ANY 0332 JOB IS SUBMITTED.  A half-overwritten tree -- some
# tiles 0331, some 0332, all the same names -- is the outcome nobody could
# untangle afterwards, and NOTHING WOULD CATCH IT: -t is unchanged (AT1), so
# the old and new tiles have identical field sizes and the field-size checker
# would pass a mixture.  The only distinguishing mark inside a tile is
# meta/input_files (0331 vs 0332 granule names) and, for 0332 tiles, the new
# meta/lineage group.
# IRREVERSIBLE: ~19 worker-hours of prelim tiles and the products made from
# them.  The deletion is a deliberate instruction (AT3), so run it once, with
# the paths in front of you, and not from inside a loop over regions.
#
# DONE 2026-09-17.  WHAT WAS THERE, listed before deleting anything:
#   local  ~/ATL14_processing/rel006/north/IS, 1.4 GB -- 28 prelim and 28
#          matched tiles with their field_sizes reports, the 41 mosaic files,
#          the five 0331 netCDFs, and input_args_IS.txt;
#   bucket s3://.../ATL14_processing/rel006/north/IS/, 112 objects,
#          1196364704 bytes (the 56 tiles and their reports).
# Both are now EMPTY: the bucket prefix lists 0 objects, and rel006/north/ has
# no IS directory until setup recreates it.
# KEPT ON PURPOSE, and none of it is in the delete paths:
#   ~/ATL14_processing/maap_ledgers/IS_*.csv   the only record of the 0331 jobs
#   ~/ATL14_processing/runs/IS_mosaic, IS_nc   the 0331 run directories and logs
#   region_files/IS_*_xy.txt                   in the repo, and still correct
#   s3://.../ATL1415/run_args/rel006/north/IS/ the published args, a different
#                                              prefix -- verified still there
# T3 WAS THEN RE-RUN, as this step requires: the region directory exists again
# and holds input_args_IS.txt (969 bytes), byte-identical to the published
# copy on the bucket.  The region directory now holds THAT FILE AND NOTHING
# ELSE, which is the clean starting point T6 and T7 need.
#
#
# ===========================================================================
# T5. [DPS] [DONE 2026-09-17 -- VERDICT: MATCH at 61a19af]
#     Rebuild, register, MATCH.
# ===========================================================================
# This carries the lineage code (plan_lineage_at_solve_time.sh L1-L5, d2960d5)
# as well.  Push, Ben registers, then:
scripts/maap/check_build_id.py \
    s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/north/IS/input_args_IS.txt \
    maap-dps-worker-16gb --expect 61a19af        # must say MATCH before any job
# RESULT 2026-09-17, Ben registered and the check was run with the IS args:
#   image stamp, live git in the image and the CWL are all
#   61a19af58e2fcdd4528edd891eea70c87e632d81, tree_state=clean,
#   algorithm_version=on_s3, maap_py=5.1.0, maap_pgt=set;
#   build_started 19:40:55Z, build_completed 19:43:57Z -- AFTER the 19:39:27Z
#   commit, so it is a real rebuild and not a stale image;
#   processID=64, job fa3e0679-e1e1-4eb8-97b0-b9b581214af0, successful in
#   ~4.5 min ([[cwl-link-is-registration-not-deploy]]: only this proves it).
# SO THE DEPLOYED IMAGE CARRIES BOTH the lineage work (L1-L5) and nothing that
# contradicts the 0332 args -- the args are an input, not part of the image.
#
#
# ===========================================================================
# T6. [DPS] [DONE 2026-09-18 -- BOTH GATES PASS]  One smoke tile, with two gates.
# ===========================================================================
# RECOMMENDATION: E1300_N-2500 -- a full-3x3 center with crossovers, the tile
# the matched smoke used.  Submit prelim for that center alone.
# GATE A, the transition: N_ATL11 > 0 and N_XO > 0 in the collector's row (a
# staging or filter mistake shows up as zero along-track data), and the tile's
# meta/input_files names 0332_007_05 granules.
# GATE B, the lineage work (that plan's L6): every name in input_files has a
# meta/lineage/<granule> group, with uuid, geoseg and orbit for along-track
# and uuid, geoseg and both rgts for crossovers; no read-failure lines in the
# job log.
# ALSO WORTH READING: wall time and peak RSS against the 0331 run (1169-4260 s,
# 4.16-9.09 GiB).  One more cycle of data and one more epoch should move both
# a little; a large jump means something else changed.
# RESULT, job 64c1c7df-b23f-46cb-ab71-d696ae106950 (submitted 2026-09-17T19:57Z,
# ledger IS_0332_smoke_jobs.csv), collected once 2026-09-18 and the tile
# fetched to rel006/north/IS/prelim/ (63689285 bytes + report):
#   successful on 61a19af, 3181 s, peak 9.32 GiB of 16, 3 iterations.
# GATE A PASSES.  N_ATL11 306563, N_AT 306475, N_XO 88, N_fit 248749.  The
#   /meta attribute input_files (an ATTRIBUTE, not a group -- "meta/input_files"
#   above means that) names 21 distinct granules: 13 along-track, ALL
#   0332_007_05, and 8 crossovers, all c01/c02 _007_03 (by design, AT5).
#   check_field_sizes.py: 1 of 1 passed, dz/dz [61,61,32].
# GATE B PASSES.  /meta/lineage holds exactly those 21 groups, none missing
#   and none extra; every uuid non-empty; along-track groups carry uuid,
#   start/end_geoseg, start/end_orbit (and start/end_rgt as well -- the
#   granule provides them, which L1 allows), crossover groups uuid, geoseg
#   and start/end_rgt; all int32 as read.  0 of 8 crossovers have
#   start_rgt == end_rgt (the 66e4a35 fix holds).  The job log has no
#   'could not read lineage' and no 'INVALID' line.
# AGAINST THE 0331 JOB FOR THE SAME TILE (c72c7cb1, b5fe447, re-collected):
#   N_ATL11 300004 -> 306563 (+2.2%), N_XO 88 -> 88, N_fit 245361 -> 248749
#   (+1.4%), sigma_hat 4.23 -> 4.22: what one more cycle should do.
#   Peak 8.80 -> 9.32 GiB: now the prelim high-water mark (was 9.09), still
#   well inside 16 GiB.  Error step 1835 -> 1952 s, in line.
#   THE FIT STEP HALVED, 2404 -> 1213 s, and it is the QR solves: ~712-735 s
#   per iteration then, ~325-354 s now, same 4 threads, same problem size.
#   NOT EXPLAINED, and harmless to the result (the numbers above agree).  The
#   error step, also compute-bound, did NOT speed up, so a faster node alone
#   does not fit; one HYPOTHESIS is contention -- the 0331 job ran inside the
#   29-job fan-out, this one ran alone.  T7's fan-out will show whether fit
#   times go back up.
# ALSO DUE NOW, AT4 ("Delete once T6 has solved"): the two deletions in T1.
#   DONE 2026-09-18 by Ben (the session's auto-mode permission check refused
#   the bucket delete) and verified -- see T1.
#
#
# ===========================================================================
# T7. [DPS+ADE] [DONE 2026-09-18 -- prelim 29/29, 28 tiles; matched 28/28 after one retry]  The IS re-run.
# ===========================================================================
# Only after T6 passes both gates.  Arctic howto steps 6-9, unchanged:
# prelim over the 29 centers, collect, fetch, check_field_sizes.py, then the
# matched list from the tiles that exist, matched, fetch, check again.
# EXPECT the counts to differ from the 0331 run: E1020_N-2580 wrote no tile
# then, and with a cycle more data it may.  Build the matched list from what
# EXISTS, as I6 says, rather than assuming 28.
# PRELIM SUBMITTED 2026-09-18T01:15:50Z (Ben: "Go ahead with the 29 centers"):
# all 29 of region_files/IS_prelim_xy.txt, E1300_N-2500 included (the smoke
# tile re-runs, as howto step 6 says), on maap-dps-worker-16gb, tag
# IS_rel006_0332_prelim, against the published args (diffed identical to
# the local copy first) and --tile_prefix .../rel006/north/IS.  29/29
# accepted.  LEDGER: ~/ATL14_processing/maap_ledgers/IS_0332_prelim_jobs.csv
# -- a NEW name; IS_prelim_jobs.csv is the 0331 run's ledger and is kept.
# DO NOT RE-REGISTER until every one of the 29 has finished (I3's split-build
# trap); collect_jobs.py's commit column must read 61a19af on all 29.
# MATCHED IS PRE-AUTHORIZED (Ben 2026-09-18): once all 29 prelim jobs have
# finished, fetch, check_field_sizes.py, build region_files/IS_0332_matched_xy.txt
# from the tiles that EXIST (bucket and local listed, must agree), commit and
# push it, then submit matched (ledger IS_0332_matched_jobs.csv) WITHOUT asking
# -- UNLESS any prelim job ended failed: then hold matched and report ("Holding
# matched on a failure is fine").  A successful job that wrote no tile is not a
# failure.  T8 is not authorized by this.
# PRELIM DONE 2026-09-18 (a background collect every 15 min, last at 02:30Z):
#   29/29 successful, ALL on 61a19af (no split build).  E1020_N-2580 again
#   took the no-data path (error step 4 s at 0.23 GiB, N_fit 349) and wrote
#   no tile, so 28 tiles -- the same 28 centers as 0331.
#   Wall 1332-3322 s (plus E1020's 490 s); 0331 was 1169-4260 s.
#   Peak 4.42-9.52 GiB; NEW HIGH-WATER MARK 9.52 GiB (E1340_N-2460, N_fit
#   273382), still well inside 16.  Memory tracks N_fit, as before.
#   THE FIT-TIME QUESTION FROM T6: E1300_N-2500's fit took 1819 s inside the
#   fan-out against 1213 s alone and 2404 s in the 0331 fan-out -- partway
#   back up, which fits contention without proving it.  Not pursued.
#   Fetched 27 + the re-run smoke tile (re-copied so the local file is the
#   fan-out's); check_field_sizes.py 28 of 28 OK.  Lineage scanned on all
#   28: every input_files name has a complete group with a uuid, every
#   along-track granule is 0332_007_05; 69 along-track + 10 crossover
#   granules in all (the ~79 T8 expects).
# MATCHED LIST: region_files/IS_0332_matched_xy.txt (2601047), 28 centers,
#   from the tiles that exist; bucket and local listings agree on names and
#   sizes; identical to the 0331 IS_matched_xy.txt.
# MATCHED SUBMITTED 2026-09-18T02:35:59Z: 28/28 accepted, maap-dps-worker-16gb,
#   tag IS_rel006_0332_matched, matched/ prefix verified empty first.
#   LEDGER ~/ATL14_processing/maap_ledgers/IS_0332_matched_jobs.csv.
#   Still DO NOT RE-REGISTER until these finish.
# MATCHED DONE 2026-09-18 (Ben: "Collect, fetch, check", one collect at 04:07Z):
#   27/28 successful on 61a19af; wall 145-830 s; peak 3.58-9.36 GiB
#   (E1340_N-2460, N_fit 273382 -- the same tile that set the prelim high).
#   ONE FAILURE, E1220_N-2460 (e48a0dd3), 136 s, and it was the PLATFORM:
#   the CWL runner's own /app/create_inputs.py -> stage_in.py timed out
#   connecting to api.maap-project.org, before our image was launched.
#   Logs in ~/tmp/triaged_job-...20260918T023520.389304Z_task-c0c6f9a3...;
#   nothing for that tile on the bucket.  The 0331 run solved it fine.
#   RETRIED (Ben: "Resubmit it"): region_files/IS_0332_matched_retry_xy.txt,
#   ledger IS_0332_matched_retry_jobs.csv, job c3e8e289 on 61a19af --
#   successful, 201 s, 4.90 GiB, N_fit 42127.
#   Fetched all 28 (0.48 GiB); local sizes == bucket listing, tile for tile;
#   same 28 names as prelim/; check_field_sizes.py 28 of 28 OK (dz/dz
#   [61,61,32], sigma_dz null).
#
#
# ===========================================================================
# T8. [ADE] [DONE 2026-09-18 -- all checks pass; rel005 comparison below]  Mosaic, netCDF, and the checks.
# ===========================================================================
# plan_IS_run.sh I9 as written, plus:
#   - the lineage must now be COMPLETE: no INVALID warning from either writer
#     (plan_lineage_at_solve_time.sh L7), a real uuid on all ~79 rows, and
#     crossover end_rgt different from start_rgt where the granule says so;
#   - the netCDFs are named ATL1[45]_IS_0332_..._006_02.nc (AT2), and the
#     0331 products are gone by then (T4, AT3);
#   - EXPECT 31 EPOCHS in ATL15, not the 30 the 0331 products had (AT1's
#     wider --t_crop), with delta_h (31, 301, 421) at 1 km and the dhdt groups
#     one longer each.  The mosaic and the tiles still have 32.
#   - I9g6 is worth more here than it was: differencing h against the rel005
#     product is the only check that the extra cycle did not shift anything.
# AUTHORIZED 2026-09-18 (Ben: "Yes, run T8").  Run dirs are NEW --
# ~/ATL14_processing/runs/IS_0332_mosaic and IS_0332_nc; runs/IS_mosaic and
# IS_nc are the 0331 run's and were left alone.
# MOSAIC: make_mosaic_jobs.py as I9c with --run_name IS_0332_mosaic -> 41
#   tasks; all 41 at -P 12 in 59 s, exit 0, error_logs/ empty, done/ 41.
#   check_mosaic_outputs.py --values: 41 files, 131 fields, PROBLEMS 0.
#   STATEMENT, measured: at 1 km every sigma_dzdt_lag* is finite on 16.9% of
#   the box against 41.4% for the values; sigma_dz is 40.2%, and the 10/20/40 km
#   sigmas match their values exactly.  The gap is IN THE PRELIM TILES, not the
#   mosaic: sigma_dzdt covers 44% of the finite dzdt cells over all 28 tiles,
#   ranging 8.5% (E1180_N-2380, E1180_N-2420) to 99.7% (E1340_N-2500) and
#   tracking data density, while sigma_dz is 96.9% on every tile.  Whether
#   0331 was the same cannot be checked -- its tiles are deleted (T4).
#   NOT PURSUED; flagged to Ben.
# netCDF: ATL14_write2nc.py 12 s, ATL15_write2nc.py 21 s, both exit 0, NO
#   INVALID warning (logs hold only GDAL's FutureWarning).  Five files:
#   ATL14_IS_0332_100m_006_02.nc (9915264) and
#   ATL15_IS_0332_3mo_{1,10,20,40}km_006_02.nc.
# CHECKS, all pass (scratch script, not committed):
#   - lineage identical in all 5 files: 79 rows = 69 ATL11 + 10 ATL11XO, the
#     same 79 names as the prelim tiles' input_files; every attribute a
#     string; 79 distinct uuids, none NOT_SET; no NOT_SET on any along-track
#     row; along-track all 0332_007_05, cycles 03-32, end_rgt == start_rgt;
#     XO end_rgt != start_rgt on 10 of 10.  XO rows are NOT_SET for
#     start/end_orbit and start/end_region ONLY -- the XO granules have no
#     such datasets (plan_lineage_at_solve_time.sh, the probe), so that is
#     the granule's absence, not a gap.  L7's gate passes.
#   - ATL14 h (3001, 4201), x/y == z0.h5; finite(h) == finite(z0) & ice_area>0
#     exactly; max |h - z0| 1.2e-4 m.
#   - ATL15 delta_h (31, 301|30|14|7, 421|42|20|10), 31 epochs 2019.00..2026.50
#     (AT1 as expected); dhdt groups 30,29,27,23,19,15,11,7,3 epochs, one
#     longer each than 0331.  1 km: mosaic epochs kept 31 of 32; finite
#     product == finite dz & ice_area>0 exactly; max |diff| 3.1e-5 m.
# I9g6 DONE -- AGAINST rel005 (ATL14_IS_0329_100m_005_02.nc and
#   ATL15_IS_0329_01km_005_02.nc, found by find_previous_product_files'
#   CMR search and read by fs.open -- fs.get does a ListBucket, which NSIDC
#   denies).  Ben's bar (2026-09-18): no >10 m errors, no major gaps.
#   GAPS -- NONE MAJOR.  ATL14: 1304 rel005 cells (0.11%) have no 0332 value,
#     82 the reverse; 24 of 12390 1-km blocks lose over half their cells.
#     ATL15 1 km: 292 cells (0.09%) over 29 common epochs.
#   ATL15 1 km delta_h -- CLOSE.  Median -0.025 m, p5/p95 -0.86/+0.68 m,
#     |d|>10 m on 166 cells of 327796 (0.05%), max 19.0 m.
#   ATL14 h -- LARGE DIFFERENCES, ALMOST ALL WHERE THERE ARE NO DATA.
#     Median +0.01 m, p5/p95 -4.42/+4.38 m, but |d|>10 m on 43129 of 1133949
#     cells (3.8%), max 330 m.  42702 of those have data_count 0 (4.2% of the
#     no-data cells); where data_count > 0 it is 427 of 118777 (0.36%).  Median
#     h_sigma is 10.9 m on the >10 m cells against 3.1 m elsewhere; 8673 of
#     them still exceed 3x the two sigmas combined.
#     RECOMMENDATION: read as interpolation differences, not errors -- but it
#     is Ben's bar, so it is reported to him rather than decided here.
# FOUND, NOT FIXED: time_coverage_duration is wrong in all five files
#   (54385 in ATL14, against ~2.37e8 s for 2019-01-01..2026-07-02).
#   ATL1415_attrs_meta.py:314 computes int((datetime_start-datetime_end).seconds)
#   -- operands reversed, and .seconds is the within-day part of a negative
#   timedelta, not total_seconds().  From 198ffee; pre_rel006 has the same
#   line; main computes it another way.  rel005 carries 2.17e8.
#
#
# ===========================================================================
# T9. [DONE 2026-09-19 for the docs; the other regions are their own howtos]  Docs, and the other regions.
# ===========================================================================
#   - plan_IS_run.sh: the run it describes becomes the 0331 run; record that
#     the products were re-made at 0332 and what the counts were.
#   - plan_lineage_at_solve_time.sh: L6/L7 stop being blocked once T6 runs.
#   - howto_MAAP_arctic.sh / _GL / _AA: nothing structural -- they read
#     default_args/latest_release.txt, which T2 repoints.
#     DONE 2026-09-18 for howto_MAAP_arctic.sh, via plan_monthly_on_maap.sh
#     M11: status banner now 0332 with complete lineage; step 10's stale
#     "LINEAGE IS INVALID" note corrected; monthly added as step 10b.
#     _GL and _AA NOT touched.
#     The other bullets here are still open.
#   DONE 2026-09-19: plan_IS_run.sh's banner records the 0332 re-run; its I9g2
#     points at the finished lineage work; howto_MAAP_GL.sh and _AA.sh were
#     rewritten from the IS run, and howto_MAAP_arctic.sh again in discover's
#     order.  The last bullet -- check GL's and AA's queues and masks against
#     the new generation before their first fan-out -- is now those howtos'
#     smoke steps (GL 3, AA 4); every mask they name is staged.
#   - GL and AA have never run, so they simply start at 0332; nothing to
#     re-do, but their queues and masks should be checked against the new
#     generation before their first fan-out.
#
#
# ===========================================================================
# QUESTIONS FOR BEN -- ALL ANSWERED 2026-09-17, in his words (AT*)
# ===========================================================================
# QT1. -t and --t_crop for 0332 (T2).  I recommend -t=2018.75,2026.75 and
#      --t_crop=2019,2026.5, from the granule's end_delta_time (2026.489) and
#      the margin the 0331 pair used.  Confirm or give the values: -t decides
#      the epoch count in every tile, so it cannot be adjusted afterwards
#      without re-solving.
# AT1: "Use -t=2018.75,2026.5 and --t_crop=2019,2026.5"
#      So -t is UNCHANGED from 0331 and the recommendation is not taken: the
#      tiles keep 32 epochs, and the crop now reaches the top of the solved
#      span instead of stopping one dt below it.  See T2 for what that means.
# QT2. --version stays 01, or does re-issuing rel006 with a longer cycle range
#      bump it (ATL14_IS_0332_100m_006_02.nc)?  --Release stays 006 either way.
# AT2: "Use release 02"
#      READ AS --version=02 (the '02' in ATL14_IS_0332_100m_006_02.nc), since
#      --Release stays 006 and 02 is the only field the question offered.
#      Say if you meant something else -- it renames every product.
# QT3. The 0331 outputs (T4): move them aside to rel006/north/IS_0331, locally
#      and on the bucket, or delete them?  They are ~19 worker-hours of prelim
#      tiles, 1.1 GB local plus the bucket copies, and the products they made
#      have invalid lineage, so they are not releasable.
# AT3: "Delete the old files"
# QT4. The 0331 index on the bucket (2.1 GB) and the ADE staging directory
#      (2.5 GB): keep, or delete once T6 has solved a tile on 0332?
# AT4: "Delete once T6 has solved"
# QT5. Does anything else need to move with the cycle range -- the crossover
#      generation (I have left it at 007_cycle_03_30_v03, cycles 1-2), or the
#      previous product (005_0329)?
# AT5: "Leave these as they are."
