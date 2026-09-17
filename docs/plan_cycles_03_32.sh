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
# [[maap-bucket-is-mountpoint-s3]]: this was `aws s3 cp`, not the mount.
#
#
# ===========================================================================
# T2. [ADE] [NEEDS EDIT: default_args]  The release args file.   QT1, QT2
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
#
# ===========================================================================
# T3. [ADE] [NOT STARTED]  Recompose and republish the args.
# ===========================================================================
# As arctic howto steps 3 and 4.  The composed file is what DPS reads, and
# nothing reconciles it with default_args:
setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
    default_args/$reg.txt default_args/quarterly.txt --Hemisphere=1
aws s3 cp $region_dir/input_args_$reg.txt $s3_run/
# CHECK, not just run: the new input_args_IS.txt must show --cycles=0332,
# --ATL11_release=007_cycle_03_32_v05, the new -t and --t_crop, and an
# UNCHANGED --ATL11_index (the root; the generation subdirectory comes from
# --ATL11_release) and --ATL11xo_version.
#
#
# ===========================================================================
# T4. [ADE+bucket] [DECIDED 2026-09-17, AT3: "Delete the old files"]
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
# then let T3 recreate the region directory.
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
#
# ===========================================================================
# T5. [DPS] [NOT STARTED]  Rebuild, register, MATCH.
# ===========================================================================
# This carries the lineage code (plan_lineage_at_solve_time.sh L1-L5, d2960d5)
# as well.  Push, Ben registers, then:
scripts/maap/check_build_id.py     # must say MATCH before any job
#
#
# ===========================================================================
# T6. [DPS] [NOT STARTED]  One smoke tile, with two gates.
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
#
#
# ===========================================================================
# T7. [DPS+ADE] [NOT STARTED]  The IS re-run.
# ===========================================================================
# Only after T6 passes both gates.  Arctic howto steps 6-9, unchanged:
# prelim over the 29 centers, collect, fetch, check_field_sizes.py, then the
# matched list from the tiles that exist, matched, fetch, check again.
# EXPECT the counts to differ from the 0331 run: E1020_N-2580 wrote no tile
# then, and with a cycle more data it may.  Build the matched list from what
# EXISTS, as I6 says, rather than assuming 28.
#
#
# ===========================================================================
# T8. [ADE] [NOT STARTED]  Mosaic, netCDF, and the checks.
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
#
#
# ===========================================================================
# T9. [NOT STARTED]  Docs, and the other regions.
# ===========================================================================
#   - plan_IS_run.sh: the run it describes becomes the 0331 run; record that
#     the products were re-made at 0332 and what the counts were.
#   - plan_lineage_at_solve_time.sh: L6/L7 stop being blocked once T6 runs.
#   - howto_MAAP_arctic.sh / _GL / _AA: nothing structural -- they read
#     default_args/latest_release.txt, which T2 repoints.
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
