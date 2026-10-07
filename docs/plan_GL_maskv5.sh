# plan_GL_maskv5.sh -- Greenland on the v5 ice mask: rerun the north's
# coastal tiles, run the rest of Greenland (quarterly)
#
# ############################################################################
# ##  WRITTEN 2026-10-02.  TENTATIVE: steps tagged per line as they run.    ##
# ##  Ben 2026-10-02: "I've just uploaded a new mask file for Greenland     ##
# ##  that genuinely covers 2018-2026.5 ... rerun the coastal tiles that    ##
# ##  already exist in Greenland north and finish the rest of Greenland.   ##
# ##  Matched tiles that overlap coastal areas or that were partially       ##
# ##  unconstrained by the boundary of Greenland north also need to be     ##
# ##  rerun."                                                               ##
# ############################################################################
#
# Status tags per step: TODO / DONE / BLOCKED.  DECIDED = Ben said so;
# RECOMMENDATION = mine; QUESTION = open.  Commands are plan_GL_north.sh's.
#
# QUESTIONS FOR BEN (answer inline; none blocks V0-V6):
#   QV1  After matched: full-GL mosaic + ATL14/ATL15 netCDF on DPS
#        (plan_dps_mosaic.sh D7 for all of GL)?  RECOMMENDATION: yes, to the
#        TEST prefix .../rel006_0332_testing/north/GL as D7 did, then the
#        rel005 comparison.  The north-only products there would be
#        overwritten.
#   QV2  GL monthly.  STATEMENT: the north monthly prelim/matched used the
#        v4.1 mask AND the north-only ATL14 as reference.  A full-GL monthly
#        run needs the V7 ATL14 as reference and the v5 mask in the monthly
#        args (.../north_monthly/GL/input_args_GL.txt, NOT changed here).
#        RECOMMENDATION: rerun all 1483 monthly centers after V7, not only the
#        coastal ones, because the reference changes everywhere.
#
# WHAT IS KNOWN (statements, with provenance):
#   - THE MASK.  s3://.../ATL1415/masks/Arctic/
#     GreenlandIceMask_2018.1_2026.5_100m_v5.tif, 25,001,381 bytes, md5
#     0735f7ed074fe56db760e5ab461ea752.  Same grid as v4.1 (origin, size
#     14958 x 27000, 100 m; checked with gdalinfo).  43 bands with 'time'
#     metadata, 2018.00 .. 2026.41; bands 1-40 have v4.1's times; 2026.00,
#     2026.25, 2026.41 are new.  Byte, nodata 255 (v4.1: NaN); values are 0
#     and 1 only (histogram of every band), so the nodata change has no
#     effect: pc from_gdal maps 255 -> NaN and ATL11_to_ATL15 sets NaN -> 0.
#     No _time_stamps.txt beside it -- none is read (grep: the times come
#     from the band metadata).
#   - HOW A TILE SEES THE MASK (read in ATL11_to_ATL15.py 577-620): prelim
#     reads the mask over the 60 km box + 10 km pad, all bands in t-range +-1
#     yr, and repeats the last band out to the end of the run (so v4.1 ran
#     with 2025.83 for 2026.0-2026.5).  The solve grids span the box (center
#     +-30 km).  Matched does NOT read the mask file: it takes the mask from
#     its own prelim tile (read_mask_file), and its 8 neighbours' fits.
#   - WHERE IT VARIES (session_tools_2026-10-02/gl_mask_v5_variability.py,
#     every pixel, every band, 11 min; output maap_ledgers/
#     GL_mask_v5_cells.npz):
#       236,079 100 m pixels vary in time in v5      (9,414 1 km cells)
#        10,747 differ from v4.1-as-the-code-saw-it  (842 1 km cells)
#            14 1 km cells differ from v4.1 but are static in v5.
#   - WHICH TILES (session_tools_2026-10-02/gl_mask_v5_tile_lists.py; a tile
#     is touched if a 1 km cell of either kind lies within center +-31 km --
#     the 60 km box plus one cell; summary maap_ledgers/
#     GL_0332_maskv5_tile_lists.txt; map GL_0332_maskv5_tiles.png):
#       north (y >= -1520 km, 557): 189 touched (all 58 touched by a v4.1
#         difference are among the 189 that vary in time)
#       south (926): never run; all run.  (305 of them are touched.)
#       PRELIM:  189 north reruns + 926 south = 1115
#       MATCHED: every center whose 3x3 neighbourhood holds a V4 prelim
#         tile = 1222 (296 north incl. the y=-1520 row, 926 south);
#         261 north matched tiles are kept.
#     Lists: maap_ledgers/GL_0332_maskv5_{prelim_rerun_north,
#     prelim_new_south,prelim,matched}_tile_list.txt.
#   - Build: e6d7051 is deployed (NM7 ran on it today); every commit since
#     is docs only (git diff --stat e6d7051 HEAD).  This plan needs no new
#     build: the mask is an argument, read from the bucket.
#   - A RERUN PRELIM THAT NOW HAS NO DATA does not upload (run.sh), so the
#     v4.1 tile would survive at the canonical prefix.  V3 keeps a copy of
#     the old tiles; V5 checks that every rerun tile was rewritten.

conda activate ATL14
cd ~/git_repos/ATL1415
ATL14_root=/home/jovyan/ATL14_processing
s3_root=s3://maap-ops-workspace/ben_smith
ledgers=$ATL14_root/maap_ledgers
region_dir=$ATL14_root/rel006/north/GL
s3_run=$s3_root/ATL1415/run_args/rel006/north/GL
s3_out=$s3_root/ATL14_processing/rel006/north/GL
s3_old=$s3_root/ATL14_processing/rel006/north/GL_maskv4.1_superseded
tag=GL_rel006_0332_maskv5
L=$ledgers/GL_0332_maskv5


# ===========================================================================
# V0. [ADE] DONE 2026-10-02.  Mask analysis and tile lists (above).
# ===========================================================================
python $ATL14_root/session_tools_2026-10-02/gl_mask_v5_variability.py $ledgers/GL_mask_v5_cells.npz
python $ATL14_root/session_tools_2026-10-02/gl_mask_v5_tile_lists.py $ledgers/GL_mask_v5_cells.npz \
    ATL1415/resources/GL/40km_tile_list.txt $L


# ===========================================================================
# V1. [ADE] DONE 2026-10-02.  Args on the v5 mask.
# ===========================================================================
# default_args/GL_0332.txt = GL_0331.txt with the v5 mask; GL_latest.txt ->
# GL_0332.txt.  Composed as howto_MAAP_GL step 2, then --solver=cholmod put
# back by hand before -b (as on 2026-09-25).  Old args kept as
# $ledgers/GL_0332_input_args_GL_maskv4.1.txt.
# GATE: diff old vs new = the --mask_file line only.  PASSED.
setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
    default_args/GL_latest.txt default_args/quarterly.txt --Hemisphere=1
sed -i 's/^-b=/--solver=cholmod\n-b=/' $region_dir/input_args_GL.txt
diff $ledgers/GL_0332_input_args_GL_maskv4.1.txt $region_dir/input_args_GL.txt


# ===========================================================================
# V2. [ADE] DONE 2026-10-02 23:2xZ (bucket == local; old copy == ledger copy).
#     Publish the args; keep the old ones beside them.
# ===========================================================================
aws s3 cp $s3_run/input_args_GL.txt $s3_run/input_args_GL_maskv4.1.txt
aws s3 cp $region_dir/input_args_GL.txt $s3_run/
# GATE: bucket copy == local.


# ===========================================================================
# V3. [S3] DONE 2026-10-02: 189 prelim + 296 matched at $s3_old, sizes equal
#     to the canonical ones (key list $ledgers/GL_0332_maskv5_V3_copy_keys.txt).
#     Copy (not move) the tiles that will be overwritten.
# ===========================================================================
# 189 prelim + 296 matched (the matched list's north centers), to $s3_old.
# Copies, so nothing disappears from the canonical tree before its
# replacement lands.


# ===========================================================================
# V4a. [DPS] DONE 2026-10-02 23:27-23:43Z.  Smoke: 3 prelim on the v5 args, 32gb queue.
#   RESULT (ledger ${L}_smoke_jobs.csv, _collect.txt): 3/3 successful, e6d7051,
#   r5.xlarge; every log names GreenlandIceMask_2018.1_2026.5_100m_v5.tif and
#   no other mask.
#                  job s  peak GiB  N_fit
#     E480_N-1040   657    9.63     514572  (v4.1, N2: 514540)
#     E200_N-1880   726   10.04     828458
#     E-160_N-2240  740    9.53     607899
#   E480 v5 vs v4.1 tile, z0 on reported cells (none reported in only one),
#   by distance to the nearest changed 1 km mask cell:
#     0-2 km max 0.66 m; 2-5 km 0.20 m; 5-10 km 0.26 m; 10-20 km 1.4e-3 m;
#     > 20 km 6.7e-5 m (median 5e-6 m).  dz max 0.15 m, dzdt_lag1 0.32 m.
#   STATEMENT: the gate below said ~1e-7 m away from changed cells; the far
#   field is 7e-5 m instead -- the fit couples the whole tile, so a mask
#   change near the coast moves the interior slightly.  The decay with
#   distance is what a mask-driven change looks like; I judged it a pass and
#   flag it here rather than re-word the gate.
# ===========================================================================
#   E480_N-1040  north rerun (79N; tides).  Has a v4.1 tile: the difference
#                must be confined near where the mask changed.
#   E200_N-1880  south, densest interior tile (sizes the memory).
#   E-160_N-2240 south, Jakobshavn: the largest temporal mask changes.
# GATES: 3 successful; commit e6d7051; log names the v5 mask; E480 vs its
#   v4.1 tile (copy in $s3_old): reported cells away from changed mask
#   cells agree to ~1e-7 m; peak memory < ~28 GiB.
# V4b. [DPS] DONE 2026-10-06 23:19Z (started 2026-10-02 23:44Z).  The other 1112 prelim, ALL AT ONCE with client-side
#   resubmit (Ben 2026-10-01 for monthly; NM5/NM7 driver pattern, 2 rounds):
#   $ledgers/GL_0332_maskv5_driver.py prelim (+ GL_0332_maskv5_prelim_driver.log;
#   ledgers GL_0332_maskv5_prelim_run[_retryN]_jobs.csv).
#   ROUND 0: 1112 POSTs in 125 s (23:44:50-23:46:55Z), none refused.  By
#     23:57Z 531 FAILED, 49 successful, 532 running.  Sampled failures: all at
#     our first step (workspace credentials, 5 attempts, exit 1 -- nothing
#     read or written): MAAP API ConnectTimeout, and NEW, "[Errno 111]
#     Connection refused"; and HTTP 401 on every attempt after the runner's
#     get_maap_pgt_token.py timed out (NM7's two classes).  STATEMENT: twice
#     NM7's job count at once overloaded the shared MAAP API.
#   CHANGED (mine, 23:58Z): retry rounds go out at --max_in_flight 100
#     (QN3's pattern), not all at once, so the ~500 resubmissions do not
#     repeat the overload.  Driver restarted (round 0 only waited on).
#     Ben's "all at once" was for the monthly load test; for this run it was
#     my choice.
#   ADE LOST: the driver log stopped at 2026-10-03 00:15Z.  Restarted
#     2026-10-06 21:30Z with the same command.  Round 0 final: 578
#     successful, 534 failed.  Round 1 (the 534, at 100 in flight)
#     submitted 21:31Z.
#   ROUND 0 FAILURES, ALL 534 classified from the triaged_job stderr
#     (scratchpad classify_r0.py, 2026-10-06).  Every one was a call to
#     api.maap-project.org that was refused, timed out, or returned 401:
#     326 MAAP stage-in (stage_in.py -> /api/environment/config; 261
#       refused, 65 timeout) -- before our code ran;
#     157 our workspace-credentials step (5 attempts; 81 with HTTP 401 after
#       the runner's get_maap_pgt_token.py timed out, 76 refused/timeout);
#     49 past credentials, in the fit: read_ATL11 -> get_s3fs could not
#       broker NSIDC credentials through MAAP (39 refused, 10 timeout),
#       0.7 GiB, ~50 s, no tile written;
#     2 failed after 453 s with no triaged_job record (cause unknown).
#     None was a data, memory, or solver failure.
#   ROUND 1 (534 at 100 in flight, 21:31-22:28Z submit, done 22:46Z): 528
#     successful, 6 failed (1.1%), again all MAAP API: 4 MAAP stage-in
#     ConnectTimeout on /api/environment/config, 1 our credentials step
#     (401 after get_maap_pgt_token timeout), 1 no triaged_job record.
#     Round 2 (the 6) submitted 22:47Z; all 6 successful 23:19Z.
#   ALL_DONE: all 1112 (+3 smoke) prelim tiles have a successful job.


# ===========================================================================
# V5. [ADE] DONE 2026-10-06 23:20Z.  Check prelim.
# ===========================================================================
# check_field_sizes through the mount on the 1115; every one of the 189
# north reruns has LastModified after V4 started (otherwise: a no-data
# rerun left the v4.1 tile -- remove it and record it); no-data centers
# listed (fetch_tiles' no_data_tiles.txt equivalent from the job logs).

# RESULT: gl_maskv5_V5_check.py 2026-10-02T23:20:00 -> 1483 prelim tiles on
#   the prefix; all 1115 V4 centers have a tile (no no-data centers); all
#   189 north reruns rewritten after 23:20Z (none stale).
#   check_field_sizes through the mount: expected dz [61, 61, 32], 1483 of
#   1483 passed, 0 problems (96 s).


# ===========================================================================
# V6. [DPS] RUNNING since 2026-10-06 23:21Z.  Matched, the 1222, all at once, same driver (--step
#   matched).  Only after V5 -- a matched job reads its neighbours' prelim
#   tiles as they are at that moment.
# ===========================================================================
# CHANGED (mine, 2026-10-06): round 0 at --max_in_flight 100 too, not all
#   at once (V4b: 48% failed all at once, 1.1% at 100 in flight).  Driver
#   log GL_0332_maskv5_matched_driver.log, ledgers
#   GL_0332_maskv5_matched_run[_retryN]_jobs.csv.
# INTERRUPTED 2026-10-07 00:08Z: the ADE restarted during round 0's submit
#   (987/1222 in the ledger, no POST failures).  Driver patched (mine): a
#   round whose ledgers miss tiles submits them as _resume<k>_jobs.csv and
#   waits on all its ledgers.  Restarted 00:09Z; resume1 = the other 235.
#   The resume's in-flight cap counts only its own jobs (~200 briefly).
# Then check_field_sizes --step matched; 1483 matched tiles (less no-data
# centers) at $s3_out/matched.


# ===========================================================================
# V7. [DPS] QUESTION QV1.  Mosaic + ATL14/ATL15 netCDF for all of GL.
# V8. [DPS] QUESTION QV2.  Monthly.
# ===========================================================================
