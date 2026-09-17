#! /usr/bin/env bash
# ===========================================================================
# PLAN: record ATL11 lineage in the prelim tiles at solve time.
# Written 2026-09-17.  NO CODE YET.  Revise as steps land.
# QL1-QL5 were ANSWERED by Ben on 2026-09-17 (AL1-AL5, at the foot of this
# file) and are folded into the steps: the layout is settled, the netCDF
# forces strings, granules are opened a second time, pre_rel006 is untouched,
# and the cycles 03-32 transition follows this work under its own plan.
# ===========================================================================
# THE DECISION (Ben, 2026-09-17; docs/plan_IS_run.sh I9g2, QI9/QI10):
#   the prelim step records the lineage attributes of every granule it reads
#   in the tile's metadata; the netCDF step reads them from the tiles and
#   NEVER opens ATL11; an attribute the tiles do not carry is marked invalid.
#   The temporary half is in 28b4f72: today every file-only attribute is
#   'NOT_SET'.  This plan is the permanent half.  Applies to GL and AA too.
#   IS will be re-run completely once more issues are fixed.
#
# Provenance on every claim: STATEMENT = verified 2026-09-17, with how;
# DECIDED = Ben said so; RECOMMENDATION = mine, overridable; QUESTION = open,
# for Ben, not guessed.  Tags: [NOT STARTED] [NEEDS CODE: x] [BLOCKED: x]
# [DONE] -- and [ADE] / [DPS] for where a step runs.
#
#
# ===========================================================================
# READ THIS FIRST: A RE-RUN OF IS IS BLOCKED ON SOMETHING ELSE ENTIRELY
# ===========================================================================
# STATEMENT, found while probing for this plan (CMR, read-only, 2026-09-17):
#   CMR NO LONGER LISTS THE ATL11 GENERATION IS RAN ON.
#   - IS ran on --ATL11_release=007_cycle_03_31_v04 (69 granules, 2026-09-15).
#   - earthaccess.search_data(short_name='ATL11') over Iceland
#     (-20, 63.5, -14, 66.5) returns 55 granules, ALL 0332_007_05
#     (e.g. ATL11_023003_0332_007_05.h5).  None is 0331_007_04, and an exact
#     granule_name search for ATL11_023003_0331_007_04.h5 returns nothing.
#     The collection list shows one ATL11 collection, version 007.
#   - The only staged index is s3://maap-ops-workspace/ben_smith/ATL11_index/
#     ATL11_index_0331_007_04/ (listed).
#   - read_ATL11_at() filters the search to --ATL11_release, so with today's
#     args every tile would find ZERO along-track granules.
#   - Crossovers are NOT affected: cycles 1 and 2 are still 007_03 (14 of
#     them in the box); the 18 granules at 007_04 are cycles 30 and 31 only.
# CONSEQUENCE: no prelim tile -- IS, GL or AA -- can be solved on the current
# args until a 0332_007_05 index is staged and --ATL11_release (and probably
# --cycles, which names the products) change.  That also blocks L6's DPS
# smoke tile.  Recorded so nobody discovers it by submitting 28 jobs that all
# come back empty.
# DECIDED (Ben 2026-09-17, AL1): "After the current round of code changes are
# complete we will transision to cycles 03-32."  So the transition follows
# this work rather than gating it, and gets its OWN plan document, written
# when Ben asks -- index staging, the args changes, and what re-runs.
# THEREFORE: L1-L5 (all the code) proceed now; L6 and L7, which need a tile
# that can actually be solved, wait for that transition.
#
#
# ===========================================================================
# BACKGROUND.  All STATEMENT, 2026-09-17.
# ===========================================================================
# WHAT THE TILE HOLDS TODAY: /meta.attrs['input_files'], a comma-joined list
# of granule BASENAMES (ATL11_to_ATL15.py:878, save_fit_to_file), one entry
# per beam pair read, so each name repeats ~3x.  Prelim tiles only; matched
# tiles write '' (they read prelim tiles, not ATL11).  The netCDF step reads
# lineage from the PRELIM tiles (write2nc --tiles_dir defaults to prelim).
#
# WHAT THE UNCERTAINTY STEP DOES TO /meta: save_errors_to_file() writes grids
# with to_h5 and adds errors_build_* attributes; it does not touch
# input_files.  A lineage group written by the fit survives it.
#
# WHAT THE GRANULES CARRY, read from one of each over NSIDC S3:
#   ATL11_023003_0332_007_05.h5                     open+read 1.18 s
#     METADATA/DatasetIdentification.attrs['uuid']  b'a805a367-...'
#     ancillary_data/ start_geoseg 354552  end_geoseg 446050
#                     start_orbit 3205     end_orbit 43428
#                     start_rgt 230 end_rgt 230, start/end_cycle 3/32,
#                     start/end_region 3/3                       (all int32)
#   ATL11XO_AR_E1200_N-2600_c01_007_03.h5           open+read 0.14 s
#     uuid b'11a11d42-...'
#     ancillary_data/ start_geoseg 352080  end_geoseg 649385
#                     start_rgt 238        end_rgt 1381          (int32)
#                     NO orbit, cycle or region datasets.
#   So every attribute the lineage needs exists in the granule.
#   AND THE XO end_rgt BUG IS REAL: attributes_for_ATL11_file() sets
#   end_rgt = start_rgt after reading, which would have written 238 for
#   1381.  The same line is on main (ATL14_attrs_meta.py:140) and pre_rel006
#   (ATL1415_attrs_meta.py:173), so products those branches wrote with XO
#   lineage carry it.  FIXED ON on_s3 2026-09-17 (Ben) -- see L4d.
#
# WHERE THE SOLVE HAS EACH GRANULE:
#   crossovers  read_ATL11_xovers() (read_ATL11.py) holds an OPEN h5py
#               handle for every tile it reads: `with pc.io_utils.open_h5(
#               xover_file, fs=fs) as h5f`.  Reading the attributes there
#               costs no extra open.  Local and cloud alike.
#   along-track read_ATL11_at() does NOT: cloud reads go through
#               pointCollection's read_ATL11_granule_cloud_items ->
#               geoIndex.query_xy_box(remote_file=s3_url), which opens and
#               closes the granule inside pointCollection.  Local reads go
#               through geoIndex().from_file(index).query_xy_box().  What
#               comes back is D11.filename per item, which is the granule's
#               path or URL (a ':pairN' suffix is possible from geoIndex;
#               pc.io_utils.strip_pair_suffix removes it).
#
#
# ===========================================================================
# L1. [DECIDED 2026-09-17 by Ben, AL2]  The storage layout.
# ===========================================================================
# DECIDED (Ben, AL2): "per-granule groups under /meta/lineage with only the
# file-only attributes".  As recommended, so:  one HDF5 group per granule,
#     /meta/lineage/<granule basename>          e.g. .../ATL11_023003_0332_007_05.h5
# holding, as attributes, exactly what was read FROM THE GRANULE, in the
# granule's own types:
#     uuid (str), start_geoseg, end_geoseg            both formats
#     start_orbit, end_orbit                          along-track
#     start_rgt, end_rgt                              crossover
# An attribute the granule did not provide is ABSENT, not a sentinel -- the
# netCDF step turns absence into 'NOT_SET' (the DECIDED rule), so there is
# one place that decides what invalid looks like.
# THE TILE KEEPS THE GRANULE'S OWN TYPES (Ben 2026-09-17, on AL4): int32
# geoseg/orbit/rgt stay int32 here, and only the netCDF writer converts to
# strings (L4).  The tile is then a faithful copy of what was read.
# WHY per-granule groups: the name is the key the reader already has
# (input_files), duplicates collapse for free, and h5ls shows it plainly.
# WHY only the file-only attributes: everything else (shortName, cycles,
# release, version, along-track rgt/region) already comes from the name in
# attributes_for_ATL11_file(), and storing it twice invites disagreement.
# KEEP input_files as it is: the reader still needs the list, and old tiles
# have only that.
#
#
# ===========================================================================
# L2. [NEEDS CODE: ATL1415/read_ATL11.py]  Read the attributes at solve time.
# ===========================================================================
# RECOMMENDATION:
#   a. one helper, lineage_attributes(h5f) -> dict, reading the L1 set from an
#      open h5py file: uuid from METADATA/DatasetIdentification, and each
#      ancillary_data/<name>[0] that exists.  Pure h5py, so it is testable
#      against a small synthetic file and works local or cloud.
#   b. read_ATL11_xovers(): call it on the handle already open, for each file
#      that is appended to xover_files_used.  No extra I/O.
#   c. read_ATL11_at(): after the in-bounds filter, for each UNIQUE
#      strip_pair_suffix(D11.filename) that contributed data, open it once
#      with pc.io_utils.open_h5(name, fs=fs) and call the helper.
#      COST, from the probe: ~0.1-1.2 s per granule.  Counted over all 28 IS
#      prelim tiles: 13-15 unique along-track granules each (the only ones
#      that need an extra open) and 2-8 crossover granules (no extra open),
#      so ~2-18 s against a 1169-4260 s prelim job (<2%).  L6
#      measures it for real.
#      DECIDED (Ben, AL3): "Open the files a second time" -- so the attributes
#      are NOT captured inside pointCollection's query, and no second repo
#      changes.
#   d. return the lineage alongside the file lists.  read_ATL11() currently
#      returns (data, file_list) and ATL11_to_ATL15 keeps S['file_list'];
#      RECOMMENDATION: a dict {basename: attrs} carried as S['lineage'],
#      leaving file_list and its callers unchanged.
# A READ FAILURE DOES NOT FAIL THE TILE.  DECIDED rule: unavailable lineage is
# marked invalid in the product.  So: catch per granule, print ONE line naming
# the granule and the exception, store what was read (possibly nothing), and
# go on.  Failing a ~2000 s solve over metadata would be the wrong trade.
# The guard against a SYSTEMATIC failure (e.g. credentials) silently emptying
# every tile is L6's smoke gate plus the netCDF INVALID warning -- not a crash.
#
#
# ===========================================================================
# L3. [NEEDS CODE: ATL1415/ATL11_to_ATL15.py]  Write it into the tile.
# ===========================================================================
# In save_fit_to_file(), beside input_files: for each granule in S['lineage'],
# require_group('/meta/lineage/<basename>') and set its attributes.
# Prelim only by construction -- a matched or error run has no lineage.
# STATEMENT: save_fit_to_file removes and recreates the file, and the
# uncertainty step appends without touching /meta/lineage (Background), so
# nothing else needs changing.  save_field_size_report is unaffected.
#
#
# ===========================================================================
# L4. [NEEDS CODE: ATL1415/ATL1415_attrs_meta.py]  Read it in the netCDF step.
# ===========================================================================
# Replaces the temporary half of 28b4f72:
#   a. set_lineage() collects, per granule name, the /meta/lineage/<name>
#      attributes from every prelim tile that lists it.
#   b. the SAME granule in two tiles must agree.  A disagreement (two uuids
#      for one name) is an ERROR naming both tiles -- it means mixed
#      generations, and a product must not paper over that.
#   c. attributes_for_ATL11_file() takes the stored attributes and fills the
#      file-only ones from them; anything absent stays 'NOT_SET'.
#   d. [DONE 2026-09-17, Ben: "fix the start_rgt=end_rgt bug"]  end_rgt =
#      start_rgt now applies to ALONG-TRACK ONLY (rgt from the name; ATL11 has
#      one rgt); a crossover keeps whatever start_rgt/end_rgt it has -- today
#      NOT_SET, after this plan the values the granule recorded.
#   e. FORCE EVERY LINEAGE ATTRIBUTE TO A STRING before setncattr.  DECIDED
#      (Ben, AL4): "Force to strings", in the netCDF only.  So an int32
#      start_geoseg is written as '354552', and the attribute's type no longer
#      depends on whether any row is invalid.  NOTE, so it is not a surprise:
#      this CHANGES the product from what the pre-28b4f72 code wrote whenever
#      every row was valid -- those attributes came out as int arrays.  The
#      values are the same, the type is now always text.
#   f. the INVALID warning counts only files still missing something, so a
#      fully recorded run prints none -- that is the check in L7.
#   g. old tiles (no /meta/lineage) keep today's behaviour: name attributes
#      only, file-only attributes invalid.
# STATEMENT, tested 2026-09-17 with netCDF4 in the ATL14 env: setncattr on
# [int32 354552, 'NOT_SET'] reads back as the STRINGS ['354552', 'NOT_SET'];
# on [int32, int32] as an int32 array.  That is the behaviour (e) removes by
# converting everything to strings first, so L7 should see text attributes
# whether or not a row is invalid.
#
#
# ===========================================================================
# L5. [NEEDS CODE: tests/]  Tests, no network.
# ===========================================================================
#   - lineage_attributes() against synthetic ATL11- and ATL11XO-shaped .h5
#     files: all present; XO without orbit; a missing dataset is absent, not
#     an exception; uuid as bytes and as str.
#   - read_ATL11_at's per-granule collection with open_h5 and the query
#     mocked: pair suffixes collapse, a granule that raises is logged and
#     skipped, others still recorded.
#   - save_fit_to_file writes /meta/lineage/<name> with the right attributes
#     (a small S dict; see tests/test_build_provenance.py for the pattern).
#   - set_lineage() over synthetic tiles: stored values used; absence ->
#     NOT_SET; conflicting uuids for one name raise, naming both tiles; XO
#     end_rgt preserved; along-track end_rgt == start_rgt; old tiles still
#     work; EVERY attribute written is a string, with and without an invalid
#     row (AL4).
#   - OPT-IN NETWORK TEST, like ATL1415_TIDE_NETWORK_TESTS: read the two
#     probe granules above through the real helper.  RECOMMENDATION: env
#     ATL1415_NSIDC_NETWORK_TESTS=1.
#
#
# ===========================================================================
# L6. [DPS] [BLOCKED: the cycles 03-32 transition (AL1), see the top]
#     Deploy and smoke ONE tile.
# ===========================================================================
# Commit and push; Ben registers; scripts/maap/check_build_id.py must say
# MATCH (image == CWL; see howto_MAAP_ogc O5) before any job.
# Smoke one NAMED prelim tile -- RECOMMENDATION: E1300_N-2500, a full-3x3
# center with crossovers.  GATE, before any fan-out:
#   - the job is successful, and its wall time vs the 2026-09-15 run says
#     what the lineage reads cost (L2c's estimate);
#   - fetch the tile; every name in meta/input_files has a /meta/lineage
#     group; every along-track group has uuid, geoseg and orbit, every
#     crossover group uuid, geoseg and rgt; no uuid is empty;
#   - the job log has no lineage read-failure lines.
# If the gate fails, stop: a fan-out would spend worker-hours writing tiles
# whose lineage is invalid.
#
#
# ===========================================================================
# L7. [ADE] [BLOCKED on L6]  The netCDF check, on the smoke tile.
# ===========================================================================
# write2nc over a --tiles_dir holding just the smoke tile (scratch, not the
# region dir): no INVALID warning; METADATA/Lineage/ATL11 has a real uuid on
# every row; XO rows carry end_rgt != start_rgt where the granule says so.
#
#
# ===========================================================================
# L8. [NOT STARTED]  Docs.
# ===========================================================================
#   - plan_IS_run.sh I9g2: LONG TERM -> points here, then DONE.
#   - howto_MAAP_arctic.sh step 10: drop the "lineage is invalid for now" note
#     once L6-L7 pass on a real run.
#   - howto_MAAP_ogc / staging: nothing -- no new CWL input, no new args.
#   - the cycles 03-32 transition gets its own plan document (AL1); when it is
#     written, link it from the blocker note at the top of this file.
#   - pre_rel006 is NOT touched (AL5): "Do not change pre_rel006.  This will
#     come in with on_s3" -- including the XO end_rgt fix already on on_s3.
#
#
# ===========================================================================
# QUESTIONS FOR BEN -- ALL ANSWERED 2026-09-17, in his words (AL*)
# ===========================================================================
# QL1. THE ATL11 GENERATION (top of file).  CMR now serves ATL11 0332_007_05
#      only, and nothing can be solved on the current 0331_007_04 args.  Is
#      the IS re-run meant to use 0332_007_05 -- which needs its index staged
#      (the per-granule tree, ~8100 files for the 0331 generation) and new
#      --ATL11_release / --cycles in the args -- and should that be planned
#      as its own step before L6?  L1-L5 do not wait on this; L6 does.
# AL1: "After the current round of code changes are complete we will
#      transision to cycles 03-32."
#      FOLLOW-UP, answered: that transition gets its OWN plan document, to be
#      written when Ben asks -- not steps in this one.  L1-L5 proceed now;
#      L6-L7 wait for it.
#
# QL2. Layout (L1): per-granule groups under /meta/lineage with only the
#      file-only attributes -- or do you want every granule attribute stored,
#      or a different location?
# AL2: "per-granule groups under /meta/lineage with only the file-only
#      attributes"
#
# QL3. Cost (L2c): an extra open per along-track granule, estimated <2% of a
#      prelim job.  Acceptable, or should the attributes be captured inside
#      pointCollection's query instead (a second repo)?
# AL3: "Open the files a second time"
#
# QL4. Attribute types (L4): a lineage attribute becomes a string array if
#      any row is 'NOT_SET', an int array otherwise.  Leave it, or force one
#      type (e.g. strings always) so the product format does not depend on
#      validity?
# AL4: "Force to strings."
#      FOLLOW-UP, answered: in the netCDF ONLY.  The tile keeps the granule's
#      own types, so it stays a faithful copy of what was read; the writer
#      converts.
#
# QL5. Discover (pre_rel006) keeps opening ATL11 at write time.  Should the
#      same change go to pre_rel006 later, or is it on_s3/MAAP only until the
#      merge (docs/plan_merge_to_main.sh)?
# AL5: "Do not change pre_rel006.  This will come in with on_s3"
