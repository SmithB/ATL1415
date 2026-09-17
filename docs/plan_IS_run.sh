# plan_IS_run.sh -- THE ICELAND (IS) RUN: prelim and matched, end to end
#
# ############################################################################
# ##  STATUS 2026-09-17: RUN THROUGH netCDF.  prelim 28/29 and matched      ##
# ##  28/28 on DPS, mosaic 41/41 and all five netCDFs in the ADE -- with   ##
# ##  lineage INVALID by design (I9g2).  IS WILL BE RE-RUN COMPLETELY once ##
# ##  more issues are fixed.  Open: I9g6, I9h.  The banner below is the  ##
# ##  2026-09-12 original, kept for the record.                            ##
# ##                                                                        ##
# ##  TENTATIVE.  Written 2026-09-12.  No IS tile has been solved on DPS.   ##
# ##  This is the sequence, the status of each step, and the decisions      ##
# ##  that were NOT recorded anywhere -- not a record of a run.             ##
# ##                                                                        ##
# ##  ALL EIGHT QUESTIONS (QI1-QI8) WERE ANSWERED BY BEN ON 2026-09-12 and  ##
# ##  are recorded at the foot of this file.  Two pieces of code remain     ##
# ##  unwritten -- I2's submit_MAAP_jobs.py and I7's run.sh/config change,  ##
# ##  the second of which needs a rebuild and a registration.               ##
# ############################################################################
#
# WHY IS: Q21's answer.  Iceland is the bounded first target, chosen because it
# is small enough to fan out under the ~10 jobs/hr public-queue throttle.
#
# THIS FILE DOES NOT RESTATE PROCEDURE.  howto_MAAP_arctic.sh steps 0-10 are
# the procedure and howto_MAAP_ogc.sh O1-O12 the job mechanics; O9's rule is
# that a procedure is replaced, never duplicated, so two copies cannot drift.
# What is here: the order, what is ready TODAY, what is blocked, and the open
# questions -- numbered QI1..QI8 -- that have to be answered by Ben before the
# blocked steps can be written.  Steps are I1..I10, matching the arctic
# howto's numbering where they correspond.
#
# Tags: [READY] runnable as written today.  [BLOCKED: x] x does not exist.
#       [DECISION] waits on a QI below.  [ADE] / [DPS] where it runs.
#
# EVERY CLAIM BELOW IS LABELLED.  STATEMENT = verified, with how.
# RECOMMENDATION = my suggestion, yours to take or drop.  QUESTION = I do not
# know and did not guess.


# ---------------------------------------------------------------------------
# THE SIZE OF THE RUN -- so the shape of the decisions is visible
# ---------------------------------------------------------------------------
# STATEMENT, measured 2026-09-12 by reading the 40 km mask off the bucket
# (/vsis3/.../RGI_reduced/06_rgi60_Iceland_reduced_40km.tif, 7x10 cells at
# 40 km) and applying the same selection make_ATL1415_queue.py's .db branch
# uses (mask_G.z == 1, no min/max filters):
#     IS HAS 29 TILE CENTERS.  x 1020..1380 km, y -2620..-2380 km.
# The "13 tiles" in Transition_to_maap.md is an Iceland BOX used for the
# crossover check, not the region.
#
# STATEMENT: all 29 lie 2774-2961 km from the north pole -- past the far end
# of the AA transect (2420 km), where the fit+error curve is flat.
# RECOMMENDATION, rough: ~0.8-0.9 h and well under 16 GiB per prelim tile
# (staging S6's curve), so IS prelim is ~25 worker-hours and ~3 h of
# submission at the public throttle.  CAVEAT: that curve is Antarctic 60 km
# tiles; Arctic data density is not the same, and IS is the first region to
# measure it.  Matched is NOT measured anywhere (S6 says so).
# STATEMENT: a prelim tile is ~240 MB (the AA transect tile E-580_N-980 is
# 239429355 bytes), so IS prelim is ~7 GB to move.


# ---------------------------------------------------------------------------
# WHAT IS ALREADY TRUE -- verified 2026-09-12, so it is not re-litigated
# ---------------------------------------------------------------------------
# S1. The ATL14 env EXISTS at /home/jovyan/.conda/envs/ATL14 -- in the
#     persistent home, not the /srv overlay that ate the 09-06 one.  Verified
#     by importing pointCollection from it and reading the mask off the bucket.
# S2/S3. Masks and the ATL11 index are on the bucket (staging S2, S3).  The
#     IS .db and its _40km.tif and _80km.tif siblings are there (959 and 591
#     bytes), listed 2026-09-12.
# The args file is composed AND published, and STILL CURRENT: 970 bytes at
#     s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/north/IS/input_args_IS.txt
#     Verified 2026-09-12 that nothing it depends on has changed since it was
#     written on 09-07: `git log` on MAAP_dps.txt, rel_006_0331.txt, IS.txt,
#     quarterly.txt and setup_ATL1415_region.py since 09-06 shows b293807
#     (before the composition) and e231ba4, whose MAAP_dps.txt diff is
#     COMMENT-ONLY -- no value line changed.  So arctic steps 3 and 4 are done.
# The tree is clean and pushed: HEAD == origin/on_s3 == 38f5ef9.


# ===========================================================================
# I0. [ADE] [DONE 2026-09-15 -- MATCH at b5fe447, recorded at I2; again at
#     6b8a2ca for I7]  Confirm the image before spending worker hours.
# ===========================================================================
conda activate ATL14
cd ~/git_repos/ATL1415
/srv/conda/envs/notebook/bin/python scripts/maap/check_build_id.py
# STATEMENT: the last registration was ab84687 (O6 run 3).  The diff
# ab84687..38f5ef9 touches run.sh in COMMENTS ONLY, plus docs and two
# ADE-side scripts -- verified with `git diff`.  So the deployed image is
# functionally current for a prelim run and NO REBUILD IS NEEDED FOR I1-I6.
# EXPECT "VERDICT: MATCH" at ab84687: O5's rule is that image == CWL is the
# test and origin is only reported, so origin being 6 commits ahead is not a
# mismatch.  A rebuild IS needed for matched-on-DPS -- see I7 and QI5.


# ===========================================================================
# I1. [ADE] [DONE 2026-09-12]  The 29 tile centers, as a frozen xy list.
# ===========================================================================
#   region_files/IS_prelim_xy.txt        29 lines, "<x0> <y0>" in meters
# WRITTEN 2026-09-12 per QI1, by reading
# /vsis3/maap-ops-workspace/ben_smith/ATL1415/masks/RGI_reduced/
#   06_rgi60_Iceland_reduced_40km.tif  (7x10 cells at 40 km)
# and taking every cell with z == 1 -- exactly the selection
# make_ATL1415_queue.py's .db branch makes (lines 216-219), with no min/max
# filter.  Re-derive it with that same read if the mask is ever restaged.
#
# TWO THINGS TO KNOW ABOUT THIS FILE:
#   - IT CARRIES NO COMMENT HEADER, deliberately.  submit_AA_queue.py parses
#     EVERY non-blank line as floats, so a '#' line would crash the submitter;
#     make_ATL1415_queue.py's --xy_list_file would merely warn.  Provenance
#     therefore lives here, in this step, and not in the file.
#   - region_files/ HELD ONLY XR/YR BOUNDS FILES until now (E_ant_test.txt),
#     and the existing xy lists live in scripts/maap/ (AA_queue_xy.txt).  This
#     is a new kind of file in that directory; it is where you asked for it.
#   - AND region_files/* WAS GITIGNORED (.gitignore:157, "mostly obsolete, and
#     were mostly for test runs"), so a list frozen there would not have been
#     versioned -- which is most of what "frozen" means.  Added an exception,
#     !region_files/*_xy.txt, rather than moving the file.  Say if you would
#     rather it sat in scripts/maap/ with the AA lists, which are tracked.
#
# WHY THE DOCUMENTED COMMAND WAS NOT USED --
#   make_ATL1415_queue.py prelim $region_dir/input_args_IS.txt --xy_out ...
# STATEMENT: it cannot run today, for three reasons read out of
# ATL1415/scripts/make_ATL1415_queue.py on 2026-09-12:
#   - there is no --xy_out.  The argparse block (lines 42-57) has
#     --xy_list_file (an INPUT) and --queue_file; the script emits shell
#     command lines, not centers.
#   - line 134, os.path.isfile(defaults['--ATL11_index']) on
#     's3://maap-ops-workspace/ben_smith/ATL11_index/' is False, so it re-joins
#     the URI onto --ATL14_root, is False again, and sys.exit(1)s.  Blocks
#     outright, before any mask is read.
#   - line 214, os.path.isfile(mask_base+'_40km.tif') on the s3 URI is False,
#     so it raises OSError on a file that is sitting on the bucket.
# STATEMENT: the two OTHER cloud bugs in the "Code that has to exist" list do
# NOT bite for IS.  Line 109's --tide_mask_file test is inside
# `if '--tide_mask_file' in defaults`, and input_args_IS.txt has no such key.
# Line 86's greedy regex does drop the bare flags (--ATL11_earthaccess,
# --previous_product_earthaccess), but nothing in the tile-center path reads
# them -- it changes no center.
# ALL THREE STAY OPEN, and GL cannot dodge them: it has no frozen list, and
# its centers come from a 1 km mask that does not exist (Q6/Q16).


# I2. [DPS] [DONE 2026-09-15 -- smoke tile, then all 29 submitted]
#     Fan out prelim.
# ===========================================================================
# ONE TILE FIRST.  -16gb is a queue nothing has run on, and algorithm_config's
# ram_min is 16: O5 notes a floor at or above what a queue offers risks a job
# that never schedules.  One job says, and costs one job.
#
# THE SMOKE TILE IS NAMED, NOT THE FIRST LINE OF THE LIST.  region_files/
# IS_smoke_xy.txt holds 1260000 -2620000 and nothing else -- the tile Ben gave
# on 2026-09-05 (43K points, an ice-sheet edge tile), which is the one arctic
# step 2 rasterizes (601x601, 72638 ice cells, 20.1%) and the one staging S7
# smoke-tested.  It IS in IS_prelim_xy.txt, on line 12.
# CORRECTED 2026-09-12: the first draft of this step said --limit 1, which
# takes the FIRST line of the sorted list -- E1020_N-2580, a tile nobody has
# ever looked at.  --limit stays in the submitter for GL-scale spot checks,
# but a smoke test should name its tile.
scripts/maap/submit_MAAP_jobs.py --xy_file region_files/IS_smoke_xy.txt \
    --step prelim --args_url $s3_run/input_args_IS.txt \
    --queue maap-dps-worker-16gb \
    --ledger ~/ATL14_processing/maap_ledgers/IS_smoke_jobs.csv
scripts/maap/collect_jobs.py ~/ATL14_processing/maap_ledgers/IS_smoke_jobs.csv
# Wait for 'successful'.  If it will not schedule or dies on memory, -32gb is
# the fallback and nothing else changes.
# DONE 2026-09-15.  STATEMENT, read with collect_jobs.py at 18:15 UTC and
# `aws s3 ls` on the canonical prefix:
#   job c5a37da1 (submitted 16:18, WITH --tile_prefix), image 6978a8a:
#     successful, 2336 s wall (fit 1277 s, error 1038 s), peak 6.55 GiB on
#     maap-dps-worker-16gb -- it schedules, and fits with room.
#     N_ATL11 239700, N_AT 239613, N_XO 87 (> 0, so crossovers are read),
#     N_fit 43789, 3 iterations.
#   s3://.../ATL14_processing/rel006/north/IS/prelim/E1260_N-2620.h5
#     20317493 bytes, and prelim/field_sizes/E1260_N-2620_report.json
#     (dz/dz 61x61x32), both written 17:01:57 -- after the error step, as
#     run.sh intends.  SO THE PRELIM UPLOAD I7 COULD NOT TEST IS VERIFIED.
# THE SIZE ESTIMATES ABOVE WERE HIGH.  RECOMMENDATION, one tile only: ~0.65 h
# per tile (not 0.8-0.9) and ~20 MB (not ~240 MB), so IS prelim is nearer
# ~19 worker-hours and ~0.6 GB.  An edge tile; interior tiles may differ.
# THEN the other 28:
scripts/maap/submit_MAAP_jobs.py --xy_file region_files/IS_prelim_xy.txt \
    --step prelim --args_url $s3_run/input_args_IS.txt \
    --queue maap-dps-worker-16gb \
    --ledger ~/ATL14_processing/maap_ledgers/IS_prelim_jobs.csv
# (the smoke tile is submitted again as part of the 29; dedup=False, so it
# really re-runs.  Simpler than excising one line, and it costs ~1 h.)
# SUBMITTED 2026-09-15 ~18:30 UTC, all 29 WITH --tile_prefix
#   s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north/IS
# (Ben chose all 29 over excising the smoke tile).  29/29 accepted, no failed
# submit; ledger ~/ATL14_processing/maap_ledgers/IS_prelim_jobs.csv.
# GATE, same day: the process had been re-registered at b5fe447 (modified
# 16:56) after the smoke job went out, so check_build_id.py was re-run first:
# VERDICT: MATCH, job 29688742 -- stamp, live git and CWL all b5fe447,
# tree_state=clean, maap_py=5.1.0, maap_pgt=set; the CWL keeps default '-'
# on both tile_prefix inputs.  b5fe447 differs from 6978a8a (the smoke
# tile's image) only in submit_MAAP_jobs.py and this file -- nothing a
# worker runs.
#
# ONCE THE REBUILD OF I7 IS DEPLOYED, add --tile_prefix to both commands so the
# tiles land in the canonical tree and the matched step can find them:
#   --tile_prefix s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north/IS
# Without it the prelim tiles are only in dps_output, and I4's fetcher is the
# only way to reach them -- which is fine for prelim and fatal for matched.
#
# DECIDED 2026-09-12 (Ben) per QI2: generalize rather than special-case.
# submit_AA_queue.py STAYS AS IT IS -- it carries Antarctica's two-width
# routing (60 km north of the 400 km line, 44 km south, deliberately
# overlapping so a tile in the band is submitted twice), which is real for AA
# and for nothing else.  The new script is for one geometry and one args file.
# Shared through ogc_jobs.py, not copied: the process lookup by name+version,
# submit_job(pid, inputs, queue, dedup=False, tag=...), the ledger columns.
#
# WHAT IT REFUSES, and why each one is there:
#   - an existing ledger, without --replace.  The ledger is the ONLY record of
#     what was submitted; overwriting it strands the worker-hours it names.
#   - --step matched without --tile_prefix (I7), and --tile_prefix at all
#     while algorithm_config.yml does not declare it.  THE SECOND CHECK IS
#     OFFLINE ON PURPOSE: the first version asked the deployed CWL and treated
#     an unreadable one as yes, and repo.maap-project.org duly timed out
#     during testing -- which would have turned the guard into a pass and cost
#     one failed job per tile.  The config is consulted first and needs no
#     network; the CWL is then a cross-check that only ever adds a refusal.
#   - a queue name that is not maap-dps-*, an args_url that is neither an
#     s3:// URI nor a file, and any xy line that does not parse (an error, not
#     a skip: a silently dropped line submits a region short by a tile).
# A submit that fails is RECORDED in the ledger and the run continues (Q11).
# Rows are flushed as they are written, so an interrupt still leaves a usable
# ledger.  --rate (default 2 s) and --max_in_flight are there for GL.
# ADDED 2026-09-15 (Ben): a ninth ledger column, tile_prefix, appended last --
# the prefix the job was sent, or "-" for none -- so QI5a's "the path is
# visible in every ledger row" is true of the ledger and not only of the job
# record.  Readers go by column name, so older ledgers still read.
# IS_smoke_jobs.csv was written before the column and was BACKFILLED with the
# prefix its submit command passed (the original is kept beside it as
# IS_smoke_jobs.csv.pre_tile_prefix).  Tested against a stubbed MAAP: the
# column holds the prefix with --tile_prefix, "-" without it or with "-",
# and the inputs sent are unchanged.
#
# TESTED 2026-09-12, all against the real deployed process (processID 64):
# the 29-center dry-run, --limit 1, and every refusal above.  RUN FOR REAL
# 2026-09-15: the smoke tile, then the 29-tile fan-out, both above.


# I3. [ADE] [DONE 2026-09-16 -- 28 successful, 1 failed (I7a)]  Watch the 29 jobs.   THE COLLECTOR.
# ===========================================================================
scripts/maap/collect_jobs.py ~/ATL14_processing/maap_ledgers/IS_prelim_jobs.csv
# THE COLLECTOR READS ABOUT TILES; IT DOES NOT MOVE THEM.  Moving them is the
# TILE FETCHER, I4.  Both were called "the collector" until 2026-09-12.
# STATEMENT: this was collect_AA_queue.py until 2026-09-12; the name was the
# only AA-specific thing about it -- it
# reads identifier/job_id/queue out of any ledger with the submitter's columns
# and reports status, wall clock, self-measured peak RSS, N_ATL11/N_AT/N_XO,
# iterations and the build each tile ran.  Read from the code, 2026-09-12.
# RECOMMENDATION: for 29 jobs this is enough, and check_MAAP_jobs.py's
# --requeue can wait for a region where hand-resubmitting is impractical.
# RECOMMENDATION: keep the ledger at ~/ATL14_processing/maap_ledgers/, OUTSIDE
# the checkout -- an untracked file in the repo makes register_algorithm.py
# refuse (O8).
# EXPECT N_XO > 0 on every tile.  If N_XO == 0 across the board, crossovers
# are not being read for this region and I7 should not start.
#
# RUN 2026-09-16 02:20 UTC, ~8 h after the 18:24 submission.  ALL 29 FINISHED.
# STATEMENT: 28 successful, 1 failed -- IS_prelim_E1020_N-2580, which is I7a.
# STATEMENT: THE N_XO GATE PASSES.  N_XO is 48..116 across the 29, never 0, so
# crossovers ARE being read for IS and I7 is not blocked on this.  N_ATL11 is
# 233k..326k, N_AT within ~100 of it every time.
# STATEMENT: cost, for sizing the production queue -- wall clock 1169..4260 s
# (median ~2000), peak RSS 4.16..9.09 GiB self-measured.  The 16 GiB queue
# (QI7) is right: the worst tile used 9.09 GiB, and nothing came near 16.
#
# STATEMENT, AND IT LOOKED LIKE A PROBLEM BUT IS NOT: the collector reports
# TWO builds ran this set -- b5fe447 for 24 tiles, a46ad52 for 5.  The 29 were
# submitted at 18:24 against b5fe447, and a46ad52 was built at 18:48, WHILE
# THEY WERE STILL QUEUED, so the later starters picked it up.
# `git diff --stat b5fe447..a46ad52` is docs/plan_IS_run.sh ONLY (+34 -3).
# THE CODE IS IDENTICAL; only this file differs.  The 29 tiles are one run.
# RECOMMENDATION: do not re-run anything over this.  But note the mechanism --
# committing and rebuilding while a fan-out is queued silently splits which
# build runs it, and only the collector's per-tile commit column shows it.
# A docs-only commit made it harmless HERE; a code commit would not have been.


# ===========================================================================
# I4. [ADE] [DONE 2026-09-16 -- 28 fetched, 652 MB]
#     THE TILE FETCHER: bring the prelim tiles down.
# ===========================================================================
scripts/maap/fetch_tiles.py ~/ATL14_processing/maap_ledgers/IS_prelim_jobs.csv \
    $region_dir --step prelim --dry-run
scripts/maap/fetch_tiles.py ~/ATL14_processing/maap_ledgers/IS_prelim_jobs.csv \
    $region_dir --step prelim
#
# NOT `aws s3 sync $s3_out/prelim/ $region_dir/prelim/`, which every howto
# still says and which cannot work: Q9's deterministic output prefix is
# unimplemented, so each job's products sit under its own timestamped
# dps_output prefix, reachable only through the job id in the ledger.
#
# STATEMENT, the finding this was built on -- resolved 2026-09-12 for AA job
# 55c01ec3 and listed on the bucket:
#   s3://maap-ops-workspace/ben_smith/dps_output/atl1415_tile_solve_1786/
#     on_s3/2026/09/11/16/43/24/959688/
#       _stdout.txt  _stderr.txt  outputs_result-*.{context,dataset,met}.json
#       prelim/E-580_N-980.h5                       (239 MB)
#       prelim/field_sizes/E-580_N-980_report.json
# run.sh's output/prelim/ SUBDIRECTORY SURVIVES THE UPLOAD, and the tile
# really is a product.  So the fetch needs no rebuild and no new CWL input.
#
# STATEMENT: THE TWO STEPS UPLOAD DIFFERENTLY, read out of run.sh 2026-09-12.
# prelim passes --base_directory $PWD/output and the solver appends '/prelim'
# (run.sh:373), so a prelim tile is at <prefix>/prelim/E<x>_N<y>.h5.  matched
# passes --out_name $PWD/output/E<x>_N<y>.h5 (run.sh:419), so a matched tile
# is at the TOP of the prefix.  The fetcher tries both for every row, so it
# does not depend on which is right -- the prelim layout is VERIFIED, the
# matched one is read from run.sh and UNVERIFIED, no matched job having run.
#
# STATEMENT: a successful job with no tile is NORMAL.  ATL11_to_ATL15 returns
# 0 without writing when a tile has too little data and run.sh then exits 0
# (run.sh:379).  Those rows report 'no tile' and are counted apart from
# failures.  Expect some among IS's 29: it is a small, coastal region.
#
# WHAT QI4 CHANGES HERE, decided 2026-09-12: once run.sh writes the canonical
# key itself, the fetcher's job shrinks to bringing the bucket tree down for
# the mosaic step -- which then really is `aws s3 sync`, as the howtos always
# said.  Keep the fetcher regardless: it is the only way to reach the tiles
# already solved on the current image, and the only one that works if a job's
# canonical write ever fails.
#
# RUN 2026-09-16, dry-run first, then for real.  28 FETCHED, 0.64 GiB, one row
# NOT FETCHED: IS_prelim_E1020_N-2580, reported as FAILED (I7a).  Local tree
# /home/jovyan/ATL14_processing/rel006/north/IS/prelim/ now holds 28 .h5
# (652 MB on disk) and field_sizes/ holds 28 _report.json.  The prelim layout
# <prefix>/prelim/E<x>_N<y>.h5 held for all 28, as VERIFIED at I2.
#
# STATEMENT, AND IT CONTRADICTS WHAT THIS STEP EXPECTED: there were ZERO
# 'no tile' rows.  "Expect some among IS's 29" was wrong -- every job that
# exited 0 wrote a tile.  The one center that produced nothing did so by
# FAILING, not by the tidy exit-0 path, which is exactly the I7a bug: IS does
# have a center with too little data, and the current code cannot express that
# outcome once the fit has already written a file.  After the I7a fix the
# 'no tile' row this step was written for is the one E1020_N-2580 will produce.


# I5. [ADE] [DONE 2026-09-17 -- scripts/check_field_sizes.py, 28/28 prelim and
#     28/28 matched pass]  Look at the tile sizes.   (arctic step 8's ADE half)
# ===========================================================================
# ATL11_to_ATL15 writes the field-size report itself; the JSON above is it.
# Check the reports with the field-size checker (docs/plan_check_field_sizes.sh):
scripts/check_field_sizes.py $region_dir/prelim @$region_dir/input_args_IS.txt
scripts/check_field_sizes.py $region_dir/matched @$region_dir/input_args_IS.txt
# RESULT 2026-09-17: both OK -- expected dz/dz [61, 61, 32] derived from -W,
# -g and -t; prelim sigma == dz/dz, matched sigma null; 28 reports, 28 tiles,
# 28 of 28 passed, 0 problems, each.
# check_tiles.ipynb, which this step first named, was NEVER FOUND anywhere
# under /home/jovyan, and Ben said on 2026-09-16 it is not needed: the checker
# replaces it.  The question below is kept for the record and is closed.
#
# QUESTION (2026-09-16), CLOSED 2026-09-17 (see above): check_tiles.ipynb DOES NOT EXIST.  Not in the repo,
# and `find /home/jovyan -name 'check_tiles*'` returns nothing at all.  Is it
# unwritten, or is it on discover / another machine?  Nothing else in this
# plan depends on it -- the mechanical half is done below -- so this is not
# blocking, but the step as written cannot be followed.
#
# THE MECHANICAL HALF, RUN 2026-09-16 over the 28 fetched reports:
# every one carries both fields and both are [61, 61, 32] -- 61x61 grid cells,
# 32 epochs -- with no zero, empty or missing value anywhere.  The tiles are
# uniform, which is what a 40 km tile at this spacing should give.
# STATEMENT: this checks SHAPE, not content.  A tile can be the right shape
# and still be wrong; that is what the notebook was presumably for.


# ===========================================================================
# I6. [ADE] [DONE 2026-09-16 -- region_files/IS_matched_xy.txt, 28]  The matched tile list.
# ===========================================================================
# The same 29 centers -- region_files/IS_prelim_xy.txt serves both steps.
# STATEMENT: expect the matched list to be SHORTER than 29, and that is not a
# failure.  ATL11_to_ATL15 returns 0 without writing when a tile has too
# little data (run.sh:379 exits 0 on that path), so some centers will have no
# prelim tile to match.  Drop those rows from the matched ledger after I4
# reports them as 'no tile'.
#
# STATEMENT 2026-09-16: the list is 28, not 29, and the dropped center is
# E1020_N-2580 -- but it was dropped for the WRONG REASON, a failed job rather
# than a clean 'no tile' (I7a).  The count is the same either way: that center
# has no data inside the mask, so after the I7a fix it will still write no
# tile and still be dropped.  I6 therefore does NOT wait on I7a.
# RECOMMENDATION: build the matched ledger from what I4 actually fetched --
# the 28 .h5 on disk -- rather than from IS_prelim_xy.txt minus a hand-kept
# exclusion, so the list cannot disagree with the tiles that exist.
# NOTE for I7: E1020_N-2580's 8 neighbours will each be missing one neighbour.
# That is the normal case QI5b already decided -- named in the log, solve
# proceeds -- and needs no action here.
#
# CORRECTION 2026-09-16, measured: THAT NOTE IS WRONG IN PRACTICE.  At the
# tile_spacing the args file actually carries (40000; -W is 60000 and run.sh
# prefers --tile_spacing), E1020_N-2580 has NO neighbour in the 28-tile set at
# all -- its 8 neighbour positions are E980/E1020/E1060 x N-2620/N-2580/N-2540,
# none of which is a center in this region.  So losing it costs NO other tile
# a neighbour, and nothing in I7 is affected by it.
#
# DONE 2026-09-16: region_files/IS_matched_xy.txt, 28 centers, built I6's
# recommended way -- from the tiles that EXIST rather than from the prelim list
# minus a hand-kept exclusion.  The S3 prefix and the local tree were listed
# separately and are IDENTICAL, 28 each; the generated list is exactly the 29
# prelim centers minus (1020000, -2580000), checked both directions.
#
# NEIGHBOUR COVERAGE OF THE 28, measured 2026-09-16 at 40 km, worth knowing
# before reading I7's results: the region is sparse and many tiles are edge
# tiles.  Only THREE centers have a full 3x3 (E1300_N-2500, E1340_N-2500,
# E1340_N-2460).  At the other end, E1020_N-2420 has ZERO neighbours present
# and E1180_N-2380 has one, so their matched solves will fetch 1/9 and 2/9.
# THAT IS EXPECTED, not a fault -- QI5b: missing neighbours are named in the
# log and the solve proceeds, and the tile's OWN prelim file is the only hard
# requirement.  But it means a 'k/9 localized' line well under 9 is the NORM
# for this region, and is not evidence of a broken fetch.


# I7. [DPS] [DONE 2026-09-16 -- 28/28 SUCCESSFUL]
#     The matched solve.
# ===========================================================================
# DECIDED 2026-09-12 (Ben) per QI6: ON DPS, implementing Q8 option (a).  The
# ADE shortcut (option c) was declined -- proving the mechanism on the
# smallest region is worth the rebuild, and GL cannot use the shortcut.
#
# WHAT LANDED, and it is NOT prelim_prefix: the input is tile_prefix, ONE
# input serving both roles -- where a prelim job WRITES its tile, and where a
# matched job READS its neighbours.  Two inputs that must always hold the same
# string would eventually hold different ones.  Say if you want them split.
#   algorithm_config.yml  a fifth input, tile_prefix, default "-" (QI5a).
#                         Empty rather than required so every prelim job that
#                         ran before this still describes a valid call.
#   run.sh                parses it; prelim uploads its tile AFTER the error
#                         step (--calc_error_for_xy writes back into the same
#                         file, so an upload in between publishes a
#                         half-finished tile); matched fetches the 3x3 first
#                         and uploads its result after.
#   scripts/s3_tiles.py   NEW, worker-side, s3fs only -- `put` one tile (and
#                         its field-size report), `get` the 3x3.  Not in
#                         scripts/maap/: that directory is the ADE's maap-py
#                         scripts, and a worker has no business talking to the
#                         job API.
#
# THE SPACING IS READ, NOT ASSUMED: run.sh takes --tile_spacing from the
# composed args file, falling back to -W, which is make_ATL1415_queue.py's own
# precedence.  AA solves two halves at different widths, and a wrong spacing
# would fetch eight tiles that exist but are not this tile's neighbours --
# which no later step could detect.
#
# MISSING NEIGHBOURS (QI5b): named in the log, and the solve proceeds.  The
# tile's OWN prelim file is still required, by the guard run.sh already had --
# one check in one place rather than two that can disagree.
#
# TESTED 2026-09-12, as far as it can be without a rebuild:
#   - s3_tiles get against a real prefix holding one real tile: 1/9 localized,
#     the 239 MB tile really downloaded, the 8 absent neighbours named at the
#     right 60 km offsets, exit 0.
#   - run.sh --step matched --tile_prefix ...: parses, reads spacing 40000 out
#     of input_args_IS.txt, fetches (0/9 for a center with nothing there),
#     then fails the own-tile guard.  exit 2.
#   - no --tile_prefix and no input/prelim: exit 2 with the new message.
#     --tile_prefix with no spacing in the args file: exit 2.
#     --step build_id with a tile_prefix, and with an EMPTY one: exit 0,
#     unchanged -- which matters, because the CWL will now bind
#     --tile_prefix "" on EVERY job, build_id and check_build_id included.
#   NOT TESTED, and untestable here: a real matched solve, and the prelim
#   upload -- both need the rebuilt image.
#   2026-09-15: THE PRELIM UPLOAD IS NOW VERIFIED by the I2 smoke tile.  The
#   matched solve is still untested.
#
# THE EMPTY DEFAULT DID NOT SURVIVE REGISTRATION.  Registered 2026-09-14 with
# default: "", and /api/build accepted it -- but the deployed CWL (deployment
# 151, s:commitHash 52e27bd) declares tile_prefix as `type: string` with NO
# default, i.e. REQUIRED.  The build form drops empty fields (howto_MAAP_ogc
# F6).  check_build_id.py and any prelim job without --tile_prefix omit the
# input, so they would have been submitted without a required value.
# FIXED (Ben chose it, 2026-09-14): default "-", read as "none" by run.sh and
# by submit_MAAP_jobs.py.  A non-empty sentinel cannot be dropped.
# AFTER RE-REGISTERING, CHECK THE CWL SAYS `default: '-'` under tile_prefix
# before trusting any job that omits it.
#
# THEN, in order:
#   1. commit and push -- register_algorithm.py refuses on unpushed work
#   2. Ben registers; wait for the build and the deploy
#      2026-09-14: registered from the WRONG hub image -- its notebook env had
#      maap-py 4.2.0 -- and the rebuild failed (Ben).  register_algorithm.py
#      and every job script now refuse maap-py < 5.0.  Re-register from the
#      right image before step 3.
#      DONE 2026-09-14: re-registered at 6978a8a.  VERIFIED 2026-09-15 by
#      reading the deployed CWL (processID 64, modified 2026-09-14T20:34):
#      s:commitHash 6978a8a, and tile_prefix carries default: '-' on both the
#      workflow input and the CommandLineTool input.
#   3. DONE 2026-09-15: scripts/maap/check_build_id.py said VERDICT: MATCH --
#      image stamp, live git and CWL all 6978a8a, tree_state=clean,
#      maap_py=5.1.0, maap_pgt=set.  Job 39a8ec02-3940-42b8-8d46-1181558a82b8,
#      submitted WITHOUT tile_prefix, so the '-' default really binds.
#      STATEMENT, same day: no IS ledger exists yet, and
#      s3://.../ATL14_processing/rel006/north/IS/ is empty -- no IS job has run.
#   4. ONE matched job before the other 28 -- but I2's prelim, WITH
#      --tile_prefix, has to fill the canonical tree first.
#      2026-09-15: I2's smoke tile is in the canonical tree (see I2), and
#      all 29 prelim jobs are submitted (image b5fe447, MATCH).  A matched job
#      waits for them: it needs its neighbours in the tree.
#      DONE 2026-09-16.  THE MATCHED SOLVE NOW WORKS ON DPS.
#      Smoke: job 38440138, IS_matched_E1300_N-2500, on 6b8a2ca.
#        successful, 402 s (step matched 376 s), peak 8.90 GiB on
#        maap-dps-worker-16gb.  N_fit 245361, 1 iteration.  N_ATL11/N_AT/N_XO
#        are blank, as they should be: matched reads prelim tiles, not ATL11.
#        Wrote .../IS/matched/E1300_N-2500.h5 (53643251 bytes) and its
#        field_sizes report at 18:13:16.  The matched prefix was VERIFIED
#        EMPTY before submitting.
#      CHOSEN DELIBERATELY: E1300_N-2500 is one of only THREE centers with a
#      full 3x3, and all 9 were confirmed on S3 first, so the smoke exercised
#      the complete s3_tiles get fetch -- the part that had never run.
#      MEMORY: 8.90 GiB of 16 is the heaviest-neighbour case; prelim peaked at
#      9.09 GiB on the same queue.  Headroom is adequate but not vast; -32gb
#      is the fallback if any tile OOMs, and nothing else changes.
#
#      NO SIGMA IN A MATCHED TILE, AND THAT IS CORRECT.  The smoke's report
#      reads {"dz/dz": [61,61,32], "dz/sigma_dz": null} where the prelim
#      report for the same center has sigma_dz [61,61,32].  CONFIRMED TWICE:
#      (1) make_ATL1415_queue.py adds the --calc_error_for_xy companion only
#      in the `if not args.step=='matched'` branch -- the matched branch emits
#      a single command with no error pass, so run.sh is faithful to it;
#      (2) Ben, 2026-09-16: "There should be no sigma in a matched result.  We
#      use the uncertainties calculated in the prelim step."
#      DO NOT read a null sigma_dz in matched/field_sizes as a fault, and do
#      not add an error pass to the matched branch of run.sh.
#
#      THEN ALL 28, submitted 2026-09-16 to IS_matched_jobs.csv, queue
#      maap-dps-worker-16gb, WITH --tile_prefix.  28/28 accepted.  The smoke
#      center is re-run as part of the 28 (dedup=False, so it really re-runs)
#      -- the same choice I2 made for prelim, and it costs ~400 s.
#
#      STATEMENT 2026-09-16, collect_jobs.py at 19:16 UTC: 28/28 SUCCESSFUL,
#      NO FAILURES.  All 28 on ONE build, 6b8a2ca -- no repeat of the split
#      build that the prelim fan-out saw (I3).  Every tile 1 iteration.
#        wall     153 .. 724 s   (E1300_N-2620 fastest, E1340_N-2460 slowest)
#        peak RSS 3.64 .. 8.99 GiB on the 16 GiB queue
#        N_fit    158 .. 271556
#      N_ATL11/N_AT/N_XO are blank on every row, correctly: matched reads
#      prelim tiles, not ATL11.
#      THE MEMORY HEADROOM HELD: the worst tile, E1340_N-2460, peaked at 8.99
#      GiB of 16 -- close to the 9.09 GiB prelim high-water mark, and under.
#      -32gb was NOT needed.  RECOMMENDATION for GL, which is far bigger: the
#      top three tiles (8.15-8.99 GiB) are all high-N_fit full-3x3 centers, so
#      memory tracks N_fit, not neighbour count alone; size the queue off the
#      largest expected N_fit rather than off this region's comfortable fit.
#
#      A COINCIDENCE THAT LOOKS LIKE A BUG, CHECKED AND CLEARED: E1180_N-2380
#      and E1180_N-2420 both report N_fit 158 exactly.  They are NOT duplicate
#      output -- the tiles differ in size (4273533 vs 4244558 bytes), as do
#      their prelim inputs (8189115 vs 8184013).  Two adjacent sparse edge
#      tiles genuinely fitting the same number of points.
#
#      28 .h5 and 28 field_sizes reports on S3 under .../IS/matched/.
#      Sizes 4244558 .. 58624350 bytes.


# ===========================================================================
# I7a. [DPS] [DONE 2026-09-16 -- DECIDED, WRITTEN, DEPLOYED AND PROVEN ON DPS
#            BY RE-RUNNING THE TILE THAT FAILED]
#      An uncertainty step with no data must clean up and exit 0.
#      RIDES I7'S REBUILD -- one rebuild and one registration cover both.
# ===========================================================================
# THE SYMPTOM: IS_prelim_E1020_N-2580 was the one failure of the 29 (I3).  Its
# FIT step succeeded -- 753 s, 4.16 GiB, 3 QR iterations, sigma_hat 2.90, tile
# written.  Its UNCERTAINTY step then died in 6 s at 0.23 GiB with
#     READING MASK DATA
#         smooth_fit.py: after masking, no data found
# and the job went permanentFail.  N_fit was 327, against thousands to
# hundreds of thousands on every other tile: this is the far-corner center,
# with enough data to fit and essentially nothing inside the ice mask.
#
# DECIDED 2026-09-16 (Ben): THE FAILED TILE IS NOT RECOVERABLE, and the fix is
# a code change -- if the uncertainty calculation has too few data to complete,
# DELETE THE FIT'S RESULTS AND EXIT 0.  Not a re-run, not a retry: that center
# legitimately has no tile, and the run must be able to SAY so.
#
# WHY IT IS UNRECOVERABLE, verified 2026-09-16: a failed OGC job gets no
# dps_output prefix at all.  get_job_result returns only the triaged_job tree
# (logs and JSON, no .h5 -- listed it), so the fit's 6th-of-an-hour of work is
# gone even though it succeeded.  Logs copied to
# ~/tmp/triaged_job-job-atl1415_tile_solve_1786__on_s3-20260915T182431.921562Z_task-3f76222d-dae9-4c21-8a99-1de1ea4da8ae/
#
# ROOT CAUSE, read from the source 2026-09-16.  THREE FACTS, and the third is
# the one that makes this a bug rather than bad luck:
#   1. LSsurf/LSsurf/smooth_fit.py:513-516 -- on data.size == 0 after masking
#      it PRINTS the message and RETURNS A NORMAL-LOOKING DICT, with 'data'
#      empty and TOC/R/RMS empty.  It does not raise and does not exit.
#   2. ATL1415/ATL11_to_ATL15.py:1334-1350 -- status defaults to 1 and is set
#      to 0 by exactly two branches: the fit branch (len(S['m']) > 0) and the
#      error branch (len(S['E']) > 0).  With no data BOTH are empty, so
#      neither fires, status stays 1, and :1356 sys.exit(1)s.  The "done with"
#      at :1349 prints BEFORE the return, which is why the log shows a tidy
#      "done with ...h5" immediately followed by failure AND NO TRACEBACK.
#      The exit-1 is silent by construction; nothing in it says what is wrong.
#   3. run.sh:425-428 guards only the case where THE FIT WROTE NOTHING:
#         if [ ! -f "${base_directory}/prelim/${tile_name}" ]; then exit 0
#      Here the fit DID write, so the guard passed and the error step ran.
#      The guard's own comment states the intent this defeats -- "Running the
#      error calculation on it would then exit 1 and mark the whole DPS job
#      failed, which at a fan-out of thousands of tiles would bury the real
#      failures."  That is precisely what happened, one tile in 29.
# AND WHY ONLY THE ERROR STEP READS THE MASK: ATL11_to_ATL15.py:582-586 sets
# read_mask_file = calc_error_file when data_file is None, so "READING MASK
# DATA" (:596) happens on the --calc_error_for_xy pass and not on the fit.
#
# THE CHANGE, TWO PLACES:
#   A. ATL11_to_ATL15.py, the status block -- a third branch: an
#      error-calculation run with NO DATA removes args.out_name, removes its
#      field-size report, says so, and sets status = 0.
#      THE TEST COVERS BOTH OF smooth_fit'S NO-DATA EXITS, because per Ben
#      (2026-09-16, below) THE CAUSE DOES NOT MATTER:
#          S.get('data') is None  or  S['data'].size == 0
#      catching smooth_fit.py:485 (`not np.any(valid_data)`, returns
#      data=None) and smooth_fit.py:513 (`data.size == 0` after masking).
#      IT IS STILL NOT "E came back empty".  Treating every empty E as success
#      would swallow a genuine error-propagation failure -- the thing LSsurf
#      49f55db's error handling exists to surface -- and turn it into a silent
#      missing tile, which is worse than the bug being fixed: a failed job is
#      at least visible in the collector.  A solver error, an OOM or anything
#      that raises never reaches this branch and still exits 1.
#      DELETE THE REPORT TOO, not just the tile: the fit step wrote
#      <dir>/field_sizes/<tile>_report.json (:913, :932), fetch_tiles.py pulls
#      prelim/field_sizes/*_report.json, and an orphan report would describe a
#      tile that does not exist.  I5 counts reports against tiles.
#   B. run.sh -- the `s3_tiles put` after the error step is UNCONDITIONAL.
#      With the tile deleted it would fail on a missing file and re-fail the
#      job in a NEW way, so it must be guarded on the file still existing.
#
# WHAT IT MUST NOT DO: make a real failure exit 0.  Only the no-data-after-
# masking case is silenced, and it is silenced LOUDLY -- it prints why, and
# the job shows up in I4 as a 'no tile' row, which is a reported outcome.
#
# ===========================================================================
# THE SCOPE DECISION, 2026-09-16 (Ben) -- THIS CLOSED THE OPEN QUESTION.
# ===========================================================================
# Ben: "Assume that tiles that fail on the uncertainty step are not critical."
# That settles what three options were circling.  The earlier instruction was
# to exit 0 for too-few-data BUT NOT for the coarse-resolution mask edge cases,
# and both causes reach the SAME line with the SAME message, so no test could
# separate them without new machinery.  Declaring the tiles not critical makes
# the separation unnecessary: BOTH causes get the same treatment.
# WHAT THIS BOUGHT: all three options are dropped.  No LSsurf change, no
# pre/post-mask count threshold, no uncoarsened retry pass.  The fix is one
# widened condition in ATL11_to_ATL15.py and one guard in run.sh.
# WHAT IT DOES NOT MEAN: "any uncertainty failure is fine".  The branch is
# still bounded by NO DATA.  Anything that raises -- solver error, OOM -- does
# not reach it and still fails the job loudly.  If those should be silenced too
# that is a separate, bigger decision and it has NOT been made.
#
# HOW IT WAS TESTED, 2026-09-16, and the honest limit.  The failing tile is
# gone and cannot be fetched, so the end-to-end path was NOT exercised.
# WHAT WAS ACTUALLY RUN:
#   1. remove_tile_and_report against a real tile + report on disk: both
#      removed, and a second call on the already-deleted pair does not raise
#      (the FileNotFoundError path), so a partial cleanup cannot fail a job.
#   2. The status elif-chain, evaluated over six cases:
#        err run, :513 no data after masking      -> TILE REMOVED, exit 0
#        err run, :485 no valid data (data=None)  -> TILE REMOVED, exit 0
#        err run, data present but E empty        -> exit 1 (FAILED)  <-- kept
#        err run, normal success                  -> errors saved, exit 0
#        fit run, normal success                  -> fit saved, exit 0
#        fit run, no data                         -> exit 1, as before,
#                                                    guarded by run.sh:425
# STILL UNTESTED, say so plainly: the real DPS round trip, i.e. that a tile
# deleted on the worker makes run.sh skip the upload and the job report a
# success.  That cannot be checked until the rebuild lands.
#
# ---------------------------------------------------------------------------
# TWO TRACE ITEMS RESOLVED 2026-09-16 (read from source + measured locally).
# Both informed the scope decision below; the question they were raised
# against is now CLOSED (see THE SCOPE DECISION).
#
# STATEMENT (traced, smooth_fit.py:495 + fd_grid.py:224-253): on the IS
# uncertainty pass it is validate_by_dz_mask THAT RUNS, NOT setup_mask.  The
# branch is `if args['mask_file'] is not None and grids['dz'].mask_3d is None`.
# IS does pass --mask_file (the RGI .db), so the first half is true -- but the
# error pass reads mask_data out of the prelim tile (:582-596), the tile's dz
# mask is 3-D (VERIFIED: dz/mask is (61,61,32) in every fetched tile, z0/mask
# is (601,601)), and a 3-D mask_data makes fd_grid set mask_3d, so mask_3d is
# NOT None and the branch goes to the else.  Note fd_grid.setup_mask also sets
# self.mask_file=None whenever mask_data is given: on this pass the RGI file is
# read for nothing.  CONSEQUENCE: option (a) below costs a change in
# validate_by_dz_mask (grid_functions.py:344-372), not in setup_mask.
#
# STATEMENT (measured 2026-09-16 over all 28 fetched tiles): COARSENING THE dz
# MASK DOES NOT BY ITSELF EMPTY THE MASK, so it is not sufficient on its own to
# explain data.size == 0.  Re-running fd_grid's own recipe -- interpolate the
# tile's dz mask onto 2x-coarser centers, threshold > 0.5 -- against the 1/4 of
# native cells you would expect: worst loss 31.6% (E1140_N-2500), median ~4%,
# and the two smallest-mask tiles (E1180_N-2380 / N-2420, 192 native cells)
# LOSE NOTHING, gaining 33% over the naive expectation.  NO TILE GOES TO ZERO,
# and a tile with only ~2 surviving mask cells per epoch still completed.
# READING: the scarce quantity is DATA POINTS, not mask cells.  E1020_N-2580
# had N_fit = 327 against thousands-to-hundreds-of-thousands elsewhere, and
# validate_by_dz_mask culls DATA where the interpolated mask is <= 0.5, so a
# handful of points against a few coarse mask cells can plausibly cull to zero
# while the mask itself stays populated.
# HONEST LIMIT: this is measured on the 28 tiles THAT SUCCEEDED; the failing
# tile's mask cannot be measured because the tile is gone.  It shows coarsening
# is not sufficient, NOT that coarsening is irrelevant -- it may still be what
# tips a 327-point tile over.  NOT SETTLED, and deliberately not settled: the
# scope decision below makes the distinction unnecessary.
#
# ===========================================================================
# DEPLOY VERIFIED 2026-09-16.  VERDICT: MATCH, and I7a rode I7's rebuild.
# ===========================================================================
# Ben registered; the image stamp, the live git in the image and the CWL are
# all 6b8a2ca8703be0f20680975d3126b68dc6198684, tree_state=clean,
# algorithm_version=on_s3, maap_py=5.1.0, maap_pgt=set.  Build ran
# 14:44:29Z -> 14:47:32Z, AFTER the 14:41:53Z commit, so it is a real rebuild
# and not an image reused under the tag.  processID=64, modified
# 2026-09-16T14:51:50.  build_id job
# job-atl1415_tile_solve_1786__on_s3-20260916T151658.131427Z, successful.
# origin/on_s3 is level at the same commit, so nothing is past the build.
# This is the MATCH the submission gate wants: registration alone would not
# have been enough (docs/howto_MAAP_ogc.sh, cwlLink is not a deploy).
#
# TWO check_build_id.py BUGS FOUND DOING IT.  One fixed, one only recorded:
#   FIXED: `--expect 6b8a2ca` (a SHORT sha) against the 40-char build stamp
#   reported VERDICT: MISMATCH, whose text says to RE-REGISTER and then to
#   take it to MAAP support -- a false alarm that would have cost a second
#   rebuild of a perfectly good image.  The compare was a bare `!=`.  It now
#   goes through same_commit(), which accepts an abbreviation of >= 7 hex
#   chars as a prefix, the way git does.  Tested over 8 cases including
#   6-char (rejected, too short to trust), a wrong short sha, and a non-hex
#   value.  The origin-is-past-this-build note uses the same compare now.
#   NOT FIXED, and it BLOCKED the re-read: `--job <id>` polls get_job_status,
#   which returns 404 forever for an OGC job id that list_jobs reports as
#   `successful` (tried the id above; 404 every 15 s until killed).  So after
#   a timeout, or to re-score a job with a different --expect, THERE IS NO
#   WORKING RE-READ -- the only paths are a fresh submission or scoring the
#   BUILD_ID line by hand through verdict(), which is what was done here.
#   This is the unfinished half of O5 in docs/howto_MAAP_ogc.sh.  [ADE]
#
# ===========================================================================
# PROVEN ON DPS 2026-09-16.  The smoke tile IS the tile that failed.
# ===========================================================================
# Ben: "submit one tile first".  The tile chosen was E1020_N-2580 itself, not
# a healthy one: re-running it reproduces the exact no-data condition, so it
# exercises the NEW branch instead of only showing nothing regressed.  Every
# parameter was copied from the failed job's own ledger row, so it is a true
# re-run.  region_files/IS_i7a_xy.txt (new) names the single center rather
# than relying on --limit 1, per this file's own rule at I2.
# Job 5cda49fa-41ed-4564-9e23-fba6db0b64cb, ledger IS_i7a_jobs.csv.
#
# STATEMENT, read with collect_jobs.py at 16:49 UTC:
#   IS_prelim_E1020_N-2580  successful  719 s  4.15 GiB  maap-dps-worker-16gb
#     step fit    697 s  4.15 GiB
#     step error    6 s  0.23 GiB
#   N_ATL11 195575, N_AT 195511, N_XO 64, N_fit 327, 3 iterations
#   commit 6b8a2ca, built 2026-09-16T14:47:32Z
# THE SAME JOB WAS permanentFail ON 2026-09-15.  It is now successful, and
# THE UNCERTAINTY STEP STILL DIES IN 6 s AT 0.23 GiB -- identical to the
# failing run, so it reached the SAME no-data path and simply exits 0 now.
# N_fit is 327, the same 327 as the failure: the condition reproduced exactly.
#
# AND IT WROTE NO TILE, which is the unusual half of the result.  Both
#   .../rel006/north/IS/prelim/E1020_N-2580.h5
#   .../rel006/north/IS/prelim/field_sizes/E1020_N-2580_report.json
# were VERIFIED ABSENT ON S3 BEFORE SUBMITTING and are STILL ABSENT after a
# successful job.  So the delete fired, run.sh's new -f guard skipped the
# upload rather than failing on the missing file, and no orphan report was
# left behind.  That is the whole DPS round trip the local tests could not
# reach, and it is now closed.
#
# THE MATCHED COUNT IS STILL 28.  This tile deliberately produces no prelim
# tile, so it does not join the matched list and I6/I8 are unaffected.
# COST: 719 s of worker time to prove it, on the 16gb queue with 4.15 GiB peak.
#
# THEN IT NEEDS A REBUILD AND A REGISTRATION, like any worker-side change, and
# check_build_id must MATCH the new commit before anything is submitted --
# a cwlLink commit does not prove a deploy (docs/howto_MAAP_ogc.sh).  Fold it
# into I7's rebuild rather than spending a second one.
# THE 28 GOOD TILES ARE UNAFFECTED and are not re-run.


# ===========================================================================
# I8. [ADE] [DONE 2026-09-16 -- 28 fetched, verified against S3]
# ===========================================================================
scripts/maap/fetch_tiles.py ~/ATL14_processing/maap_ledgers/IS_matched_jobs.csv \
    $region_dir --step matched
# The fetcher already tries the top-of-prefix layout matched uses.  If I7 runs
# in the ADE this step disappears.
#
# DONE 2026-09-16 into /home/jovyan/ATL14_processing/rel006/north/IS/matched/:
# 28 fetched, 0.48 GiB transferred, 490 MB on disk, plus 28 field_sizes
# reports.  Every tile came down at the top-of-prefix layout, as the comment
# above predicts -- no fallback path was needed.
# VERIFIED, not just counted: the 28 local sizes were diffed against the S3
# listing and are IDENTICAL, tile for tile.  All 28 reports read
# dz/dz [61,61,32] and dz/sigma_dz null -- one shape, no anomalies, and the
# null sigma is correct (see I7, and Ben 2026-09-16).
# THIS COMPLETED THE IS RUN AS FIRST SCOPED (prelim + matched).  I9 is next.


# ===========================================================================
# I9. [ADE] [MOSAIC DONE 2026-09-16, netCDF DONE 2026-09-17 with invalid lineage;
#     browse open.  IN SCOPE -- Ben 2026-09-16; was OUT OF SCOPE 2026-09-12]  Mosaic, netCDF.
# ===========================================================================
# TENTATIVE.  Written 2026-09-16 before any of it ran; revise as steps land.
# Scope: MOSAIC (I9a-f, done) and netCDF (I9g1-5 and 7 done; I9g6 open).  Browse (I9h) is listed, not planned.
# Runs in the ADE, not DPS: no submission, no registration, nothing polled.
# Ben restarted the instance for memory: 16 cores, 124 GiB (101 available).
region_dir=/home/jovyan/ATL14_processing/rel006/north/IS
mosaic_run=/home/jovyan/ATL14_processing/runs/IS_mosaic

# I9a. [OK 2026-09-16]  Tiles on disk -- collected at I8, rechecked after the
#      restart: 28 prelim + 28 matched, the SAME 28 names, 652 + 490 MB.
#      Matched tiles give the values, prelim tiles give every sigma_* (the
#      queue globs 'prelim/*.h5' for those), so both directories are inputs.

# I9b. [DONE 2026-09-16 -- fixed in 2ed36c0, tests/test_lazy_imports.py]
#      A REGRESSION blocked I9c and I9g.
#      STATEMENT: make_mosaic_jobs.py dies writing slurm_run.sh with
#        TypeError: 'module' object is not callable
#      at ATL1415.make_slurm_file(...).  Cause: 0c6ea35 (lazy imports, on_s3
#      only).  The old __init__ star-imported each submodule, which bound the
#      FUNCTION over the same-named MODULE; the lazy __getattr__ returns the
#      module.  Seven names are affected -- make_slurm_file,
#      make_nc_projection_variable, make_tile_stats_group, read_ATL11,
#      assign_firn_variable, SMB_corr_from_grid, ATL11_to_ATL15.
#      Callers that break: make_mosaic_jobs.py, make_200km_tiles.py,
#      make_200km_to_mosaic_jobs.py, setup_slurm_run.py, setup_ATL1415_run.py
#      (all ATL1415.make_slurm_file), and ATL14_write2nc.py / ATL15_write2nc.py
#      (from ATL1415 import make_nc_projection_variable, make_tile_stats_group).
#      The DPS solve is unaffected: it imports submodules by full path.
#      FIX: resolve those seven names to the function, as before 0c6ea35, with
#      a test in tests/ so it cannot regress silently again.

# I9c. [DONE 2026-09-16 -- 41 tasks; I9e's done/ holds all 41]  Build the queue.
mkdir -p $(dirname $mosaic_run) && cd $(dirname $mosaic_run)
make_mosaic_jobs.py -b $region_dir -rr IS -t 2018.75,2026.5 -e ATL14 \
    --run_name IS_mosaic @/home/jovyan/git_repos/ATL1415/default_args/quarterly.txt
#      -e ATL14: the default, IS2, is discover's env and does not exist here.
#      Dry run in scratch (before I9b) gave 41 tasks: z0; for each of 40/20/10 km
#      the avg_dz grid + 9 avg_dzdt lags; dz; 9 dzdt lags.  Every task writes a
#      DIFFERENT output file, so tasks are safe to run in parallel.
#      Lags inferred from -t and -g = 1,2,4,8,12,16,20,24,28, identical to
#      --dzdt_lags in input_args_IS.txt.  No bounds.txt in $region_dir, so no
#      crop -- expected, cropping moved to the to_nc step in bacb2ef.

# I9d. [DONE 2026-09-16]  Smoke one task: task_1 (z0, the 100 m grid, the biggest).
cd $mosaic_run
SLURM_ARRAY_TASK_ID=1 /home/jovyan/git_repos/ATL1415/scripts/run_with_rusage.py \
    IS_mosaic_task_1 bash slurm_run.sh
#      slurm_run.sh IS plain bash -- the #SBATCH lines are comments -- so run
#      locally it does the same queue -> running -> done moves and writes
#      logs/ and error_logs/ exactly as on discover.  THIS REPLACES
#      run_queue_local.sh (arctic step 10): no new runner is needed.
#      Record wall time and peak RSS; they set -P for I9e.
#      DONE 2026-09-16: exit 0, 41.7 s, peak RSS 0.25 GiB.  z0.h5 is 105 MB,
#      grid (3001, 4201) at 100 m, x 990..1410 km, y -2650..-2350 km.  All 7
#      fields present (z0, misfit_rms, misfit_scaled_rms, mask, cell_area,
#      count, sigma_z0); z0/mask/cell_area/sigma_z0 finite over 49.6% of the
#      box, count and misfit_* over 0.9% (data cells only).
#      MEMORY IS NOT THE CONSTRAINT for IS -- the restart was not needed for it.

# I9e. [DONE 2026-09-16]  Run the other 40.
seq 2 41 | xargs -P 12 -I{} env SLURM_ARRAY_TASK_ID={} bash slurm_run.sh
#      -P 12 of 16 cores: memory is no limit at 0.25 GiB per task (I9d).
#      20:16:55Z -> 20:17:54Z, under a minute for all 40.  40/40 exit code 0,
#      error_logs/ empty, done/ holds 41.  41 files, 178 MiB, in $region_dir.

# I9f. [DONE 2026-09-16]  Verify.  EXIT CODES ARE NOT ENOUGH:
#      - make_mosaic.py returns 0 when pc.grid.mosaic fails (it prints the
#        message only under -v, and the queue does not pass -v);
#      - a task file has no `set -e`, so only its LAST line's status counts.
#      So: error_logs/ empty AND done/ holds 41 AND every -O file named in the
#      queue exists and holds every -F field under its --in_group.
#      Also, sigma fields must be present (from prelim) -- a missing sigma_*
#      is the likeliest silent failure, because its line is never the last.
#      RESULT: 41 -O files, 131 -F fields, every one present, none all-NaN,
#      PROBLEMS: 0.  Beyond the checker: within each file every field --
#      matched values and prelim sigmas alike -- has ONE shape, and the time
#      axes are consistent: dz (301,421,32) at 1 km; 10/20/40 km grids
#      (30,42) (14,20) (7,10); dzdt_lagL has 32-L epochs.
check_mosaic_outputs.py $mosaic_run --values
#      The checker parses make_mosaic.py lines out of the task files.  It is
#      now ATL1415/scripts/check_mosaic_outputs.py, an installed entry point
#      (re-run `pip install -e .` for the command to appear), and is in every
#      howto's mosaic step, discover and MAAP.  METADATA ONLY BY DEFAULT
#      (1.5 s on IS); --values reads the data to flag all-NaN fields (7.4 s on
#      IS, minutes for AA) -- Ben 2026-09-16.  The result above is --values.

# ---------------------------------------------------------------------------
# I9g. [ADE] [I9g1-5, I9g7 DONE 2026-09-17; I9g6 OPEN]  netCDF.
# ---------------------------------------------------------------------------
# Ben 2026-09-16: "plan the netCDF step".  Planned from a PROBE, not from
# reading alone: both writers were run on the real IS mosaics with -b pointing
# at a scratch directory of symlinks, so nothing was written into $region_dir.
# The writers read, besides the mosaics: the prelim tiles (tile_stats and
# lineage), the package's attrs CSVs and metadata templates, and
# region_extent_polygons.json.  NOT the DEM, masks, geoid or previous product.

# I9g1. [OK -- PROBED 2026-09-16]  Everything except lineage already works on IS.
#      STATEMENT.  With set_lineage() stubbed out by a probe script (scratch
#      only), both writers exit 0:
#        ATL14_IS_0331_100m_006_01.nc      h (3001, 4201)            12 s
#        ATL15_IS_0331_3mo_{1,10,20,40}km_006_01.nc                  18 s total
#          delta_h (30, 301, 421) at 1 km; 30 epochs = t_crop 2019..2026.25
#          out of the mosaic's 32; groups delta_h + 9 dhdt_NNNmo.
#      tile_stats reads data/three_sigma_edit, RMS/*, E_RMS/*, bias/* from the
#      28 prelim tiles -- all present.  The ATL15 attrs CSV names exactly the
#      files I9e wrote (dz.h5, dz_{res}km.h5, dzdt_lag{lag}.h5,
#      dzdt_{res}km_lag{lag}.h5).  I9b's fix is what lets both scripts import.
#      Memory and time are no concern for IS.

# I9g2. [DONE 2026-09-17 -- TEMPORARY FIX; lineage INVALID]  THE BLOCKER.
#      STATEMENT, from the unstubbed probe:
#        FileNotFoundError ... 's3://maap-ops-workspace/ben_smith/ATL11_053503_0331_007_04.h5'
#      set_lineage() collects the ATL11 file names from each prelim tile's
#      meta.input_files (basenames only -- ATL11_to_ATL15.py:878), then
#      attributes_for_ATL11_file() OPENS EVERY ONE with h5py by local path.
#      IS needs 79 files: 69 ATL11 and 10 ATL11XO (5 tiles x cycles 1, 2).
#      STATEMENT, checked 2026-09-17: the tiles hold the NAMES and nothing
#      else.  No tile has /meta/sensors (the sensor_N attributes belong to
#      the old multi-sensor code: LSsurf two_mission_dhdt.py, main's
#      reread_data_from_fits.py); data/sensor is 0 everywhere.  Matched tiles
#      carry input_files = '' -- the names are in the PRELIM tiles only.
#      From the name: shortName, start_rgt/start_region/cycles (along-track),
#      cycle (XO), release, version.  Only by opening the granule: uuid,
#      start/end_geoseg, start/end_orbit (along-track), start/end_rgt (XO).
#
#      DECISION, Ben 2026-09-17 (closes QI9 and QI10):
#        - LONG TERM: the prelim step records the granule attributes in the
#          tile metadata at solve time, when it has every granule open.  The
#          netCDF step reads them from the tiles and never opens ATL11.  An
#          attribute the tiles do not carry is marked INVALID in the output.
#          PLANNED 2026-09-17: docs/plan_lineage_at_solve_time.sh.
#          Lands before the IS re-run -- which is ALSO blocked: CMR now lists
#          ATL11 0332_007_05 only, not the 0331_007_04 IS ran on (that plan,
#          top, and its QL1).
#        - TODAY: the tiles carry none of them, so those attributes are
#          invalid.  IS will be RE-RUN COMPLETELY once more issues are fixed,
#          so these netCDFs are for getting the step working, not for release.
#
#      TEMPORARY FIX:
#        a. attributes_for_ATL11_file() stops opening files: it fills what the
#           name gives and leaves every file-only attribute at 'NOT_SET', the
#           marker the function already initializes with.  The local-path
#           and XO-schema lookups go, and with them the dead
#           start/end_delta_time reads.
#        b. a name matching neither pattern raises, naming the file and the
#           tile, instead of an AttributeError on None; an empty input_files
#           (a matched tile) contributes no names rather than a '' name.
#        c. set_lineage() prints ONE warning line saying how many files have
#           invalid lineage and which attributes -- visible, not silent.
#        d. unchanged: end_rgt = start_rgt, including for XO (see below).
#      TEST: tests/test_lineage.py -- along-track and XO names give the name
#      attributes and 'NOT_SET' for the rest without touching the filesystem;
#      a bad name raises; set_lineage on synthetic tiles dedupes, skips '',
#      and warns.  Then I9g3-4 for real.
#      NO DPS REBUILD: write2nc runs in the ADE, and the tiles are unchanged.
#      pre_rel006 (discover) keeps the file-opening code; on_s3 diverges here.
#      RESULT: tests/test_lineage.py adds 5; the suite is 59 passed, 2 skipped
#      (conda run -n ATL14 python -m pytest tests -- pytest had to be
#      reinstalled in ATL14 after the restart).
#
#      NOTED, NOT FIXED: attributes_for_ATL11_file() sets end_rgt = start_rgt
#      AFTER the file read, overwriting the end_rgt it read from each ATL11XO
#      granule.  Moot while XO rgts are invalid; fix it with the long-term
#      change.

# I9g3. [DONE 2026-09-17]  ATL14.
ATL14_write2nc.py @$region_dir/input_args_IS.txt
#      Output $region_dir/ATL14_IS_0331_100m_006_01.nc.  -b, cycles, release,
#      version, t_crop, region and the earthaccess flags all come from the args.
#      RESULT: exit 0, 10 s, 9895823 bytes.  Run as
#        cd ~/ATL14_processing/runs/IS_nc; conda run -n ATL14 ATL14_write2nc.py @... > ATL14.log 2>&1
#      Log has the two set_lineage INVALID warnings (69 along-track, 10 xo).

# I9g4. [DONE 2026-09-17]  ATL15 -- one call writes all four resolutions.
ATL15_write2nc.py @$region_dir/input_args_IS.txt
#      --avg_scales=40000,20000,10000 in the args makes it loop over
#      [None, 40000, 20000, 10000]: ATL15_IS_0331_3mo_{1,40,20,10}km_006_01.nc.
#      No slurm runner: the two commands run directly (30 s together, I9g1).
#      RESULT: exit 0, 15 s, 4 files (1km 17282097, 10km 879394, 20km 696733,
#      40km 636317 bytes); ATL15.log has the INVALID warning pair once per file.

# I9g5. [DONE 2026-09-17, by a scratch script]  Verify.
#      - 5 files, exit 0 each;
#      - METADATA/Lineage/ATL11 lists 79 files, 69 ATL11 + 10 ATL11XO; the
#        name attributes are set on every row, and -- TODAY, per I9g2 --
#        uuid, geoseg, orbit and XO rgt are 'NOT_SET' on every row;
#      - h, delta_h and each dhdt group have the shapes of I9g1, and time has
#        30 epochs inside t_crop;
#      - finite cells of h/delta_h match the finite cells of the mosaic fields
#        they came from (ice_area masks them, so this checks the masking);
#      - the files open with netCDF4 and with GDAL (NETCDF:"file":h).
#      RESULT, all five files:
#      - lineage 79 = 69 ATL11 + 10 ATL11XO; shortName, cycles, release,
#        version set on every row; uuid, start/end_geoseg, start_orbit all
#        'NOT_SET'; start_rgt has 70 values = 69 rgts + 'NOT_SET' for XO.
#      - ATL14 h (3001, 4201); ATL15 delta_h (30, 301|30|14|7, 421|42|20|10),
#        time 30 epochs, 365.25..3013.31 days = 2019.0..2026.25; groups
#        delta_h + dhdt_{003,006,012,024,036,048,060,072,084}mo with 29..2 epochs.
#      - ATL14 h vs z0.h5 and ATL15 1 km delta_h vs dz.h5: same x/y; finite in
#        the product == finite in the mosaic AND ice_area > 0, EXACTLY (no
#        product-only cells); max |difference| 1.2e-4 m (h), 3.1e-5 m (delta_h).
#      - netCDF4 opens all.  GDAL opens h (4201x3001) and delta_h (421x301x30)
#        from the NOTEBOOK env only: ATL14's GDAL lacks the netCDF and HDF5
#        plugins (libgdal-netcdf, libgdal-hdf5 not installed).  Not a file problem.

# I9g6. [RECOMMENDATION, NOT PLANNED IN DETAIL]  Compare against rel005.
#      The run used --previous_product=005_0329.  Differencing h against
#      ATL14_IS_0329_100m_005_* (median and spread over common ice cells) is
#      the cheapest science-level check that nothing is shifted or flipped --
#      the I9f and I9g5 checks are all structural.

# I9g7. [DONE 2026-09-17]  Howtos.  Once I9g3-5 pass, record in
#      docs/howto_MAAP_arctic.sh step 10 (it already names both commands)
#      that they run directly in the ADE, and which lineage flags they need.
#      RESULT: arctic step 10 now runs slurm_run.sh directly (no
#      run_queue_local.sh) and the writers with NO lineage flags, noting the
#      invalid lineage.  Steps 6-9 were rewritten to the commands IS ran;
#      step 5 stays NEEDS CODE.  howto_MAAP_GL.sh and howto_MAAP_AA.sh still
#      cite run_queue_local.sh -- NOT updated, see the commit.
#
# QI9.  CLOSED -- Ben 2026-09-17: lineage at solve time, see I9g2 DECISION.
# QI10. CLOSED -- same decision covers GL/AA.
# I9h. [NOT STARTED, NOT PLANNED]  Browse.


# ===========================================================================
# I10. [ADE] [SUGGESTION, NO SOFTWARE]  Annotate the run's build history.
# ===========================================================================
# howto_MAAP_ogc O12b, unchanged.


# ###########################################################################
# QUESTIONS -- undocumented decisions.  I did not guess at any of these.
# ###########################################################################
#
# QI1. ANSWERED 2026-09-12 (Ben): FREEZE THE LIST in region_files/, from the
#      40 km tif -- the 29 centers already read off the bucket.  IS is then not
#      gated on make_ATL1415_queue.py.  Its three cloud bugs (I1) stay open and
#      must be fixed before GL, which has no frozen list and cannot dodge them.
#
# QI2. ANSWERED 2026-09-12 (Ben): generalize to submit_MAAP_jobs.py now.  I2.
#
# QI3. ANSWERED 2026-09-12 (Ben): build the tile fetcher, ADE-side.
#      Q9's answer had approved a deterministic prefix "written by the job
#      itself", but not the route to it, and the I4 finding changed the
#      balance: because a job's output/prelim/ survives the upload intact, an
#      ADE-side fetcher driven by the ledger works TODAY -- no rebuild, no new
#      CWL input, nothing re-registered.
#      DONE: scripts/maap/fetch_tiles.py, and collect_AA_queue.py renamed to
#      scripts/maap/collect_jobs.py so the two halves of "collect" have
#      different names -- the collector reads the jobs' logs, the fetcher
#      moves the .h5 tiles.
#      STILL OPEN, and NOT answered by this: the bucket copy stays
#      non-addressable, so this does not by itself unblock matched-on-DPS (a
#      matched job still has no key to fetch its 8 neighbours from), and Q9's
#      deterministic prefix is still wanted before production.  See QI4, QI5.
#
# QI4. ANSWERED 2026-09-12 (Ben), both halves:
#      WHO WRITES THE TREE: run.sh, at solve time -- Q9 taken literally.  The
#      ADE-sync alternative was declined: every tile would cross the wire
#      twice, and the tree would only exist after a manual step.  Since the
#      rebuild for prelim_prefix is happening anyway, the marginal cost is
#      small.  The tile fetcher is NOT made redundant -- the mosaic step still
#      needs a local copy -- but it stops being the only route.
#      THE LAYOUT: the hemi suffix, as Q17.
#        quarterly  s3://maap-ops-workspace/ben_smith/ATL14_processing/
#                     rel006/north/IS/{prelim,matched}/E<x>_N<y>.h5
#        monthly    .../rel006/north_monthly/IS/{prelim,matched}/E<x>_N<y>.h5
#      Mirrors the local tree and discover exactly, so run_arctic_*.sh's
#      hemi_suffix grep keeps working unchanged.
#
# QI5. ANSWERED 2026-09-12 (Ben), both mechanics:
#      (a) A FIFTH CWL INPUT, prelim_prefix -- not derived inside run.sh.  The
#          submitter passes it per job, so the path is visible in every job
#          record and every ledger row rather than hidden in a convention.
#          It costs a re-registration, which this rebuild needs anyway.
#      (b) COPY WHAT EXISTS AND PROCEED.  A missing neighbour is the normal
#          case on IS -- a small coastal region where most of the 29 tiles
#          have fewer than 8 neighbours, and where a tile with too little data
#          legitimately writes nothing at all.  Failing on it would fail most
#          of the region.  run.sh stays STRICT ON THE TILE'S OWN prelim file,
#          which it already is (run.sh:403).
#
#      ONE DETAIL THIS LEAVES, and it is mine to propose rather than yours to
#      decide unless you disagree: QI4 makes a PRELIM job write to the bucket
#      too, so it also needs to be told a prefix.
#      RECOMMENDATION: ONE input serves both -- for a prelim job it is where
#      to WRITE the tile, for a matched job it is where to READ the
#      neighbours' -- and it should therefore be named tile_prefix, not
#      prelim_prefix.  Two inputs that are always the same string invite them
#      being different by accident.  Say if you would rather have two.
#
# QI6. ANSWERED 2026-09-12 (Ben): DPS, not the ADE.  I7, and it is the item
#      that costs a rebuild and a registration.
#
# QI7. ANSWERED 2026-09-12 (Ben): maap-dps-worker-16gb, not the -32gb used so
#      far.  The ram_min=16 scheduling risk is noted at I2, with the one-tile
#      test that settles it cheaply.
#
# QI8. ANSWERED 2026-09-12 (Ben): MATCHED TILES ON DISK.  Mosaic, netCDF and
#      browse are the next pass, so run_queue_local.sh is NOT on this critical
#      path.
#
#
# ###########################################################################
# STILL OPEN -- both newly unblocked by the QI6 answer
# ###########################################################################
# QI4 and QI5 below are what matched-on-DPS needs before any code is written.
