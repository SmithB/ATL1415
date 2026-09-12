# plan_IS_run.sh -- THE ICELAND (IS) RUN: prelim and matched, end to end
#
# ############################################################################
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
# I0. [ADE] [READY]  Confirm the image before spending worker hours.
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


# I2. [DPS] [READY -- scripts/maap/submit_MAAP_jobs.py, new 2026-09-12]
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
# the fallback and nothing else changes.  THEN the other 28:
scripts/maap/submit_MAAP_jobs.py --xy_file region_files/IS_prelim_xy.txt \
    --step prelim --args_url $s3_run/input_args_IS.txt \
    --queue maap-dps-worker-16gb \
    --ledger ~/ATL14_processing/maap_ledgers/IS_prelim_jobs.csv
# (the smoke tile is submitted again as part of the 29; dedup=False, so it
# really re-runs.  Simpler than excising one line, and it costs ~1 h.)
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
#
# TESTED 2026-09-12, all against the real deployed process (processID 64):
# the 29-center dry-run, --limit 1, and every refusal above.  Not yet run for
# real -- no IS job has been submitted.


# I3. [ADE] [READY]  Watch the 29 jobs.   THE COLLECTOR.
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


# ===========================================================================
# I4. [ADE] [READY -- scripts/maap/fetch_tiles.py, new 2026-09-12]
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


# I5. [ADE] [READY]  Look at the tile sizes.   (arctic step 8's ADE half)
# ===========================================================================
# ATL11_to_ATL15 writes the field-size report itself; the JSON above is it.
# Inspect with check_tiles.ipynb.


# ===========================================================================
# I6. [ADE] [READY]  The matched tile list.
# ===========================================================================
# The same 29 centers -- region_files/IS_prelim_xy.txt serves both steps.
# STATEMENT: expect the matched list to be SHORTER than 29, and that is not a
# failure.  ATL11_to_ATL15 returns 0 without writing when a tile has too
# little data (run.sh:379 exits 0 on that path), so some centers will have no
# prelim tile to match.  Drop those rows from the matched ledger after I4
# reports them as 'no tile'.


# I7. [DPS] [CODE WRITTEN 2026-09-12; BLOCKED ON A REBUILD + REGISTRATION]
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
#   algorithm_config.yml  a fifth input, tile_prefix, default "" (QI5a).
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
#
# ONE THING TO WATCH AT REGISTRATION: default: "" on a string input is not
# something this config has used before, and whether /api/build accepts an
# empty default is UNVERIFIED.  If it is rejected, the fallbacks are to drop
# the default (making tile_prefix required, which run.sh already tolerates
# since it treats an empty value as absent) or to default it to "-" and have
# run.sh read that as empty.  The build log will say.
#
# THEN, in order:
#   1. commit and push -- register_algorithm.py refuses on unpushed work
#   2. Ben registers; wait for the build and the deploy
#   3. scripts/maap/check_build_id.py -- must say MATCH at the new commit
#   4. ONE matched job before the other 28


# I8. [ADE] [READY]  Bring the matched tiles down.
# ===========================================================================
scripts/maap/fetch_tiles.py ~/ATL14_processing/maap_ledgers/IS_matched_jobs.csv \
    $region_dir --step matched
# The fetcher already tries the top-of-prefix layout matched uses.  If I7 runs
# in the ADE this step disappears.


# ===========================================================================
# I9. [ADE] [OUT OF SCOPE -- CONFIRMED 2026-09-12 (Ben)]  Mosaic, netCDF, browse.
# ===========================================================================
# You asked for prelim and matched.  Recording the dependency only: arctic
# step 10 is still [NEEDS CODE: run_queue_local.sh] -- make_mosaic_jobs.py
# emits queue/task_N plus a slurm_run.sh and there is no sbatch in the ADE.
# -> QI8 asks whether "IS end to end" (Q21) stops at matched or goes through
# netCDF, because that decides whether run_queue_local.sh is on this critical
# path or the next one.


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
