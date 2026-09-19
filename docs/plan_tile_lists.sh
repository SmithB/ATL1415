#! /usr/bin/env bash
# ===========================================================================
# PLAN: per-region TILE LISTS drive MAAP submissions, and a no-data prelim fit
# exits cleanly.  From docs/plan_monthly_on_maap.sh AM5, AM7, AM8.
# Written 2026-09-18, before the code.  STATUS 2026-09-18: TL0-TL5 DONE in
# the ADE and tested; TL6 (Ben registers) is next; TL7 (DPS proof) waits on it.
# Every step carries its own status tag.
# ===========================================================================
# WHAT BEN DECIDED (plan_monthly_on_maap.sh, in his words):
#   AM5: "I have pushed lists of tiles for each region into
#        ATL1415/resources/{region} in the on_s3 branch.  If any tiles fail
#        for lack of data on the prelim step, they should be deleted from this
#        list."
#   AM7: "B. Both.  The solver change AND the lists."
#   AM8: "A. Yes, for prelim AND matched"
#
# Provenance per claim: STATEMENT = verified 2026-09-18, with how;
# DECIDED = Ben said so; READING = my interpretation of what he said,
# overridable; RECOMMENDATION = mine; QUESTION = open.
# Tags: [ADE] / [DPS] / [BEN]; [NOT STARTED] [NEEDS CODE: x] [DONE].
#
#
# ===========================================================================
# BACKGROUND.  STATEMENT, 2026-09-18, by reading the code named.
# ===========================================================================
# THE LISTS (738bbd2): ATL1415/resources/<region>/40km_tile_list.txt for AA,
#   CN, CS, GL, IS, RA; AA also has 200km_tile_list.txt; there is NO SV list.
#   One tile file name per line, E<x km>_N<y km>.h5.  The IS list is 28
#   names, set-equal to the 28 prelim tiles that exist, E1020_N-2580 absent.
# THE FAILURE THIS FIXES: the monthly E1020_N-2580 prelim job (9266c3d7).
#   smooth_fit.py:485 found no valid data and returned normally with data=None;
#   ATL11_to_ATL15.main() saves nothing and falls through with status 1,
#   because its no-data branch (I7a, plan_IS_run.sh) fires only on the error
#   step (calc_error_file set).  run.sh runs under `set -euo pipefail`
#   (run.sh:39), so that exit 1 ends the job before run.sh:424 -- which
#   ALREADY handles "the fit wrote no tile" by skipping the error step and
#   exiting 0.
#   THEREFORE THE SOLVER IS THE ONLY CHANGE NEEDED; run.sh is untouched.
# SUBMISSIONS TODAY: submit_MAAP_jobs.py --xy_file reads "x0 y0" per line
#   (read_centers, :85).  The IS fan-outs used region_files/IS_prelim_xy.txt
#   (29, from the 40 km mask, E1020 included) and matched lists built from the
#   tiles that exist.  On discover, make_ATL1415_queue.py --tile_list_file
#   filters ONLY the matched step (:272-274).
#
#
# ===========================================================================
# TL0. [DONE 2026-09-18]  Pull the lists.
# ===========================================================================
# 738bbd2 fast-forwarded onto on_s3 (1b5fbaa's parent chain intact).  Its
# .gitattributes commit d5ef8f5 changes nothing net: the file is identical on
# both sides.
#
#
# ===========================================================================
# TL1. [ADE] [DONE 2026-09-18 -- tested locally; DPS proof is TL7]  A no-data prelim FIT exits 0.
# ===========================================================================
# A fourth branch in main()'s status chain, beside the I7a one:
#   args.prelim and calc_error_file is None and
#   (S.get('data') is None or S['data'].size == 0)
#     -> remove_tile_and_report(out_name, 'no data for the prelim fit');
#        status 0.
# Both smooth_fit no-data exits land there, as in I7a (:485 data None, :513
# data empty after masking).  remove_tile_and_report is idempotent -- the fit
# writes nothing on this path, so it normally finds nothing to remove; it is
# called anyway so a stale tile or report can never survive a no-data exit.
# READING: PRELIM ONLY.  AM5 says "on the prelim step", and a MATCHED fit
#   with no data is unexpected -- its own prelim tile had data -- so it keeps
#   exiting 1, loudly.
# BOUNDARY UNCHANGED: anything that raises -- solver error, OOM, a missing
#   input -- never reaches the branch and still exits 1 (Ben, 2026-09-16).
# NO run.sh CHANGE: run.sh:424 then skips the error step and exits 0.
# DONE: the branch sits after I7a's in main(); remove_tile_and_report's
#   docstring now names both callers.  tests/test_no_data_exit.py, 10 cases.
#   MUTATION-CHECKED: against the solver WITHOUT the branch, exactly the 3
#   new-branch cases fail (exit 1 -- E1020's monthly behaviour) and the 7
#   boundary and unchanged-path cases pass; with it, all 10 pass.
#
#
# ===========================================================================
# TL2. [ADE] [DONE 2026-09-18 -- tested; dry-run verified on real data]  --tile_list.
# ===========================================================================
# --tile_list <file>, exactly one of it and --xy_file required.  Each line
# E<x>_N<y>.h5 -> (x*1000, y*1000) meters.  Same rule as read_centers: every
# non-blank line must parse, and one that does not is an error, not a skip.
# --xy_file stays, for the single-center smoke and retry files.
# MATCHED PRE-FLIGHT, with --step matched and --tile_prefix: every listed
#   center's prelim tile must exist at <tile_prefix>/prelim/<name>.  If any
#   is missing, REFUSE -- exit 2, name them, submit nothing.  A missing tile
#   is either a no-data center not yet pruned (run TL3) or a failure under
#   investigation; either way submitting its matched job only fails on DPS.
#   RECOMMENDATION, from the standing rule to stop and name the problem
#   rather than work around it.
# DONE: read_tile_list, s3_names (one listing per region; a FAILED listing
#   exits 2 rather than returning an empty set, which would read as "every
#   tile missing"), missing_prelim_tiles.  --tile_list and --xy_file are a
#   required mutually exclusive pair.  tests/test_submit_tile_list.py, 17.
# VERIFIED LIVE, --dry-run, nothing submitted, no ledger written:
#   - matched from ATL1415/resources/IS/40km_tile_list.txt against the
#     monthly prefix: 28 centers, pre-flight passes;
#   - matched for E1020 (region_files/IS_i7a_xy.txt): REFUSED, exit 2,
#     naming E1020_N-2580.h5;
#   - prelim from AA/40km_tile_list.txt: REFUSED, exit 2, at line 8945,
#     'field_sizes' -- see OPEN.
#   The real IS list parses to exactly the 28 centers of
#   region_files/IS_0332_monthly_matched_xy.txt (a test asserts it).
#
#
# ===========================================================================
# TL3. [ADE] [DONE 2026-09-18 -- tested; dry run matches both IS ledgers]  Prune.
# ===========================================================================
scripts/maap/prune_tile_list.py <prelim_ledger> ATL1415/resources/<region>/40km_tile_list.txt [--write]
# For each PRELIM row of the ledger:
#   successful, and no tile at <tile_prefix>/prelim/<name>  -> PRUNE
#   failed                                                   -> INVESTIGATE, kept
#   still running/accepted                                    -> REFUSE, exit 2
# Dry run by default: prints what it would remove.  --write rewrites the list
# in place, same format and order, and the change is committed for Ben to
# see, never edited by hand.
# WHY "successful, no tile" IS THE WHOLE RULE: after TL1, both no-data exits
#   (fit and error step) are successful jobs that leave no tile, and a failed
#   job again means a real fault.  So the script never reads logs.
#   Consequence: a no-data job from a build BEFORE TL1 (like 9266c3d7) shows
#   as failed and is NOT pruned -- correct, since it is the old build's
#   verdict, and the IS list is already pruned by hand.
# READING, recorded in plan_monthly_on_maap.sh: one list per region serves
#   quarterly and monthly -- a center with no quarterly tile has no monthly
#   reference coverage, so quarterly prunes it first.  The script prunes on
#   any prelim ledger, as AM5 says, and prints which ledger caused each
#   removal.
# DONE: scripts/maap/prune_tile_list.py; the rule is classify(), a pure
#   function.  Also refuses a non-prelim ledger, a list with non-tile lines,
#   and rows with no tile_prefix.  Exit 0 clean / 1 investigate / 2 refused.
#   tests/test_prune_tile_list.py, 10.
# VERIFIED LIVE, dry run (list unchanged), against the IS list:
#   - IS_0332_prelim_jobs.csv: 28 kept; E1020_N-2580 "no data, already
#     absent" -- its quarterly job SUCCEEDED with no tile (I7a).  Exit 0.
#   - IS_0332_monthly_prelim_jobs.csv: 28 kept; E1020_N-2580 INVESTIGATE --
#     job 9266c3d7 FAILED on the pre-TL1 build, so it is not pruned.  Exit 1.
#   Both are the answers the rule predicts.  --write not yet used on a real
#   list: the IS list needs no change.
#
#
# ===========================================================================
# TL4. [ADE] [DONE 2026-09-18 -- suite 156 passed, 2 skipped (was 119)]  Tests.
# ===========================================================================
#   - TL1: main()'s status chain with the fit monkeypatched -- no-data prelim
#     fit -> 0 and no tile; the same for matched -> 1; I7a's error-step branch
#     unchanged; a normal fit -> 0 with a tile.
#   - TL2: tile-name parsing (negative coordinates, a bad line exits 2),
#     --tile_list and --xy_file mutually exclusive, the matched pre-flight
#     refusal -- all without MAAP or S3.
#   - TL3: the prune rule over successful/failed/running rows, dry run vs
#     --write, order and format preserved.
#   Full suite green before anything is pushed.
#
#
# ===========================================================================
# TL5. [ADE] [DONE 2026-09-18]  Docs.
# ===========================================================================
#   - howto_MAAP_arctic.sh steps 6, 9 and 10b: --tile_list in place of the
#     region_files center lists; prune after each prelim fan-out.
#   - region_files/*_prelim_xy.txt and *_matched_xy.txt are RETIRED for
#     fan-outs (AM8), NOT deleted: the IS plans cite them as what ran.
# DONE in howto_MAAP_arctic.sh: step 6 submits --tile_list; step 8 prunes
#   (dry run, then --write, then commit); step 9 submits matched from the
#   same list and explains the pre-flight; step 10b uses the same list for
#   monthly and replaces its "expect a failed job" note with the TL1
#   behaviour.  Tagged UNTESTED on DPS where they are.
#   howto_MAAP_GL.sh and _AA.sh NOT touched -- GL and AA have not run.
#
#
# ===========================================================================
# TL6. [BEN] [READY 2026-09-18 -- Ben told to register]  Register, then check_build_id MATCH.
# ===========================================================================
# TL1 is solver code, so it reaches DPS only through a rebuild.  The checkout
# must be clean and pushed first (register_algorithm.py refuses otherwise).
# Never while a fan-out is queued.  Then check_build_id.py --expect <sha>.
#
#
# ===========================================================================
# TL7. [DPS] [NOT STARTED -- needs Ben's go]  Prove TL1 on DPS.
# ===========================================================================
# Resubmit the monthly E1020_N-2580 prelim with --xy_file
# region_files/IS_i7a_xy.txt (one center, 1020000 -2580000; the file the I7a
# proof used) and the monthly args -- the job that failed as 9266c3d7.  EXPECT: SUCCESSFUL, "no fit written ... skipping error
# calculation" in the log, no tile at the prefix.  The same proof I7a had.
# Then prune_tile_list.py on that ledger should propose removing E1020 --
# which the IS list has already dropped, so it must report "not in list".
#
#
# ===========================================================================
# OPEN
# ===========================================================================
#   - SV has no list.  Not needed until SV runs.
#   - FIXED 2026-09-19 at Ben's request: AA/40km_tile_list.txt line 8945 was
#     'field_sizes' -- the report subdirectory's name, so the list was likely
#     made from a listing of a prelim/ directory.  Removed; nothing else
#     changed.  The list is now 8944 names, all tile names, all distinct, and
#     read_tile_list parses it.  All six 40 km lists are now clean.
#   - AA has two lists (200km, 40km); which one gates which step is not
#     settled here.  Not needed until AA runs.
