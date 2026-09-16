#! /usr/bin/env bash
# ===========================================================================
# PLAN: a checker for the per-tile field-size reports.
# Written 2026-09-16 as an ASSIGNMENT BRIEF for an agent.  NO CODE YET.
# ===========================================================================
# THE ASK (Ben, 2026-09-16): "there should be a script that checks the field
# size reports written by each tile job and makes sure that the outputs are
# the right size and that the prelim step produces a sigma_dz that's the right
# size."  check_tiles.ipynb is explicitly NOT needed (same message), so this
# replaces it rather than restoring it.
#
# READ THE STATUS TAG ON EVERY STEP.  Nothing below is written yet.
# Provenance is marked per claim: VERIFIED = checked in this repo or against
# the real IS outputs on 2026-09-16; DECIDED = Ben said so; RECOMMENDATION =
# mine, and overridable; QUESTION = genuinely open, do not guess.
#
#
# ===========================================================================
# BACKGROUND THE AGENT NEEDS.  All VERIFIED 2026-09-16.
# ===========================================================================
# WHO WRITES THE REPORTS.  Two places, same per-file logic:
#   ATL1415/ATL11_to_ATL15.py  save_field_size_report()  -- each tile job
#     self-reports as it writes its tile.  This is the live path.
#   ATL1415/scripts/make_field_size_report.py -- the older batch pass over a
#     whole step directory.  STILL EXISTS.  docs/howto_MAAP_GL.sh:170 says it
#     "is no longer needed"; docs/workflow_overview.md:104 still documents it.
#     THIS PLAN DOES NOT TOUCH EITHER.  Both WRITE reports; the new script
#     READS them.  Do not "consolidate" them into the checker.
#
# WHERE THEY LIVE: <step_dir>/field_sizes/<tile>_report.json, beside the
# tiles, e.g. .../rel006/north/IS/prelim/field_sizes/E1300_N-2500_report.json
#
# WHAT ONE LOOKS LIKE.  Exactly three keys, no more:
#   prelim : {"file": "/PlKweE/output/prelim/E1300_N-2500.h5",
#             "dz/dz": [61,61,32], "dz/sigma_dz": [61,61,32]}
#   matched: {"file": "/dZiINo/output/E1300_N-2500.h5",
#             "dz/dz": [61,61,32], "dz/sigma_dz": null}
# TRAP: "file" is the path INSIDE THE DPS WORKER ("/PlKweE/output/..."), not
# anywhere on the ADE.  It is useless for locating anything.  MATCH TILES TO
# REPORTS BY BASENAME ONLY.  Do not os.path.exists() the "file" value.
#
# A MISSING FIELD IS null, NOT an absent key: save_field_size_report catches
# the per-field exception and stores None, so `"dz/sigma_dz": null` means "the
# field was not in the file", and a KeyError means something else went wrong.
#
#
# ===========================================================================
# C1. [NOT WRITTEN]  What "the right size" is.  DERIVE IT, DO NOT HARD-CODE.
# ===========================================================================
# VERIFIED 2026-09-16, by deriving it and comparing against all 56 real IS
# reports (28 prelim + 28 matched): the expected dz shape follows from the
# run's own args file, and for IS it comes out at exactly [61,61,32], which is
# what every one of the 56 reports carries.
#
#   nx = ny = W / dz_spacing + 1        W from -W, dz_spacing from -g's 2nd
#   nt      = (t1 - t0) / dt + 1        t from -t, dt from -g's 3rd
#
# For IS (~/ATL14_processing/rel006/north/IS/input_args_IS.txt):
#   -W=60000   -g=100,1000,0.25   -t=2018.75,2026.5
#   -> 60000/1000+1 = 61,  (2026.5-2018.75)/0.25+1 = 32   -> [61,61,32]
#
# HARD-CODING 61x61x32 WOULD BE WRONG, and silently so:
#   - VERIFIED: -t is the span that sets nt, NOT --t_crop.  IS carries
#     --t_crop=2019,2026.25, which would give 30 epochs, not the 32 that the
#     real files have.  Use -t.  A checker built on t_crop would fail all 56.
#   - VERIFIED (docs/plan_IS_run.sh, I7): AA solves two halves at DIFFERENT
#     widths, so one region can legitimately have two expected shapes.
#   - GL and AA are far larger runs than IS and are the whole point of having
#     this script; a constant tuned to IS is worthless there.
# RECOMMENDATION: take the args file as an argument, parse the three flags,
# derive the shape, and PRINT the derivation in the output so a human can see
# what it checked against rather than trusting it.
# QUESTION FOR BEN, do not guess: should a tile whose shape merely DISAGREES
# with its neighbours (but matches no derived expectation, e.g. because no
# args file was passed) be an error, or is the derived check the only one?
#
#
# ===========================================================================
# C2. [NOT WRITTEN]  The four checks.
# ===========================================================================
# 1. SHAPE.  Every report's "dz/dz" equals the derived [nx,ny,nt].
# 2. SIGMA IN PRELIM.  For step=prelim, "dz/sigma_dz" MUST be present and
#    equal to "dz/dz".  This is the half Ben named explicitly.
#    VERIFIED: all 28 IS prelim reports satisfy it.
# 3. SIGMA IN MATCHED.  For step=matched, "dz/sigma_dz" MUST be null.
#    DECIDED (Ben, 2026-09-16): "There should be no sigma in a matched result.
#    We use the uncertainties calculated in the prelim step."
#    VERIFIED in the code: make_ATL1415_queue.py adds the --calc_error_for_xy
#    companion only in its `if not args.step=='matched'` branch, and run.sh's
#    matched branch calls run_solve once.  VERIFIED in the data: all 28 IS
#    matched reports have null.
#    SO THE CHECK IS STEP-DEPENDENT, AND INVERTED BETWEEN THE TWO STEPS.  A
#    checker that just demands sigma everywhere would fail all 28 matched
#    tiles; one that ignores sigma would miss the bug Ben actually asked for.
#    A matched tile WITH a sigma is as much a fault as a prelim tile without.
# 4. PAIRING.  Every tile .h5 in the step dir has a report, and every report
#    has a tile.  ORPHAN REPORTS ARE A REAL FAILURE MODE, not a hypothetical:
#    ATL11_to_ATL15.remove_tile_and_report() deletes a tile and its report
#    together precisely so a report cannot outlive its tile, and if one ever
#    does, this is the check that catches it.
#
#
# ===========================================================================
# C3. [NOT WRITTEN]  Shape of the tool.
# ===========================================================================
# RECOMMENDATION, and the agent may argue with it:
#   PATH   scripts/check_field_sizes.py, executable, python3.
#          NOT scripts/maap/ -- VERIFIED, docs/plan_IS_run.sh I7 states that
#          directory is the ADE's maap-py job-API scripts and nothing else
#          belongs there.  NOT ATL1415/scripts/ -- that is the installed
#          package's console scripts; this is an operator tool like
#          scripts/s3_tiles.py.
#   USAGE  check_field_sizes.py <step_dir> [--args_file F] [--step prelim|matched]
#          Infer --step from the trailing directory name when it is literally
#          'prelim' or 'matched'; REQUIRE it explicitly otherwise rather than
#          guessing, because guessing wrong inverts check 3.
#   OUTPUT One line per problem, then a summary: how many reports, how many
#          tiles, how many passed, and the derived shape it checked against.
#          Silence on success is wrong here -- print the counts, so an
#          operator can tell "all fine" from "found nothing to check".
#   EXIT   0 all good; 1 at least one check failed; 2 could not run at all
#          (no such directory, no reports, unparseable args file).  Matches
#          check_build_id.py's convention of 2 for "the check did not happen",
#          which is VERIFIED in that file's header.
# NO NEW DEPENDENCIES: json, os, glob, argparse, sys.  h5py is NOT needed --
# this reads the reports, not the tiles.  If the agent finds itself opening a
# .h5, it has misread the assignment.
#
#
# ===========================================================================
# C4. [NOT WRITTEN]  Tests, and the free gift of real data.
# ===========================================================================
# THERE ARE 56 REAL REPORTS ON DISK RIGHT NOW, and they are the best fixture
# available.  VERIFIED present 2026-09-16:
#   ~/ATL14_processing/rel006/north/IS/prelim/field_sizes/   28, sigma [61,61,32]
#   ~/ATL14_processing/rel006/north/IS/matched/field_sizes/  28, sigma null
#   ~/ATL14_processing/rel006/north/IS/input_args_IS.txt     the args file
# BOTH DIRECTORIES MUST PASS CLEAN.  That is the acceptance test: if the new
# script reports a problem with the IS run, the script is wrong, because that
# run is known good (28/28 prelim usable, 28/28 matched successful, all
# verified against S3 -- docs/plan_IS_run.sh I7/I8).
# THEN THE FAILURE CASES, on COPIES in a temp dir, never on the real tree:
#   - a prelim report with "dz/sigma_dz": null            -> must fail (2)
#   - a matched report with a sigma shape filled in       -> must fail (3)
#   - a report with dz/dz [61,61,30]                      -> must fail (1)
#   - a report whose tile has been deleted                -> must fail (4)
#   - a tile with no report                               -> must fail (4)
#   - an empty field_sizes dir                            -> exit 2, not 0
# DO NOT WRITE INTO ~/ATL14_processing.  Those 56 tiles cost ~19 worker-hours
# and a failed tile is unrecoverable (VERIFIED: a failed OGC job gets no
# dps_output prefix at all).  Copy what you need somewhere else.
#
#
# ===========================================================================
# C5. [NOT WRITTEN]  Wiring it into the docs, LAST.
# ===========================================================================
# Once it passes C4, and not before:
#   - docs/plan_IS_run.sh I5 currently points at check_tiles.ipynb, which
#     VERIFIED does not exist anywhere under /home/jovyan.  Repoint I5 at this
#     script and say plainly that the notebook was never found.
#   - docs/workflow_overview.md section 5 documents make_field_size_report.py
#     (the WRITER).  Add the checker beside it; do not conflate them.
# RECOMMENDATION: one commit for the script + tests, a second for the doc
# rewiring, so the doc change is revertible on its own.
#
#
# ===========================================================================
# HOW TO HAND THIS TO AN AGENT
# ===========================================================================
# The brief is self-contained: point the agent at this file and the four
# numbered steps.  It needs no MAAP access, no S3 access and no credentials --
# everything it must read is local.  It must not submit jobs, and it has no
# reason to.
# SCOPE FENCE, state it explicitly when assigning: the agent writes ONE new
# script plus its tests, and touches the two docs in C5.  It does NOT modify
# ATL11_to_ATL15.py, make_field_size_report.py, run.sh, or anything under
# scripts/maap/.  If a check it is asked for seems to require changing the
# WRITERS, that is a QUESTION FOR BEN, not a licence to edit them.
