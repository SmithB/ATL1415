#! /usr/bin/env bash
# ===========================================================================
# PLAN: solve the fit's least-squares problem with CHOLMOD on the normal
# equations (plus one refinement step) instead of SPQR, in LSsurf's
# iterate_fit.
# Written 2026-09-24, before any code.  TENTATIVE.  Every step carries its
# own status tag.  Ben: "Go ahead with the bench jobs and the CHOLMOD plan."
# QC1-QC3 answered 2026-09-24.  C0-C4 DONE, ALL FIVE TILES PASS.  NEXT: Ben merges
# LSsurf cholmod_fit (00eef0c), then registers (C5-C6); C7 needs Ben's go.
# UPDATE 2026-09-25: C5-C8 DONE.  --solver=cholmod is IN the IS quarterly
# args file (local + run_args on S3); and, at Ben's request, in the IS
# MONTHLY args file too (local + run_args on S3), 2026-09-25.
# ===========================================================================
# Provenance per claim: STATEMENT = verified 2026-09-24, with how;
# DECIDED = Ben said so; RECOMMENDATION = mine; QUESTION = open.
# Tags: [ADE] / [DPS] / [BEN]; [NOT STARTED] [NEEDS CODE: x] [DONE].
#
# QUESTIONS FOR BEN (answer on these lines):
#   QC1  The pass bar for the FIT (z0, dz, dzdt, biases).  Ben's 5% rule is
#        for the sigma fields only ([[error-estimates-5pct-tolerance]]).
#        RECOMMENDATION: identical editing (same in_TSE on every iteration)
#        and max |d model| <= 1e-4 m on every tile tested; anything else is
#        reported tile by tile for Ben to judge.
#        QC1 answer (Ben, 2026-09-24): "10^-3 m is the standard."  DECIDED:
#          pass = max |d model| <= 1e-3 m on every model field; editing
#          (in_TSE) identity is reported, not required.
#   QC2  Default solver once validated: cholmod-with-SPQR-fallback for every
#        caller, or opt-in (ATL1415 passes it; discover and the deprecated/
#        callers keep SPQR)?  RECOMMENDATION: an LSsurf argument
#        solver='spqr'|'cholmod', default 'spqr'; ATL1415 passes 'cholmod'.
#        Nothing changes for anyone who does not ask.
#        QC2 answer (Ben): "Cholmod should be opt-in."  DECIDED: LSsurf
#          solver='spqr' default; ATL11_to_ATL15 --solver {spqr,cholmod}, default
#          spqr -- a region opts in with --solver=cholmod in its args file.
#   QC3  Which tiles to validate on (C4).  RECOMMENDATION: the list in C4.
#        QC3 answer (Ben): "go ahead with those tiles."  DECIDED: the C4 list.
#
#
# ===========================================================================
# BACKGROUND.  STATEMENT, 2026-09-24, local ADE (r5.4xlarge), E1340_N-2420
# quarterly prelim, iteration-0 system saved from a real run:
# A = Ip_r TCinv G, 1,812,872 x 477,008, nnz 7.7M.
# ===========================================================================
# WHERE THE FIT'S TIME GOES.  cProfile of the whole fit (463 s, 4 threads):
#   3 x sparseqr.solve  280 s (84 / 97 / 98)   <- the target
#   read_ATL11 (S3)     126 s
#   setup_averaging_ops  16 s
#   rest of iterate_fit  ~3 s total: residuals, sigma_extra, editing
#
# THE ALTERNATIVE, timed on the same saved system, uncontended:
#                              1 thread   4 threads
#   SPQR solve (as now)         167 s      84 s
#   CHOLMOD on A'A: analyze      --        7.1 s
#                   factorize    --       10.8 s
#                   total (A'A, analyze, factorize, solve)
#                                30 s      18 s
#   One refinement step  x += F.solve(A'(b - A x))  0.5 s.
#   max |x - x_SPQR|: 3.4e-5 m without refinement, 1.1e-9 m with one step.
#   CHOLMOD warns "nearly singular (rcond=1.10e-12)" -- cond(A'A) is about
#   1e12, so the unrefined normal-equation solve loses ~12 digits; the
#   refinement step restores them (the corrected semi-normal equations).
#   Gain per solve: 4.7x at 4 threads, 5.5x at 1 -- and the DPS workers are
#   2 physical cores (r5.xlarge, 2026-09-24 build_id job), where the
#   1-thread ratio is the one that matters most.
#   Prototype: session scratchpad, fit/solvers.py and fit/refac.py;
#   scikit-sparse 0.5.0 built into a scratch --target dir (the ATL14 env is
#   untouched).
#
# ONE UNEXPLAINED OBSERVATION.  A factorize() on an already-used
#   CholeskyFactor once took 471 s, versus 11 s fresh; it has not recurred
#   uncontended (two fresh cho_factor calls: 18.2, 17.6 s).  The plan builds
#   a fresh factor every call, and C1's timing line would show a recurrence.
#
# WHAT IT DOES NOT TOUCH.  The error step (calc_and_parse_errors) keeps
#   sparseqr.rz: its R feeds inv_tr_upper, rz is ~10 s of it, and the reach
#   kernel is validated against it.  (R from Cholesky equals SPQR's R up to
#   row signs for the same column order -- a possible later step, not this.)
#
#
# ===========================================================================
# C0. [ADE] [DONE 2026-09-24, ATL1415 8f3456b]  The dependency.
# ===========================================================================
# STATEMENT: conda-forge ships scikit-sparse 0.5.0 for py313 (conda search,
#   2026-09-24), linked against conda-forge suitesparse -- the same stack
#   sparseqr uses.  Add `scikit-sparse` to environment.yml (conda, not pip:
#   no compile at build time).  0.5.0 has a NEW API (cho_factor /
#   CholeskyFactor(A, order=...)); 0.4.x code examples do not apply.
# LSsurf declares it as an OPTIONAL import: solver='cholmod' without it
#   raises, loudly, naming the package (fail loudly, no silent SPQR).
#
#
# ===========================================================================
# C1. [ADE] [DONE 2026-09-24, LSsurf cholmod_fit 00eef0c; ATL1415 --solver 8f3456b]  The solver.
# AS BUILT: LSsurf/ls_solvers.py; up to 3 refinement steps (stop at 1e-10);
#   fallback also when rcond(A'A) < 1e-14 (a repeated column gave 2.9e-16 and
#   NO CHOLMOD error).  ATL11_to_ATL15 --solver checks LSsurf.ls_solvers and
#   sksparse at PARSE time: smooth_fit ignores unknown keywords, so an old
#   LSsurf would otherwise run SPQR silently.
# ===========================================================================
# One function, used by iterate_fit where sparseqr.solve is called now:
#   solve_ls(A, b, threads, solver) -> x
#   'spqr':    exactly today's loop (orderings 6 then 5).
#   'cholmod': N = A'A (CSC), cho_factor(N, order='metis'), x = solve(A'b),
#              then ONE refinement step.  Report per call: times for A'A,
#              factor, solve, refine, and the refinement's relative update
#              |dx|/|x| -- the measure of how much the raw solve lost.
#   FALLBACK to SPQR, with a printed reason, when: CHOLMOD raises (not
#              positive definite, out of memory), OR the refinement update
#              is not small (|dx|/|x| > 1e-6, a threshold to confirm in C4).
#              A fallback is visible in the log; C4 counts them.
# threadpool_limits stays as it is (CHOLMOD's supernodal factor threads
#   through BLAS, like SPQR).
# MEMORY: the factor of A'A is the SAME SIZE as SPQR's R (same fill, same
#   ordering family) and there is no Q; A'A itself is ~nnz(R)-ish at most.
#   Measure in C4 against the fit peaks (5.15 / 9.52 GiB on the two E1340s).
#
#
# ===========================================================================
# C2. [ADE] [DONE 2026-09-24, LSsurf tests/test_ls_solvers.py, 22 passed]  Unit tests.
# Note: column scaling did NOT make a test the unrefined solve fails;
#   near-collinear columns (eps 1e-4, 1e-5) do.  Mutants caught: refinement
#   disabled (2 fail), rcond check removed (rank-deficient case fails).
# ===========================================================================
# Random sparse least-squares problems, well and badly conditioned: cholmod
#   matches spqr within 1e-9 relative; a rank-deficient A falls back to SPQR
#   and says so; solver='cholmod' without scikit-sparse raises.  Mutation:
#   drop the refinement step and confirm the badly conditioned case fails.
#
#
# ===========================================================================
# C3. [ADE] [DONE 2026-09-24]  Install into the ATL14 env (reversible).
# Dry run: scikit-sparse 0.5.0 the ONLY new package, nothing changed.
#   LSsurf cholmod_fit installed --no-deps from the checkout.
# ===========================================================================
# conda install -n ATL14 -c conda-forge scikit-sparse=0.5.0, dry-run first:
#   it must NOT change suitesparse, numpy or scipy (if it would, stop and
#   ask).  Then the LSsurf branch, pip install --no-deps as in RK5.
#
#
# ===========================================================================
# C4. [ADE] [DONE 2026-09-24: ALL PASS]  Validate on real tiles.
# AS RUN: --solver=cholmod prelim fit + error step, --THREADS=4, locally,
#   compared with the DPS SPQR tiles (0332 quarterly / monthly, 61a19af) --
#   earlier the local SPQR rerun of E1340_N-2420 matched DPS's edits exactly.
#   STATEMENT:
#   tile               solve s/iter  rcond      edit flips   worst model |d|  fit peak
#   q E1340_N-2420     18            2e-13..1e-12  0/65207      1.2e-9 m      5.16 GiB
#   q E1340_N-2460     28-29         2e-12..7e-12  0/273383     4.2e-9 m      7.89 GiB (DPS SPQR 9.52)
#   q E1180_N-2420     14            2e-12..4e-12  0/159        3.4e-9 m      4.24 GiB
#   q E1020_N-2420     15            3e-12..1e-11  0/6819       1.1e-9 m      4.35 GiB
#   m E1340_N-2460     5-5.5         8e-13..4e-12  0/280956     5.1e-11 m     2.41 GiB
#   Worst sigma |d|/sigma 5e-11.  No fallback.  Refinement converged in 1-2
#   steps.  SPQR for comparison: E1340_N-2420 84-98 s/iter locally; DPS
#   185-461 (q 2420), 311-726 (q 2460), ~210 (monthly 2460) s/iter.
#   Whole local fit, E1340_N-2420: 463 s SPQR -> 220 s cholmod (the rest is
#   mostly the S3 ATL11 read).  Logs: session scratchpad c4/.
# ===========================================================================
# Full prelim fit + error step, spqr vs cholmod, --THREADS=4 (and 2, the DPS
#   physical-core count), scratch output only:
#     E1340_N-2420  moderate, the tile everything above was measured on
#     E1340_N-2460  densest IS prelim (N_fit 273k, 9.52 GiB fit peak)
#     E1180_N-2420  sparse (158 matched points; 65k-ish prelim)
#     E1020_N-2420  isolated, no neighbours
#     E1340_N-2460 monthly  different grid ([25,25,94])
#   (QC3: add a GL or AA tile if Ben wants a different region's geometry.)
# Report per tile: fit time per iteration, peak RSS, fallbacks, in_TSE
#   identical per iteration (Y/N, count of flips), max |d| per model field
#   (z0, dz, dzdt lags, biases), and every sigma field's |dsigma|/sigma
#   (should be ~0: the error step is unchanged, but its input edits are
#   the fit's).  Pass = QC1.
#
#
# ===========================================================================
# C5. [BEN] [DONE 2026-09-25]  LSsurf main 04abfbb (merged 02:50:21Z).
# C6. [BEN] [DONE 2026-09-25]  Registered; check_build_id --expect d59b140 MATCH,
#     built 02:52:57-02:55:55Z (after the merge).  That job: c5.4xlarge, lifecycle=spot.
# C7. [DPS] [DONE 2026-09-25, PASS]  Ben: add --solver=cholmod to input_args_IS.txt.
#     Job bb01eeea (ledger IS_cholmod_smoke_jobs.csv), NO --tile_prefix (the
#     canonical tile was not overwritten; the smoke tile is in dps_output).
#     t3.xlarge spot, threads 2 (physical cores): fit 263 s (was 1214),
#     error 254 s (was 1270), job 530 s (was 2510), peak 5.13 GiB (was 6.43).
#     cholmod 25-26 s/iter, rcond 1.9e-13..1.2e-12, no fallback; vs the SPQR
#     tile: 0/65207 edit flips, worst model 1.2e-9 m, sigma 4.2e-11.
#   (original C7 text:) One smoke prelim (E1340_N-2420):
#     the fit's QR lines replaced by cholmod timing lines, the WORKER and cpu
#     lines (plan_dps_speed D1) for the cores it got, and the tile compared
#     with the SPQR one as in C4.
# C8. [ADE] [DONE 2026-09-25]  maap_resource_estimate.txt: dated UPDATE appended
#     (one-tile numbers; GL/AA totals deliberately NOT rescaled).
# ===========================================================================
