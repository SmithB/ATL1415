#! /usr/bin/env bash
# ===========================================================================
# PLAN: a faster error-propagation kernel in LSsurf -- replace inv_tr_upper
# with a sparse-reach, multithreaded solve, and stop recomputing Ip_c.Rinv
# once per averaging operator.
# Written 2026-09-24, before the code.  TENTATIVE.  Every step carries its
# own status tag.  STATUS 2026-09-24 (late): RK0-RK4 DONE.  Ben merged
# reach_kernel into LSsurf main (87ed699); the ADE ATL14 env is reinstalled
# from it.  RK3 tests on LSsurf branch reach_kernel_tests (65779e8) for Ben
# to merge.  NEXT: RK6 (Ben registers), then RK7 (needs Ben's go).
# ===========================================================================
# WHAT BEN DECIDED (2026-09-24, in his words):
#   "Error estimates that differ by 5% are functionally identical."
#   "Let's keep 1e-5 as the tolerance.  Relaxing that doesn't buy much time.
#    We can replace inv_tr_upper."
#   On the uncommitted np.float -> np.float64 edits: "B" (delete the two
#   unused modules), "then write the plan".
#   On the plan: "Once LSsurf is committed, I'll merge into main, then LSsurf
#   should get reinstalled on MAAP.  I'd trust that reinstall, so I'd say
#   'no' to QR2 -- the build ID can just check the ATL1415 commit.  Before
#   RK3/4 I'll need to restart the MAAP instance to get a larger memory
#   allocation."
#
# Provenance per claim: STATEMENT = verified 2026-09-24, with how;
# DECIDED = Ben said so; RECOMMENDATION = mine; QUESTION = open.
# Tags: [ADE] / [DPS] / [BEN]; [NOT STARTED] [NEEDS CODE: x] [DONE].
#
#
# ===========================================================================
# BACKGROUND.  STATEMENT, 2026-09-24, measured on this ADE (4 vCPU).
# ===========================================================================
# WHERE THE TIME GOES.  Error step = LSsurf smooth_fit.calc_and_parse_errors
#   (smooth_fit.py:231-306): sparseqr.rz -> R; inv_tr_upper(R, nnz_max, 1e-5)
#   -> thresholded Rinv; sigma = row norms of Rinv; then, for each of the 40
#   averaging ops, op.grid_error(Ip_c.dot(Rinv)).
#   IS 0332 prelim on DPS (collect_jobs.py on IS_0332_prelim_jobs.csv): the
#   error step is 782-1599 s on every data-bearing tile, 40-70% of each job,
#   and barely depends on N_fit (158 points -> 782 s).  The cost is set by
#   the grid (~44.5k unknowns at -s 5,2), not by the data.
#   cProfile, E1180_N-2420 error step (641 s local):
#     inv_tr_upper                        443 s  (timed directly)
#     csr_tocsc + csr_matmat in the op loop ~155 s (41 calls each: Ip_c.dot
#                                          (Rinv) is rebuilt inside the loop)
#     sparseqr.rz                          10 s
#   E1340_N-2420 (21% dz mask; 1270 s error step on DPS): inv_tr_upper 677 s,
#   N=44853, nnz(R)=29.3M, nnz(Rinv)=113.2M.
#
# THE NEW KERNEL.  Same algebra (back substitution, R x = e_col, one column
#   at a time), walked by COLUMNS of R (CSC) instead of rows: once x[j] is
#   final it is scattered up column j, and a column j whose x[j] never
#   received a contribution is skipped.  Only the "reach" of col is touched;
#   the old kernel walks every row i <= col for every col.  Explained in full
#   in the 2026-09-24 session; prototype in the session scratchpad
#   (errprof/kern/inv_tr_fast.pyx: _reach_block; reach_threads.py).
#   Columns are independent, so blocks of columns run on threads (nogil).
#
# MEASURED (prototype vs the current kernel, tol 1e-5, stored entries
#   IDENTICAL -- same rows, cols, order; sigma differs by <= 1.2e-14 rel):
#                       current   reach x1   x2     x4
#     E1180_N-2420       443 s      87 s     44 s   40 s
#     E1340_N-2420       677 s     119 s     62 s   54 s
#   Op loop with Ip_c.dot(Rinv) hoisted: E1180 1.2 s, E1340 10 s (1 thread) /
#   5 s (4 threads), sigma bit-identical to what the real run wrote (all 40
#   ops).  2 -> 4 threads gains little here; this VM may be 2 physical cores.
#   DPS worker core count: UNKNOWN (no saved log records it).
#
# DPS GETS LSsurf UNPINNED: pyproject.toml:47 is
#   "LSsurf @ git+https://github.com/smithb/LSsurf.git" -> whatever main is
#   at build time.  Nothing records WHICH LSsurf commit a build used:
#   check_build_id.py stamps the ATL1415 commit only.
#
#
# ===========================================================================
# RK0. [DONE 2026-09-24]  Delete the dead Cython modules (Ben: "B").
# ===========================================================================
# LSsurf branch reach_kernel (off main b5b93fc), commit ec79a3f, NOT pushed:
# propagate_qz_errors.pyx and spsolve_tr_upper.pyx removed with their
# setup.py Extension lines.  Nothing imported them; both failed at import
# under NumPy 2.4.6 (np.float).  Verified: a clean scratch build_ext produces
# only inv_tr_upper.so, and LSsurf.smooth_fit imports.  The two commented-out
# imports in LSsurf/deprecated/smooth_xyt_fit.py:22-23 are left as they are.
#
#
# ===========================================================================
# RK1. [ADE] [DONE 2026-09-24, LSsurf 556bb8a -- formal tests are RK3]  Replace the kernel.
# ===========================================================================
# Keep the module, function name and return contract so the four callers
# (smooth_fit + three in deprecated/) need nothing beyond RK2:
#     inv_tr_upper(R, nnz, tol, threads=1) -> (rows, cols, vals, status)
#   rows/cols int32, vals float64, columns in descending order, rows within a
#   column descending -- the order the old kernel emits (verified identical).
#   status=1 when more than `nnz` entries would be stored, as today, so the
#   nnz_max retry loop in smooth_fit.py:276-282 still works.
#   `tol` stays a C float, as in the old signature, so the threshold is the
#   same float32-rounded 1e-5 and the stored entries stay identical.
# Inside: one CSC copy of R (sorted indices; diagonal = last entry per
#   column); columns cut into fixed 256-column blocks, highest first, so the
#   blocks and the output do not depend on the thread count; each thread
#   owns its x / mark work vectors; each block solves with the GIL released
#   into a buffer bounded by min(block cols x col_hi, nnz) entries.
# MEMORY: blocks are copied, in order, into ONE preallocated output
#   (np.empty(nnz), as today) and dropped; at most 2 x threads blocks are in
#   flight, so the peak is about the result plus a few blocks, not 2x the
#   result.  Rinv is 113M entries = 1.8 GB on E1340_N-2420, and the error
#   step already peaks at 8.3 GiB of 16 on the densest tiles.  The CSC copy
#   of R adds nnz(R) x 12 B (~350 MB on E1340_N-2420).
# DEFAULT threads=1 keeps the deprecated/ callers unchanged.
# AS BUILT, two things the prototype did not tell us (STATEMENT, measured on
#   the saved E1180_N-2420 R, same machine, back to back):
#   - comparing the C-float tol against a double in the inner loop cost 14%;
#     the kernel takes a double (the public argument stays a float, so the
#     float32-rounded threshold and the stored entries are unchanged);
#   - LSsurf builds at -O2 by default; setup.py now passes -O3 for this
#     extension (111 -> 100 s before the tol fix).
#   Final: 86 s (1 thread), 47 s (2), 41 s (4) vs 443 s; entries identical to
#   the old kernel's; block width (32-1024 cols) makes no difference.
#   Pre-commit check only: dense-inverse reference on random R (N=1..1000,
#   1 and 3-4 threads), overflow -> status=1.
#
#
# ===========================================================================
# RK2. [ADE] [DONE 2026-09-24, LSsurf 556bb8a]  Threads in; hoist Ip_c.Rinv.
# ===========================================================================
# a. smooth_fit.py:279: pass threads=args['THREADS'] (already set from
#    --THREADS / nproc by run.sh and used for sparseqr at :244).
# b. smooth_fit.py:302-306: compute Ip_c.dot(Rinv) ONCE before the loop,
#    `del Rinv` after E0 and that product, and pass the product to every
#    op.grid_error.  Same arithmetic, same result; removes 40 rebuilds.
# c. (added while writing) `del RR, CC, VV` once Rinv is built: the triplets
#    were held to the end of the function, ~16 B per Rinv entry.
#
#
# ===========================================================================
# GATE: BEN RESTARTS THE MAAP INSTANCE (more memory) BEFORE RK3/RK4.
# The restart likely clears /tmp, including the 2026-09-24 scratch captures;
# RK4 re-runs the real error step, so nothing there is needed.
# ===========================================================================
# RK3. [ADE] [DONE 2026-09-24, LSsurf 65779e8 on branch reach_kernel_tests]  Unit tests.
# ===========================================================================
# Reference = dense inverse of small random upper-triangular sparse R
#   (scipy.linalg.solve_triangular), thresholded at tol, in the kernel's
#   output order.  Check: same (rows, cols) entries; vals rtol 1e-12;
#   threads=1 and threads=4 identical; a too-small nnz returns status=1; a
#   1x1 and a diagonal R.  Mutation check: break the reach skip and confirm a
#   test fails.
# AS DONE (STATEMENT): tests/test_inv_tr_upper.py, 12 tests; LSsurf suite 13
#   passed.  Mutants built in scratch: dropping the reach bound (`lo`) fails 6;
#   keeping the diagonal only when above tol failed NOTHING at first (every
#   test diagonal of Rinv was > tol) -- test_diagonal now has 1/1e6 and
#   catches it.  TRAP: the LSsurf checkout holds a stale, git-ignored in-place
#   inv_tr_upper.cpython-313 .so (Sep 3, OLD kernel); `python -m pytest` or
#   `python -c` run from the checkout imports it.  test_new_kernel_is_installed
#   fails loudly on it (TypeError: 'threads') -- verified.  Run the tests with
#   the `pytest` executable, or delete that .so (Ben's call; not deleted).
#
#
# ===========================================================================
# RK4. [ADE] [DONE 2026-09-24, PASS]  End to end on real tiles, locally.
# ===========================================================================
# pip install the branch into the ATL14 env (reversible: reinstall main).
# Run the ERROR step only (--calc_error_for_xy, on scratch copies, as on
# 2026-09-24) for E1340_N-2420 and the densest prelim tile, E1340_N-2460
# (N_fit 273382, 9.52 GiB fit peak), under scripts/run_with_rusage.py.
# Pass: every sigma field within 1e-10 relative of the existing tile (Ben's
#   bar is 5%; anything above round-off means a bug, not an approximation);
#   error-step time and peak RSS recorded against the DPS numbers above.
# AS DONE (STATEMENT): LSsurf main 87ed699 installed in ATL14 (git+https,
#   --no-deps, the RK5 command).  Error step on scratch copies of the fetched
#   0332 prelim tiles, --THREADS=4, both tiles at once on the 16-core ADE,
#   under run_with_rusage.py.  Compared every dataset against the original:
#                   DPS old kernel (step error)   ADE new kernel
#     E1340_N-2420     1270 s   6.43 GiB            168 s   4.83 GiB
#     E1340_N-2460     1553 s   8.30 GiB            257 s   5.82 GiB
#   Worst |dsigma|/sigma over all 45 sigma fields: 6.8e-15 / 8.3e-15 (pass
#   bar 1e-10); no NaN-pattern change; every non-sigma field identical except
#   bias/sigma/val (also a sigma, 2.3e-15).  Caveat: different hardware from
#   DPS, whose core count is still unknown -- RK7 gives the DPS number.
#   Peak RSS FELL ~1.6-2.5 GiB (triplets freed, Ip_c.Rinv built once).
#
#
# ===========================================================================
# RK5. [BEN] [MERGED 2026-09-24 (87ed699); ADE reinstalled; RK3 tests branch still to merge]  Merge LSsurf; reinstall it on MAAP.
# ===========================================================================
# RECOMMENDATION: merge after RK4 passes -- RK3/RK4 are the real tests.
# QR1 DECIDED (Ben): I commit on LSsurf branch reach_kernel; BEN merges it
#   into LSsurf main.  No pin -- ATL1415 stays on the unpinned
#   git+https://github.com/smithb/LSsurf.git, so the next DPS build takes
#   main.  (Every discover install takes it too.)
# Then reinstall LSsurf on MAAP:
#   ADE:  conda run -n ATL14 python -m pip install --force-reinstall \
#           --no-deps git+https://github.com/smithb/LSsurf.git
#         (--no-deps: only LSsurf changes; the rest of the env stays put).
#   DPS:  the rebuild at RK6 reinstalls it.
# QR2 DECIDED (Ben): NO -- the build ID checks the ATL1415 commit only; the
#   reinstall is trusted.
#
#
# ===========================================================================
# RK6. [BEN] [NOT STARTED]  Register; check_build_id --expect <new sha> MATCH.
# ===========================================================================
# ONLY AFTER RK5's merge: the build pulls LSsurf main at build time, so a
# registration before the merge rebuilds the OLD kernel and still MATCHes.
# No DPS jobs may be in flight (split-build trap, plan_IS_run.sh).
#
#
# ===========================================================================
# RK7. [DPS] [NOT STARTED -- needs Ben's go]  One smoke prelim job.
# ===========================================================================
# E1340_N-2420 prelim on the new build.  Pass: successful; step error well
# under its 1270 s (expect ~100-200 s including setup and rz); sigma fields
# within 1e-10 of the existing tile; peak RSS no higher than 8.30 GiB + ~2
# GiB.  Also read run.sh's "threads :" line -- the first record of the DPS
# worker's core count.
#
#
# ===========================================================================
# RK8. [ADE] [NOT STARTED]  Records.
# ===========================================================================
# Update ~/ATL14_processing/maap_resource_estimate.txt (the error step was
# about half of every prelim job) and this plan's status tags.
