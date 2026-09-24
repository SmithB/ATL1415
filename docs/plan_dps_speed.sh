#! /usr/bin/env bash
# ===========================================================================
# PLAN: find out why the SuiteSparse solves run 2-5x slower on DPS than on
# the ADE.
# Written 2026-09-24, before any code.  TENTATIVE.  Every step carries its
# own status tag.  Ben asked: "track down why the DPS cholmod runs are slow"
# (the SPQR solves in LSsurf iterate_fit; SPQR sits on CHOLMOD).
# ===========================================================================
# Provenance per claim: STATEMENT = verified 2026-09-24, with how;
# RECOMMENDATION = mine; QUESTION = open.
# Tags: [ADE] / [DPS] / [BEN]; [NOT STARTED] [NEEDS CODE: x] [DONE].
#
#
# ===========================================================================
# WHAT THE EXISTING LOGS SAY.  STATEMENT, from the _stderr.txt of all 116
# IS 0332 jobs (prelim, smoke, matched, monthly; fetched to the session
# scratchpad) and get_job_metrics().
# ===========================================================================
# 1. DPS reports nothing about the hardware: get_job_metrics() machine_type,
#    architecture, machine_memory_size and operating_system are all null.
#    run.sh logs "threads : 4" (nproc) on every job.
# 2. cwltool logs, on every job: "Skipping Docker software container
#    '--cpus' limit despite presence of ResourceRequirement with coresMin
#    and/or coresMax".  The container is not CPU-limited by cwltool.
# 3. SAME TILE ALONE vs IN A FAN-OUT (QR seconds, iterations 0/1/2; jobs
#    running at the time from the job start/end times):
#      E1300_N-2500 quarterly smoke, 1 job running      353 / 334 / 325
#      E1300_N-2500 quarterly prelim, 21-28 running     855 / 329 / 325
#      E1340_N-2460 monthly smoke, 1 running            216 / 207 / 229
#      E1340_N-2460 monthly prelim, 10-27 running       219 / 210 / 234
#    -> iteration 0 of a burst start is up to 2.6x slower; once the burst
#       settles, a job in a fan-out runs as fast as a job alone.  Something
#       is shared during the burst (CPU on a shared node, or I/O).
# 4. THE STEADY STATE IS STILL SLOW.  E1340_N-2420 quarterly prelim, same
#    system, same TSE counts (57310, 56042), 4 threads:
#      DPS  461 / 277 / 185 s        ADE  84 / 97 / 98 s (1 thread: 167 s)
#    DPS's settled 185 s is slower than the ADE on ONE thread.  The single-
#    threaded old-kernel error step shows a smaller gap (~1250 s DPS vs
#    ~840 s ADE), so the multithreaded QR loses more than a per-core speed
#    difference explains.
# 5. THE ADE: r5.4xlarge (IMDS), Xeon Platinum 8259CL 2.5 GHz, 8 cores x 2
#    hyperthreads, cgroup quota 15.4 CPUs.  IMDS answers from the ADE.
#
# HYPOTHESES, none tested (RECOMMENDATION: test all three at once, D1-D2):
#   H1  slower or older CPUs on the 16gb queue's instances.
#   H2  "4 threads" is 2 physical cores x 2 hyperthreads, or a CPU quota
#       below 4, so 4 BLAS threads fight over fewer real cores.
#   H3  several jobs share one node; with no --cpus limit, each runs 4
#       threads against the same cores.  Fits the burst effect (3); by
#       itself does not explain the steady gap (4).
#   (H4, less likely: the image's OpenBLAS chose a generic kernel on the
#    worker CPU.  environment.yml does not pin BLAS; conda-forge gives
#    OpenBLAS by default, as on the ADE, but no job has reported it.)
#
#
# ===========================================================================
# D1. [ADE] [NEEDS CODE: run.sh, collect_jobs.py]  Every job reports its worker.
# ===========================================================================
# One "WORKER:" line at job start and a load line after each solve, cheap
# and always on (like BUILD_ID, so every tile records it):
#   instance_type, instance_id, az   -- IMDSv2, then v1, 2 s timeout;
#                                       "unreachable" if blocked (the
#                                       container's hop limit may stop it)
#   cpu_model, nproc, threads_per_core, cpu.max, cpuset  -- lscpu, cgroup
#   /proc/loadavg                    -- host-wide: shows neighbours
#   OpenBLAS version / architecture / num_threads -- threadpoolctl
# collect_jobs.py gains instance, cpu and load columns.  instance_id tells
# which jobs shared a node (H3).
# No change to the solve.
#
#
# ===========================================================================
# D2. [ADE] [NEEDS CODE: run.sh step=bench]  One identical system, both places.
# ===========================================================================
# step=bench: fetch the saved E1340_N-2420 iteration-0 system (A0.npz,
#   b0.npy, 17+15 MB, staged to the bucket), time sparseqr.solve at 1, 2
#   and 4 threads, and a fixed single-thread CPU loop, then exit.  ~6 min.
#   The same script runs on the ADE (numbers above).  Separates H1 (all
#   times slower by one factor) from H2 (1 thread fine, 4 threads do not
#   scale) directly, with no ATL11 read or S3 noise.
#
#
# ===========================================================================
# D3. [BEN] [NOT STARTED]  Register; check_build_id MATCH.
# D4. [DPS] [NOT STARTED -- needs Ben's go]  Run bench x2 alone, then bench
#     x8 at once (H3), on the 16gb queue; optionally one on -32gb to compare
#     instance types.
# D5. [ADE] [NOT STARTED]  Read, decide.  If the cause is the instance type or
#     node packing, the fix is a queue choice or a MAAP admin question, not
#     code; ~/ATL14_processing/maap_resource_estimate.txt gets the numbers.
# ===========================================================================
