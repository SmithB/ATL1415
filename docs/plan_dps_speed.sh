#! /usr/bin/env bash
# ===========================================================================
# PLAN: find out why the SuiteSparse solves run 2-5x slower on DPS than on
# the ADE.
# Written 2026-09-24, before any code.  TENTATIVE.  Every step carries its
# own status tag.  STATUS 2026-09-24 late: Ben said "Go ahead with D1-D3";
# D1-D2 WRITTEN (tests/test_worker_facts.py); D3 = Ben registers.  Ben asked: "track down why the DPS cholmod runs are slow"
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
# D1. [ADE] [DONE 2026-09-24, needs D3 to deploy]  Every job reports its worker.
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
# AS BUILT: scripts/worker_facts.py prints the WORKER: line (build_id report,
#   every tile job's header, the bench job).  scripts/run_with_rusage.py adds
#   one "=== cpu [step]:" line per solve: CPU time / wall = cores actually
#   got, load1 min/med/max sampled every 10 s, /proc/stat steal, cgroup
#   cpu.stat throttled time -- the direct test of H2 (a quota or too few real
#   cores shows as cores << threads, or throttled > 0).  ogc_jobs.py
#   parse_worker / parse_cpu / parse_bench; collect_jobs.py: instance column,
#   cores/load/steal/throttled on each step line, a Workers summary (one line
#   per EC2 instance with its job count).  Old ledgers read unchanged ('-').
# ADE, for reference (STATEMENT): r5.4xlarge, 8259CL, vcpus 16,
#   threads_per_core 2, cpu_quota 15.39, blas openblas-0.3.34-SkylakeX-
#   pthreads, AND mem_limit_gib 14.9 -- the ADE's cgroup memory limit is
#   16.0 GB although MemTotal is 124 GiB.
#
#
# ===========================================================================
# D2. [ADE] [DONE 2026-09-24, needs D3 to deploy]  One identical system, both places.
# ===========================================================================
# step=bench: fetch the saved E1340_N-2420 iteration-0 system (A0.npz,
#   b0.npy, 17+15 MB, staged to the bucket), time sparseqr.solve at 1, 2
#   and 4 threads, and a fixed single-thread CPU loop, then exit.  ~6 min.
#   The same script runs on the ADE (numbers above).  Separates H1 (all
#   times slower by one factor) from H2 (1 thread fine, 4 threads do not
#   scale) directly, with no ATL11 read or S3 noise.
# AS BUILT: run.sh step=bench (args_file = the system's directory; default
#   s3://maap-ops-workspace/ben_smith/ATL1415/bench/E1340_N-2420_it0, staged
#   2026-09-24, byte-identical to the local copy); scripts/maap/bench_solve.py
#   (every timed solve is checked against x0 -- a wrong answer fails the job);
#   scripts/maap/submit_bench.py <n> [--queue q] writes a ledger for
#   collect_jobs.py.
# ADE BASELINE (STATEMENT: run.sh step=bench on this ADE, 2026-09-24):
#   BENCH: qr_t1=173.2 qr_t2=114.6 qr_t4=83.0 dgemm1=0.261 spmv=2.4 (seconds);
#   every solve matches x0 to <1e-9; peak RSS 3.49 GiB, 377 s total.
#
#
# ===========================================================================
# D3. [BEN] [DONE 2026-09-24]  Ben registered; check_build_id --expect 4d2c9ae
#     VERDICT MATCH.  That build_id job (default queue -32gb) ran on an
#     r5.xlarge: 8259CL, vcpus 4, threads_per_core 2, no cpu quota, 31 GiB.
# D4. [DPS] [DONE 2026-09-24, Ben: "Go ahead with the bench jobs"]  bench x2
#     alone, then x8 at once on -16gb, plus x1 on -32gb.  Ledgers
#     ~/ATL14_processing/maap_ledgers/bench_16gb_2x_20260924T234029_jobs.csv,
#     bench_16gb_8x_20260924T235345_jobs.csv, bench_32gb_1x_20260924T235348_jobs.csv.
#     STATEMENT (collect_jobs.py + job metrics), seconds, SPQR on the saved system:
#       machine                  qr_t1  qr_t2  qr_t4  steal    n
#       ADE r5.4xlarge            173    115     83    0       1
#       -16gb t3.xlarge, steal<1s 181-189 121-125 129-135 <1 s  5
#       -16gb t3.xlarge, steal>40s 212-342 118-316 139-175 45-397 s 5
#       -32gb c5.4xlarge (8275CL 3.0 GHz, 16 vcpus)
#                                 148    103     76    0.1     1
#     - EVERY -16gb job ran on a t3.xlarge: BURSTABLE, 4 vcpus = 2 physical
#       cores, same 8259CL CPU as the ADE.  -32gb is NOT one type: r5.xlarge
#       (build_id) and c5.4xlarge (bench) on the same night.
#     - H1 (slower CPU): NO -- same model; with no steal, 1 thread is within
#       5-9% of the ADE.
#     - H2 (fewer real cores): YES -- 2 physical cores.  4 threads is SLOWER
#       than 2 on every no-steal t3 (129-135 vs 121-125 s); run.sh uses nproc=4.
#     - STEAL: half the t3 jobs lost 45-397 s of CPU to the hypervisor
#       (13% of all vcpu time at worst), and those ran 1.2-2.5x slower.  NOT
#       simple credit exhaustion: instance ...cf1642 lost 334 s on its second
#       job and 0.7 s on its third, straight after.  Cause unknown from here
#       (burstable hosts are shared; credit mode not visible to the job).
#     - H3 (jobs sharing a node): NO for these -- every instance ran one job
#       at a time (start/end times; instances are reused back to back).
#     - The prelim fan-out's slow iteration 0 (up to 855 s) is therefore more
#       likely steal on t3s than co-located jobs; the cpu line on the next
#       real tiles will say (steal per step).
# D5. [ADE] [IN PROGRESS]  Read, decide.  RECOMMENDATIONS (for Ben):
#     R1  threads = physical cores (nproc / threads_per_core) in run.sh: 2 on
#         a t3.xlarge.  Small (~6% on the QR), free, and it stops 4 threads
#         fighting over 2 cores in the error kernel too.
#     R2  run CPU-bound steps on a non-burstable queue.  The best t3 is ~1.7x
#         slower than the c5.4xlarge -32gb gave, and the worst ~3x -- but
#         -32gb is a mix of types.  QUESTION for the MAAP admins: which
#         instance types does each queue map to, can -16gb avoid t3, and are
#         the t3s in standard or unlimited credit mode?
#     R3  a real prelim tile on -32gb vs -16gb (E1340_N-2420) to measure the
#         whole job, not just the bench.  Needs Ben's go.
# Original D5 text:  If the cause is the instance type or
#     node packing, the fix is a queue choice or a MAAP admin question, not
#     code; ~/ATL14_processing/maap_resource_estimate.txt gets the numbers.
# ===========================================================================
