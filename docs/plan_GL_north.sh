# plan_GL_north.sh -- the northern third of Greenland on MAAP (quarterly)
#
# ############################################################################
# ##  WRITTEN 2026-09-29.  TENTATIVE: nothing below has been run.          ##
# ##  Purpose: the "few hundred jobs" MAAP asked for, to size a full-scale ##
# ##  run.  The commands are docs/howto_MAAP_GL.sh's, on a subset list.     ##
# ############################################################################
#
# Status tags per step: TODO / DONE / BLOCKED.  DECIDED = Ben said so;
# RECOMMENDATION = mine; QUESTION = open.
#
# QUESTIONS FOR BEN (answer inline):
#   QN1  Where to cut.  RECOMMENDATION: the northern third by EXTENT,
#        tile centers y >= -1520 km: 557 prelim + <=557 matched jobs.
#        (By COUNT it would be y >= -1440 km, 498 tiles.)  GL centers run
#        y = -3320 .. -640 km; 1483 in all.
#        QN1 answer (Ben 2026-09-29): cut at -1520.  DECIDED.
#   QN2  How far.  RECOMMENDATION: prelim + matched + mosaic (N3-N6); no
#        netCDF, no monthly.  The DPS load MAAP asked about is prelim +
#        matched; the mosaic is ADE-only but is the one GL-scale step never
#        timed.  A netCDF of a third of GL would be mostly empty.
#        QN2 answer (Ben 2026-09-29): go as far as netCDF.  DECIDED:
#          N3-N7 (N7 = ATL14 + ATL15 netCDF, local only); no monthly.
#   QN3  Concurrency.  RECOMMENDATION: --max_in_flight 100 and ask MAAP
#        whether they want it higher -- how the cluster behaves at 100+
#        concurrent jobs is part of what they want to learn.  IS ran 29 at
#        once; nothing larger has run.
#        QN3 answer (Ben 2026-09-29): try 100 at a time.  DECIDED.
#   QN4  Output prefix.  RECOMMENDATION: the CANONICAL GL prefix
#        (.../ATL14_processing/rel006/north/GL), so these tiles count toward
#        the full run, which then submits only the remaining centers.
#        Caveat: the southernmost matched row (y = -1520 km) is matched
#        without its southern neighbours and must be re-run as part of the
#        full run.  The alternative is a dated test prefix (as the transects
#        used), which throws the tiles away afterwards.
#        QN4 answer (Ben 2026-09-29): standard GL output location.  DECIDED.
#
# WHAT IS KNOWN (statements, with provenance):
#   - GL transect, 19 tiles, quarterly prelim, on d00568c (2026-09-25,
#     maap_ledgers/GL_transect_cholmod*_collect.txt): 19/19 OK, mean 12 min,
#     max 17.5 min, peak 11.9 GiB (E80_N-920).  -> -32gb queue (howto step 3
#     rule: -16gb only under ~10 GiB).
#   - d00568c predates the ATL11 read fix (pointCollection PR #57 blockcache,
#     ATL1415 cb8e674).  Measured locally on E200_N-1880: read 252 -> 109 s,
#     4.37 -> 0.58 GB, identical solve inputs (plan_rerun_timing.sh E).
#   - input_args_GL.txt on the bucket is 0332 + --solver=cholmod (composed
#     2026-09-25).  Nothing at the canonical GL prefix yet.
#   - ADE ATL14 env: pointCollection reinstalled from main d36570b on
#     2026-09-29 (blockcache + mosaic spacing fix, PR #58); pC 269 tests and
#     ATL1415 181 tests pass.
#   - ADE disk: 374 GB free.  GL prelim tile ~0.14 GB -> ~80 GB prelim,
#     similar again for matched (matched dirs include fetched neighbours).
#   - The GL mask ends 2026.0, t_crop 2026.5 (plan_rerun_timing QT4: OK for
#     timing).  NOT a publishable product for that reason alone.
#
# ESTIMATE (from the transect, before the read fix, so an upper bound):
#   prelim 557 x ~12 min ~= 110 job-hours; matched ~1-3 min est ~= 20.
#   At 100 in flight, prelim is ~1.5 h of wall clock if the queue keeps up.

conda activate ATL14
cd ~/git_repos/ATL1415
repo=$PWD
rel=006; cyc=0332; tspan=2018.75,2026.5
ATL14_root=/home/jovyan/ATL14_processing
s3_root=s3://maap-ops-workspace/ben_smith
ledgers=$ATL14_root/maap_ledgers
runs=$ATL14_root/runs
region_dir=$ATL14_root/rel$rel/north/GL
s3_run=$s3_root/ATL1415/run_args/rel$rel/north/GL
s3_out=$s3_root/ATL14_processing/rel$rel/north/GL       # QN4
tag=GL_rel${rel}_${cyc}_north
L=$ledgers/GL_${cyc}_north
ymin=-1520                                              # QN1, km


# ===========================================================================
# N0. [ADE] DONE 2026-09-29: Ben registered; check_build_id MATCH at
#     18e3936 (built 17:14Z, maap_pgt=set; output saved as
#     maap_ledgers/GL_0332_north_check_build_id.txt).
#     Register, then prove the build.
# ===========================================================================
# Ben registers (checkout must be clean and pushed).  The new image installs
# pointCollection and LSsurf main at build time (pyproject.toml is unpinned),
# which is what brings in the read fix.
/srv/conda/envs/notebook/bin/python register_algorithm.py --dry-run
scripts/maap/check_build_id.py $s3_run/input_args_GL.txt maap-dps-worker-32gb
# GATE: VERDICT MATCH, maap_pgt=set, image commit >= this plan's commit.
# Run check_build_id WITHOUT a short timeout -- it waits on the queue.


# ===========================================================================
# N1. [ADE] DONE 2026-09-29: 557 centers, both N2 smoke tiles included.
#     The subset list (in ledgers, NOT in the checkout).
# ===========================================================================
sed 's/.*_N\(-\{0,1\}[0-9]*\).*/\1 &/' ATL1415/resources/GL/40km_tile_list.txt \
    | awk -v ymin=$ymin '$1 >= ymin {print $2}' > ${L}_tile_list.txt   # mawk: no match() arrays
wc -l ${L}_tile_list.txt                      # expect 557 for ymin=-1520


# ===========================================================================
# N2. [DPS] DONE 2026-09-29 on 18e3936, canonical prefix, ALL GATES PASS
#     (ledger maap_ledgers/GL_0332_north_smoke_jobs.csv, _collect.txt):
#                  job s (transect)   fit s (transect)   peak GiB   instance
#     E80_N-920     889  (1049)        418  (584)         12.07      r5.xlarge
#     E480_N-1040   621  (778)         319  (534)          9.54      r5.xlarge
#     N_ATL11/N_AT/N_XO/N_fit identical to the transect.  Tiles vs transect,
#     reported cells: E80 identical (0 m), E480 4.3e-12 m.  check_field_sizes
#     2/2 OK.  The read fix is live on DPS: fit step -28% / -40%; the error
#     step is unchanged, so whole jobs are -15% / -20%.  Both tiles fetched
#     to $region_dir/prelim; N3 uses ${L}_N3_tile_list.txt (555 centers,
#     these two removed) so they are not re-run.
#     Smoke on the new build: two transect tiles re-run.
# ===========================================================================
# E80_N-920 (the transect's memory peak) and E480_N-1040 (79N, Gr1km-v2
# tides).  Both are in the subset.  Their d00568c results are the baseline,
# so this measures the read fix ON DPS -- the number the estimate lacks.
printf '80000 -920000\n480000 -1040000\n' > ${L}_smoke_xy.txt
scripts/maap/submit_MAAP_jobs.py --xy_file ${L}_smoke_xy.txt \
    --step prelim --args_url $s3_run/input_args_GL.txt \
    --tile_prefix $s3_out --queue maap-dps-worker-32gb \
    --tag ${tag}_smoke --ledger ${L}_smoke_jobs.csv
scripts/maap/collect_jobs.py ${L}_smoke_jobs.csv > ${L}_smoke_collect.txt
# GATES: both successful, commit = N0 build, N_AT/N_XO/N_fit equal to the
# transect's (maap_ledgers/GL_transect_cholmod*_collect.txt); fit step
# shorter.  Compare tiles on reported cells only
# (session_tools_2026-09-25/compare_reported.py) -- expect ~1e-7 m.
# If QN4 = canonical, these two tiles are already production tiles, and N3
# will re-run them unless they are dropped from the list.


# ===========================================================================
# N3. [DPS] SUBMITTING since 2026-09-29 (Ben's go; start time in
#     ${L}_prelim_start.txt): 555 centers, 100 in flight, on 18e3936.
#     Ledger ${L}_prelim_jobs.csv, log ${L}_prelim_submit.log.
#     RESULT (2026-09-29 ~23:10Z): 522 successful, 32 FAILED, 1 STUCK
#     "running" (E280_N-1240, 47fdc904, submitted 18:20:51; ~5 h).  Failures
#     spread over 17:57-19:02Z, interleaved with successes -- NOT only the
#     start-up burst (an earlier note here said so; corrected).  By log:
#       5  MAAP runner /app/create_inputs.py -> stage_in.py: TCP connect
#          timeout (Errno 110, ~130 s) to api.maap-project.org
#          /api/environment/config; job ends before our container.
#       8  MAAP runner /app/get_maap_pgt_token.py: same connect timeout on
#          /api/members/ben_smith; runner STILL starts our container with
#          MAAP_PGT set; our credential step then fails in ~6 s.
#       10 our credential step (pC _s3fs_from_maap) fails after 135-140 s at
#          ~0 CPU -> earthaccess fallback -> AttributeError.  INFERRED (timing
#          only) to be the same connect timeout: the reason is not logged
#          because pointCollection/ps_scale_for_lat.py line 3 calls
#          warnings.filterwarnings("ignore") on import, silencing
#          _s3fs_from_maap's warning.
#       1  E-40_N-760: botocore NoCredentialsError at 8.6 s.  Undiagnosed.
#       8  no logs (get_job_result {} / one HTTP 500), all submitted
#          18:12-18:21Z.  Unknown.
#     NOT KNOWN: where connections failed (API, load balancer, network
#     path) or whether load caused it.  Draft note to MAAP with job IDs:
#     ~/ATL14_processing/maap_note_api_timeouts_2026-09-29.txt (Ben sends).
#     QUESTIONS (not now): pC ps_scale_for_lat global warnings filter;
#     credential failure should stop with a clear error (memory: fail
#     loudly) instead of falling back to earthaccess on a worker.
#     RETRY 1 (Ben's go): the 32 failed, submitted ~23:10Z, ledger
#     ${L}_prelim_retry1_jobs.csv.  E280_N-1240 not included (still running).
#     RETRY 1 DONE 23:30Z: 32/32 successful (mean 557 s, peak 12.32 GiB)
#     -- consistent with transient failures.  Prelim now 556/557 (incl. N2);
#     E280_N-1240 still 'running' at 23:32Z (~5 h 10 min) -- AWAITS BEN
#     (dismiss + retry?).  Not yet fetched: N4 fetch/check still to do.
#     Prelim fan-out.
# ===========================================================================
grep -vxE "E80_N-920\.h5|E480_N-1040\.h5" ${L}_tile_list.txt > ${L}_N3_tile_list.txt   # N2 tiles done (555)
nohup scripts/maap/submit_MAAP_jobs.py --tile_list ${L}_N3_tile_list.txt \
    --step prelim --args_url $s3_run/input_args_GL.txt \
    --tile_prefix $s3_out --queue maap-dps-worker-32gb \
    --tag ${tag}_prelim --ledger ${L}_prelim_jobs.csv \
    --max_in_flight 100 > ${L}_prelim_submit.log 2>&1 &      # QN3
# NEVER register while this runs.


# ===========================================================================
# N4. [ADE] DONE 2026-09-30 ~01:00Z.  Collect, fetch, check.
#     E280_N-1240: cancel_job -> HTTP 202 'dismissed', but status stayed
#     'running' (checked 3 min); it never wrote a tile.  Retry 2 (531a4e12)
#     successful, 743 s, 11.34 GiB.  PRELIM COMPLETE: 557/557 tiles, 82 GB
#     local; check_field_sizes 557/557 OK (3 s); no no-data centers.
#     All 557 successful jobs (main+retries+smoke): median 702 s, mean 644,
#     p90 811, max 937 s; 99.6 job-hours.  Peak memory median 10.8 GiB, p90
#     11.8, max 12.92 (none > 14).  Instances: 554 r5.xlarge, 2 m5.2xlarge.
#     ADE TIME: collect_jobs 555 jobs 1962 s (3.5 s/job); fetch_tiles 555
#     rows 2593 s (4.7 s/job, 522 tiles); retry1 fetch 32 tiles.
#     Files: ${L}_prelim{,_retry1,_retry2}_{collect,fetch}.txt,
#     ${L}_prelim_field_sizes.txt.
# ===========================================================================
scripts/maap/collect_jobs.py ${L}_prelim_jobs.csv > ${L}_prelim_collect.txt
scripts/maap/fetch_tiles.py  ${L}_prelim_jobs.csv $region_dir --step prelim
scripts/check_field_sizes.py $region_dir/prelim @$region_dir/input_args_GL.txt
# Time these three passes too: ADE-side cost per job is part of the answer
# for MAAP (they make API/S3 calls per job and have only run at 29).
# Failed jobs: retry per arctic howto step 5, NEW ledger.


# ===========================================================================
# N5. [DPS+ADE] 479/557; retry 1 78/78; MATCHED COMPLETE 557/557.  Submitted from 2026-09-30
#     01:26:50Z on 615d5ef (Ben's go): 557 centers, dry run skipped none,
#     100 in flight; ledger ${L}_matched_jobs.csv, log _matched_submit.log
#     (stdout buffered -- watch the ledger), start _matched_start.txt.
#     RESULT (03:10Z): all submitted by 01:55:41Z, 0 submit failures; nothing
#     in flight.  479 successful (479 matched tiles on the bucket), 78 FAILED.
#     Failures come in bursts: runs of 12-20 consecutive submissions
#     (01:39:32-01:54:29), with successes between the bursts.  Evidence:
#       74 no logs, no metrics, get_job_result {}, no triaged_job dir, no
#          tile -- look like they never reached a worker (NOT KNOWN why)
#        2 MAAP API ConnectTimeout at maap.MAAP() start-up (as N3)
#        2 'cannot kill container: No such container' (docker daemon)
#     Failed list: scratchpad only; the ledger + status re-derives it.
#     RETRY 1 (Ben's go): the 78, submitted from 03:14:17Z, --max_in_flight 40
#     (lower than 100 in case the bursts are load-related); list
#     ${L}_matched_retry1_tile_list.txt, ledger _matched_retry1_jobs.csv.
#     RETRY 1 DONE 03:30:43Z: 78/78 successful.  557/557 matched tiles on
#     the bucket.  Collect/fetch/check: NM0.
# ===========================================================================
nohup scripts/maap/submit_MAAP_jobs.py --tile_list ${L}_tile_list.txt \
    --step matched --args_url $s3_run/input_args_GL.txt \
    --tile_prefix $s3_out --queue maap-dps-worker-32gb \
    --tag ${tag}_matched --ledger ${L}_matched_jobs.csv \
    --max_in_flight 100 > ${L}_matched_submit.log 2>&1 &
# The y = -1520 km row lacks southern neighbours (QN4 caveat).


# ===========================================================================
# N6. [ADE] SUPERSEDED -> DONE ON DPS 2026-10-01 (plan_dps_mosaic.sh D7).  Mosaic, timed.
# ===========================================================================
cd $runs
make_mosaic_jobs.py -b $region_dir -rr GL -t $tspan -e ATL14 \
    --run_name GL_${cyc}_north_mosaic @$repo/default_args/quarterly.txt
cd GL_${cyc}_north_mosaic
seq 1 $(ls queue | wc -l) | xargs -P 4 -I{} env SLURM_ARRAY_TASK_ID={} \
    python $repo/scripts/run_with_rusage.py task{} bash slurm_run.sh
check_mosaic_outputs.py $runs/GL_${cyc}_north_mosaic --values
cd $repo


# ===========================================================================
# N7. [ADE] SUPERSEDED -> DONE ON DPS 2026-10-01 (plan_dps_mosaic.sh D7; products at .../rel006_0332_testing/north/GL).  netCDF, timed.  LOCAL ONLY -- do not publish.
# ===========================================================================
mkdir -p $runs/GL_${cyc}_north_nc && cd $runs/GL_${cyc}_north_nc
python $repo/scripts/run_with_rusage.py ATL14 \
    ATL14_write2nc.py @$region_dir/input_args_GL.txt > ATL14.log 2>&1
python $repo/scripts/run_with_rusage.py ATL15 \
    ATL15_write2nc.py @$region_dir/input_args_GL.txt > ATL15.log 2>&1
cd $repo
# No INVALID line in either log; XO rows NOT_SET in four attributes only
# (as IS).  The files take the CANONICAL names
# (ATL14_GL_0332_100m_006_02.nc ...) in $region_dir, but hold one third of
# GL and a mask ending 2026.0: NOT products.  Do NOT copy them to $s3_out
# (howto step 11 does, and the monthly args would then read this partial
# ATL14 as the reference).  The full run overwrites them.
# NEW: a netCDF of a partial region has never been written; the rest of
# the grid should come out empty (fill), not as an error.  Look at a quick
# plot of h and delta_h before calling it done.  Compare with rel005 over
# the northern third only (Ben's bar: no >10 m errors, no major gaps).


# ===========================================================================
# NM. MONTHLY PRELIM, GL north.  WRITTEN 2026-09-30; QM-A..C answered.
# ===========================================================================
# DECIDED (Ben 2026-09-30): after the current batch, take GL north through
#   MONTHLY PRELIM (revises QN2's "no monthly"), and submit the ~500 jobs
#   ALL AT ONCE -- no --max_in_flight -- because that is what most MAAP
#   users do and it is the case our N3/N5 numbers do not cover.
#
# STATEMENT: monthly is NOT one more submit after N5.  Its args carry
#   --ATL14_reference_file, the quarterly ATL14 of the same release, read by
#   URI from the bucket (plan_monthly_on_maap.sh 2, M0-M1; howto_MAAP_GL
#   11-18).  So N5 collect/fetch, N6 mosaic and N7 ATL14 netCDF come first.
#   N7 as written says the partial ATL14 must NOT go to $s3_out, because the
#   canonical name would then be read as the reference for the full run.
#
# QUESTIONS FOR BEN:
#   QM-A Reference DEM.  RECOMMENDATION: the N7 north-only ATL14, after it
#        passes your rel005 bar over the northern third (as QM1 for IS),
#        copied to a NON-canonical key, e.g.
#          $s3_root/ATL1415/run_args/rel006/north_monthly/GL/ref/
#            ATL14_GL_0332_100m_006_02_north_partial.nc
#        so nothing canonical exists until the full run.  Monthly reads the
#        reference only at data points, and every GL-north center has a
#        quarterly tile (557/557), so a north-only reference covers them.
#        Other options: B. canonical $s3_out name (overwritten by the full
#        run -- risky); C. rel005 ATL14 (not on the bucket; departs from IS).
#        QM-A answer (Ben 2026-09-30): north ATL14 at the side key.  DECIDED.
#   QM-B Smoke first?  RECOMMENDATION: yes, ONE center (E80_N-920, as N2)
#        before the ~556.  A bad reference fails every job the same way
#        (IS W5: a missing reference edits away every point), which would
#        spoil the all-at-once test as well as the run.  ~10-20 min.
#        QM-B answer (Ben 2026-09-30): yes, one center.  DECIDED.
#   QM-C What "all at once" means.  RECOMMENDATION: no --max_in_flight and
#        --rate 0 -- the submitter then POSTs back to back, as a user's
#        plain loop would.  Queue maap-dps-worker-32gb, as N3/N5, so the
#        only change from N3 is the submission pattern.  (IS monthly
#        needed ~2.3x less memory than quarterly, so 16gb would fit -- but
#        that queue is burstable t3 (memory: DPS queue hardware) and would
#        confound the comparison.)
#        QM-C answer (Ben 2026-09-30): no cap, --rate 0, 32gb.  DECIDED.
#
# REVISED 2026-09-30 ~17:00Z -- HOME QUOTA.  STATEMENT: /home/jovyan has a
#   150 GB quota (MAAP admin, via Ben); df does not show it.  NM0's fetch
#   filled it: 395 then 403 of 557 matched local, "[Errno 28] No space left
#   on device" (visible once fetch_tiles.py printed the whole error).
#   DECIDED (Ben): no DPS output is copied to /home.  The authoritative tiles
#   are on $s3_out (557 prelim + 557 matched, verified); a local file of a
#   different size is an error.  DONE: all 1920 local prelim/matched files
#   size-matched S3 and were deleted (list: $ledgers/
#   GL_0332_north_local_delete_2026-09-30.txt); home 9.3 GiB after.
#   DECIDED (Ben): mosaic (N6) and netCDF (N7) run as DPS JOBS, not on the
#   ADE ("instances are too unpredictable").  PLAN: docs/plan_dps_mosaic.sh
#   (TENTATIVE, QD1-QD7 open for Ben).  TODO: agree a plan for that
#   (entry points, CWL, registration, how a job reads the 557+557 tiles and
#   where it writes) -- written and agreed before any code.  UNTIL THEN
#   NM0-NM6 ARE BLOCKED and both drivers stay stopped; the drivers below
#   (NM0 fetch, local mosaic/nc) are superseded.
# NM0-NM1 RUN UNATTENDED since 2026-09-30 03:42Z by a RESTARTABLE driver on
#   the NFS home (survives the ADE instance closing; /tmp does not):
#     ~/ATL14_processing/maap_ledgers/GL_0332_north_NM_driver.sh
#     log  ..._NM_driver.log;  markers ..._NM_state/<step>.done
#   IF THE INSTANCE CLOSED: re-run it (command in its header); finished
#   steps are skipped, interrupted mosaic tasks requeued.  It stops at NM2.
#   Matched failure evidence (N5): maap_ledgers/GL_0332_north_matched_evidence/.
# NM0. [ADE] SUPERSEDED (no fetch; tiles verified on $s3_out).  N5 finish: collect/fetch/check matched (main + retry1).
# NM1. DONE ON DPS 2026-10-01 (plan_dps_mosaic.sh D7).  N6 mosaic, N7 ATL14 (+ATL15) netCDF, as written above.
# NM2-NM6 RUN UNATTENDED by driver 2 (launched 2026-09-30 ~04:00Z; waits for
#   driver 1's ALL_DONE): ~/ATL14_processing/maap_ledgers/GL_0332_north_NM_driver2.sh,
#   log _NM_driver2.log, markers GL_0332_monthly_NM_state/.  Helpers (tested
#   on existing jobs/files): ~/ATL14_processing/session_tools_2026-09-30/.
#   Ben 2026-09-30: "Go ahead regardless and we'll evaluate the ATL14 after
#   the fact" -- NM2 is REPORT ONLY, not a gate.  Remaining hard stops:
#   ref upload/readback, args upload/readback, smoke gate, field sizes.
#   MONTHLY ARGS composed 03:55Z (setup_ATL1415_region.py as NM3); sorted diff
#   vs quarterly = exactly IS M4's four (+ref, dzdt_lags, -b, -g) after I
#   APPENDED --solver=cholmod by hand, as the quarterly GL args were
#   (2026-09-25).  CHOICE MINE, flagged to Ben: monthly on cholmod is
#   untested for accuracy (C4 was quarterly); reversible by re-running on QR.
# NM2. [ADE] DONE 2026-10-01 (result below).  Compare ATL14 with rel005 over the north; REPORT ONLY.
#   python ~/ATL14_processing/session_tools_2026-09-30/compare_rel005.py \
#       $region_dir/ATL14_GL_0332_100m_006_02.nc $region_dir/ATL15_GL_0332_3mo_1km_006_02.nc \
#       --region GL --ymin=-1540000 > ${L}_rel005_compare.txt
#   (script VALIDATED: reproduces T8 I9g6 for IS exactly; rel005 GL =
#   ATL14_GL_0329_100m_005_02.nc + ATL15_GL_0329_01km_005_02.nc, CMR.)
#   DONE 2026-10-01, on the DPS-made products (plan_dps_mosaic.sh D7), read
#   through the mount at .../rel006_0332_testing/north/GL.  The script above
#   was OOM-killed at the ADE's 16 GB; run instead with
#   ~/ATL14_processing/session_tools_2026-10-01/compare_rel005_lowmem.py
#   (float32 grids, views not copies; its IS output is IDENTICAL to the
#   original's).  Output: ${L}_rel005_compare.txt.  STATEMENTS:
#     ATL14, (9301, 14601) common cells, y >= -1540 km:
#       gaps: 27,586 rel005 cells have no new value (0.04%); 453,916 the
#         reverse; 529 of 646,573 1 km blocks lose over half their cells.
#       h new - rel005, N = 62,563,420: median +0.000 m, p5/p95 -0.80/+0.82
#         m; |d| > 10 m on 615,099 (0.98%), max 1972.6 m.  Of those, 573,547
#         have data_count 0 (1.17% of the no-data cells); where data_count
#         > 0, 41,552 of 13,621,763 (0.31%).  Median h_sigma (new) on the
#         > 10 m cells 5.9 m, elsewhere 0.1 m; 264,404 exceed 3x the
#         combined sigma.
#     ATL15 1 km, 29 common epochs, (931, 1461) common cells:
#       gaps: 312 cells with a rel005 value and no new one at some common
#         epoch (0.05%).
#       delta_h new - rel005, N = 18,082,003 cell-epochs: median +0.000 m,
#         p5/p95 -0.08/+0.07 m; |d| > 10 m on 824, max 54.0 m.
#   FOR IS (T8) the same figures were: ATL14 |d| > 10 m on 3.8%, 99% of
#   them in data_count 0 cells.  BEN'S BAR (no > 10 m errors, no major
#   gaps) is his to apply: REPORT ONLY.
#   ALSO 2026-10-01 (Ben asked for two maps; script
#   session_tools_2026-10-01/gl_north_maps.py, output
#   ~/ATL14_processing/runs/GL_0332_north_D7_maps/ with counts.txt):
#     ATL15 delta_h finite -> not finite along time: 1 km, 850 of 632,441
#       cells (735 not finite at the last epoch, 115 finite again; 482 with
#       more than one such step); 10 km, 17 of 7,304.  All at the margins.
#       The finite pattern of delta_h equals that of ice_area in every
#       cell-epoch.  Largest steps (1 km): 385 cells lost at 2021.25 and
#       296 regained at 2021.50; 288 lost at 2022.00; 174 at 2020.25.
#       Median ice_area in the epoch before a loss: 0.66 of a cell.
#     ATL14 h below the EGM2008 geoid (geoid_h, bilinear): 3,247 of
#       63,510,449 cells (0.005%), in 403 1 km blocks, all at the margins;
#       2,105 by more than 5 m, 1,468 by more than 20 m, 646 by more than
#       100 m, lowest -846 m.  h_sigma at each block's lowest cell: median
#       15.2 m, <= 5 m in 79 of the 403 blocks.
# NM3. [ADE] DONE 2026-10-01 (QM-A; result below NM6).  Copy the ATL14 to the reference key; compose the
#      monthly args:
#   ref=<QM-A key>
#   setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
#       default_args/GL_latest.txt default_args/monthly.txt --Hemisphere=1 \
#       --ATL14_reference_file=$ref
#   (paths _monthly: region_dir/s3_run/s3_out/tag/L carry north_monthly)
#   aws s3 cp $region_dir/input_args_GL.txt $s3_run/
#   GATE: the args name the s3:// ref; --hemi_suffix=_monthly; -g 1/12.
# NM4. [DPS] DONE 2026-10-01 (QM-B; result below NM6).  Smoke E80_N-920, --step prelim, 32gb queue.
#   GATE: successful; log names the s3:// reference; N_fit same order as
#   quarterly E80_N-920; record wall time and peak memory.
# NM5. [DPS] DONE 2026-10-01 (QM-C; result below NM6).  Fan-out ALL AT ONCE -- the other 556:
#   grep -vx E80_N-920.h5 $ledgers/GL_0332_north_tile_list.txt > ${L}_north_NM5_tile_list.txt
#   nohup scripts/maap/submit_MAAP_jobs.py --tile_list ${L}_north_NM5_tile_list.txt \
#       --step prelim --args_url $s3_run/input_args_GL.txt \
#       --tile_prefix $s3_out --queue maap-dps-worker-32gb \
#       --tag ${tag}_north_prelim --ledger ${L}_north_prelim_jobs.csv \
#       --rate 0 > ${L}_north_prelim_submit.log 2>&1 &
#   RECORD for MAAP (the point of the test): how long the 556 POSTs take;
#   any submit failures; queue wait (accepted -> started) and the number
#   running over time, from job metrics -- this answers the open question of
#   whether the 32gb queue's own cap is above 100; failures by class and
#   whether they cluster when instances come up.
# NM6. [ADE] DONE 2026-10-01 (result below).  Collect/fetch/check as N4; retry failures (Ben's go).
#   Monthly matched, mosaic and netCDF are NOT in scope.
#
# NM3-NM6 RESULT, 2026-10-01 (Ben: NM2 "passes my bar.  Go ahead with
#   monthly.  The monthly run should submit all tiles at once with retry
#   enabled.").  All on build 0b29127, maap-dps-worker-32gb (r5.xlarge).
#   NM3 DONE.  ref = aws s3 cp (S3 to S3) of the DPS-made north ATL14 to
#     $s3_root/ATL1415/run_args/rel006/north_monthly/GL/ref/
#       ATL14_GL_0332_100m_006_02_north_partial.nc  (582,455,883 bytes ==
#     source); nm_tools.py readback at E80_N-920: h and h_sigma identical
#     to the source, 410,881 finite.  Monthly args (composed 2026-09-30)
#     uploaded to $s3_run (north_monthly); bucket copy identical; sorted
#     diff vs the quarterly args = ref, dzdt_lags, -b, -g only.
#     NOTE: the gate text above says --hemi_suffix=_monthly; the args carry
#     none, as IS's did not (plan_monthly_on_maap M4).
#   NM4 DONE.  Smoke E80_N-920, job cc89aae7: successful, 404 s (fit 319 s,
#     error 70 s), peak 5.42 GiB; nm_tools.py smoke_gate PASS (log names the
#     s3:// ref; N_fit 1,199,471 == quarterly).  Quarterly was 889 s,
#     12.07 GiB.
#   "RETRY ENABLED".  STATEMENT: the submission API has no retry option --
#     maap-py 5.1.0a2 submit_job sends inputs/queue/dedup/tag only; the live
#     https://api.maap-project.org/api/swagger.json names none; maap-api-nasa
#     main (api/endpoints/ogc.py) reads only those four fields.  DECIDED
#     (Ben 2026-10-01): client-side resubmit.  DRIVER:
#     maap_ledgers/GL_0332_monthly_NM5_driver.py (+ .log): round 0 all
#     tiles, --rate 0, no cap; each later round resubmits the tiles without
#     a successful job, new ledger, at most 2 rounds; restartable.
#   NM5 DONE.  556 POSTs in 69 s, none refused (21:53-21:54Z).  555 were
#     running at once by 22:01Z, so the 32gb queue's cap is above 100.
#       round 0  556 jobs: 390 successful, 166 failed (30%); all terminal
#                by 22:33Z
#       round 1  166 jobs (23:00Z): 157 successful, 9 failed
#       round 2    9 jobs (23:22Z):   9 successful
#     ALL_DONE 23:39Z: every tile has a successful job, 1 h 45 min after
#     the first POST (27 min of that is collect_jobs on 556 jobs).
#     Successful jobs: 237-1383 s, median 369 s (round 0), peak RSS median
#     4.8, max 5.74 GiB; 57 job-hours in the 557 successful jobs.  Queue
#     wait (submitted -> job start, 635 jobs with metrics): median 353 s,
#     p95 410 s, max 2017 s.  Metrics: ${L}_north_prelim_metrics.json
#     (L here = maap_ledgers/GL_0332_monthly).
#     FAILURES BY CLASS (get_job_result + triaged _stderr.txt; file
#     GL_0332_monthly_north_prelim_failure_classes.txt):
#       round 0:  95 no logs (empty result record)
#                 61 NSIDC broker call failed -> earthaccess fallback ->
#                    AttributeError 'NoneType' ... 'get_s3_filesystem'
#                    (fit step median 136 s at ~0 CPU before failing)
#                  9 MAAP runner ConnectTimeout to api.maap-project.org
#                    /api/environment/config, before our container
#                  1 result record unreadable (HTTP 500 from the API)
#       round 1:   6 NSIDC broker call failed; 3 runner ConnectTimeout
#       NoCredentialsError (workspace bucket): 0 in 731 jobs.
#     INFERRED (timing only): the NSIDC failures are connect timeouts to
#       the MAAP API when hundreds of jobs start together; the reason is
#       not logged (pointCollection/ps_scale_for_lat.py line 3 silences
#       the warning that carries it).  Same class as N3's 10.
#   NM6 DONE.  557 tiles and 557 field-size reports at
#     $s3_root/ATL14_processing/rel006/north_monthly/GL/prelim (87.1 GB);
#     check_field_sizes.py through the mount: expected dz [25, 25, 94],
#     557 of 557 passed, 0 problems.  Nothing fetched to /home.

# NM7. [DPS] IN PROGRESS 2026-10-02.  MONTHLY MATCHED, GL north, all at once.
#   DECIDED (Ben 2026-10-02): "Run matched - see if the updates to the
#   credential passing has reduced the failure rate.  Submit all jobs at
#   once."  This revises NM6's "monthly matched ... NOT in scope".
#   Build e6d7051 (workspace keys from MAAP's broker, worker role off;
#   pointCollection 2371978: NSIDC broker retry).  32gb queue, as NM5, so
#   the submission pattern and queue are the same as the round it is
#   compared with (NM5 round 0: 556 jobs, 166 failed = 30%).
#   STATEMENT: a matched job reads no ATL11, so it makes NO NSIDC broker
#     call -- NM5's 61 NSIDC failures cannot recur here whatever the code
#     does.  What this run tests is the NEW workspace-credentials call
#     (one per job, the same MAAP API, 5 tries) and the two classes that
#     are not ours: MAAP runner timeouts (9 in NM5) and jobs with no logs
#     (95).  The NSIDC retry gets its test at the next prelim fan-out.
#   a. smoke E80_N-920 (full 3x3 of neighbours), as QM-B: a fault common to
#      every job would spoil the all-at-once test.  GATE: successful, tile
#      at <monthly prefix>/matched/, field size [25, 25, 94].
#   b. the other 556 at once: maap_ledgers/GL_0332_monthly_NM7_driver.py
#      (the NM5 driver with --step matched; automatic resubmit, 2 rounds),
#      ledgers GL_0332_monthly_north_matched[_retryN]_jobs.csv.
#   c. RECORD: first-round failures by class vs NM5; "h left at exit";
#      time and memory; check_field_sizes --step matched through the mount.

# ===========================================================================
# N8. [ADE] DRAFT WRITTEN 2026-10-01, AWAITS BEN'S REVIEW.  Numbers for MAAP.
# ===========================================================================
# From the collect files: per-job time, peak memory, instance mix, queue
# wait, failures; ADE time for N4/N5 passes; mosaic and netCDF time/memory.  Replace
# the GL rows of ~/ATL14_processing/maap_resource_estimate.txt (prelim
# measured at 557, matched measured instead of scaled from IS) -- Ben
# reviews before anything is sent.
# MOSAIC AND NETCDF, MEASURED 2026-10-01 on DPS (plan_dps_mosaic.sh D7; all
# on maap-dps-worker-16gb, t3/t3a.xlarge; collect files
# maap_ledgers/GL_0332_north_D7_*_collect.txt).  GL north = 557 tiles:
#   200 km tiles  32 jobs, 186-818 s each (median ~590 s), peak 0.80 GiB;
#                 4.9 job-hours for the 32 that succeeded, plus 1.9
#                 job-hours lost in 14 jobs that failed on credentials
#                 (before the run.sh retry settings)
#   mosaics       41 jobs, 33-504 s (median 60 s; z0 504 s), peak 1.65 GiB
#                 (z0); 0.87 job-hours
#   netCDF        ATL14 1155 s, peak 4.37 GiB; ATL15 3752 s, peak 2.48 GiB;
#                 1.36 job-hours
#   output        11.2 GB at the out_prefix (1358 objects), of which the
#                 five netCDFs are 1.28 GB
# DRAFT 2026-10-01: ~/ATL14_processing/maap_resource_estimate.txt updated
# (previous text kept as maap_resource_estimate_2026-09-25.txt).  Changed
# lines carry [GL north]: GL job-hours (quarterly ~320, monthly ~200 with
# matched still estimated, mosaic + netCDF ~20 per product), GL memory and
# times measured over 557 tiles, GL tile sizes, products now built on DPS,
# and a new RELIABILITY section (6% / 14% / 30% first-round failures at
# ~100 in flight / ~100 in flight / all at once).  Scaling to all of GL is
# x 1,483/557 -- an ASSUMPTION that the north is representative.
# NOT DONE: monthly matched for GL has never been run (out of scope, NM6);
# its row stays an estimate.  Ben reviews before anything is sent.
