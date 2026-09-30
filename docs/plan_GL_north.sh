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
# N6. [ADE] TODO.  Mosaic, timed.
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
# N7. [ADE] TODO (QN2).  netCDF, timed.  LOCAL ONLY -- do not publish.
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
# NM0. [ADE] TODO.  N5 finish: collect/fetch/check matched (main + retry1).
# NM1. [ADE] TODO.  N6 mosaic, N7 ATL14 (+ATL15) netCDF, as written above.
# NM2. [ADE] TODO (QM-A).  Compare ATL14 with rel005 over the north; Ben's
#      bar.  GATE: Ben passes it.
# NM3. [ADE] TODO (QM-A).  Copy the ATL14 to the reference key; compose the
#      monthly args:
#   ref=<QM-A key>
#   setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
#       default_args/GL_latest.txt default_args/monthly.txt --Hemisphere=1 \
#       --ATL14_reference_file=$ref
#   (paths _monthly: region_dir/s3_run/s3_out/tag/L carry north_monthly)
#   aws s3 cp $region_dir/input_args_GL.txt $s3_run/
#   GATE: the args name the s3:// ref; --hemi_suffix=_monthly; -g 1/12.
# NM4. [DPS] TODO (QM-B).  Smoke E80_N-920, --step prelim, 32gb queue.
#   GATE: successful; log names the s3:// reference; N_fit same order as
#   quarterly E80_N-920; record wall time and peak memory.
# NM5. [DPS] TODO (QM-C).  Fan-out ALL AT ONCE -- the other 556:
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
# NM6. [ADE] TODO.  Collect/fetch/check as N4; retry failures (Ben's go).
#   Monthly matched, mosaic and netCDF are NOT in scope.


# ===========================================================================
# N8. [ADE] TODO.  Numbers for MAAP.
# ===========================================================================
# From the collect files: per-job time, peak memory, instance mix, queue
# wait, failures; ADE time for N4/N5 passes; mosaic and netCDF time/memory.  Replace
# the GL rows of ~/ATL14_processing/maap_resource_estimate.txt (prelim
# measured at 557, matched measured instead of scaled from IS) -- Ben
# reviews before anything is sent.
