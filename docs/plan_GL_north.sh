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
#        QN1 answer:
#   QN2  How far.  RECOMMENDATION: prelim + matched + mosaic (N3-N6); no
#        netCDF, no monthly.  The DPS load MAAP asked about is prelim +
#        matched; the mosaic is ADE-only but is the one GL-scale step never
#        timed.  A netCDF of a third of GL would be mostly empty.
#        QN2 answer:
#   QN3  Concurrency.  RECOMMENDATION: --max_in_flight 100 and ask MAAP
#        whether they want it higher -- how the cluster behaves at 100+
#        concurrent jobs is part of what they want to learn.  IS ran 29 at
#        once; nothing larger has run.
#        QN3 answer:
#   QN4  Output prefix.  RECOMMENDATION: the CANONICAL GL prefix
#        (.../ATL14_processing/rel006/north/GL), so these tiles count toward
#        the full run, which then submits only the remaining centers.
#        Caveat: the southernmost matched row (y = -1520 km) is matched
#        without its southern neighbours and must be re-run as part of the
#        full run.  The alternative is a dated test prefix (as the transects
#        used), which throws the tiles away afterwards.
#        QN4 answer:
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
# N0. [ADE] BLOCKED on Ben.  Register, then prove the build.
# ===========================================================================
# Ben registers (checkout must be clean and pushed).  The new image installs
# pointCollection and LSsurf main at build time (pyproject.toml is unpinned),
# which is what brings in the read fix.
/srv/conda/envs/notebook/bin/python register_algorithm.py --dry-run
scripts/maap/check_build_id.py $s3_run/input_args_GL.txt maap-dps-worker-32gb
# GATE: VERDICT MATCH, maap_pgt=set, image commit >= this plan's commit.
# Run check_build_id WITHOUT a short timeout -- it waits on the queue.


# ===========================================================================
# N1. [ADE] TODO.  The subset list (in ledgers, NOT in the checkout).
# ===========================================================================
sed 's/.*_N\(-\{0,1\}[0-9]*\).*/\1 &/' ATL1415/resources/GL/40km_tile_list.txt \
    | awk -v ymin=$ymin '$1 >= ymin {print $2}' > ${L}_tile_list.txt   # mawk: no match() arrays
wc -l ${L}_tile_list.txt                      # expect 557 for ymin=-1520


# ===========================================================================
# N2. [DPS] TODO.  Smoke on the new build: two transect tiles re-run.
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
# N3. [DPS] TODO.  Prelim fan-out.
# ===========================================================================
nohup scripts/maap/submit_MAAP_jobs.py --tile_list ${L}_tile_list.txt \
    --step prelim --args_url $s3_run/input_args_GL.txt \
    --tile_prefix $s3_out --queue maap-dps-worker-32gb \
    --tag ${tag}_prelim --ledger ${L}_prelim_jobs.csv \
    --max_in_flight 100 > ${L}_prelim_submit.log 2>&1 &      # QN3
# NEVER register while this runs.


# ===========================================================================
# N4. [ADE] TODO.  Collect, fetch, check.
# ===========================================================================
scripts/maap/collect_jobs.py ${L}_prelim_jobs.csv > ${L}_prelim_collect.txt
scripts/maap/fetch_tiles.py  ${L}_prelim_jobs.csv $region_dir --step prelim
scripts/check_field_sizes.py $region_dir/prelim @$region_dir/input_args_GL.txt
# Time these three passes too: ADE-side cost per job is part of the answer
# for MAAP (they make API/S3 calls per job and have only run at 29).
# Failed jobs: retry per arctic howto step 5, NEW ledger.


# ===========================================================================
# N5. [DPS+ADE] TODO.  Matched, then collect/fetch/check as N4.
# ===========================================================================
nohup scripts/maap/submit_MAAP_jobs.py --tile_list ${L}_tile_list.txt \
    --step matched --args_url $s3_run/input_args_GL.txt \
    --tile_prefix $s3_out --queue maap-dps-worker-32gb \
    --tag ${tag}_matched --ledger ${L}_matched_jobs.csv \
    --max_in_flight 100 > ${L}_matched_submit.log 2>&1 &
# The y = -1520 km row lacks southern neighbours (QN4 caveat).


# ===========================================================================
# N6. [ADE] TODO (QN2).  Mosaic, timed.
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
# N7. [ADE] TODO.  Numbers for MAAP.
# ===========================================================================
# From the collect files: per-job time, peak memory, instance mix, queue
# wait, failures; ADE time for N4/N5 passes; mosaic time/memory.  Replace
# the GL rows of ~/ATL14_processing/maap_resource_estimate.txt (prelim
# measured at 557, matched measured instead of scaled from IS) -- Ben
# reviews before anything is sent.
