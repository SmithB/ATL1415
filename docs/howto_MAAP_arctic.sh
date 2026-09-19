# howto_MAAP_arctic.sh -- the arctic regions (RA IS CN CS SV) on MAAP
#
# ############################################################################
# ##  REWRITTEN 2026-09-19 from what the Iceland run learned.               ##
# ##  IS HAS RUN EVERY STEP, quarterly and monthly, at cycles 03-32         ##
# ##  (plan_cycles_03_32.sh, plan_monthly_on_maap.sh, plan_tile_lists.sh).  ##
# ##  RA, CN, CS AND SV HAVE NEVER RUN ON MAAP: for them every step is      ##
# ##  TENTATIVE until it has, and the five-region loops below have not been ##
# ##  run as loops.  Revise the tags as regions go through.                 ##
# ############################################################################
#
# THE ORDER IS docs/howto_arctic.sh's (the discover/SLURM variant, which is
# still the production path and is NOT replaced by this file): every region's
# prelim, then every region's matched, then every region's mosaic and netCDF,
# then all of it again for the monthly product.  discover runs each phase with
# scripts/run_arctic_*.sh; here each phase is a loop over existing tools, so no
# MAAP variant of those scripts is needed.
#
# Tags, per step:  [ADE] runs in the ADE, [DPS] submits jobs.
#   [OK on IS]    has run on Iceland as written (cited)
#   [UNTESTED]    has not run in this form
# Steps are numbered 0-20 so they can be cited ("arctic step 7").
#
# Prerequisite: docs/howto_MAAP_staging.sh S1-S5 (ADE env, masks and ATL11
# index on the bucket, the algorithm registered).  All ten arctic masks -- the
# RGI .db files and their _40km.tif siblings -- are staged (checked
# 2026-09-19).
#
# WHAT IS DIFFERENT FROM DISCOVER, in one place:
#   a. TILES ARE SOLVED ON DPS, one job per tile, submitted from the ADE by
#      scripts/maap/submit_MAAP_jobs.py.  Every job is written to a LEDGER
#      (csv, job ids); collect_jobs.py and fetch_tiles.py read it back.
#      KEEP LEDGERS OUTSIDE THE CHECKOUT -- an untracked file there makes
#      register_algorithm.py refuse.
#   b. THE TILE CENTERS ARE ATL1415/resources/<REG>/40km_tile_list.txt, not
#      make_ATL1415_queue.py (plan_tile_lists.sh).  Centers with no data come
#      out of the list by hand after a run (step 20).
#   c. EACH JOB PUTS ITS TILE AT <s3_out>/prelim/ (--tile_prefix), which is
#      how a matched job finds its neighbours.  The ADE copy is fetched.
#   d. MOSAIC AND netCDF RUN IN THE ADE.  make_mosaic_jobs.py's slurm_run.sh
#      is plain bash (the #SBATCH lines are comments), so it runs locally with
#      SLURM_ARRAY_TASK_ID set.  No sbatch, no run_queue_local.sh.
#   e. NOTHING DELETES THE REGION DIRECTORY.  discover's run_arctic_prelim.sh
#      does `rm -r $base` first; here the ADE tree holds fetched tiles and
#      products that exist nowhere else locally.  Setup only writes the args
#      file.
#   f. DPS RUNS THE REGISTERED BUILD, NOT THIS CHECKOUT.  Never register while
#      a fan-out is queued: IS's first 29 silently split across two builds
#      that way (plan_IS_run.sh I3).  collect_jobs.py's "Builds that ran these
#      tiles" line is the check.

conda activate ATL14
cd ~/git_repos/ATL1415
repo=$PWD
regions="RA IS CN CS SV"      # at 0332 IS is DONE: use regions="RA CN CS SV"
rel_file=default_args/latest_release.txt
rel=$(grep '^--Release=' $rel_file | cut -d= -f2)     # 006
cyc=$(grep '^--cycles=' $rel_file | cut -d= -f2)      # 0332
ver=$(grep '^--version=' $rel_file | cut -d= -f2)     # 02
tspan=$(grep '^-t=' $rel_file | cut -d= -f2)          # 2018.75,2026.5
ATL14_root=/home/jovyan/ATL14_processing
s3_root=s3://maap-ops-workspace/ben_smith
ledgers=$ATL14_root/maap_ledgers
runs=$ATL14_root/runs
mkdir -p $ledgers $runs

# paths <REG> [_monthly] -- every per-region path, from one place
paths () {
    region_dir=$ATL14_root/rel$rel/north$2/$1
    s3_run=$s3_root/ATL1415/run_args/rel$rel/north$2/$1
    s3_out=$s3_root/ATL14_processing/rel$rel/north$2/$1
    tile_list=$repo/ATL1415/resources/$1/40km_tile_list.txt
    tag=$1_rel${rel}_${cyc}$2
    L=$ledgers/$1_${cyc}$2          # ledger stem: ${L}_prelim_jobs.csv ...
}


# ###########################################################################
# SETUP
# ###########################################################################

# ===========================================================================
# 0. [ADE] [OK on IS]  The DPS build is the commit you mean to run.
# ===========================================================================
# Register only if solver code (ATL1415/, run.sh, the image) changed since the
# last build; docs, tile lists and scripts/maap/ are ADE-side and need none.
# Either way, CHECK the image before spending jobs on it.  The checkout must be
# clean and pushed first, or register_algorithm.py refuses (by design).
/srv/conda/envs/notebook/bin/python register_algorithm.py --dry-run   # "push check: OK"
/srv/conda/envs/notebook/bin/python register_algorithm.py
scripts/maap/check_build_id.py $s3_root/ATL1415/run_args/rel006/north/IS/input_args_IS.txt \
    maap-dps-worker-16gb
# MUST say "VERDICT: MATCH" with maap_pgt=set (without it a worker cannot read
# ATL11 at all).  With no --expect it tests the image against the commit the
# build service recorded -- the right test.  Pass --expect <sha> only for the
# sha you REGISTERED, never HEAD: after a docs-only push HEAD is past the build
# ("origin is now past this build" is normal) and would read as a MISMATCH.
# check_build_id.py SUBMITS a small job; it has no --help (it refuses stray
# options since 7d64824 -- once, --help cost a real job).


# ===========================================================================
# 1. [ADE] [OK on IS]  Point the release symlink at this release.
# ===========================================================================
ln -sf rel_006_0332.txt default_args/latest_release.txt
# The ATL11 generation must be the one CMR serves: 0331_007_04 is gone
# (plan_cycles_03_32.sh), and its release file cannot be solved.


# ===========================================================================
# 2. [ADE] [OK on IS; UNTESTED for RA CN CS SV]  Compose and publish the args.
# ===========================================================================
for reg in $regions; do
    paths $reg
    setup_ATL1415_region.py default_args/MAAP_dps.txt $rel_file \
        default_args/$reg.txt default_args/quarterly.txt --Hemisphere=1
    aws s3 cp $region_dir/input_args_$reg.txt $s3_run/
done
# CHECK: every cloud input in the composed file is an s3:// URI or a CMR
# search (--ATL11_earthaccess, --previous_product_earthaccess,
# --previous_product=005_0329); -b is the one local path, and run.sh
# overrides it on the worker.  SV's region file sets --DEM_tol=200 (50
# elsewhere).  RECOMPOSE AND REPUBLISH after any change to default_args/ --
# DPS reads the bucket copy, and nothing reconciles the two.
for reg in $regions; do paths $reg
    aws s3 cp $s3_run/input_args_$reg.txt - | diff - $region_dir/input_args_$reg.txt \
        && echo "$reg: bucket copy identical"
done


# ###########################################################################
# PRELIM -- every region
# ###########################################################################

# ===========================================================================
# 3. [DPS] [OK on IS; UNTESTED for RA CN CS SV]  One smoke tile per region.
# ===========================================================================
# A region that has never run gets ONE job first: it proves the args, the
# mask, the ATL11 read and the crossovers on a worker for one job's cost.
# --limit 1 takes the list's first center.  ALWAYS --dry-run first.
for reg in $regions; do paths $reg
    scripts/maap/submit_MAAP_jobs.py --tile_list $tile_list --limit 1 \
        --step prelim --args_url $s3_run/input_args_$reg.txt \
        --tile_prefix $s3_out --queue maap-dps-worker-16gb \
        --tag ${tag}_smoke --ledger ${L}_smoke_jobs.csv --dry-run
done
# (the same loop without --dry-run), then, when they finish:
for reg in $regions; do paths $reg; scripts/maap/collect_jobs.py ${L}_smoke_jobs.csv; done
# GATES: successful; N_XO > 0 (crossovers are read); commit = the step 0 build.
# A tile with no data is ALSO successful, with no tile -- a legitimate result,
# but then it proves nothing: smoke the next center (--xy_file with that one
# center, in a file under $ledgers, not in the checkout).


# ===========================================================================
# 4. [DPS] [OK on IS; UNTESTED for RA CN CS SV]  Fan out, every region.
# ===========================================================================
for reg in $regions; do paths $reg
    scripts/maap/submit_MAAP_jobs.py --tile_list $tile_list \
        --step prelim --args_url $s3_run/input_args_$reg.txt \
        --tile_prefix $s3_out --queue maap-dps-worker-16gb \
        --tag ${tag}_prelim --ledger ${L}_prelim_jobs.csv
done
# The smoke centers re-run; the solve is deterministic (IS: bit-identical).
# An existing ledger is never overwritten -- a region already submitted stops
# with an error and the loop goes on to the next.
# QUEUE: IS quarterly prelim peaked at 9.52 GiB of 16 (memory tracks N_fit,
# not tile count).  The other regions are unmeasured; an OOM is a failed job
# -- resubmit those centers on maap-dps-worker-32gb (step 5).
# THROUGHPUT, measured on IS 2026-09-18: 29 jobs submitted within one minute
# all finished within ~45 min, each running 8-18 min -- they ran largely in
# parallel, and the "~10 jobs/hr" limit older notes feared was not seen.
# Unmeasured at hundreds of jobs; --max_in_flight N holds the queue at N.


# ===========================================================================
# 5. [ADE] [OK on IS]  Watch.
# ===========================================================================
for reg in $regions; do paths $reg; echo "== $reg"; scripts/maap/collect_jobs.py ${L}_prelim_jobs.csv; done
# Per tile: status, wall, peak RSS, N_ATL11 / N_AT / N_XO, N_fit, iterations,
# and the build that ran it.  IS quarterly: 1332-3322 s, 4.4-9.5 GiB.
# A FAILED JOB'S TILE IS UNRECOVERABLE -- it gets no dps_output prefix, only a
# triaged_job log (s3://maap-ops-workspace/dataset/triaged_job/...).  Read
# it before resubmitting: a deterministic failure fails again.  To retry
# centers, write them ("x0 y0" per line) to a file under $ledgers and submit
# with --xy_file and a NEW --ledger (IS: one matched job lost to MAAP's own
# input-staging timeout, rerun cleanly).


# ===========================================================================
# 6. [ADE] [OK on IS]  Fetch and check the prelim tiles.
# ===========================================================================
for reg in $regions; do paths $reg; echo "== $reg"
    scripts/maap/fetch_tiles.py ${L}_prelim_jobs.csv $region_dir --step prelim
    scripts/check_field_sizes.py $region_dir/prelim @$region_dir/input_args_$reg.txt
done
# fetch_tiles: "fetched" / "have it" per tile; FAILED jobs under NOT FETCHED;
# and a NO DATA block -- successful jobs that left no tile -- whose names it
# adds to $region_dir/prelim/no_data_tiles.txt for step 20.  Nothing else is
# done with them during the run.
# check_field_sizes: dz/dz [61, 61, 32] for -W 60000 at quarterly -g; prelim
# sigma_dz == dz/dz; a report for every tile.  Exit 0 OK, 1 problems, 2 the
# check did not happen.


# ###########################################################################
# MATCHED -- every region
# ###########################################################################

# ===========================================================================
# 7. [DPS] [OK on IS]  Submit matched, every region.
# ===========================================================================
for reg in $regions; do paths $reg
    scripts/maap/submit_MAAP_jobs.py --tile_list $tile_list \
        --step matched --args_url $s3_run/input_args_$reg.txt \
        --tile_prefix $s3_out --queue maap-dps-worker-16gb \
        --tag ${tag}_matched --ledger ${L}_matched_jobs.csv
done
# The SAME tile list.  The submitter lists $s3_out/prelim/ once and SKIPS, by
# name, every center with no prelim tile -- no-data centers and failed
# prelims -- so matched runs exactly where a prelim tile exists.  A matched
# job fetches its own and its 8 neighbours' prelim tiles from --tile_prefix;
# a missing neighbour is logged, not fatal (IS is sparse: few centers have a
# full 3x3, and "k/9 localized" well under 9 is normal there).
# IS quarterly matched: 145-830 s, 3.6-9.4 GiB, 1 iteration.


# ===========================================================================
# 8. [ADE] [OK on IS]  Watch, fetch and check matched.
# ===========================================================================
for reg in $regions; do paths $reg; echo "== $reg"
    scripts/maap/collect_jobs.py ${L}_matched_jobs.csv
    scripts/maap/fetch_tiles.py ${L}_matched_jobs.csv $region_dir --step matched
    scripts/check_field_sizes.py $region_dir/matched @$region_dir/input_args_$reg.txt
done
# Matched tiles have NO sigma_dz, by design: the uncertainties come from the
# prelim tiles.  A matched job with no tile is NOT normal (fetch_tiles lists
# it under NOT FETCHED).


# ###########################################################################
# MOSAIC AND netCDF -- every region
# ###########################################################################

# ===========================================================================
# 9. [ADE] [OK on IS]  Mosaic, every region.
# ===========================================================================
# Run directories go under $runs, OUTSIDE the checkout.  -e ATL14: the
# default, IS2, is discover's env.  IS quarterly: 41 tasks, ~0.25 GiB each,
# 59 s at -P 12 on the ADE's 16 cores.  Larger regions are untimed.
cd $runs
for reg in $regions; do paths $reg
    make_mosaic_jobs.py -b $region_dir -rr $reg -t $tspan -e ATL14 \
        --run_name ${reg}_${cyc}_mosaic @$repo/default_args/quarterly.txt
    ( cd ${reg}_${cyc}_mosaic
      seq 1 $(ls queue | wc -l) | xargs -P 12 -I{} env SLURM_ARRAY_TASK_ID={} bash slurm_run.sh )
    check_mosaic_outputs.py $runs/${reg}_${cyc}_mosaic --values
done
cd $repo
# EXIT CODES ARE NOT ENOUGH: check_mosaic_outputs.py --values reads every
# field and flags all-NaN ones.  Then error_logs/ must be empty and done/
# must hold every task.  Values come from matched/, sigmas from prelim/.
# KNOWN, NOT FIXED: on IS quarterly, sigma_dzdt is finite on only ~44% of the
# cells dzdt is (plan_cycles_03_32.sh T8).  Flagged, not pursued.


# ===========================================================================
# 10. [ADE] [OK on IS]  netCDF, every region: ATL14 and ATL15.
# ===========================================================================
for reg in $regions; do paths $reg
    mkdir -p $runs/${reg}_${cyc}_nc
    ( cd $runs/${reg}_${cyc}_nc
      ATL14_write2nc.py @$region_dir/input_args_$reg.txt > ATL14.log 2>&1
      ATL15_write2nc.py @$region_dir/input_args_$reg.txt > ATL15.log 2>&1 )
done
# -> ATL14_<REG>_<cyc>_100m_<rel>_<ver>.nc and ATL15_<REG>_<cyc>_3mo_{1,10,20,40}km_...
# in $region_dir.  IS: 12 s + 18 s.  EXPECT NO "INVALID" LINE in either log:
# lineage comes from the prelim tiles (plan_lineage_at_solve_time.sh), and an
# INVALID warning means a tile solved on a build without it.  The only
# NOT_SET lineage values are on the ATL11XO rows, in start/end_orbit and
# start/end_region, which the XO granules do not carry -- correct.
# ATL15 has 31 epochs, 2019.0-2026.5 (--t_crop); the tiles have 32.
# TO LOOK AT THEM with GDAL use the notebook env: ATL14's GDAL has no netCDF
# or HDF5 plugin.
#
# COMPARE WITH THE PREVIOUS PRODUCT before building anything on it -- Ben's
# bar: no >10 m errors, no major gaps.  [NO SOFTWARE: IS's was a scratch
# script -- plan_cycles_03_32.sh T8 "I9g6" has the method and numbers.]


# ###########################################################################
# MONTHLY -- every region, again  (plan_monthly_on_maap.sh)
# ###########################################################################
# The SAME pipeline with default_args/monthly.txt (-g=1250,2500,1/12, the
# monthly --dzdt_lags, hemi suffix _monthly) and the region's QUARTERLY ATL14
# as a reference DEM, subtracted from every point.  Only ATL15 is written.
# PREREQUISITE: steps 0-10 done and the quarterly product accepted, per region.
# IS monthly: prelim 496-1071 s at 1.3-4.1 GiB, matched 85-307 s at
# 1.1-3.9 GiB -- ~3x faster and ~2.3x lighter than quarterly.

# ===========================================================================
# 11. [ADE] [OK on IS]  Publish each region's quarterly products.
# ===========================================================================
for reg in $regions; do paths $reg
    for f in $region_dir/ATL1[45]_${reg}_${cyc}_*_${rel}_${ver}.nc; do aws s3 cp $f $s3_out/; done
done
# Beside the quarterly tiles.  The ATL14 is the reference the solver reads by
# URI; the ATL15s go up so the bucket holds the whole product.  Check each
# bucket size against its local file.  RE-WRITE AND REPUBLISH BEFORE MONTHLY
# if anything in the quarterly products changes: the reference must be final.


# ===========================================================================
# 12. [ADE] [OK on IS]  Compose and publish the monthly args.
# ===========================================================================
for reg in $regions; do paths $reg
    ref=$s3_out/ATL14_${reg}_${cyc}_100m_${rel}_${ver}.nc     # quarterly s3_out
    setup_ATL1415_region.py default_args/MAAP_dps.txt $rel_file \
        default_args/$reg.txt default_args/monthly.txt --Hemisphere=1 \
        --ATL14_reference_file=$ref
    paths $reg _monthly
    aws s3 cp $region_dir/input_args_$reg.txt $s3_run/
done
# The monthly args differ from the quarterly in exactly four lines: -g,
# --dzdt_lags, --ATL14_reference_file and -b (-> rel006/north_monthly/<REG>).
# One URI, no wildcard: a URI with a wildcard raises.


# ===========================================================================
# 13. [DPS] [OK on IS]  Monthly smoke, one tile per region.
# ===========================================================================
# As step 3 with "paths $reg _monthly" (tag and ledger then carry _monthly).
# THE GATE THAT MATTERS: N_fit about the same as that center's QUARTERLY
# N_fit (IS E1340_N-2460: 280955 vs 273382).  An unreadable reference does
# not raise -- every point becomes invalid -- so a collapsed N_fit is the
# only sign.  Also: dz/dz [25, 25, 94] (check_field_sizes reads 1/12).


# ===========================================================================
# 14. [DPS] [OK on IS]  Monthly prelim, every region -- steps 4-6.
# ===========================================================================
for reg in $regions; do paths $reg _monthly
    scripts/maap/submit_MAAP_jobs.py --tile_list $tile_list \
        --step prelim --args_url $s3_run/input_args_$reg.txt \
        --tile_prefix $s3_out --queue maap-dps-worker-16gb \
        --tag ${tag}_prelim --ledger ${L}_prelim_jobs.csv
done
# then the loops of steps 5 and 6 with "paths $reg _monthly".
# The SAME tile list serves both periods.  A center the quarterly ATL14 does
# not cover has an all-NaN reference and no data: a successful job, no tile
# (since plan_tile_lists.sh TL1), saved to no_data_tiles.txt like any other.


# ===========================================================================
# 15. [DPS] [OK on IS]  Monthly matched, every region -- steps 7-8.
# ===========================================================================
# The loops of steps 7 and 8 with "paths $reg _monthly".


# ===========================================================================
# 16. [ADE] [OK on IS]  Monthly mosaic, every region.
# ===========================================================================
cd $runs
for reg in $regions; do paths $reg _monthly
    make_mosaic_jobs.py -b $region_dir -rr $reg -t $tspan -e ATL14 \
        --run_name ${reg}_${cyc}_monthly_mosaic @$repo/default_args/monthly.txt
    cat ${reg}_${cyc}_monthly_mosaic/queue/* | grep -oE "lag[0-9]+" | sort -t g -k2 -n -u | tr '\n' ' '; echo
done
# INSPECT THAT OUTPUT BEFORE RUNNING: make_mosaic_jobs.py INFERS the dzdt lags
# from -t and -g rather than reading --dzdt_lags; each region must print
# lag1 lag3 lag6 lag12 lag24 lag36 lag48 lag60 lag72 lag84.  No z0 task
# (skip_z0: z0 spacing over 1000 m).  IS: 44 tasks, 74 s.  Then run and check
# each run directory exactly as step 9.
cd $repo


# ===========================================================================
# 17. [ADE] [OK on IS]  Monthly netCDF, every region: ATL15 ONLY.
# ===========================================================================
for reg in $regions; do paths $reg _monthly
    mkdir -p $runs/${reg}_${cyc}_monthly_nc
    ( cd $runs/${reg}_${cyc}_monthly_nc
      ATL15_write2nc.py @$region_dir/input_args_$reg.txt > ATL15.log 2>&1 )
done
# -> ATL15_<REG>_<cyc>_1mo_{2.5,10,20,40}km_<rel>_<ver>.nc, 91 epochs.
# COMPARING WITH QUARTERLY: raw delta_h IS comparable -- both are 0 at the
# 2020 reference with sigma 0 (IS: median difference -0.05 m).  Do NOT
# difference against the first epoch; it spreads that epoch's errors across
# every other.  Expect a SEASONAL difference of a few tenths of a metre (IS:
# Apr -0.37 m, Oct +0.31 m), so many values exceed 3x sigma without anything
# being wrong; >10 m differences on IS sat in 6 weakly-constrained cells at
# the ends of the record.


# ===========================================================================
# 18. [ADE] [OK on IS, at Ben's request]  Publish the monthly products.
# ===========================================================================
for reg in $regions; do paths $reg _monthly
    for f in $region_dir/ATL15_${reg}_${cyc}_1mo_*_${rel}_${ver}.nc; do aws s3 cp $f $s3_out/; done
done


# ###########################################################################
# AFTER THE RUN
# ###########################################################################

# ===========================================================================
# 19. [ADE] [SUGGESTION, NO SOFTWARE]  Annotate the run's build history.
# ===========================================================================
# The git history across the builds the run's tiles ran, with what changed
# and which tiles ran on which build: howto_MAAP_ogc.sh O12b.


# ===========================================================================
# 20. [ADE] [OK on IS]  Take the no-data centers out of the tile lists.
# ===========================================================================
# Both periods' lists, merged; then commit AND PUSH -- an unpushed commit
# blocks register_algorithm.py.  Never mid-run: matched already skips them.
for reg in $regions; do
    tile_list=$repo/ATL1415/resources/$reg/40km_tile_list.txt
    cat $ATL14_root/rel$rel/north{,_monthly}/$reg/prelim/no_data_tiles.txt 2>/dev/null \
        | sort -u > $ledgers/${reg}_no_data.txt
    [ -s $ledgers/${reg}_no_data.txt ] || continue
    grep -vxFf $ledgers/${reg}_no_data.txt $tile_list > $ledgers/t && mv $ledgers/t $tile_list
done
git diff --stat ATL1415/resources/
git commit -m "Drop no-data centers from the arctic tile lists" ATL1415/resources/ && git push
