# howto_MAAP_AA.sh -- Antarctica on MAAP (per-tile solves on DPS)
#
# ############################################################################
# ##  REWRITTEN 2026-09-19 from what the Iceland run learned.               ##
# ##  AA HAS NOT RUN IN PRODUCTION ON MAAP.  EVERY STEP IS TENTATIVE.       ##
# ##  What HAS run: a 17-job cost transect on DPS (2026-09-11, build        ##
# ##  ab84687) -- tides, crossovers and the pole hole all worked; numbers   ##
# ##  in scripts/maap/AA_cost_results.csv.  The per-tile procedure is the   ##
# ##  one IS ran end to end (docs/howto_MAAP_arctic.sh explains each step   ##
# ##  and cites the IS record); this file carries it for AA's two halves    ##
# ##  and adds the stages only AA has.                                      ##
# ############################################################################
#
# The discover/SLURM variant is docs/howto_AA.sh, which is still the
# production path and is NOT replaced by this file.
#
# Tags:  [ADE] / [DPS];  [OK on IS] the same command ran for Iceland;
#   [OK, transect] ran in the cost transect;  [UNTESTED];  [NEEDS CODE: x].
# Steps are numbered 0-17 so they can be cited ("AA step 9").
#
# WHAT MAKES AA DIFFERENT, in one place:
#   a. TWO HALVES WITH DIFFERENT TILE WIDTHS: 60 km north of the 400 km line,
#      44 km south of it, overlapping ON PURPOSE (Ben, 2026-09-08): a center
#      with max(|x|,|y|) >= 360 km belongs to the 60 km half, one with
#      max(|x|,|y|) <= 440 km to the 44 km half, so the 240 centers between
#      are solved at BOTH widths.  DO NOT "FIX" THE LIMITS.
#      Each half has its own region directory (AA, AA_44km), args file and
#      bucket prefix.  Because a matched job reads its neighbours from its
#      own --tile_prefix, keeping the prefixes separate is what keeps a
#      matched neighbourhood from mixing 44 km and 60 km fits under identical
#      file names -- which would give a solved tile, not an error.
#   b. THE 44 km ARGS FILE KEEPS --region=AA.  It is DERIVED from the 60 km
#      file by scripts/maap/make_AA_44km_args.py (-W and -b only).  Composing
#      it with --region=AA_44km broke all four 44 km jobs on 2026-09-08:
#      ATL11_to_ATL15.py loads the gridded mask only for region AA or GL.
#   c. SCALE: 8944 list centers -> 9184 jobs (8724 at 60 km, 460 at 44 km,
#      the 240 overlap centers twice).  Tools built and tested at 29 jobs.
#   d. MEMORY: the transect's worst tile, 44km E220_N20 at the pole-hole
#      edge, peaked at 21.47 GiB; every other tile <= 18.92 GiB.  Use the
#      32 GiB queue.  Roughly an hour per tile (median fit+error 0.88 h).
#   e. AFTER THE TILES, THREE EXTRA ADE STAGES: 200 km tiles per half, the
#      four sectors A1-A4, and a mosaic and netCDF per sector.  WHETHER THE
#      ADE CAN MOSAIC AA IN REASONABLE TIME IS UNMEASURED (IS's mosaic was
#      trivial; AA is ~300x the tiles).
#   f. MONTHLY IS BLOCKED on a reference DEM that spans four sector files
#      (step 16, NEEDS CODE).

conda activate ATL14
cd ~/git_repos/ATL1415
repo=$PWD
rel_file=default_args/latest_release.txt
rel=$(grep '^--Release=' $rel_file | cut -d= -f2)     # 006
cyc=$(grep '^--cycles=' $rel_file | cut -d= -f2)      # 0332
ver=$(grep '^--version=' $rel_file | cut -d= -f2)     # 02
tspan=$(grep '^-t=' $rel_file | cut -d= -f2)          # 2018.75,2026.5
ATL14_root=/home/jovyan/ATL14_processing
s3_root=s3://maap-ops-workspace/ben_smith
ledgers=$ATL14_root/maap_ledgers
runs=$ATL14_root/runs
south=$ATL14_root/rel$rel/south
tile_list=$repo/ATL1415/resources/AA/40km_tile_list.txt     # 8944 centers
# half <60km|44km> -- every per-half path, from one place
half () {
    if [ $1 = 60km ]; then name=AA; else name=AA_44km; fi
    region_dir=$south/$name
    args=input_args_$name.txt
    s3_run=$s3_root/ATL1415/run_args/rel$rel/south/AA      # both args files live here
    s3_out=$s3_root/ATL14_processing/rel$rel/south/$name
    half_list=$ledgers/AA_$1_tile_list.txt
    tag=AA_$1_rel${rel}_${cyc}
    L=$ledgers/AA_$1_${cyc}
}


# ===========================================================================
# 0. [ADE] [OK on IS]  The DPS build is the commit you mean to run.  (arctic 0)
# ===========================================================================
/srv/conda/envs/notebook/bin/python register_algorithm.py --dry-run   # "push check: OK"
scripts/maap/check_build_id.py $s3_root/ATL1415/run_args/rel006/north/IS/input_args_IS.txt \
    maap-dps-worker-16gb                  # VERDICT: MATCH, maap_pgt=set


# ===========================================================================
# 1. [ADE] [OK]  Point the release symlinks at this release.
# ===========================================================================
ln -sf rel_006_0332.txt default_args/latest_release.txt
ln -sf AA_0331.txt      default_args/AA_latest.txt
# AA_0331.txt names only masks, the tide model and adjustment, and the z0
# scaling map -- nothing tied to the ATL11 generation -- so it serves cycles
# 03-32.  All of them are staged (checked 2026-09-19).


# ===========================================================================
# 2. [ADE] [OK, transect; UNTESTED at 0332]  Compose, derive and publish the args.
# ===========================================================================
setup_ATL1415_region.py default_args/MAAP_dps.txt $rel_file \
    default_args/AA_latest.txt default_args/quarterly.txt --Hemisphere=-1
scripts/maap/make_AA_44km_args.py $south/AA/input_args_AA.txt \
    $south/AA_44km/input_args_AA_44km.txt
half 60km; aws s3 cp $south/AA/input_args_AA.txt $s3_run/
aws s3 cp $south/AA_44km/input_args_AA_44km.txt $s3_run/
# CHECK: --tide_adjustment survived (the greedy-regex bug once dropped it),
# the previous product is a CMR search, and the two files differ in exactly
# -W and -b:
grep -E '^(--tide_adjustment|--tide_model|--mask_file|--previous_product)' $south/AA/input_args_AA.txt
diff $south/AA/input_args_AA.txt $south/AA_44km/input_args_AA_44km.txt


# ===========================================================================
# 3. [ADE] [UNTESTED]  Split the tile list into the two halves.
# ===========================================================================
# By the halves' own rule, into tile-list files under $ledgers (not the
# checkout).  Checked 2026-09-19 on the list: 8724 + 460, 240 in both.
python - <<EOF
import re
cs = [l.strip() for l in open('$tile_list') if l.strip()]
ext = lambda n: max(abs(int(v)) for v in re.match(r'E(-?\d+)_N(-?\d+)\.h5', n).groups()) * 1000
open('$ledgers/AA_60km_tile_list.txt', 'w').writelines(n + '\n' for n in cs if ext(n) >= 360000)
open('$ledgers/AA_44km_tile_list.txt', 'w').writelines(n + '\n' for n in cs if ext(n) <= 440000)
EOF
wc -l $ledgers/AA_60km_tile_list.txt $ledgers/AA_44km_tile_list.txt


# ===========================================================================
# 4. [DPS] [UNTESTED]  Smoke two known tiles on the current build.
# ===========================================================================
# The transect ran on ab84687; the solver has changed since (lineage in the
# tiles, the no-data exit).  Re-run two of its tiles and compare:
#   44km  220000 20000     E220_N20, the pole-hole edge: the worst memory
#                          (21.47 GiB) and the most crossovers (N_XO 171488)
#   60km  -580000 -980000  E-580_N-980, a grounding-line tide tile (0.25
#                          tide adjustment scale)
for h in 44km 60km; do half $h
    if [ $h = 44km ]; then xy='220000 20000'; else xy='-580000 -980000'; fi
    echo "$xy" > $ledgers/AA_${h}_smoke_xy.txt
    scripts/maap/submit_MAAP_jobs.py --xy_file $ledgers/AA_${h}_smoke_xy.txt \
        --step prelim --args_url $s3_run/$args --tile_prefix $s3_out \
        --queue maap-dps-worker-32gb --tag ${tag}_smoke --ledger ${L}_smoke_jobs.csv --dry-run
done
# (then without --dry-run; collect_jobs.py on each ledger)
# GATES: successful; N_XO the transect's order; N_AT HIGHER than the
# transect's (0332 has one more cycle than the 0331 it ran on); peak memory
# under 32 GiB.  /meta/lineage present in the fetched tile.


# ===========================================================================
# 5. [DPS] [UNTESTED]  Fan out, both halves.
# ===========================================================================
for h in 60km 44km; do half $h
    nohup scripts/maap/submit_MAAP_jobs.py --tile_list $half_list \
        --step prelim --args_url $s3_run/$args --tile_prefix $s3_out \
        --queue maap-dps-worker-32gb --tag ${tag}_prelim --ledger ${L}_prelim_jobs.csv \
        --max_in_flight 200 > ${L}_prelim_submit.log 2>&1 &
done
# --max_in_flight N=200 is a RECOMMENDATION, not a measured limit (GL step 4
# explains).  Never register while these run.


# ===========================================================================
# 6. [ADE] [OK on IS; UNTESTED at AA scale]  Watch, fetch, check -- both halves.
# ===========================================================================
for h in 60km 44km; do half $h
    scripts/maap/collect_jobs.py ${L}_prelim_jobs.csv > ${L}_prelim_collect.txt
    scripts/maap/fetch_tiles.py  ${L}_prelim_jobs.csv $region_dir --step prelim
    scripts/check_field_sizes.py $region_dir/prelim @$region_dir/$args
done
# check_field_sizes derives the shape from each half's own -W: [61, 61, 32]
# at 60 km, [45, 45, 32] at 44 km.  Each half's no-data centers go to its own
# $region_dir/prelim/no_data_tiles.txt (step 17).  The pole-hole tiles (e.g.
# E100_N20, 8 s in the transect) are no-data by construction.


# ===========================================================================
# 7. [DPS] [OK on IS; UNTESTED for AA]  Matched, both halves, each on its own prefix.
# ===========================================================================
for h in 60km 44km; do half $h
    nohup scripts/maap/submit_MAAP_jobs.py --tile_list $half_list \
        --step matched --args_url $s3_run/$args --tile_prefix $s3_out \
        --queue maap-dps-worker-32gb --tag ${tag}_matched --ledger ${L}_matched_jobs.csv \
        --max_in_flight 200 > ${L}_matched_submit.log 2>&1 &
done
# EACH HALF'S OWN --tile_prefix: see (a) at the top.  Centers without a prelim
# tile in THAT half are skipped by name.  discover packs 4 matched tiles per
# task (--lines_per_task 4); on DPS a job is one tile.


# ===========================================================================
# 8. [ADE] [OK on IS; UNTESTED for AA]  Watch, fetch, check matched.
# ===========================================================================
for h in 60km 44km; do half $h
    scripts/maap/collect_jobs.py ${L}_matched_jobs.csv > ${L}_matched_collect.txt
    scripts/maap/fetch_tiles.py  ${L}_matched_jobs.csv $region_dir --step matched
    scripts/check_field_sizes.py $region_dir/matched @$region_dir/$args
done


# ===========================================================================
# 9. [ADE] [UNTESTED]  200 km tiles, each half.
# ===========================================================================
# make_200km_tiles.py writes tile_run_<name>/ in the current directory, with
# queue/task_N and a runner named slurm_mos_run (not slurm_run.sh).  Plain
# bash, like the IS mosaic's, so it runs locally.
cd $runs
make_200km_tiles.py $south/AA AA -t $tspan
make_200km_tiles.py $south/AA_44km AA --name AA_south --W 44000 --spacing 40000 -t $tspan
for d in tile_run_AA tile_run_AA_south; do
    ( cd $d; seq 1 $(ls queue | wc -l) | xargs -P 8 -I{} env SLURM_ARRAY_TASK_ID={} bash slurm_mos_run )
done
cd $repo
# Each half's 200 km tile centers come from <region_dir>/200km_tile_list.txt
# if it exists, else are derived from that half's prelim tiles and written
# there.  QUESTION for Ben: ATL1415/resources/AA/200km_tile_list.txt (413
# centers) is in that format but matches NEITHER half: it contains every cell
# derivable from the 40 km list (411 for both halves together, 407 for the
# 60 km half) plus 2 more.  Which region directory is it for?  Until
# answered, each half derives its own.
# -P 8 is a guess: UNTIMED on AA.


# ===========================================================================
# 10. [ADE] [UNTESTED]  The four sectors.
# ===========================================================================
setup_AA_sectors.py $south
# Symlinks each half's 200 km tiles and per-tile outputs into A1-A4 and
# writes each sector's bounds.txt and input_args_A{n}.txt (from
# input_args_AA.txt).  Defaults: --north_name AA --south_name AA_44km.


# ===========================================================================
# 11. [ADE] [UNTESTED]  Mosaic, per sector.
# ===========================================================================
cd $runs
for s in A1 A2 A3 A4; do
    make_200km_to_mosaic_jobs.py -b $south/$s -rr $s -t $tspan    # NO --run: it would sbatch
    ( cd mosaic_run_$s; seq 1 $(ls queue | wc -l) | xargs -P 4 -I{} env SLURM_ARRAY_TASK_ID={} bash slurm_run.sh )
    check_mosaic_outputs.py $runs/mosaic_run_$s --values
done
cd $repo
# Its run directory is mosaic_run_<sector>; it has no --run_name and no -e
# (the template's environment defaults to ATL14, which is right here).
# THIS IS THE STEP THE ADE MAY NOT MANAGE: time one sector before the others.


# ===========================================================================
# 12. [ADE] [UNTESTED]  netCDF, per sector.
# ===========================================================================
# In place of scripts/run_antarctic_tonc.sh, which queues the same eight
# commands for sbatch (and with discover's env, IS2):
for s in A1 A2 A3 A4; do
    mkdir -p $runs/AA_${cyc}_nc_$s
    ( cd $runs/AA_${cyc}_nc_$s
      ATL14_write2nc.py @$south/$s/input_args_$s.txt > ATL14.log 2>&1
      ATL15_write2nc.py @$south/$s/input_args_$s.txt > ATL15.log 2>&1 )
done
# No INVALID line; XO rows NOT_SET in four attributes only (arctic 10).


# ===========================================================================
# 13. [ADE] [NO SOFTWARE]  Compare with the previous product.
# ===========================================================================
# Ben's bar: no >10 m errors, no major gaps.  Method: plan_cycles_03_32.sh T8
# "I9g6" (IS, against rel005).


# ===========================================================================
# 14. [ADE] [UNTESTED]  Publish the sector products.
# ===========================================================================
for s in A1 A2 A3 A4; do
    for f in $south/$s/ATL1[45]_${s}_${cyc}_*_${rel}_${ver}.nc; do
        aws s3 cp $f $s3_root/ATL14_processing/rel$rel/south/$s/
    done
done


# ===========================================================================
# 15. [ADE] [SUGGESTION, NO SOFTWARE]  Annotate the build history.  (ogc O12b)
# ===========================================================================


# ===========================================================================
# 16. [ADE+DPS] [NEEDS CODE: a multi-file --ATL14_reference_file]  Monthly.
# ===========================================================================
# As IS monthly (arctic 11-18) for both halves, with south_monthly in every
# path and setup_AA_sectors.py --near_pole_radius 0 -- EXCEPT the reference
# DEM.  AA's quarterly ATL14 is FOUR sector files, and the solver takes ONE:
# discover passes a glob ("rel005_0329/south/A*/ATL14_*_0329_100m_005_02.nc"),
# and a URI with a wildcard RAISES (by design: a glob over s3:// used to
# return nothing and silently edit away every point).  The error says "Name
# the granules explicitly, one --ATL14_reference_file each", but the argument
# is not action='append' (ATL11_to_ATL15.py), so a second one overwrites the
# first.  Either the message or the argument is wrong; deciding which is the
# first job of AA monthly.  Nothing else in this file blocks on it.


# ===========================================================================
# 17. [ADE] [OK on IS]  Take the no-data centers out of the list.
# ===========================================================================
# ONE list serves both halves, so a center in the overlap band comes out only
# if it had no data in EVERY half that solved it.
python - <<EOF
import os, re
root = '$south'
nd = {}
for h, name in (('60km', 'AA'), ('44km', 'AA_44km')):
    f = os.path.join(root, name, 'prelim', 'no_data_tiles.txt')
    nd[h] = {l.strip() for l in open(f)} if os.path.isfile(f) else set()
ext = lambda n: max(abs(int(v)) for v in re.match(r'E(-?\d+)_N(-?\d+)\.h5', n).groups()) * 1000
halves = lambda n: [h for h, ok in (('60km', ext(n) >= 360000), ('44km', ext(n) <= 440000)) if ok]
drop = {n for n in nd['60km'] | nd['44km'] if all(n in nd[h] for h in halves(n))}
lines = open('$tile_list').readlines()
open('$tile_list', 'w').writelines(l for l in lines if l.strip() not in drop)
print(f'removed {len(drop)} of {len(lines)}')
EOF
git commit -m "Drop no-data centers from the AA tile list" $tile_list && git push
