# howto_MAAP_AA.sh -- Antarctica on MAAP.
# The discover/SLURM variant is docs/howto_AA.sh.  Reasoning, the transect
# record and open design work: ~/ATL14_processing/dev/ (howto_long/howto_MAAP_AA.sh,
# plans/plan_pack_tiles.sh, plan_rerun_timing.sh).
# STATUS: rel006 0332 quarterly ran through step 9 on MAAP (2026-10-15..18;
# dev/plans/plan_full_run_iceCaps.sh A0-A6, plan_AA_mosaic_dps.sh).  Steps 1-6
# ran with dev drivers (drive_packed.py), not the commands below, which keep
# their [UNTESTED] tags.  Step 10 (monthly) has not run.
# Before a full run: plan_pack_tiles.sh K4-K8 (several tiles per job,
# start-up cap) -- AA's job count makes start-up failures the main risk.

conda activate ATL14
cd ~/git_repos/ATL1415
repo=$PWD
rel_file=default_args/latest_release.txt
rel=$(grep '^--Release=' $rel_file | cut -d= -f2)     # 006
cyc=$(grep '^--cycles=' $rel_file | cut -d= -f2)      # 0332
ATL14_root=/home/jovyan/ATL14_processing
s3_root=s3://maap-ops-workspace/ben_smith
mnt=~/my-private-bucket
ledgers=$ATL14_root/maap_ledgers
south=$ATL14_root/rel$rel/south
tile_list=$repo/ATL1415/resources/AA/40km_tile_list.txt
half () {     # half 60km|44km
    if [ $1 = 60km ]; then name=AA; else name=AA_44km; fi
    region_dir=$south/$name
    args=input_args_$name.txt
    s3_run=$s3_root/ATL1415/run_args/rel$rel/south/AA      # both args files live here
    s3_out=$s3_root/ATL14_processing/rel$rel/south/$name
    half_list=$ledgers/AA_$1_tile_list.txt
    tag=AA_$1_rel${rel}_${cyc}
    L=$ledgers/AA_$1_${cyc}
    sub="scripts/maap/submit_MAAP_jobs.py --args_url $s3_run/$args --tile_prefix $s3_out"
}

# 1. Build check: howto_MAAP_ogc.sh step 3 (VERDICT: MATCH).

# 2. [UNTESTED at 0332] Args for both halves.  AA_latest.txt still points at
#    AA_0331.txt, and 0331 is gone from CMR: update it first.
setup_ATL1415_region.py default_args/MAAP_dps.txt $rel_file \
    default_args/AA_latest.txt default_args/quarterly.txt --Hemisphere=-1
scripts/maap/make_AA_44km_args.py $south/AA/input_args_AA.txt $south/AA_44km/input_args_AA_44km.txt
half 60km; aws s3 cp $south/AA/input_args_AA.txt $s3_run/
aws s3 cp $south/AA_44km/input_args_AA_44km.txt $s3_run/

# 3. [UNTESTED] Split the tile list into the halves (60 km: |x| or |y| >= 360 km;
#    44 km: <= 440 km; the band in between is in both).
python - <<EOF
import re
cs = [l.strip() for l in open('$tile_list') if l.strip()]
ext = lambda n: max(abs(int(v)) for v in re.match(r'E(-?\d+)_N(-?\d+)\.h5', n).groups()) * 1000
open('$ledgers/AA_60km_tile_list.txt', 'w').writelines(n + '\n' for n in cs if ext(n) >= 360000)
open('$ledgers/AA_44km_tile_list.txt', 'w').writelines(n + '\n' for n in cs if ext(n) <= 440000)
EOF

# 4. [UNTESTED] Smoke two transect tiles, one per half.
for h in 44km 60km; do half $h
    if [ $h = 44km ]; then xy='220000 20000'; else xy='-580000 -980000'; fi
    echo "$xy" > ${L}_smoke_xy.txt
    $sub --xy_file ${L}_smoke_xy.txt --step prelim --queue maap-dps-worker-32gb \
        --tag ${tag}_smoke --ledger ${L}_smoke_jobs.csv
done

# 5. [UNTESTED at scale] Prelim, both halves; resubmit failures; field sizes.
for h in 60km 44km; do half $h
    nohup $sub --tile_list $half_list --step prelim --queue maap-dps-worker-32gb \
        --tag ${tag}_prelim --ledger ${L}_prelim_jobs.csv \
        --max_in_flight 100 > ${L}_prelim_submit.log 2>&1 &
done
for h in 60km 44km; do half $h
    scripts/check_field_sizes.py ${s3_out/$s3_root/$mnt}/prelim @$region_dir/$args
done

# 6. [UNTESTED at scale] Matched: step 5 with --step matched.

# 7. [OK] The four sectors (A1-A4), once both halves' prelim tiles are complete:
#    bounds, args, tile manifest (prelim_tiles.txt: a tile belongs to the 200 km
#    cell holding its center) and each half's 200 km cells, then publish them.
s3_sec=$s3_root/ATL14_processing/rel$rel/south                  # sector files
s3_prod=$s3_root/ATL14_processing/rel${rel}_${cyc}_testing/south  # 200 km tiles, mosaics, nc
scripts/maap/setup_AA_sectors_s3.py $south --tiles_60km $s3_sec/AA --tiles_44km $s3_sec/AA_44km
for n in 1 2 3 4; do
    aws s3 cp $south/A$n/bounds.txt $s3_sec/A$n/
    aws s3 cp $south/A$n/prelim_tiles.txt $s3_sec/A$n/
    aws s3 cp $south/A$n/input_args_A$n.txt $s3_root/ATL1415/run_args/rel$rel/south/A$n/
done

# 8. [OK] 200 km tiles, each half's cells written into their sector; then the
#    mosaics per sector (A4 z0 fits 16gb).  422 + 164 jobs.
for n in 1 2 3 4; do for h in 60km 44km; do half $h
    $sub --out_prefix $s3_prod/A$n --step mosaic200 --queue maap-dps-worker-16gb \
        --tag ${tag}_mosaic200_A$n --ledger ${L}_mosaic200_A${n}_jobs.csv --max_in_flight 60 \
        $(sed 's/^/--task=/' $south/A$n/mosaic200_tasks_$h.txt)
done; done
sector () {   # sector <n>: submit for sector n
    scripts/maap/submit_MAAP_jobs.py --args_url $s3_root/ATL1415/run_args/rel$rel/south/A$1/input_args_A$1.txt \
        --tile_prefix $s3_sec/A$1 --out_prefix $s3_prod/A$1 "${@:2}"
}
for n in 1 2 3 4; do
    sector $n --step mosaic --queue maap-dps-worker-16gb --tag AA_rel${rel}_${cyc}_mosaic_A$n \
        --ledger $ledgers/AA_${cyc}_mosaic_A${n}_jobs.csv
done

# 9. [OK] netCDF per sector (ATL14 peaks 17 GiB): needs a build with the sector
#    manifest (a0e157b or later).  Checks as howto_MAAP_GL.sh step 9.
for n in 1 2 3 4; do
    sector $n --step nc --queue maap-dps-worker-32gb --tag AA_rel${rel}_${cyc}_nc_A$n \
        --ledger $ledgers/AA_${cyc}_nc_A${n}_jobs.csv
done

# 10. [UNTESTED] Monthly: one partition (all centers, 60 km), reference = the
#     quarterly sectors' ATL14s, read with from_list (an s3 glob is allowed):
#     --ATL14_reference_file=$s3_prod/A?/ATL14_A?_${cyc}_100m_${rel}_02.nc
