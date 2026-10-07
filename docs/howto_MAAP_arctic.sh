# howto_MAAP_arctic.sh -- the smaller Arctic regions (IS RA CN CS SV) on MAAP, all on DPS.
# The discover/SLURM variant is docs/howto_arctic.sh.  Greenland: howto_MAAP_GL.sh.
# Reasoning, records and results: ~/ATL14_processing/dev/ (plans/plan_IS_run.sh,
# plan_dps_mosaic.sh D5-D6m; howto_long/howto_MAAP_arctic.sh).
# Tags: [OK on IS] run end to end on Iceland; [UNTESTED] for RA CN CS SV.
# These regions mosaic directly (no 200 km step).  Every DPS step: watch with
# collect_jobs.py, resubmit what did not succeed under a NEW ledger
# (`failed`, as in howto_MAAP_GL.sh) until nothing is left.

conda activate ATL14
cd ~/git_repos/ATL1415
repo=$PWD
regions="RA CN CS SV"        # IS is done at 0332
rel_file=default_args/latest_release.txt
rel=$(grep '^--Release=' $rel_file | cut -d= -f2)     # 006
cyc=$(grep '^--cycles=' $rel_file | cut -d= -f2)      # 0332
ver=$(grep '^--version=' $rel_file | cut -d= -f2)     # 02
ATL14_root=/home/jovyan/ATL14_processing
s3_root=s3://maap-ops-workspace/ben_smith
mnt=~/my-private-bucket                               # = $s3_root, mounted
ledgers=$ATL14_root/maap_ledgers
paths () {    # paths <region> [_monthly]
    region_dir=$ATL14_root/rel$rel/north$2/$1
    s3_run=$s3_root/ATL1415/run_args/rel$rel/north$2/$1
    s3_out=$s3_root/ATL14_processing/rel$rel/north$2/$1
    s3_prod=$s3_root/ATL14_processing/rel${rel}_${cyc}_testing/north$2/$1
    tile_list=$repo/ATL1415/resources/$1/40km_tile_list.txt
    tag=$1_rel${rel}_${cyc}$2
    L=$ledgers/$1_${cyc}$2
    sub="scripts/maap/submit_MAAP_jobs.py --args_url $s3_run/input_args_$1.txt --tile_prefix $s3_out"
}

# 1. [OK on IS] Build check: howto_MAAP_ogc.sh step 3 (VERDICT: MATCH).

# 2. [OK on IS] Compose and publish the quarterly args.
ln -sf rel_006_0332.txt default_args/latest_release.txt
for reg in $regions; do paths $reg
    setup_ATL1415_region.py default_args/MAAP_dps.txt $rel_file \
        default_args/$reg.txt default_args/quarterly.txt --Hemisphere=1
    aws s3 cp $region_dir/input_args_$reg.txt $s3_run/
done

# 3. [OK on IS] Smoke: one tile per region; check it before fanning out.
for reg in $regions; do paths $reg
    $sub --tile_list $tile_list --limit 1 --step prelim --queue maap-dps-worker-16gb \
        --tag ${tag}_smoke --ledger ${L}_smoke_jobs.csv
done

# 4. [OK on IS] Prelim, then field sizes through the mount.
for reg in $regions; do paths $reg
    $sub --tile_list $tile_list --step prelim --queue maap-dps-worker-16gb \
        --tag ${tag}_prelim --ledger ${L}_prelim_jobs.csv --max_in_flight 100
done
for reg in $regions; do paths $reg
    scripts/check_field_sizes.py ${s3_out/$s3_root/$mnt}/prelim @$region_dir/input_args_$reg.txt
done

# 5. [OK on IS] Matched (after ALL prelim tiles exist), then field sizes.
for reg in $regions; do paths $reg
    $sub --tile_list $tile_list --step matched --queue maap-dps-worker-16gb \
        --tag ${tag}_matched --ledger ${L}_matched_jobs.csv --max_in_flight 100
done
for reg in $regions; do paths $reg
    scripts/check_field_sizes.py ${s3_out/$s3_root/$mnt}/matched @$region_dir/input_args_$reg.txt --step matched
done

# 6. [OK on IS] Mosaics (one job per field group), then netCDF (ATL14 + ATL15).
for reg in $regions; do paths $reg
    $sub --out_prefix $s3_prod --step mosaic --queue maap-dps-worker-16gb \
        --tag ${tag}_mosaic --ledger ${L}_mosaic_jobs.csv
done
for reg in $regions; do paths $reg
    $sub --out_prefix $s3_prod --step nc --queue maap-dps-worker-16gb \
        --tag ${tag}_nc --ledger ${L}_nc_jobs.csv
done

# 7. [OK on IS] Checks: as howto_MAAP_GL.sh step 9.

# 8. [OK on IS] Monthly: args with the quarterly ATL14 as reference, then
#    steps 3-6 after `paths $reg _monthly`; the nc step with --task ATL15 only.
for reg in $regions; do paths $reg
    ref=$s3_prod/ATL14_${reg}_${cyc}_100m_${rel}_${ver}.nc
    setup_ATL1415_region.py default_args/MAAP_dps.txt $rel_file \
        default_args/$reg.txt default_args/monthly.txt --Hemisphere=1 \
        --ATL14_reference_file=$ref
    paths $reg _monthly
    aws s3 cp $region_dir/input_args_$reg.txt $s3_run/
done

# 9. [OK on IS] Take no-data centers out of each tile list (from both periods'
#    prelim/no_data_tiles.txt); commit AND push.
