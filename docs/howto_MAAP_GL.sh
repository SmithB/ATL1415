# howto_MAAP_GL.sh -- Greenland on MAAP: tiles, mosaics and netCDFs all solved as DPS jobs.
# The discover/SLURM variant is docs/howto_GL.sh.  Reasoning, records and
# results: ~/ATL14_processing/dev/ (plans/plan_GL_north.sh, plan_GL_maskv5.sh,
# plan_dps_mosaic.sh; howto_long/howto_MAAP_GL.sh).
# Tags: [OK] run end to end on GL (2026-09-29 .. 2026-10-07); [UNTESTED].
# Every DPS step: watch with collect_jobs.py, then resubmit what did not
# succeed under a NEW ledger (step 0's `failed`) until nothing is left.
# ~1-2% of jobs fail at start-up on MAAP API calls (MAAP Community #1334);
# resubmitting fixes them.  Never register while jobs are in flight.

conda activate ATL14
cd ~/git_repos/ATL1415
repo=$PWD
rel_file=default_args/latest_release.txt
rel=$(grep '^--Release=' $rel_file | cut -d= -f2)     # 006
cyc=$(grep '^--cycles=' $rel_file | cut -d= -f2)      # 0332
ver=$(grep '^--version=' $rel_file | cut -d= -f2)     # 02
ATL14_root=/home/jovyan/ATL14_processing
s3_root=s3://maap-ops-workspace/ben_smith
mnt=~/my-private-bucket                               # = $s3_root, mounted
ledgers=$ATL14_root/maap_ledgers
tile_list=$repo/ATL1415/resources/GL/40km_tile_list.txt     # 1483 centers
paths () {    # paths [_monthly]
    region_dir=$ATL14_root/rel$rel/north$1/GL                       # args only
    s3_run=$s3_root/ATL1415/run_args/rel$rel/north$1/GL             # args the jobs read
    s3_out=$s3_root/ATL14_processing/rel$rel/north$1/GL             # tiles
    s3_prod=$s3_root/ATL14_processing/rel${rel}_${cyc}_testing/north$1/GL   # mosaics + nc
    tag=GL_rel${rel}_${cyc}$1
    L=$ledgers/GL_${cyc}$1
    sub="scripts/maap/submit_MAAP_jobs.py --args_url $s3_run/input_args_GL.txt --tile_prefix $s3_out"
}
failed () {   # failed <ledger>: tiles (or mosaic tasks) without a successful job
python - "$1" <<'EOF'
import csv, sys
from maap.maap import MAAP
m = MAAP()
for r in csv.DictReader(open(sys.argv[1])):
    try: ok = m.get_job_status(r['job_id']).json().get('status') == 'successful'
    except Exception: ok = False
    if not ok:
        print(r['task'] if r['task'] != '-' else
              'E%d_N%d.h5' % (int(int(r['x0']) / 1000), int(int(r['y0']) / 1000)))
EOF
}

# 1. [OK] Build check: howto_MAAP_ogc.sh step 3 (VERDICT: MATCH).

# 2. [OK] Compose and publish the quarterly args.
paths
ln -sf rel_006_0332.txt default_args/latest_release.txt
ln -sf GL_0332.txt      default_args/GL_latest.txt
setup_ATL1415_region.py default_args/MAAP_dps.txt $rel_file \
    default_args/GL_latest.txt default_args/quarterly.txt --Hemisphere=1
aws s3 cp $region_dir/input_args_GL.txt $s3_run/

# 3. [OK] Smoke: two tiles, then check them before fanning out.
printf '80000 -920000\n480000 -1040000\n' > ${L}_smoke_xy.txt
$sub --xy_file ${L}_smoke_xy.txt --step prelim --queue maap-dps-worker-32gb \
    --tag ${tag}_smoke --ledger ${L}_smoke_jobs.csv
scripts/maap/collect_jobs.py ${L}_smoke_jobs.csv

# 4. [OK] Prelim, every tile, 100 in flight (runs for hours: nohup).
nohup $sub --tile_list $tile_list --step prelim --queue maap-dps-worker-32gb \
    --tag ${tag}_prelim --ledger ${L}_prelim_jobs.csv \
    --max_in_flight 100 > ${L}_prelim_submit.log 2>&1 &
scripts/tile_dash.py --ledger "${L}_prelim*_jobs.csv" --tile_list $tile_list   # live view (Ctrl+C)
scripts/maap/collect_jobs.py ${L}_prelim_jobs.csv > ${L}_prelim_collect.txt
failed ${L}_prelim_jobs.csv > ${L}_prelim_retry1_tile_list.txt
$sub --tile_list ${L}_prelim_retry1_tile_list.txt --step prelim --queue maap-dps-worker-32gb \
    --tag ${tag}_prelim_retry1 --ledger ${L}_prelim_retry1_jobs.csv --max_in_flight 100
#    ... retry2, retry3 the same way until `failed` prints nothing.

# 5. [OK] Field sizes, read through the mount.
scripts/check_field_sizes.py ${s3_out/$s3_root/$mnt}/prelim @$region_dir/input_args_GL.txt

# 6. [OK] Matched: steps 4-5 with --step matched (after ALL prelim tiles exist).
nohup $sub --tile_list $tile_list --step matched --queue maap-dps-worker-32gb \
    --tag ${tag}_matched --ledger ${L}_matched_jobs.csv \
    --max_in_flight 100 > ${L}_matched_submit.log 2>&1 &
scripts/check_field_sizes.py ${s3_out/$s3_root/$mnt}/matched @$region_dir/input_args_GL.txt --step matched

# 7. [OK] 200 km tiles (one job per 200 km center), then the mosaics (one job per field group).
$sub --out_prefix $s3_prod --step mosaic200 --queue maap-dps-worker-16gb \
    --tag ${tag}_mosaic200 --ledger ${L}_mosaic200_jobs.csv --max_in_flight 100
#    resubmit:  $sub --out_prefix $s3_prod --step mosaic200 --task=<x_y> [--task=...] ... (new ledger;
#               "--task=", as a negative x reads as an option otherwise)
$sub --out_prefix $s3_prod --step mosaic --queue maap-dps-worker-16gb \
    --tag ${tag}_mosaic --ledger ${L}_mosaic_jobs.csv --max_in_flight 100

# 8. [OK] netCDF: ATL14 and ATL15.
$sub --out_prefix $s3_prod --step nc --queue maap-dps-worker-32gb \
    --tag ${tag}_nc --ledger ${L}_nc_jobs.csv

# 9. [OK] Checks.  Nothing under $s3_prod may predate step 7 (stale products):
aws s3 ls --recursive $s3_prod/ | sort | head
#    check_mosaic_outputs.py --values on a run dir built with -b ${s3_prod/$s3_root/$mnt}
#    (method: dev/plans/plan_dps_mosaic.sh D7); compare with the previous
#    release (method: dev/plans/plan_cycles_03_32.sh T8).

# 10. [OK on GL north] Monthly: args with the quarterly ATL14 as reference,
#     then steps 3-8 after `paths _monthly`; step 8 with --task ATL15 only.
ref=$s3_prod/ATL14_GL_${cyc}_100m_${rel}_${ver}.nc
setup_ATL1415_region.py default_args/MAAP_dps.txt $rel_file \
    default_args/GL_latest.txt default_args/monthly.txt --Hemisphere=1 \
    --ATL14_reference_file=$ref
paths _monthly
aws s3 cp $region_dir/input_args_GL.txt $s3_run/

# 11. [UNTESTED on MAAP] Take no-data centers out of the tile list; commit AND push.
