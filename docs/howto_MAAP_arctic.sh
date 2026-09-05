# howto_MAAP_arctic.sh -- the arctic regions (RA IS CN CS SV) on MAAP
#
# ############################################################################
# ##  TENTATIVE.  Written 2026-09-05 BEFORE any of it has been run end to   ##
# ##  end -- no ATL1415 tile has been solved on DPS yet.  This is the plan,  ##
# ##  not a record of a successful run.  Expect steps to move, split and    ##
# ##  change as testing advances; revise this file as that happens.         ##
# ############################################################################
#
# The discover/SLURM variant is docs/howto_arctic.sh, which is still the
# production path and is NOT replaced by this file.
#
# READ docs/howto_MAAP_GL.sh FIRST.  GL is the reference workflow; this file
# marks only what the arctic regions do differently.  Same tags:
#   [ADE] / [DPS]   [OK] [UNTESTED] [NEEDS CODE: x]
# Steps are numbered 1..11.
#
# Prerequisite: docs/howto_MAAP_staging.sh S1-S6; step 5 is gated on S7.
#
# WHY THIS IS THE ONE TO TEST FIRST (Q21):
#   ICELAND (IS) IS THE TEST REGION.  It is small enough to fan out under the
#   ~10 jobs/hr public-queue throttle, it is a real region rather than a toy,
#   and the crossover work was verified against an Iceland box (13 tiles).
#   Do the whole of this file for IS alone before running the five-region loop.
#
# WHAT MAKES THE ARCTIC DIFFERENT:
#   a. FIVE REGIONS, not one.  The discover workflow drives them through
#      scripts/run_arctic_{prelim,matched,mosaic,to_nc}.sh, which loop over
#      "RA IS CN CS SV".  Those greps already read --Release/--ATL14_root/
#      --cycles/--version out of MAAP_dps.txt unchanged; ONLY the
#      `setup_slurm_run.py ...; sbatch` tail has to change.
#      PUT THE MAAP VARIANTS IN scripts/maap/ so the SLURM originals stay put.
#   b. THE MASKS ARE VECTOR .db FILES, not geotiffs (Q12).  masks/RGI_reduced/
#      is on the bucket as real files, and GDAL's SQLite driver reads a .db
#      over /vsis3.  ATL1415/make_mask_from_vector.py called ogr.Open() on the
#      raw URI, which failed; as of 2026-09-05 it routes through
#      pc.io_utils.as_gdal_path().  THAT FIX HAS NOT BEEN EXERCISED AGAINST
#      THE BUCKET -- step 2 is where it gets checked, and it is the single
#      thing most likely to stop this file at the first hurdle.
#   c. the region args files have no cycles suffix: default_args/{RA,IS,CN,CS,SV}.txt.
#   d. all five are northern hemisphere, so rel006/north/<REG>/.

conda activate ATL14
cd ~/git_repos/ATL1415
reg=IS                       # <-- do IS alone first (Q21)
region_dir=/home/jovyan/ATL14_processing/rel006/north/$reg
s3_run=s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/north/$reg
s3_out=s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north/$reg


# ===========================================================================
# 1. [ADE] [OK]  Point the release symlink at this release.
# ===========================================================================
ln -sf rel_006_0331.txt default_args/latest_release.txt


# ===========================================================================
# 2. [ADE] [UNTESTED]  CHECK THE .db MASK READS FROM THE BUCKET.   <-- do this first
# ===========================================================================
# Q12's one-line fix, never yet exercised.  If this raises "No such file or
# directory", nothing else in this file will work, and the problem is in
# ATL1415/make_mask_from_vector.py, not in the staging.
#
# THE TEST TILE, given by Ben 2026-09-05: x0,y0 = 1260, -2620 km.  43K points --
# a medium-sized dataset, likely on the edge of the ice sheet, so it exercises
# the mask edge rather than a saturated interior tile.  On the grid: centers sit
# at odd multiples of tile_spacing/2, and 1260 = 63 x 20 km, -2620 = -131 x 20 km.
# The same tile is the DPS smoke test -- howto_MAAP_staging.sh S7.
python - <<'EOF'
import pointCollection as pc
from ATL1415.make_mask_from_vector import make_mask_from_vector
db = 's3://maap-ops-workspace/ben_smith/ATL1415/masks/RGI_reduced/06_rgi60_Iceland_reduced.db'
print(pc.io_utils.as_gdal_path(db))          # expect /vsis3/maap-ops-workspace/...
m = make_mask_from_vector(db, W={'x':6e4, 'y':6e4},
                          ctr={'x':1260000., 'y':-2620000.}, spacing=100,
                          srs_proj4=('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 '
                                     '+k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))
print(m.z.shape, m.z.sum())    # expect a nonzero count, and NOT the full 601x601 --
                               # an edge tile should be part ice, part not
EOF


# ===========================================================================
# 3. [ADE] [UNTESTED]  Compose the args file.        (as GL step 2)
# ===========================================================================
setup_ATL1415_region.py default_args/MAAP_dps.txt default_args/latest_release.txt \
    default_args/$reg.txt default_args/quarterly.txt --Hemisphere=1


# ===========================================================================
# 4. [ADE] [UNTESTED]  Publish the args file.        (as GL step 3)
# ===========================================================================
aws s3 cp $region_dir/input_args_$reg.txt $s3_run/


# ===========================================================================
# 5. [ADE] [NEEDS CODE: make_ATL1415_queue.py --xy_out]  Tile centers.
# ===========================================================================
# Same four blockers as GL step 4.  The 1 km mask question (Q6/Q16) takes a
# different shape here: the arctic regions carry a `_40km.tif` sibling
# alongside the .db, so the grid may come from that rather than from a 1 km
# decimation.  CONFIRM WHICH before assuming the GL/AA answer applies.
make_ATL1415_queue.py prelim $region_dir/input_args_$reg.txt --xy_out ${reg}_prelim_xy.txt


# ===========================================================================
# 6. [DPS] [NEEDS CODE: scripts/submit_MAAP_jobs.py]  Fan out.  (as GL step 5)
# ===========================================================================
# GATED ON THE SMOKE TEST (staging S7).  IS is small enough that even the
# ~10 jobs/hr public throttle is survivable -- which is the point of testing here.
submit_MAAP_jobs.py --xy_file ${reg}_prelim_xy.txt --step prelim \
    --args_url $s3_run/input_args_$reg.txt --out_prefix $s3_out/prelim \
    --queue maap-dps-worker-32gb \
    --tag ${reg}_rel006_prelim --ledger ${reg}_prelim_jobs.csv


# ===========================================================================
# 7. [DPS] [NEEDS CODE: scripts/check_MAAP_jobs.py]  Watch.
# ===========================================================================
# The MAAP analogue of the discover idiom
#   for j in RA IS CN CS SV; do echo $j; slurm_run_status.py $j"_prelim"; done
for j in RA IS CN CS SV; do echo $j; check_MAAP_jobs.py ${j}_prelim_jobs.csv; done


# ===========================================================================
# 8. [ADE] [NEEDS CODE: deterministic output prefix]  Collect.  (as GL step 7)
# ===========================================================================
aws s3 sync $s3_out/prelim/ $region_dir/prelim/


# ===========================================================================
# 9. [DPS] [NEEDS CODE: run.sh prelim_prefix input]  Matched.  (as GL step 9)
# ===========================================================================
make_ATL1415_queue.py matched $region_dir/input_args_$reg.txt --xy_out ${reg}_matched_xy.txt
# ... submit_MAAP_jobs.py --step matched --prelim_prefix $s3_out/prelim ...
aws s3 sync $s3_out/matched/ $region_dir/matched/


# ===========================================================================
# 10. [ADE] [NEEDS CODE: run_queue_local.sh]  Mosaic and netCDF.
# ===========================================================================
make_mosaic_jobs.py -b $region_dir -rr $reg -t 2018.75,2026.5 \
    --run_name ${reg}_mosaic @default_args/quarterly.txt
run_queue_local.sh ${reg}_mosaic -P 8
ATL14_write2nc.py @$region_dir/input_args_$reg.txt
ATL15_write2nc.py @$region_dir/input_args_$reg.txt


# ===========================================================================
# 11. [ADE+DPS] [NEEDS CODE: scripts/maap/run_arctic_*.sh]  All five regions.
# ===========================================================================
# ONLY AFTER IS HAS GONE THROUGH END TO END.  These are the MAAP counterparts
# of scripts/run_arctic_{prelim,matched,mosaic,to_nc}.sh -- same argument
# order (release file, location file, period file, optional region), same
# per-region loop, same directory bookkeeping; the tail submits to DPS or runs
# the queue locally instead of calling sbatch.
bash scripts/maap/run_arctic_prelim.sh  default_args/latest_release.txt \
     default_args/MAAP_dps.txt default_args/quarterly.txt
for j in RA IS CN CS SV; do echo $j; check_MAAP_jobs.py ${j}_prelim_jobs.csv; done

bash scripts/maap/run_arctic_matched.sh default_args/latest_release.txt \
     default_args/MAAP_dps.txt default_args/quarterly.txt
for j in RA IS CN CS SV; do echo $j; check_MAAP_jobs.py ${j}_matched_jobs.csv; done

bash scripts/maap/run_arctic_mosaic.sh  default_args/latest_release.txt \
     default_args/MAAP_dps.txt default_args/quarterly.txt
bash scripts/maap/run_arctic_to_nc.sh   default_args/latest_release.txt \
     default_args/MAAP_dps.txt default_args/quarterly.txt

# monthly: the same four with default_args/monthly.txt.  The run_arctic_*.sh
# scripts already derive hemi_suffix="_monthly" by grepping the period file,
# so that part needs no change (Q17).
