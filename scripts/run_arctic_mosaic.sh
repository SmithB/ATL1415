#! /usr/bin/env bash

# run the arctic mosaic.  The first argument should be the release file, and should specify the release number
# as Release=[NNN]
# the second argument should be the location file and should specify
# ATL14_root=[path]
# If no other arguments are specified, will run RA IS CN CS and SV, otherwise will run the specified region

release_file=$1
shift

loc_file=$1
shift

period_file=$1
shift

echo "release_file=$release_file, loc_file=$loc_file"
if [ ! -f $release_file ]; then
    echo "release file not found"; exit
fi

if $(grep -q monthly $period_file) ; then
    time_resolution="monthly"
    hemi_suffix="_monthly"
else
    time_resolution="quarterly"
    hemi_suffix=""
fi
echo $hemi_suffix

release=`grep Release $release_file | sed s/\=/\ / | awk '{print $NF}'`
root=`grep ATL14_root $loc_file | sed s/\=/\ / | awk '{print $NF}'`

if [ $# -eq 0 ]; then
    regions="RA IS CN CS SV"
else
    regions=$1
fi

echo $release
echo $root
for reg in $regions; do

    run_dir=${reg}${hemi_suffix}_mosaic
    [ -d $run_dir ] && rm -r $run_dir

    base=${root}/rel${release}/north${hemi_suffix}/${reg}/
    echo $base
    if [ $time_resolution == "monthly" ] ; then
        make_mosaic_jobs_monthly $base
    else
        make_mosaic_jobs $base
    fi
    pushd $run_dir
    sbatch slurm_mos_run
    popd
done
