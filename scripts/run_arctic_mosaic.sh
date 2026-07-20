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
    hemi_suffix="_monthly"
else
    hemi_suffix=""
fi
echo $hemi_suffix

release=`grep '^--Release=' $release_file | sed s/\=/\ / | awk '{print $NF}'`
root=`grep '^--ATL14_root=' $loc_file | sed s/\=/\ / | awk '{print $NF}'`
time_span=`grep '^-t=' $release_file | sed s/\=/\ / | awk '{print $NF}'`

if [ -z "$release" ]; then
    echo "could not find --Release= in $release_file"; exit
fi
if [ -z "$root" ]; then
    echo "could not find --ATL14_root= in $loc_file"; exit
fi
if [ -z "$time_span" ]; then
    echo "could not find -t= in $release_file"; exit
fi
version=`grep '^--version=' $release_file | sed s/\=/\ / | awk '{print $NF}'`

echo ""
echo "======================================================"
echo "  RELEASE: $release   VERSION: $version"
echo "  Release file: $(readlink -f $release_file)"
echo "======================================================"
echo ""

if [ $# -eq 0 ]; then
    regions="RA IS CN CS SV"
else
    regions=$1
fi

echo $release
echo $root
for reg in $regions; do

    base=${root}/rel${release}/north${hemi_suffix}/${reg}/
    echo $base
    run_dir=${reg}${hemi_suffix}_mosaic
    [ -d $run_dir ] && rm -r $run_dir
    make_mosaic_jobs.py -b $base -rr $reg -t $time_span --run_name $run_dir @$period_file
    slurm_script=slurm_run.sh
    pushd $run_dir
    sbatch $slurm_script
    popd
done
