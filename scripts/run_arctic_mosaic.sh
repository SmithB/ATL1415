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
echo "release_file=$release_file, loc_file=$loc_file"
if [ ! -f $release_file ]; then
    echo "release file not found"; exit
fi

release=`grep Release $release_file | sed s/\=/\ / | awk '{print $NF}'`
root=`grep ATL14_root $loc_file | sed s/\=/\ / | awk '{print $NF}'`

if [ $# -eq 0 ]; then
    regions="RA IS CN CS SV"
fi
echo $release
echo $root

for reg in $regions; do 
    [ -d mosaic_run_${reg} ] && rm -r mosaic_run_${reg} 
    [ -d ${reg}_mosaic ] && rm -r ${reg}_mosaic

    scripts/make_mosaic_jobs ${root}/rel${release}/north/${reg}
    pushd $reg"_mosaic"
    sbatch slurm_mos_run
    popd
done
