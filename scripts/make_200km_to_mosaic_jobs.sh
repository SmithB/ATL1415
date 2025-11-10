#! /usr/bin/env bash

source activate IS2

base=$1
region=`basename $base`

mosaic_run="mosaic_run_"$region

if [ -f $base/bounds.txt  ]; then
    crop="-c "$(head -1 $base/bounds.txt)
else
    crop=""
fi

echo "crop="$crop


[ -d $mosaic_run ] || mkdir $mosaic_run

for thedir in queue running done logs; do
    [ -d $mosaic_run/$thedir ] || mkdir $mosaic_run/$thedir
done

compute_sigma='True'

field=dz

glob_str="'matched/*.h5'"

task=0

task=$(($task+1))
this_replace='-R'
field_list="dz sigma_dz count misfit_rms misfit_scaled_rms mask cell_area"
echo "source activate IS2" > $mosaic_run/task_${task}
for field in $field_list; do
    echo $field
    echo "make_mosaic.py $crop $this_replace  -d ${base}/200km_tiles/dz  -g '*.h5' -O $base/dz.h5 --in_group dz/ -F  $field" >> $mosaic_run/task_${task} 
    this_replace=''
done

for lag in _lag1 _lag4 _lag8 _lag12 _lag16 _lag20 _lag24; do
    task=$(($task+1))
    TR='--t_range '$(echo $lag | sed s/_lag// | awk '{print 2019+0.25*$1/2, 2050}')
    echo "lag=$lag"
    field=dzdt$lag
    echo "source activate IS2" > $mosaic_run/task_${task}
    this_replace='-R'
    for ff in $field sigma_$field cell_area; do
        echo "make_mosaic.py $crop $TR $this_replace -d ${base}/200km_tiles/$field  -g '*.h5' -O $base/$field.h5 --in_group $field/ -F  $ff" >> $mosaic_run/task_${task} 
        this_replace=''
    done
done


for group in avg_dz_40000m avg_dz_20000m avg_dz_10000m; do
    field=$group
    echo "$group $field"    
    this_replace='-R'
    out=`echo $group | sed s/000m/km/ | sed s/avg_//`
    task=$(($task+1))
    echo "source activate IS2" > $mosaic_run/task_${task}
    this_replace='-R'
    for ff in $field sigma_$field cell_area; do
        echo "make_mosaic.py $crop $this_replace -d ${base}/200km_tiles/$group  -g '*.h5' -O $base/$out.h5 --in_group $group/ -F  $ff" >> $mosaic_run/task_${task} 
        this_replace=''
    done
     
    group=`echo $group | sed s/dz/dzdt/`
    out=`echo $out | sed s/dz/dzdt/`
    for lag in _lag1 _lag4 _lag8 _lag12 _lag16 _lag20 _lag24; do
        field=$group$lag
	TR='--t_range '$(echo $lag | sed s/_lag// | awk '{print 2019+0.25*$1/2, 2050}')
        task=$(($task+1))
        echo "source activate IS2" > $mosaic_run/task_${task}
        this_replace='-R'

        echo "lag=$lag, group=$group, task=$task"
        for ff in $field sigma_$field cell_area; do
            echo "make_mosaic.py $crop $TR $this_replace -d ${base}/200km_tiles/$field  -g '*.h5' -O $base/$out$lag.h5 --in_group $field/ -F  $ff" >> $mosaic_run/task_${task}
            this_replace=''
        done
    done
done


field_list="z0 misfit_rms misfit_scaled_rms mask cell_area count"


task=$(($task+1))

this_replace='-R'
if [ -d $base/200km_tiles/z0 ] ; then
    echo "200km" 
    field_list=$field_list" sigma_z0"
    for field in $field_list; do
        echo $field
        echo "make_mosaic.py $crop $this_replace  -d ${base}/200km_tiles/z0  -g '*.h5' -O $base/z0.h5 --in_group z0/ -F  $field" >> $mosaic_run/task_${task} 
        this_replace=''
    done
fi
mv $mosaic_run/task* $mosaic_run/queue

cat slurm_scripts/slurm_mos_run | sed s/LAST_TASK/$task/ > $mosaic_run/slurm_mos_run
pushd $mosaic_run
sbatch slurm_mos_run
