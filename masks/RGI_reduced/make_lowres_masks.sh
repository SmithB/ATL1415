#! /usr/bin/env bash


for j in `ls *.db`; do
    k=$(echo $j | sed s/.db/_10km.tif/)
    echo $j $k
    [ -f $k ] && rm $k
    gdal_rasterize -init 0 -burn 1 -tr 10000 10000 -at -tap $j $k
done
