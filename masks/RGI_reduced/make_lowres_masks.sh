#! /usr/bin/env bash


for j in `ls *.db`; do
    k=$(echo $j | sed s/.db/_40km.tif/)
    echo $j $k
    [ -f $k ] && rm $k
    gdal_rasterize -init 0 -burn 1 -tr 40000 40000 -at -tap $j $k
done

for j in `ls *.db`; do
    k=$(echo $j | sed s/.db/_80km.tif/)
    echo $j $k
    [ -f $k ] && rm $k
    gdal_rasterize -init 0 -burn 1 -tr 80000 80000 -at -tap $j $k
done



