gdal_translate -a_srs EPSG:3031 -co COMPRESS=LZW -co PREDICTOR=1 -co TILED=yes ../bedmap2_thickness_gt_50_plus_sio_shelves.tif bedmap2_thickness_gt_50_plus_sio_shelves_edited_v2.tif

for j in subtract/*.geojson;
do
    gdal_rasterize -burn 0 $j bedmap2_thickness_gt_50_plus_sio_shelves_edited_v2.tif
done

for j in add/*.geojson;
do
    gdal_rasterize -burn 1 $j bedmap2_thickness_gt_50_plus_sio_shelves_edited_v2.tif
done
