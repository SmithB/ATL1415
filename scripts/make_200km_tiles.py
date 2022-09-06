import glob
import numpy as np
import re
import sys
import pointCollection as pc
import os
import stat


def make_fields():
    fields={}
    fields['z0']="z0 sigma_z0 misfit_rms misfit_scaled_rms mask cell_area count".split(' ')
    #NOTE: This skipps the dz tiling
    return fields

    fields['dz']="dz sigma_dz count misfit_rms misfit_scaled_rms mask cell_area".split(' ')

    lags=['_lag1', '_lag4', '_lag8']
    for lag in lags:
        fields['dzdt'+lag]=["dzdt"+lag, "sigma_dzdt"+lag]

    for res in ["_40000m", "_20000m", "_10000m"]:
        fields['avg_dz'+res] = ["avg_dz"+res, "sigma_avg_dz"+res,'cell_area']
        for lag in lags:
            field_str='avg_dzdt'+res+lag
            fields[field_str]=[field_str, 'sigma_'+field_str]

    #for key, item in fields.items():
    #print(key+" : "+str(item))
    #print(fields)
    return fields

def make_200km_tiles(region_dir):
    tile_ctr_file=os.path.join(region_dir,'200km_tile_list.txt')

    if os.path.isfile(tile_ctr_file):
        with open(tile_ctr_file) as fh:
            xyc=[ [*map(float, line.rstrip().split(' '))] for line in fh]
        return xyc
                
    tile_files=[]
    for sub in ['matched']:
        tile_files += glob.glob(os.path.join(region_dir, sub, 'E*N*.h5'))

    tile_list=[]
    for tile_name in tile_files:
        try:
            tile_list += [np.array([*map(int, tile_re.search(tile_name).groups())])]
        except Exception as e:
            print(tile_name)
            print(e)

    xy0=np.c_[tile_list]*1000
    xyc=pc.unique_by_rows(np.floor(xy0/tile_W)*tile_W+tile_W/2)
    with open(tile_ctr_file,'w') as fh:
        for line in xyc:
            fh.write(str(line[0])+' '+str(line[1])+'\n')
    return xyc


tile_W=2.e5

tile_re = re.compile('E(.*)_N(.*).h5')

region_dir=sys.argv[1]
region=sys.argv[2]

print("region_dir is " +region_dir)

fields=make_fields()
xyc=make_200km_tiles(region_dir)

tile_dir_200km=os.path.join(region_dir,'200km_tiles')
if not os.path.isdir(tile_dir_200km):
    os.mkdir(tile_dir_200km)

if not os.path.isdir(f'tile_run_{region}'):
    os.mkdir(f'tile_run_{region}')

non_sigma_fields={}
sigma_fields={}
for group, field_list in fields.items():
    non_sigma_fields[group]=" ".join([field for field in field_list if 'sigma' not in field])
    sigma_fields[group]=" ".join([field for field in field_list if 'sigma' in field])

for count, xy in enumerate(xyc):
    search_bounds =[xy[0]-tile_W/2-1.e4, xy[0]+tile_W/2+1.e4, xy[1]-tile_W/2-1.e4, xy[1]+tile_W/2+1.e4]
    search_bounds_str = " ".join([str(ii) for ii in search_bounds])
    tile_bounds = [xy[0]-tile_W/2, xy[0]+tile_W/2, xy[1]-tile_W/2, xy[1]+tile_W/2]
    tile_bounds_1km = "_".join([str(int(ii/1000)) for ii in tile_bounds])
    tile_bounds_str = " ".join([str(ii) for ii in tile_bounds])
    
    task_file=f'tile_run_{region}/task_{count}'
    with open(task_file,'w') as fh:
        fh.write("source activate IS2\n") 
        for group in fields.keys():
            if "40000m" in group:
                pad=0
                feather=0
                spacing_str="-S 40000 40000"
            else:
                pad=10000
                feather=10000
                spacing_str=""

            out_dir = os.path.join(tile_dir_200km, group)
            if not os.path.isdir(out_dir):
                os.mkdir(out_dir)
            out_file = os.path.join(out_dir, f"{group}{tile_bounds_1km}.h5")

            fh.write("#\n")
            fh.write(f"make_mosaic.py -w -R -d {region_dir} -g matched/E*.h5 -r {search_bounds_str} -f {feather} -p {pad} -c {tile_bounds_str} -G {group} -F {non_sigma_fields[group]} -O {out_file} {spacing_str}\n")
            fh.write(f"make_mosaic.py -w  -d {region_dir} -g prelim/E*.h5 -r {search_bounds_str} -f {feather} -p {pad} -c {tile_bounds_str} -G {group} -F {sigma_fields[group]} -O {out_file} {spacing_str}\n")
    st=os.stat(task_file)
    os.chmod(task_file, st.st_mode | stat.S_IEXEC)   

