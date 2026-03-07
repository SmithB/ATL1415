#!/usr/bin/env python

import glob
import pointCollection as pc
import matplotlib.pyplot as plt
import numpy as np
import glob
import shutil
import re
import os

import sys

def get_200km_loc(filename, field):
    file=os.path.basename(filename)
    bounds=np.array([*map(int, file.replace(field,'').replace('.h5','').split('_'))])
    return [bounds[0:2], bounds[2:]], np.array([np.mean(bounds[0:2]), np.mean(bounds[2:])])

def get_quadrant(xy):
    if xy.ndim==1:
        quadrant=np.zeros([1])
    else:
        quadrant=np.zeros(xy.shape[0], dtype=int)
    ii = (xy[0] >= 0) & (xy[1] >= 0)
    quadrant[ii] = 1
    ii = (xy[0] <0) & (xy[1] >= 0)
    quadrant[ii] = 2
    ii = (xy[0] < 0) & (xy[1] < 0)
    quadrant[ii] = 3
    ii = (xy[0] >= 0) & (xy[1] < 0)
    quadrant[ii] = 4
    return quadrant

def mkdir(thedir):
    if not os.path.isdir(thedir):
        os.mkdir(thedir)

def main():
    release=sys.argv[1]
    top=f'/discover/nobackup/projects/icesat2/ATL14_processing/rel{release}/south_monthly/'
    assert(os.path.isdir(top))
    tile_re=re.compile('E(.*)_N(.*).h5')

    # # make symbolic links to 200km tiles for each region

    for quadrant in range(1,5):
        mkdir(top+f'/A{quadrant}')
        mkdir(top+f'/A{quadrant}/200km_tiles')
        mkdir(top+f'/A{quadrant}/prelim')
        mkdir(top+f'/A{quadrant}/matched')

    # get the subdirectories, replicate them into the quadrant directories
    subs_200=glob.glob(top+'/AA/200km_tiles/*')
    for quadrant in range(1,5):
        quad_200km=top+f'/A{quadrant}/200km_tiles'
        for sub in subs_200:
            new_dir=quad_200km+'/'+os.path.basename(sub)
            if not os.path.isdir(new_dir):
                os.mkdir(new_dir)

    src_subs=sorted(glob.glob(top+'AA/200km_tiles/*'))
    for src_sub in src_subs:
        field=os.path.basename(src_sub)
        tile_files=glob.glob(src_sub+'/*')
        for tile_file in tile_files:
            bds, ctr=get_200km_loc(tile_file, field)
            quad = get_quadrant(ctr)
            dst=top+f'A{int(quad[0])}/200km_tiles/{field}/{os.path.basename(tile_file)}'
            if not os.path.islink(dst):
                os.symlink(tile_file, dst)

    xy_re=re.compile('E(.*)_N(.*).h5')
    for step in ['prelim','matched']:
        for tile in glob.glob(top+'/AA/'+step+'/E*.h5'):
            xy=np.array([*map(int, xy_re.search(tile).groups())])
            quad=get_quadrant(xy)
            dst=top+f'A{int(quad[0])}/{step}/{os.path.basename(tile)}'
            if not os.path.islink(dst):
                os.symlink(tile, dst)

    # make the bounds for each quadrant
    files=glob.glob(top+'AA/prelim/E*')
    xy=[[*map(int, xy_re.search(file).groups())] for file in files]
    xy=np.c_[xy]
    full_bounds=[[np.min(xy[:,0])-40, np.max(xy[:,0])+40],[np.min(xy[:,1])-40, np.max(xy[:,1])+40]]
    full_bounds=[1000*np.array(jj) for jj in  full_bounds]

    with open(f'{top}/A1/bounds.txt','w') as fh:
        fh.write('%d %d %d %d\n' % (0, full_bounds[0][1], 0, full_bounds[1][1]))
    with open(f'{top}/A2/bounds.txt','w') as fh:
        fh.write('%d %d %d %d\n' % (full_bounds[0][0],0,  0, full_bounds[1][1]))
    with open(f'{top}/A3/bounds.txt','w') as fh:
        fh.write('%d %d %d %d\n' % (full_bounds[0][0], 0, full_bounds[1][0],0))
    with open(f'{top}/A4/bounds.txt','w') as fh:
        fh.write('%d %d %d %d\n' % (0, full_bounds[0][1], full_bounds[1][0],0))

    # ## Make the input arguments for each quadrant

    src=f'{top}/AA/input_args_AA.txt'

    for sub in range(1, 5):
        this_dir=f'{top}/A{sub}'
        dst=this_dir+f'/input_args_A{sub}.txt'
        with open(src,'r') as fh_in:
            with open(dst,'w') as fh_out:
                for line in fh_in:
                    if line.startswith('-b=') or line.startswith('--region='):
                        fh_out.write(line.replace('AA', f'A{sub}'))
                    else:
                        fh_out.write(line)




