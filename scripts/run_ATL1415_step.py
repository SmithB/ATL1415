#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Jun  6 21:00:02 2019

@author: ben
"""

#import matplotlib.pyplot as plt
import numpy as np
import pointCollection as pc
import scipy.ndimage as snd
import sys
import os
import re
import argparse

def pad_mask_canvas(D, N):
    dx=np.diff(D.x[0:2])
    left=np.arange(-N*dx,0, dx)
    right=np.arange(0, N*dx, dx)
    x1=np.unique(np.concatenate([left+D.x[0], D.x, D.x[-1]+right]))
    y1=np.unique(np.concatenate([left+D.y[0], D.y, D.y[-1]+right]))
    cols=np.flatnonzero(np.in1d(x1, D.x))
    rows=np.flatnonzero(np.in1d(y1, D.y))
    z1=np.zeros([y1.size, x1.size], dtype='bool')
    z1[rows[0]:rows[-1]+1,cols[0]:cols[-1]+1]=D.z.astype('bool')
    return pc.grid.data().from_dict({'x':x1, 'y':y1,'z':z1})

def get_xy_from_mask(args, Hxy, XR, YR):
    mask_base, mask_ext = os.path.splitext(args.mask_file)
    if mask_ext == '.tif':
        mask_re=re.compile('(_\d+m.tif)')
        m=mask_re.search(args.mask_file)
        if m is not None:
            tif_1km = args.mask_file.replace(m.group(1), '_1km.tif')
    # make a 1km tif if it does not exist:
    if mask_ext in ['.shp','.db','.h5']:
        tif_1km = args.mask_file.replace(mask_ext, '_1km.tif')
        if mask_ext in ['.shp','.db'] and not os.path.isfile(tif_1km):
            os.system(f'gdal_rasterize -init 0 -burn 1 -tr 1000 1000 -at -ot byte -co "COMPRESS=LZW" -co "PREDICTOR=1" {args.mask_file} {tif_1km}')
    temp=pc.grid.data().from_geotif(tif_1km)
    mask_G=pad_mask_canvas(temp, 200)
    # downsample the mask to 4 km
    mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([1, 5], dtype='bool'))
    mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([5, 1], dtype='bool'))
    mask_G=mask_G[2::4, 2::4]
    # blur the mask to two times the half tile
    mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([1, int(2*Hxy/4000)+1], dtype='bool'))
    mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([int(2*Hxy/4000)+1, 1], dtype='bool'))
    x0=np.unique(np.round(mask_G.x/Hxy)*Hxy)
    y0=np.unique(np.round(mask_G.y/Hxy)*Hxy)
    x0, y0 = np.meshgrid(x0, y0)
    xg=x0.ravel()
    yg=y0.ravel()
    good=(np.abs(mask_G.interp(xg, yg)-1)<0.1) & (np.mod(xg, Wxy)==0) & (np.mod(yg, Wxy)==0)

    if XR is not None:
        good &= (xg>=XR[0]) & (xg <= XR[1]) & (yg > YR[0]) & (yg < YR[1])

    xg=xg[good]
    yg=yg[good]

    return xg, yg

def rewrite_args_file(args_file, args):
    out_file=args_file.replace('.txt', f'_{args.step}.txt')
    arg_re=re.compile('(.*)=(.*)')
    arg_dict={}
    skip_args=['--region_file', '--submit']
    with open(args_file,'r') as fh:
        for line in fh:
            m=arg_re.search(line)
            if m is None:
                continue
            arg_dict[m.group(1)]=m.group(2)
    for key, val in vars(args).items():
        if '-'+key in arg_dict:
            arg_dict['-'+key]=val
        else:
            arg_dict['--'+key]=val
    
    with open(out_file,'w') as fh:
        for key, val in arg_dict.items():
            if val is None:
                continue
            if key in skip_args:
                continue
            fh.write(f'{key}={val}\n')
    return out_file


# define the script.  This is assumed to be in the path of the environment
# that is running 
prog = "ATL11_to_ATL15.py"
environment = "IS2"

replacement_args=''

# account for a bug in argparse that misinterprets negative agruents
argv=sys.argv
for i, arg in enumerate(argv):
    if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

defaults_file=None
for ii in sys.argv:    
    if ii.startswith('@'):        
        defaults_file=ii[1:]   
    
parser = argparse.ArgumentParser(description="generate a list of commands to run ATL11_to_ATL15", fromfile_prefix_chars='@')
parser.add_argument('step', type=str)
parser.add_argument('--region_file', '-R', type=str)
parser.add_argument('--mask_dir', type=str)
parser.add_argument('-W', type=float)
parser.add_argument('--skip_errors','-s', action='store_true')
parser.add_argument('--tile_spacing', type=int)
parser.add_argument('--ATL14_root', type=str)
parser.add_argument('--ATL11_index', type=str)
parser.add_argument('--region', type=str)
parser.add_argument('--Release', type=str)
parser.add_argument('--Hemisphere', type=str)
parser.add_argument('--mask_file', type=str)
parser.add_argument('--d2z0_file', type=str)
parser.add_argument('--name', type=str)
parser.add_argument('--submit', action='store_true')

args, _ = parser.parse_known_args()

if args.step not in ['centers', 'edges','corners', 'prelim','matched']:
    raise(ValueError('step argument not known: must be one of : centers, edges, corners'))
    sys.exit()

if args.skip_errors:
    calc_errors=False
else:
    calc_errors=True
    
if args.Hemisphere == 1 or args.Hemisphere == "1":
    hemisphere_name='north'
else:
    hemisphere_name='south'

XR=None
YR=None
if args.region_file is not None:
    line_re=re.compile('(..)\s*=\s*\[\s*(\S+),\s*(\S+)\s*]')
    temp={}
    with open(args.region_file,'r') as fh:
        for line in fh:
            m = line_re.search(line)
            temp[m.group(1)]=[np.float(m.group(2)), np.float(m.group(3))]
    XR=temp['XR']
    YR=temp['YR']

# check if enough parameters have been specified to allow a run
required_keys_present=True

if args.mask_dir is not None:
    args.mask_file=os.path.join(args.mask_dir, args.mask_file)
    if args.tide_mask_file is not None and not os.path.isfile(args.mask_file):
        args.tide_mask_file=os.path.join(args.mask_dir, args.tide_mask_file)

if not os.path.isfile(args.ATL11_index):
    original_index_file = args.ATL11_index
    args.ATL11_index = os.path.join(args.ATL14_root, args.ATL11_index)
    if not os.path.isfile(defaults['--ATL11_index']):
        print("could not find ATL11 index in " + args.ATL11_index + " or " + original_index_file)
        sys.exit(1)

E_d2z0=None
if args.d2z0_file is not None:
    if not os.path.isfile(args.d2z0_file):
        args.d2z0_file = os.path.join(args.mask_dir, args.d2z0_file)
    if os.path.isfile(args.d2z0_file):
        E_d2z0 = pc.grid.data().from_geotif(args.d2z0_file)

# figure out what directories we need to make
release_dir = os.path.join(args.ATL14_root, "rel"+args.Release)
hemi_dir=os.path.join(release_dir, hemisphere_name)
region_dir=os.path.join(hemi_dir, args.region)

for this in [release_dir, hemi_dir, region_dir]:
    if not os.path.isdir(this):
        print("missing directory: "+ this)
        sys.exit(1)

step_dir=os.path.join(region_dir, args.step)
if not os.path.isdir(step_dir):
    os.mkdir(step_dir)

# generate the center locations
if args.tile_spacing is None:
    Wxy=float(args.W)
    Hxy=Wxy/2
else:
    Wxy=args.tile_spacing
    Hxy=args.tile_spacing


print(f"Wxy={Wxy}")


xg, yg  = get_xy_from_mask(args, Hxy, XR, YR)


if args.step in ['centers', 'prelim','matched']:
    delta_x=[0]
    delta_y=[0]
#elif  args.step=='prelim' or args.step=='matched':
#    delta_x, delta_y = [ii.ravel() for ii in np.meshgrid([-1., 0., 1.], [-1., 0., 1.])]
elif args.step=='edges':
    delta_x=[-1, 0, 0, 1.]
    delta_y=[0, -1, 1, 0.]
elif args.step=='corners':
    delta_x=[-1, 1, -1, 1.]
    delta_y=[-1, -1, 1, 1.]
if args.step=='matched':
    calc_errors=False
    args.prior_edge_include=1000
    args.max_iterations=1

queued=[];

if args.name is not None:
    run_name=args.name
else:
    run_name=f"1415_run_{args.region}_{args.step}"
run_dir=run_name
queue_dir=f"{run_dir}/queue"
if not os.path.isdir(run_dir):
    os.mkdir(run_dir)
if not os.path.isdir(queue_dir):
    os.mkdir(queue_dir)

for name in ["logs", "running","done", "slurm_logs"]:
    if not os.path.isdir(run_dir+'/'+name):
        os.mkdir(run_dir+'/'+name)

count=0

step_file=rewrite_args_file(defaults_file, args)

xyE=open('xyE.txt','w')

for xy0 in zip(xg, yg):
    for dx, dy in zip(delta_x, delta_y):  
        xy1=np.array(xy0)+np.array([dx, dy])*Hxy
        if  np.mod(xy0[0], Hxy)>0 or np.mod(xy0[1], Hxy)>1:
            continue

        if args.step=='matched':
            in_file=os.path.join(region_dir, 'prelim', 'E%d_N%d.h5' % (xy1[0]/1000, xy1[1]/1000))
            if not os.path.isfile(in_file):
                continue

        if tuple(xy1) in queued:
            continue
        else:
            queued.append(tuple(xy1))
        
        out_file='%s/E%d_N%d.h5' % (step_dir, xy1[0]/1000, xy1[1]/1000)  
        if os.path.isfile(out_file):
            continue
        count +=1
        task_file=f'{queue_dir}/calc_dh_{count}'
        with open(task_file,'w') as fh_out:
            fh_out.write(f'source activate {environment}\n')
            if args.step=='matched':
                cmd = '%s --data_file %s --%s @%s ' % (prog, in_file, args.step, step_file)
            else:
                cmd = '%s --xy0 %d %d --%s @%s ' % (prog, xy1[0], xy1[1], args.step, step_file)
            fh_out.write(cmd+'\n')

            if calc_errors:
                fh_out.write(cmd+' --calc_error_for_xy'+'\n')
print("Wrote commands to "+queue_dir)

replacements={"[[JOB_NAME]]":run_name+'_dh', "[[TIME]]":"04:00:00", '[[NUM_TASKS]]':'3', "[[JOB_DIR]]":queue_dir, "[[JOB_NUMBERS]]":f"1-{count}", "[[PREFIX]]":"calc_dh"}

with open('/home/besmith4/slurm_files/templates/worker','r') as temp:
    with open(run_dir+'/slurm_tile_run','w') as out:
        for line in temp:
            for search, replace in replacements.items():
                line=line.replace(search, replace)
            out.write(line)

if args.submit:
    if args.step=='centers' or args.step=='prelim':
        os.system(f'cd {run_dir}; sbatch slurm_tile_run > dependents')
    else:
        if args.step=='edges':
            last_dir=run_dir.replace('edges','centers')
        elif args.step=='corners':
            last_dir=run_dir.replace('corners','edges')
        dep_job=os.popen(f'tail -1 {last_dir}/dependents').read().rstrip().split()[-1]
        os.system(f'cd {run_dir}; sbatch --dependency=afterany:{dep_job} slurm_tile_run > dependents')

