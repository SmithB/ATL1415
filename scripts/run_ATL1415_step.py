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
    if mask_ext in ('.tif'):
 
        tif_1km=args.mask_file.replace('100m','1km').replace('125m','1km')
        temp=pc.grid.data().from_geotif(tif_1km)
    
        mask_G=pad_mask_canvas(temp, 200)
        # downsample the mask to 4 km
        mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([1, 5], dtype='bool'))
        mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([5, 1], dtype='bool'))
        mask_G=mask_G[2::4, 2::4]
        # blur the mask to three times the half tile
        mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([1, int(3*Hxy/4000)+1], dtype='bool'))
        mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([int(3*Hxy/4000)+1, 1], dtype='bool'))
        x0=np.unique(np.round(mask_G.x/Hxy)*Hxy)
        y0=np.unique(np.round(mask_G.y/Hxy)*Hxy)
        x0, y0 = np.meshgrid(x0, y0)
        xg=x0.ravel()
        yg=y0.ravel()
        good=(np.abs(mask_G.interp(xg, yg)-1)<0.1) & (np.mod(xg, Wxy)==0) & (np.mod(yg, Wxy)==0)
    elif mask_ext in ['.shp','.db']:
        # the mask is a shape.
        # We require that an 80-km grid based on the mask exists
        if not os.path.isfile(mask_base+'_80km.tif'):
            raise(OSError(f"gridded mask file {mask_base+'_80km.tif'} not found"))
        mask_G=pc.grid.data().from_geotif(mask_base+'_80km.tif')
        xg, yg = np.meshgrid(mask_G.x, mask_G.y)
        xg=xg.ravel()[mask_G.z.ravel()==1]
        yg=yg.ravel()[mask_G.z.ravel()==1]
        good=np.ones_like(xg, dtype=bool)

    if XR is not None:
        good &= (xg>=XR[0]) & (xg <= XR[1]) & (yg > YR[0]) & (yg < YR[1])

    xg=xg[good]
    yg=yg[good]


    return xg, yg


# define the script.  This is assumed to be in the path of the environment
# that is running 
prog = "ATL11_to_ATL15.py"
environment = "IS2"

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
parser.add_argument('--region', type=str)
parser.add_argument('--Release', type=str)
parser.add_argument('--Hemisphere', type=str)
parser.add_argument('--mask_file', type=str)
parser.add_argument('--d2z0_file', type=str)
parser.add_argument('--slurm', action='store_true')
args, _ = parser.parse_known_args()


if args.step not in ['centers', 'edges','corners']:
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
else:
    Wxy=args.tile_spacing

Hxy=Wxy/2
xg, yg  = get_xy_from_mask(args, Hxy, XR, YR)

if args.step=='centers':
    delta_x=[0]
    delta_y=[0]
elif args.step=='edges':
    delta_x=[-1, 0, 0, 1.]
    delta_y=[0, -1, 1, 0.]
elif args.step=='corners':
    delta_x=[-1, 1, -1, 1.]
    delta_y=[-1, -1, 1, 1.]

queued=[];

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

xyE=open('xyE.txt','w')

for xy0 in zip(xg, yg):
    for dx, dy in zip(delta_x, delta_y):  
        xy1=np.array(xy0)+np.array([dx, dy])*Hxy
        if tuple(xy1) in queued:
            continue
        else:
            queued.append(tuple(xy1))
        
        if E_d2z0 is not None:
            c = np.argmin(np.abs(E_d2z0.x-xy1[0]))
            r = np.argmin(np.abs(E_d2z0.y-xy1[1]))
            E_d2z0dx2 = np.minimum(1.e-2, np.maximum(1.e-4, E_d2z0.z[r,c]))
            if np.isnan(E_d2z0dx2):
                E_d2z0dx2=1.e-2
            xyE.write(f'{xy1[0]} {xy1[1]} {E_d2z0dx2}\n')
        out_file='%s/E%d_N%d.h5' % (step_dir, xy1[0]/1000, xy1[1]/1000)  
        if os.path.isfile(out_file):
            continue
        count +=1
        task_file=f'{queue_dir}/calc_dh_{count}'
        with open(task_file,'w') as fh_out:
            fh_out.write(f'source activate {environment}\n')
            cmd = '%s --xy0 %d %d --%s @%s ' % (prog, xy1[0], xy1[1], args.step, defaults_file)
            if E_d2z0 is not None:
                cmd += f' --E_d2z0dx2={E_d2z0dx2:.6f}'
            fh_out.write(cmd+'\n')
        calc_sigma_file=f'{queue_dir}/calc_sigma_{count}'
        if calc_errors:
            calc_sigma_file=f'{queue_dir}/calc_sigma_{count}'
            with open(calc_sigma_file,'w') as fh_out:
                cmd += ' --calc_error_for_xy'
                fh_out.write(f'source activate {environment}\n')
                fh_out.write(cmd+'\n')
print("Wrote commands to "+queue_dir)

replacements={"[[JOB_NAME]]":run_name+'_dh', "[[TIME]]":"04:00:00", '[[NUM_TASKS]]':'3', "[[JOB_DIR]]":queue_dir, "[[JOB_NUMBERS]]":f"1-{count}", "[[PREFIX]]":"calc_dh"}

with open('/home/besmith4/slurm_files/templates/worker','r') as temp:
    with open(run_dir+'/slurm_tile_run','w') as out:
        for line in temp:
            for search, replace in replacements.items():
                line=line.replace(search, replace)
            out.write(line)

if args.slurm:
    if args.step=='centers':
        os.system(f'cd {run_dir}; sbatch slurm_tile_run > dependents')
    else:
        if args.step=='edges':
            last_dir=run_dir.replace('edges','centers')
        elif args.step=='corners':
            last_dir=run_dir.replace('corners','edges')
        dep_job=os.popen(f'tail -1 {last_dir}/dependents').read().rstrip().split()[-1]
        os.system(f'cd {run_dir}; sbatch --dependency=afterany:{dep_job} slurm_tile_run > dependents')

    replacements={"[[JOB_NAME]]":run_name+'_sigma', "[[TIME]]":"03:00:00", '[[NUM_TASKS]]':'3', "[[JOB_DIR]]":queue_dir, "[[JOB_NUMBERS]]":f"1-{count}", "[[PREFIX]]":"calc_sigma"}

    with open('/home/besmith4/slurm_files/templates/worker','r') as temp:
        with open(run_dir+'/slurm_sigma_run','w') as out:
            for line in temp:
                for search, replace in replacements.items():
                    line=line.replace(search, replace)
            out.write(line)


    dep_job=os.popen(f'tail -1 {run_dir}/dependents').read().rstrip().split()[-1]
    os.system(f'cd {run_dir}; sbatch --dependency=afterany:{dep_job} slurm_sigma_run >> dependents')


