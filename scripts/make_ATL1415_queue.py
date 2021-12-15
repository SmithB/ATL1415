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


# define the script.  This is assumed to be in the path of the environment
# that is running 
prog = "ATL11_to_ATL15.py"
environment = "IS2"

# account for a bug in argparse that misinterprets negative agruents
argv=sys.argv
for i, arg in enumerate(argv):
    if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg


parser = argparse.ArgumentParser(description="generate a list of commands to run ATL11_to_ATL15")
parser.add_argument('step', type=str)
parser.add_argument('defaults_files', nargs='+', type=str)
parser.add_argument('--region_file', '-R', type=str)
parser.add_argument('--skip_errors','-s', action='store_true')
parser.add_argument('--tile_spacing', type=int)
args = parser.parse_args()

if args.step not in ['centers', 'edges','corners']:
    raise(ValueError('step argument not known: must be one of : centers, edges, corners'))
    sys.exit()

if args.skip_errors:
    calc_errors=False
else:
    calc_errors=True
    
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

defaults_re=re.compile('(.*)\s*=\s*(.*)')

# read in all defaults files (must be of syntax --key=value or -key=value)
defaults={}

for defaults_file in args.defaults_files:
    with open(defaults_file,'r') as fh:
        for line in fh:
            m=defaults_re.search(line)
            if m is not None:
                defaults[m.group(1)]=m.group(2)

# check if enough parameters have been specified to allow a run
required_keys_present=True
for key in ['--ATL14_root', '--region', '--Release','--Hemisphere', '--mask_file']:
    if key not in defaults:
        print(f"make_1415_queue.py:\n\tError: required key {key} not in defaults files")
        required_keys_present=False
if not required_keys_present:
    sys.exit(1)

if '--mask_dir' in defaults:
    defaults['--mask_file']=os.path.join(defaults['--mask_dir'], defaults['--mask_file'])   
    if '--tide_mask_file' in defaults and not os.path.isfile(defaults['--tide_mask_file']):
        defaults['--tide_mask_file']=os.path.join(defaults['--mask_dir'], defaults['--tide_mask_file'])
    defaults.pop('--mask_dir', None)
    

if defaults['--Hemisphere']==1 or defaults['--Hemisphere']=="1":
    hemisphere_name='north'
else:
    hemisphere_name='south'

# figure out what directories we need to make
release_dir = os.path.join(defaults['--ATL14_root'], "rel"+defaults['--Release'])
hemi_dir=os.path.join(release_dir, hemisphere_name)
region_dir=os.path.join(hemi_dir, defaults['--region'])

for this in [release_dir, hemi_dir, region_dir]:
    if not os.path.isdir(this):
        print("missing directory: "+ this)
        sys.exit(1)

# write out the composite defaults file
defaults_file=os.path.join(region_dir, f'input_args_{defaults["--region"]}.txt')
with open(defaults_file, 'w') as fh:
    for key, val in defaults.items():
        fh.write(f'{key}={val}\n')
    fh.write(f"-b={region_dir}\n")

step_dir=os.path.join(region_dir, args.step)
if not os.path.isdir(step_dir):
    os.mkdir(step_dir)

# generate the center locations
if args.tile_spacing is None:
    Wxy=float(defaults['-W'])
else:
    Wxy=args.tile_spacing

Hxy=Wxy/2

mask_base, mask_ext = os.path.splitext(defaults['--mask_file'])
if mask_ext in ('.tif'):
    
    tif_1km=defaults['--mask_file'].replace('100m','1km').replace('125m','1km')
    temp=pc.grid.data().from_geotif(tif_1km)
        
    mask_G=pad_mask_canvas(temp, 200)
    mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([1, int(3*Hxy/1000)+1], dtype='bool'))
    mask_G.z=snd.binary_dilation(mask_G.z, structure=np.ones([int(3*Hxy/1000)+1, 1], dtype='bool'))
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
queue_file=f"1415_queue_{defaults['--region']}_{args.step}.txt"

with open(queue_file,'w') as qh:
    for xy0 in zip(xg, yg):
        for dx, dy in zip(delta_x, delta_y):  
            xy1=np.array(xy0)+np.array([dx, dy])*Hxy
            if tuple(xy1) in queued:
                continue
            else:
                queued.append(tuple(xy1))
            out_file='%s/E%d_N%d.h5' % (step_dir, xy1[0]/1000, xy1[1]/1000)  
            if not os.path.isfile(out_file):
                cmd='%s --xy0 %d %d --%s @%s ' % (prog, xy1[0], xy1[1], args.step, defaults_file)
                if calc_errors:
                    cmd += '; '+cmd+' --calc_error_for_xy'
                qh.write(f'source activate {environment}; '+cmd+'; echo COMPLETE\n')
print("Wrote commands to "+queue_file)


