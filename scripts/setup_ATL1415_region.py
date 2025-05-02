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


# account for a bug in argparse that misinterprets negative agruents
argv=sys.argv
for i, arg in enumerate(argv):
    if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg


parser = argparse.ArgumentParser(description="generate directories and defaults files for ATL11_to_ATL15")
parser.add_argument('defaults_files', nargs='+', type=str)

args = parser.parse_args()
 
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
    for key in ['--mask_file','--d2z0_file','--tide_mask_file', '--tide_adjustment_file', '--geoid_file', '--E_d2z0dx2_file']:
        if key in defaults and not (os.path.isabs(defaults[key]) and  os.path.isfile(defaults[key])):
            defaults[key] = os.path.join(defaults['--mask_dir'], defaults[key])   
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
        os.mkdir(this)

# ATL11 index may be specified relative to ATL1r_root
if '--ATL11_index' in defaults and not os.path.isfile(defaults['--ATL11_index']):
    temp1=os.path.join(defaults['--ATL14_root'], defaults['--ATL11_index'])
    if os.path.isfile(temp1):
        defaults['--ATL11_index']=temp1
    else:
        print(temp1 + ' not found')
        temp2=os.path.join(os.path.dirname(defaults['--ATL14_root']), defaults['--ATL11_index'])
        print(f'looking for {temp2}')
        if os.path.isfile(temp2):
            defaults['--ATL11_index']=temp2
        else:
            raise(OSError(f"ATL11 index file {defaults['--ATL11_index']} does not exist in \n\t{temp1}\n\t or\n\t {temp2}"))


# if ATL11 release is specified and ATL11 geoindex is not specified, guess the location
if '--ATL11_index' not in defaults and '--ATL11_release' in defaults:
    defaults['--ATL11_index'] = os.path.join(defaults['--ATL14_root'], 'ATL11_'+defaults['--ATL11_release'], hemisphere_name, 'index','GeoIndex.h5')
    defaults.pop('--ATL11_release')

if not os.path.isfile(defaults['--ATL11_index']):
    raise(OSError(f"ATL11 index file {defaults['--ATL11_index']} does not exist"))        
        
# write out the composite defaults file
defaults_file=os.path.join(region_dir, f'input_args_{defaults["--region"]}.txt')
with open(defaults_file, 'w') as fh:
    for key, val in defaults.items():
        fh.write(f'{key}={val}\n')
    fh.write(f"-b={region_dir}\n")

print("setup_ATL1415_region.py: wrote defaults to:")
print("\t"+defaults_file)

