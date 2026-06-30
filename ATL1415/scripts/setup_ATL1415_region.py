#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 21:00:02 2019

@author: ben
"""

import sys
import os
import re
import glob
import argparse

def main(argv=None):

    # account for a bug in argparse that misinterprets negative agruents

    if argv is None:
        argv=sys.argv
    for i, arg in enumerate(argv):
        if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg


    parser = argparse.ArgumentParser(description="generate directories and defaults files for ATL11_to_ATL15")
    parser.add_argument('defaults_files', nargs='+', type=str)
    parser.add_argument('--ATL14_reference_file', type=str)
    parser.add_argument('--Hemisphere', type=int, choices=[1, -1],
                        help='hemisphere: 1=north, -1=south')
    parser.add_argument('--previous_product_top', type=str, default=None,
                        help='top-level directory of previous ATL14/15 products; '
                             'subdirs are north/<REGION> (Arctic) or south/A1..A4 (Antarctic)')

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

    if args.ATL14_reference_file is not None:
        defaults['--ATL14_reference_file']=args.ATL14_reference_file

    if args.Hemisphere is not None:
        defaults['--Hemisphere'] = str(args.Hemisphere)

    # check if enough parameters have been specified to allow a run
    required_keys_present=True
    for key in ['--ATL14_root', '--region', '--Release','--Hemisphere', '--mask_file']:
        if key not in defaults:
            print(f"make_1415_queue.py:\n\tError: required key {key} not in defaults files")
            required_keys_present=False
    if not required_keys_present:
        sys.exit(1)

    if '--mask_dir' in defaults:
        for key in ['--mask_file','--d2z0_file','--tide_mask_file', '--tide_adjustment_file', '--geoid_file', '--E_d2z0dx2_file', '--E_d3zdx2dt_scale_file']:
            if key in defaults and not (os.path.isabs(defaults[key]) and  os.path.isfile(defaults[key])):
                defaults[key] = os.path.join(defaults['--mask_dir'], defaults[key])
        defaults.pop('--mask_dir', None)


    if defaults['--Hemisphere']==1 or defaults['--Hemisphere']=="1":
        hemisphere_base = 'north'
    else:
        hemisphere_base = 'south'

    hemisphere_name = hemisphere_base
    if '--hemi_suffix' in defaults:
        hemisphere_name += defaults['--hemi_suffix']
        
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


    # if ATL11 release is specified and ATL11 geoindex is not specified, build the location
    if '--ATL11_index' not in defaults and '--ATL11_release' in defaults:
        defaults['--ATL11_index'] = os.path.join(defaults['--ATL14_root'], defaults['--ATL11_release'], hemisphere_base, 'index','GeoIndex.h5')
        defaults.pop('--ATL11_release')

    if not os.path.isfile(defaults['--ATL11_index']):
        raise(OSError(f"ATL11 index file {defaults['--ATL11_index']} does not exist"))

    # derive xover dir from the index location unless explicitly specified
    if '--ATL11_xover_dir' not in defaults:
        defaults['--ATL11_xover_dir'] = os.path.join(
            os.path.dirname(os.path.dirname(defaults['--ATL11_index'])), 'xover')

    # resolve previous-product directories if a top-level path was given
    pp_dirs = []
    if args.previous_product_top is not None:
        if hemisphere_base == 'north':
            pp_dirs = [os.path.join(args.previous_product_top, 'north', defaults['--region'])]
        else:
            south_dir = os.path.join(args.previous_product_top, 'south')
            pp_dirs = sorted(d for d in glob.glob(os.path.join(south_dir, '*'))
                             if os.path.isdir(d))

    # write out the composite defaults file
    defaults_file=os.path.join(region_dir, f'input_args_{defaults["--region"]}.txt')
    with open(defaults_file, 'w') as fh:
        for key, val in defaults.items():
            if key in ["--hemi_suffix"]:
                continue
            fh.write(f'{key}={val}\n')
        for pp_dir in pp_dirs:
            fh.write(f'--previous_product={pp_dir}\n')
        fh.write(f"-b={region_dir}\n")

    print("setup_ATL1415_region.py: wrote defaults to:")
    print("\t"+defaults_file)


if __name__=='__main__':
    main()
