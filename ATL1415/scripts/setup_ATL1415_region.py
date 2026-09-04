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
import ATL1415
from ATL1415.paths import join_path_or_uri, exists as path_exists

# marks a bare store_true flag (e.g. --tide_adjustment) in the defaults dict,
# as distinct from a key with an empty value
FLAG = object()

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

    args = parser.parse_args()

    defaults_re=re.compile(r'(.*?)\s*=\s*(.*)')

    # read in all defaults files.  A line is either '--key=value' or a bare
    # '--flag' for a store_true option; a bare flag is recorded with the value
    # FLAG and written back out as a bare line.  Without that, an options file
    # asking for --tide_adjustment or --ATL11_earthaccess would have the
    # request silently dropped, since the key=value regex does not match it.
    defaults={}

    for defaults_file in args.defaults_files:
        with open(defaults_file,'r') as fh:
            for line in fh:
                line=line.strip()
                if not line or line.startswith('#'):
                    continue
                m=defaults_re.search(line)
                if m is not None:
                    defaults[m.group(1).strip()]=m.group(2).strip()
                elif line.startswith('-'):
                    defaults[line]=FLAG

    if args.ATL14_reference_file is not None:
        defaults['--ATL14_reference_file']=args.ATL14_reference_file

    if args.Hemisphere is not None:
        defaults['--Hemisphere'] = str(args.Hemisphere)

    # check if enough parameters have been specified to allow a run
    required_keys_present=True
    for key in ['--ATL14_root', '--region', '--Release','--Hemisphere', '--mask_file']:
        if key not in defaults:
            print(f"setup_ATL1415_region.py:\n\tError: required key {key} not in defaults files")
            required_keys_present=False
    if not required_keys_present:
        sys.exit(1)

    if '--mask_dir' in defaults:
        # join_path_or_uri leaves an absolute path or a URI alone and joins
        # only a relative one, so --mask_dir may itself be an s3:// prefix.
        # (os.path.join would have produced '<mask_dir>/s3://bucket/key',
        # since a URI does not start with '/' and so does not read as absolute.)
        for key in ['--mask_file','--d2z0_file','--tide_mask_file', '--tide_adjustment_file', '--geoid_file', '--E_d2z0dx2_file', '--E_d3zdx2dt_scale_file']:
            if key in defaults:
                defaults[key] = join_path_or_uri(defaults['--mask_dir'], defaults[key])
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

    # --ATL11_earthaccess reinterprets --ATL11_index as the ROOT holding one
    # 'ATL11_index_<cycles>_<release>_<version>/' subtree of per-granule
    # geoIndex files, rather than as a single whole-archive GeoIndex.h5.  It is
    # therefore a directory (or an s3:// prefix), it cannot be derived from
    # --ATL11_release, and it must be given explicitly.
    cloud_ATL11 = '--ATL11_earthaccess' in defaults

    # if ATL11 release is specified and ATL11 geoindex is not specified, build the location
    if not cloud_ATL11 and '--ATL11_index' not in defaults and '--ATL11_release' in defaults:
        atl11_release = defaults['--ATL11_release']
        if not atl11_release.startswith('ATL11_'):
            atl11_release = f'ATL11_{atl11_release}'
        defaults['--ATL11_index'] = os.path.join(defaults['--ATL14_root'], atl11_release, hemisphere_base, 'index','GeoIndex.h5')
        defaults.pop('--ATL11_release')

    if '--ATL11_index' not in defaults:
        if cloud_ATL11:
            raise(OSError('--ATL11_index must be given explicitly when --ATL11_earthaccess '
                          'is set: it is the root of the per-granule geoIndex tree, and '
                          'cannot be derived from --ATL11_release'))
        raise(OSError('--ATL11_index must be given, or --ATL11_release so it can be derived'))

    # os.path.isfile() is both too strict for a cloud run (the index root is a
    # directory) and always False for a URI, so check for existence instead
    if not path_exists(defaults['--ATL11_index']):
        what = 'ATL11 index root' if cloud_ATL11 else 'ATL11 index file'
        raise(OSError(f"{what} {defaults['--ATL11_index']} does not exist"))

    # derive xover dir from the index location unless explicitly specified.
    # The derivation walks up from '<...>/index/GeoIndex.h5', which only makes
    # sense for a local index file -- a cloud run must name it explicitly.
    if '--ATL11_xover_dir' not in defaults and not cloud_ATL11:
        defaults['--ATL11_xover_dir'] = os.path.join(
            os.path.dirname(os.path.dirname(defaults['--ATL11_index'])), 'xover_tiles')

    # infer dzdt lags from grid spacing / time span if not explicitly given
    if '--dzdt_lags' not in defaults and '-g' in defaults and '-t' in defaults:
        dt_str = defaults['-g'].split(',')[2]
        if '/' in dt_str:
            num, den = dt_str.split('/')
            t_res = float(num) / float(den)
        else:
            t_res = float(dt_str)
        time_span = [*map(float, defaults['-t'].split(','))]
        lags = ATL1415.infer_dzdt_lags(t_res, time_span)
        defaults['--dzdt_lags'] = ','.join(map(str, lags))

    # resolve previous-product directories if a top-level path was given
    pp_dirs = []
    if '--previous_product_top' in defaults: 
        if hemisphere_base == 'north':
            pp_dirs = [os.path.join(defaults['--previous_product_top'], 'north', defaults['--region'])]
        else:
            south_dir = os.path.join(defaults['--previous_product_top'], 'south')
            pp_dirs = sorted(d for d in glob.glob(os.path.join(south_dir, 'A?'))
                             if os.path.isdir(d))

    # write out the composite defaults file
    defaults_file=os.path.join(region_dir, f'input_args_{defaults["--region"]}.txt')
    with open(defaults_file, 'w') as fh:
        for key, val in defaults.items():
            if key in ["--hemi_suffix"]:
                continue
            if val is FLAG:
                fh.write(f'{key}\n')
            else:
                fh.write(f'{key}={val}\n')
        for pp_dir in pp_dirs:
            fh.write(f'--previous_product={pp_dir}\n')
        fh.write(f"-b={region_dir}\n")

    print("setup_ATL1415_region.py: wrote defaults to:")
    print("\t"+defaults_file)


if __name__=='__main__':
    main()
