#! /usr/bin/env python
import glob
import h5py
import argparse
import numpy as np
import re
import os
parser=argparse.ArgumentParser(\
            fromfile_prefix_chars="@")
parser.add_argument('--glob_str','-g', type=str, required=True)
parser.add_argument('--defaults_file','-d', type=str, required=True)
parser.add_argument('--step', type=str)
args, _= parser.parse_known_args()

xy_re=re.compile('E(.*)_N(.*).h5')

pad=np.array([-2.e3, 2.e3])
with open('large_sigma_queue.txt','w') as fh_large:
    with open('small_sigma_queue.txt','w') as fh_small:
        for file in glob.glob(args.glob_str):
            with h5py.File(file,'r') as h5f:
                if 'sigma_dz' in h5f['dz']:
                    continue

            xy0=np.array([*map(float, xy_re.search(file).groups())])*1000
            if xy0[0]**2 + xy0[1]**2 > (800e3)**2:
                fh_small.write(f"source activate IS2; ATL11_to_ATL15.py --xy0 {xy0[0]} {xy0[1]} --{args.step} @{args.defaults_file} --calc_error_for_xy; echo COMPLETE\n")
            else:
                fh_large.write(f"source activate IS2; ATL11_to_ATL15.py --xy0 {xy0[0]} {xy0[1]} --{args.step} @{args.defaults_file} --calc_error_for_xy; echo COMPLETE\n")



