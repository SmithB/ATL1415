#! /usr/bin/env python

import os 
import glob
import h5py
import sys

def_file=sys.argv[1]
region=sys.argv[2]

thedir=os.path.dirname(def_file)

count=0

if not os.path.isdir(region+'_sigma_calc'):
    os.mkdir(region+'_sigma_calc')
    os.mkdir(region+'_sigma_calc'+'/queue')
    os.mkdir(region+'_sigma_calc'+'/running')
    os.mkdir(region+'_sigma_calc'+'/done')
    os.mkdir(region+'_sigma_calc'+'/logs')
    
queue_dir=region+'_sigma_calc'+'/queue'

for step in ['corners', 'edges', 'centers']:
    files=glob.glob(os.path.join(thedir, step, 'E*.h5'))
    for file in files:
        with h5py.File(file, 'r') as h5f:
            if 'sigma_dz' in h5f['/dz/'].keys():
                continue
        count += 1
        with open(os.path.join(queue_dir, 'calc_sigma_'+str(count)),'w') as qh:
           qh.write('source activate IS2\n')
           qh.write('ATL11_to_ATL15.py --calc_error_file '+file+' @'+def_file+'\n')


       
