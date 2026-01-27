#! /usr/bin/env python


import sys
import h5py
import os
import glob
import numpy as np
import json

def main():

    if len(sys.argv)==3:
        tile_root=sys.argv[1]
        region=sys.argv[2]
        step=sys.argv[3]
        step_dir=os.path.join(tile_root, region, step)
    else:
        step_dir=sys.argv[1]

    #print(step_dir)
    #sys.exit()

    dst_directory=os.path.join(step_dir,  'field_sizes')
    if not os.path.isdir(dst_directory):
        os.mkdir(dst_directory)

    files=glob.glob(os.path.join(step_dir,'E*.h5'))
    for count, thefile in enumerate(files):
        if np.mod(count, 500)==0:
            print(f'{count} out of {len(files)}')
        report={'file':thefile}
        try:
            with h5py.File(thefile,'r') as h5f:
                for field in ['dz/dz','dz/sigma_dz']:
                    try:
                        report[field]=list(h5f[field].shape)
                    except Exception:
                        report[field]=None
        except Exception:
            pass
        report_file=os.path.join(dst_directory, os.path.basename(thefile).replace('.h5','_report.json'))
        with open(report_file,'w') as fh:
            json.dump(report, fh)

if __name__=='__main__':
    main()
