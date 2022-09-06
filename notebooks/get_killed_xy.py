#! /usr/bin/env python


import re
import numpy as np

#1415_run_AA_prelim_dh.63745071_3327

job_re=re.compile('dh.\d+_(\d+)')
xy_re=re.compile('working on .*E(.*)_N(.*).h5')

all_xy=[]

with open('killed_jobs','r') as fh:
    for line in fh:
        job_number=job_re.search(line).group(1)
        log_file=f'logs/calc_dh_{job_number}.log'
        with open(log_file,'r') as log:
            for log_line in log:
                try:
                    all_xy += [[*map(int, xy_re.search(log_line).groups())]]
                    break
                except Exception:
                    pass

print(all_xy)

#all_xy=np.c_[all_xy]

#for row in all_xy:
#    print(f'{row[0]}, {row[1]}')
