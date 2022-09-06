#! /usr/bin/env python


import re
import numpy as np
import glob
#1415_run_AA_prelim_dh.63745071_3327

job_re=re.compile('dh_(\d+)')
xy_re=re.compile('working on .*E(.*)_N(.*).h5')



log_re=re.compile('(\d+).log')



xyJ=[]


for log_file in glob.glob('logs/*.log'):
    #print(log_file)
    job_num=job_re.search(log_file).group(1)
    with open(log_file,'r') as log:
        for log_line in log:
            try:
                xy=xy_re.search(log_line).groups()[0:2]
                xyJ += [(xy[0], xy[1], job_num)]
                print('%s %s %s' % xyJ[-1])
                break
            except Exception:
                pass

print(all_xy)


