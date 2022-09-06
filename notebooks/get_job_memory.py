#! /usr/bin/env python


import re
import numpy as np
import glob
#1415_run_AA_prelim_dh.63745071_3327

job_re=re.compile('dh.\d+_(\d+)')

mem_ex_re = re.compile('memory limit\s+.(\d+)')
mem_re = re.compile('Memory Used.* (\d+)K')

job_mem={}


for file in glob.glob('1415_run*'):
    try:
        job_num=job_re.search(file).group(1)
        job_mem[job_num]=[0, 0]
    except Exception as e:
        print(file)
        print(e)
        continue
    with open(file,'r') as fh:
        for line in fh:
            try:
                job_mem[job_num][1]=mem_ex_re.search(line).group(1)
            except Exception:
                pass
            try:
                job_mem[job_num][0]=mem_re.search(line).group(1)
                #print(mem_re.search(line).groups())
            except Exception:
                pass


for job_num, mem in job_mem.items():
    print(f'{job_num} {mem[0]} {mem[1]}')

