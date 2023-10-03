#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 00:04:37 2023

@author: ben
"""

import glob
import ATL1415
import argparse
import os
import stat
import numpy as np
import re

def setup_directories(run_dir):
    # setup the directories
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)

    sub_list=['queue','running','done', 'logs']
    for sub in sub_list:
        thedir=os.path.join(run_dir, sub)
        if not(os.path.isdir(thedir)):
            os.mkdir(thedir)

def get_last_task(run_name):
    task_file=os.path.join(run_name, 'last_task')
    if os.path.isfile(task_file):
        with open(task_file, 'r') as fh:
            for line in fh:
                temp=line;
            last_file_num=int(temp);
    else:
        last_file_num=0;
    return last_file_num

def add_files_to_queue(run_name=None, task_list_file=None, task_glob=None, shell=None, env=None, R_range=None):
    last_file_num=get_last_task(run_name)
    xy0_re=re.compile('--xy0\s+(\S+)\s+(\S+)')
    if task_list_file is not None:
        add_count=0
        with open(task_list_file,'r') as fh:
            for line in fh:
                if R_range is not None:
                    xy=np.array([*map(float, xy0_re.search(line).groups())])
                    R2=np.sum(xy**2)
                    if (R2 < R_range[0]**2)  | (R2 >= R_range[1]**2) : 
                        continue
                last_file_num=last_file_num+1;
                this_file=os.path.join(run_name,'queue','task_%d' % last_file_num)
                with open(this_file,'w') as out_fh:
                    #print("adding %s to queue" % this_file)
                    add_count +=1
                    if shell is not None:
                        out_fh.write(f'#! /usr/bin/env {shell}\n')
                    if env is not None and "source activate" not in line:
                        out_fh.write("source activate %s\n" % env)
                    out_fh.write('%s\n'% line.rstrip());
                os.chmod(this_file, os.stat(this_file).st_mode | stat.S_IEXEC)


    if task_glob is not None:
        task_files=glob.glob(task_glob)
        for file in task_files:
            if R_range is not None:
                skip=False
                with open(file,'r') as fh:
                    for line in fh:
                        m=xy0_re.search(line).groups()
                        if m is None:
                            continue
                        xy=np.array([*map(float, m.groups())])
                        R2=np.sum(xy**2)
                        if (R2 < R_range[0]**2)  | (R2 >= R_range[1]) : 
                            skip=True
                if skip:
                    continue
            last_file_num=last_file_num+1;
            this_file=os.path.join(run_name,'queue','task_%d' % last_file_num)
            add_count += 1
            os.rename(file, this_file)
    print(f"added {add_count} files to the queue")
    with open(os.path.join(run_name,'last_task'),'w+') as last_task_fh:
        last_task_fh.write('%d\n'% last_file_num)
    return add_count

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_name','-r', type=str, default='ATL_run', help="name to assign to jobs and temporary directories")
    parser.add_argument('--queue_file', '-q', type=str, help="filename containing jobs, one per line")
    parser.add_argument('--task_glob', '-g', type=str, help='glob to match jobs')
    parser.add_argument('--environment','-e', type=str, default='ATL1415', help="environment that each job will activate")
    parser.add_argument('--shell','-s', type=str, default=None, help="shell to specify for each job (may not be needed)")
    parser.add_argument('--jobs_per_task','-j', type=int, help="number of jobs per node")
    parser.add_argument('--time','-t', type=str, default='02:00:00', help="time limit per job (hh:mm:ss)")
    parser.add_argument('--css', action='store_true', help="if set, the run will use the constraint=cssro argument, needed for ATL11")
    args=parser.parse_args()

    R_dict=None
    R_vals=[0, 1.e7]
    if args.jobs_per_task is None:
        R_dict={1.e8:{'dir':'queue_3cpu', 'ncpu':3, 'xy':[], 'src_file':[]},
        6.e5:{'dir':'queue_14cpu','ncpu':14, 'xy':[], 'src_file':[]},
        3.0e5:{'dir':'queue_20cpu','ncpu':20, 'xy':[], 'src_file':[]},
        2.e5:{'dir':'queue_28cpu','ncpu':28, 'xy':[], 'src_file':[]}}
        R_vals=list(R_dict.keys())[::-1]
    setup_directories(args.run_name)
    for ii in range(len(R_vals)-1):
        first_task=get_last_task(args.run_name)
        if R_dict is not None:
            R_range=[R_vals[ii], R_vals[ii+1]]
            slurm_file='slurm_run_'+ R_dict[R_vals[ii]]['dir']+'.sh'
            N_tasks=R_dict[R_vals[ii]]['ncpu']
        else:
            R_range=None
            slurm_file='slurm_run.sh'
            N_tasks=args.jobs_per_task
        N_added = add_files_to_queue(run_name=args.run_name, task_list_file=args.queue_file, task_glob=args.task_glob, shell=args.shell, env=args.environment, R_range=R_range)
        last_task=get_last_task(args.run_name)
        if N_added <1:
            continue
        ATL1415.make_slurm_file(os.path.join(args.run_name, slurm_file),
                                subs={'JOB_NAME':args.run_name,
                                      'TIME':args.time,
                                      'NUM_TASKS':str(N_tasks),
                                      'JOB_NUMBERS':f'{first_task+1}-{last_task}'}, css=args.css)


if __name__ == '__main__':
    __main__()
