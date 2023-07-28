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

def add_files_to_queue(run_name=None, task_list_file=None, task_glob=None, shell=None, env=None):
    last_file_num=get_last_task(run_name)
    if task_list_file is not None:
        with open(task_list_file,'r') as fh:
            add_count=0
            for line in fh:
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
            print(f"added {add_count} files to the queue")

    if task_glob is not None:
        task_files=glob.glob(task_glob)
        for file in task_files:
            last_file_num=last_file_num+1;
            this_file=os.path.join(run_name,'queue','task_%d' % last_file_num)
            os.rename(file, this_file)

    with open(os.path.join(run_name,'last_task'),'w+') as last_task_fh:
        last_task_fh.write('%d\n'% last_file_num)

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_name','-r', type=str, default='ATL_run', help="name to assign to jobs and temporary directories")
    parser.add_argument('--queue_file', '-q', type=str, help="filename containing jobs, one per line")
    parser.add_argument('--task_glob', '-g', type=str, help='glob to match jobs')
    parser.add_argument('--environment','-e', type=str, default='ATL1415', help="environment that each job will activate")
    parser.add_argument('--shell','-s', type=str, default=None, help="shell to specify for each job (may not be needed)")
    parser.add_argument('--jobs_per_task','-j', type=int, default=1, help="number of jobs per node")
    parser.add_argument('--time','-t', type=str, default='02:00:00', help="time limit per job (hh:mm:ss)")
    parser.add_argument('--css', action='store_true', help="if set, the run will use the constraint=cssro argument, needed for ATL11")
    args=parser.parse_args()
    
    first_task=get_last_task(args.run_name)+1
    setup_directories(args.run_name)
    add_files_to_queue(run_name=args.run_name, task_list_file=args.queue_file, task_glob=args.task_glob, shell=args.shell, env=args.environment)
    last_task=get_last_task(args.run_name)
    ATL1415.make_slurm_file(os.path.join(args.run_name, 'slurm_script.sh'),
                          subs={'JOB_NAME':args.run_name,
                                'TIME':args.time,
                                'NUM_TASKS':str(args.jobs_per_task),
                                'JOB_NUMBERS':f'{first_task}-{last_task}'}, css=args.css)

if __name__ == '__main__':
    __main__()
