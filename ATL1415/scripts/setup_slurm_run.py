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

    sub_list=['queue','running','done', 'logs', 'active_logs', 'error_logs']
    for sub in sub_list:
        thedir=os.path.join(run_dir, sub)
        if not(os.path.isdir(thedir)):
            os.mkdir(thedir)

def get_last_task(run_name):
    task_file=os.path.join(run_name, 'last_task')
    if os.path.isfile(task_file):
        temp='0'
        with open(task_file, 'r') as fh:
            for line in fh:
                temp=line
        last_file_num=int(temp)
    else:
        last_file_num=0
    return last_file_num

def add_files_to_queue(run_name=None, task_list_file=None,  shell=None, env=None, max_jobs=None, lines_per_task=1, line_count=None):
    last_file_num=get_last_task(run_name)

    add_count=0
    with open(task_list_file,'r') as fh:
        this_count=0
        task_lines = []
        for line_number, line in enumerate(fh):
            # fast forward to first requested line
            if line_number <= line_count[0]:
                continue

            task_lines.append(line)
            this_count += 1
            if this_count >= lines_per_task:
                # we have added enough jobs to the current file, write out
                last_file_num=last_file_num+1;
                this_file=os.path.join(run_name,'queue','task_%d' % last_file_num)
                add_count += 1
                with open(this_file,'w') as out_fh:
                    #print("adding %s to queue" % this_file)
                    if shell is not None:
                        out_fh.write(f'#! /usr/bin/env {shell}\n')
                    if env is not None and "source activate" not in line:
                        out_fh.write("source activate %s\n" % env)
                    for line in task_lines:
                        out_fh.write('%s\n'% line.rstrip());
                os.chmod(this_file, os.stat(this_file).st_mode | stat.S_IEXEC)
                # reset count and empty task list
                this_count = 0
                task_lines = []
                if add_count >= max_jobs:
                    break

    print(f"added {add_count} files to the queue")
    with open(os.path.join(run_name,'last_task'),'w+') as last_task_fh:
        last_task_fh.write('%d\n'% last_file_num)
    line_count[0]=line_number
    return add_count

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_name','-r', type=str, default='ATL_run', help="name to assign to jobs and temporary directories")
    parser.add_argument('--queue_file', '-q', type=str, help="filename containing jobs, one per line")
    parser.add_argument('--environment','-e', type=str, default='ATL1415', help="environment that each job will activate")
    parser.add_argument('--shell','-s', type=str, default=None, help="shell to specify for each job (may not be needed)")
    parser.add_argument('--lines_per_task', type=int, default=1, help="combine this number of lines from the input into each task")
    parser.add_argument('--jobs_per_task','-j', type=int, required=True, help="number of jobs per node")
    parser.add_argument('--time','-t', type=str, default='02:00:00', help="time limit per job (hh:mm:ss)")
    parser.add_argument('--css', action='store_true', help="if set, the run will use the constraint=cssro argument, needed for ATL11")
    parser.add_argument('--max_jobs', type=int, default=4000, help="maximum number of task files to add to the queue in one call")
    args=parser.parse_args()

    setup_directories(args.run_name)
    slurm_file='slurm_run.sh'
    N_tasks=args.jobs_per_task
    N_added=1
    part_count=0
    line_count=[-1]
    while N_added > 0:
        first_task=get_last_task(args.run_name)
        N_added = add_files_to_queue(run_name=args.run_name, task_list_file=args.queue_file, shell=args.shell, env=args.environment, max_jobs=args.max_jobs, lines_per_task=args.lines_per_task, line_count=line_count)
        if N_added <1:
            continue
        last_task=get_last_task(args.run_name)
        this_slurm_file=slurm_file
        if part_count > 0:
            this_slurm_file=slurm_file.replace('.sh',f'_part{part_count}.sh')
        ATL1415.make_slurm_file(os.path.join(args.run_name, this_slurm_file),
                            subs={'JOB_NAME':args.run_name,
                                  'TIME':args.time,
                                  'NUM_TASKS':str(N_tasks),
                                  'JOB_NUMBERS':f'{first_task+1}-{last_task}'}, css=args.css)
        part_count += 1

if __name__ == '__main__':
    main()
