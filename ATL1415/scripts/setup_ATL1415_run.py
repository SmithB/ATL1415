#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set up a SLURM run for an ATL1415 queue file, splitting tasks into
per-radius-tier queues so that near-pole tiles (which cost more to fit)
get more CPUs per SLURM task than far-from-origin tiles.

@author: ben
"""

import argparse
import json
import os
import re
from importlib import resources

import numpy as np

import ATL1415
from ATL1415.scripts.setup_slurm_run import setup_directories, get_last_task, add_files_to_queue

xy0_re=re.compile(r'--xy0\s+(\S+)\s+(\S+)')
xyfile_re=re.compile(r'/E(\S*)_N(\S*).h5')

DEFAULT_RADIUS_CONFIG='radius_ncpu_tiers.json'

def load_radius_tiers(radius_config=None):
    if radius_config is None:
        radius_config=os.path.join(resources.files('ATL1415'), 'resources', DEFAULT_RADIUS_CONFIG)
    with open(radius_config, 'r') as fh:
        tiers=json.load(fh)['tiers']
    return sorted(tiers, key=lambda tier: tier['R_min'])

def split_queue_by_radius(queue_file, tiers, run_name):
    """
    Read queue_file and bucket each line into a per-tier queue file under
    run_name, based on the tile location parsed from a '--xy0 X Y'
    substring or an 'E<x>_N<y>.h5' filename substring in the line.

    Returns a dict mapping tier index (into tiers) to the tier's queue-file
    path, for tiers that received at least one line.
    """
    tier_files={}
    tier_handles={}
    try:
        with open(queue_file, 'r') as fh:
            for line in fh:
                try:
                    xy=np.array([*map(float, xy0_re.search(line).groups())])
                except Exception:
                    try:
                        xy=1000*np.array([*map(float, xyfile_re.search(line).groups())])
                    except Exception:
                        continue
                R2=np.sum(xy**2)
                for ii, tier in enumerate(tiers):
                    if (R2 >= tier['R_min']**2) and (R2 < tier['R_max']**2):
                        if ii not in tier_handles:
                            tier_file=os.path.join(run_name, f"tier_{tier['ncpu']}cpu_queue.txt")
                            tier_files[ii]=tier_file
                            tier_handles[ii]=open(tier_file, 'w')
                        tier_handles[ii].write(line)
                        break
    finally:
        for fh in tier_handles.values():
            fh.close()
    return tier_files

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_name','-r', type=str, default='ATL_run', help="name to assign to jobs and temporary directories")
    parser.add_argument('--queue_file', '-q', type=str, required=True, help="filename containing jobs, one per line")
    parser.add_argument('--environment','-e', type=str, default='ATL1415', help="environment that each job will activate")
    parser.add_argument('--shell','-s', type=str, default=None, help="shell to specify for each job (may not be needed)")
    parser.add_argument('--lines_per_task', type=int, default=1, help="combine this number of lines from the input into each task")
    parser.add_argument('--time','-t', type=str, default='02:00:00', help="time limit per job (hh:mm:ss)")
    parser.add_argument('--css', action='store_true', help="if set, the run will use the constraint=cssro argument, needed for ATL11")
    parser.add_argument('--max_jobs', type=int, default=4000, help="maximum number of task files to add to the queue in one call")
    parser.add_argument('--radius_config', type=str, default=None,
                        help="JSON file mapping radius-from-origin ranges to cpu counts; "
                             f"defaults to the packaged ATL1415/resources/{DEFAULT_RADIUS_CONFIG}")
    args=parser.parse_args()

    tiers=load_radius_tiers(args.radius_config)

    setup_directories(args.run_name)
    tier_files=split_queue_by_radius(args.queue_file, tiers, args.run_name)

    for ii, tier in enumerate(tiers):
        if ii not in tier_files:
            continue
        N_tasks=tier['ncpu']
        slurm_file=f"slurm_run_queue_{N_tasks}cpu.sh"
        N_added=1
        part_count=0
        line_count=[-1]
        while N_added > 0:
            first_task=get_last_task(args.run_name)
            N_added = add_files_to_queue(run_name=args.run_name, task_list_file=tier_files[ii], shell=args.shell, max_jobs=args.max_jobs, lines_per_task=args.lines_per_task, line_count=line_count)
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
                                      'ENVIRONMENT':args.environment,
                                      'JOB_NUMBERS':f'{first_task+1}-{last_task}'}, css=args.css)
            part_count += 1

if __name__ == '__main__':
    main()
