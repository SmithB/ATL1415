#! /usr/bin/env python

"""
Summarize the status of a slurm task-queue run (queue/running/done/logs/
active_logs/error_logs, as created by setup_slurm_run.py, make_mosaic_jobs.py,
or make_200km_to_mosaic_jobs.py), and optionally requeue failed or killed
tasks.

A task that completed with a nonzero exit status has its log copied to
error_logs/ by the slurm job template. A task that was killed before it
could finish (slurm time/memory limit, node failure, scancel) never makes
it past the active_logs/ stage: its log is never moved to logs/, and its
task file is never moved out of running/. Those two states are
distinguished here directly from where the files ended up, rather than by
grepping log text for words like "killed" or "error".
"""

import argparse
import glob
import json
import os
import re
import shutil

tile_re = re.compile(r'working on .*/E(-?\d+)_N(-?\d+)\.h5')
trace_line_re = re.compile(r'(File "(.*)", line (\d+), in (\S+))')
error_re = re.compile(r'^((\S+Error): (.*))')
annon_error_re = re.compile(r'^((\S+)Error)')


def task_num_from_path(path):
    m = re.search(r'task_(\d+)', os.path.basename(path))
    return m.group(1) if m else None


def parse_error_log(log_file):
    """Return a dict describing the last traceback/exception found in a log file, or None."""
    trace_line_info = None
    this_xy = None
    with open(log_file, 'r', errors='replace') as fh:
        for line in fh:
            m = tile_re.search(line)
            if m is not None:
                this_xy = [*map(int, m.groups())]
            m = trace_line_re.search(line)
            if m is not None:
                trace_line_info = {'str': m.group(1), 'file': m.group(2),
                                    'line': m.group(3), 'function': m.group(4)}
            m = annon_error_re.search(line)
            if m is not None and trace_line_info is not None:
                trace_line_info['exception'] = m.group(1)
                m = error_re.search(line)
                if m is not None:
                    trace_line_info['exception'] = m.group(1)
                    break
    if trace_line_info is None:
        return None
    if 'exception' not in trace_line_info:
        trace_line_info['exception'] = 'unrecognized (traceback present, no Error line found)'
    trace_line_info['xy'] = this_xy
    return trace_line_info


def summarize(run_dir):
    counts = {
        'queued': len(glob.glob(os.path.join(run_dir, 'queue', 'task_*'))),
        'running': len(glob.glob(os.path.join(run_dir, 'running', 'task_*'))),
        'done': len(glob.glob(os.path.join(run_dir, 'done', 'task_*'))),
    }

    error_logs = sorted(glob.glob(os.path.join(run_dir, 'error_logs', 'task_*.log')))
    killed_logs = sorted(glob.glob(os.path.join(run_dir, 'active_logs', 'task_*.log')))

    errors = []
    by_exception = {}
    for log_file in error_logs:
        task_num = task_num_from_path(log_file)
        info = parse_error_log(log_file)
        if info is None:
            info = {'exception': 'unrecognized (nonzero exit, no traceback found)'}
        entry = {'task': task_num, 'log': log_file, **info}
        errors.append(entry)
        by_exception.setdefault(entry['exception'], []).append(task_num)

    killed = [{'task': task_num_from_path(f), 'log': f} for f in killed_logs]

    return {
        'counts': counts,
        'errors': errors,
        'killed': killed,
        'by_exception': {key: len(val) for key, val in by_exception.items()},
    }


def print_summary(summary):
    c = summary['counts']
    print(f"queued: {c['queued']}   running: {c['running']}   done: {c['done']}")
    print(f"errors (completed, nonzero exit): {len(summary['errors'])}")
    print(f"killed (started, never finished): {len(summary['killed'])}")
    if summary['by_exception']:
        print("\nerrors by exception type:")
        for exc, n in sorted(summary['by_exception'].items(), key=lambda kv: -kv[1]):
            print(f"\t{n}\t{exc}")
    if summary['killed']:
        print("\nkilled tasks:")
        for k in summary['killed']:
            print(f"\ttask_{k['task']}")


def requeue_errors(run_dir, summary, dry_run=False):
    archive_dir = os.path.join(run_dir, 'error_logs', 'requeued')
    if not dry_run:
        os.makedirs(archive_dir, exist_ok=True)
    requeued = []
    for entry in summary['errors']:
        task_num = entry['task']
        done_file = os.path.join(run_dir, 'done', f'task_{task_num}')
        queue_file = os.path.join(run_dir, 'queue', f'task_{task_num}')
        if not os.path.isfile(done_file):
            print(f"task_{task_num}: expected done file {done_file} not found, skipping")
            continue
        print(f"requeuing task_{task_num} (error: {entry.get('exception', 'unknown')})")
        requeued.append(task_num)
        if dry_run:
            continue
        shutil.move(done_file, queue_file)
        shutil.move(entry['log'], os.path.join(archive_dir, os.path.basename(entry['log'])))
    return requeued


def requeue_killed(run_dir, summary, dry_run=False):
    archive_dir = os.path.join(run_dir, 'active_logs', 'requeued')
    if not dry_run:
        os.makedirs(archive_dir, exist_ok=True)
    requeued = []
    for entry in summary['killed']:
        task_num = entry['task']
        running_file = os.path.join(run_dir, 'running', f'task_{task_num}')
        queue_file = os.path.join(run_dir, 'queue', f'task_{task_num}')
        if not os.path.isfile(running_file):
            print(f"task_{task_num}: expected running file {running_file} not found, skipping")
            continue
        print(f"requeuing killed task_{task_num}")
        requeued.append(task_num)
        if dry_run:
            continue
        shutil.move(running_file, queue_file)
        shutil.move(entry['log'], os.path.join(archive_dir, os.path.basename(entry['log'])))
    return requeued


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('run_dir', nargs='?', default='.',
                         help="run directory containing queue/running/done/logs/active_logs/error_logs (default: current directory)")
    parser.add_argument('--requeue', choices=['errors', 'killed', 'both'],
                         help="move failed ('errors') and/or stalled ('killed') tasks back into queue/")
    parser.add_argument('--dry_run', action='store_true',
                         help="report what would be requeued without moving any files")
    parser.add_argument('--json', type=str, help="also write the summary to this file as json")
    args = parser.parse_args()

    summary = summarize(args.run_dir)
    print_summary(summary)

    if args.requeue in ('errors', 'both'):
        print()
        requeue_errors(args.run_dir, summary, dry_run=args.dry_run)
    if args.requeue in ('killed', 'both'):
        print()
        requeue_killed(args.run_dir, summary, dry_run=args.dry_run)

    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(summary, fh, indent=2)


if __name__ == '__main__':
    main()
