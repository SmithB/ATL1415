#!/usr/bin/env python3
"""
Check a mosaic run's outputs against the run's own task files.

Exit codes are not enough to tell whether a mosaic run worked:
  - make_mosaic.py returns 0 when pc.grid.mosaic fails (it prints the reason
    only under -v, and the queue builders do not pass -v);
  - a task file has no `set -e`, so a task reports only its LAST line's status,
    and the sigma_* lines are never last.

So for every make_mosaic.py line in <run>/{queue,running,done}/task_*, this
checks that the output file exists and holds every -F field under the output
group, and reports each field's shape.  It also counts error_logs/ and any
task not yet in done/.  Exits 1 if anything is wrong.

That default reads HDF5 METADATA ONLY, so it is fast at any size.  --values
also reads every field's data to count finite values and flags a field that
is entirely NaN.  That costs a full read of every output (7 s for IS; minutes
for an Antarctic run, whose 1 km dz fields are several GB), so it is opt-in.
A failed make_mosaic.py writes nothing, which the metadata check already sees
as a missing field.

Works on any run built by make_mosaic_jobs.py or make_200km_to_mosaic_jobs.py,
on discover or MAAP: it reads only the run directory and the output files.
Run it after the run has finished -- a task still in queue/ or running/ counts
as a problem, and a file being written may not open.

Usage: check_mosaic_outputs.py <mosaic_run_dir> [--values] [-q]
"""
import argparse
import glob
import os
import shlex
import sys

import h5py
import numpy as np

TASK_STATES = ('queue', 'running', 'done')

# Rows along the first axis read at a time when counting finite values, so an
# Antarctic 1 km dz field (several GB as float64) is never held whole.
BLOCK_ROWS = 256


def parse_mosaic_line(line):
    """
    Parse one make_mosaic.py command line.

    Returns (output_path, output_group, fields), or None if the line is not a
    make_mosaic.py call.  The output path and group are resolved the way
    make_mosaic.py resolves them: -O is joined to -d, and the output group is
    --out_group if given, else --in_group (default '/').
    """
    try:
        words = shlex.split(line.strip().rstrip(';'))
    except ValueError:
        return None
    if not words or os.path.basename(words[0]) != 'make_mosaic.py':
        return None
    directory, output, in_group, out_group, fields = '', None, '/', None, []
    i = 1
    while i < len(words):
        word = words[i]
        if word in ('-d', '--directory'):
            directory = words[i + 1]; i += 2
        elif word in ('-O', '--output'):
            output = words[i + 1]; i += 2
        elif word in ('-G', '--in_group'):
            in_group = words[i + 1]; i += 2
        elif word == '--out_group':
            out_group = words[i + 1]; i += 2
        elif word in ('-F', '--fields'):
            i += 1
            while i < len(words) and not words[i].startswith('-'):
                fields.append(words[i]); i += 1
        else:
            i += 1
    if output is None:
        return None
    group = (out_group if out_group is not None else in_group).strip('/')
    return os.path.join(directory, output), group, fields


def finite_fraction(dataset):
    """Fraction of finite values, read BLOCK_ROWS rows at a time."""
    if dataset.dtype.kind not in 'fc':
        return 1.0
    if dataset.size == 0:
        return 0.0
    if dataset.ndim == 0:
        return float(np.isfinite(dataset[()]))
    n_finite = 0
    for start in range(0, dataset.shape[0], BLOCK_ROWS):
        n_finite += np.count_nonzero(np.isfinite(dataset[start:start + BLOCK_ROWS]))
    return n_finite / dataset.size


def task_number(path):
    try:
        return int(path.rsplit('_', 1)[1])
    except ValueError:
        return -1


def check_run(run_dir, values=False, quiet=False, out=sys.stdout):
    """Check one mosaic run directory.  Returns the number of problems found."""
    if not os.path.isdir(run_dir):
        print(f'NOT A DIRECTORY: {run_dir}', file=out)
        return 1
    tasks = sorted((path for state in TASK_STATES
                    for path in glob.glob(os.path.join(run_dir, state, 'task_*'))),
                   key=task_number)
    by_state = {state: [t for t in tasks if os.path.basename(os.path.dirname(t)) == state]
                for state in TASK_STATES}
    errors = sorted(glob.glob(os.path.join(run_dir, 'error_logs', '*')))

    print('tasks: ' + ', '.join(f'{state} {len(by_state[state])}' for state in TASK_STATES),
          file=out)
    print(f'error_logs: {len(errors)}', file=out)
    problems = 0
    if not tasks:
        print('NO TASK FILES found under queue/, running/ or done/', file=out)
        problems += 1
    for state in ('queue', 'running'):
        for task in by_state[state]:
            print(f'NOT DONE      {state}/{os.path.basename(task)}', file=out)
            problems += 1
    for error in errors:
        print(f'ERROR LOG     {os.path.relpath(error, run_dir)}', file=out)
        problems += 1

    n_fields = 0
    output_files = set()
    for task in tasks:
        name = os.path.basename(task)
        with open(task) as fh:
            lines = fh.readlines()
        for line in lines:
            parsed = parse_mosaic_line(line)
            if parsed is None:
                continue
            path, group, fields = parsed
            output_files.add(path)
            if not os.path.isfile(path):
                print(f'MISSING FILE  {name}  {path}', file=out)
                problems += len(fields) or 1
                continue
            try:
                h5 = h5py.File(path, 'r')
            except OSError as exc:
                print(f'CANNOT OPEN   {name}  {path}: {exc}', file=out)
                problems += len(fields) or 1
                continue
            with h5:
                for field in fields:
                    n_fields += 1
                    key = f'{group}/{field}' if group else field
                    if key not in h5:
                        print(f'MISSING FIELD {name}  {path}:{key}', file=out)
                        problems += 1
                        continue
                    summary = f'{os.path.basename(path):24s} {key:48s} {str(h5[key].shape):16s}'
                    if values:
                        frac = finite_fraction(h5[key])
                        if frac == 0:
                            print(f'ALL NAN       {name}  {path}:{key}', file=out)
                            problems += 1
                            continue
                        summary += f' finite {frac:6.1%}'
                    if not quiet:
                        print(f'ok  {summary}', file=out)

    print(f'{len(output_files)} output files, {n_fields} fields checked'
          f'{" (metadata and values)" if values else " (metadata only; --values also reads the data)"}',
          file=out)
    print(f'PROBLEMS: {problems}', file=out)
    return problems


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('run_dir', help='mosaic run directory (holds queue/, done/, error_logs/)')
    parser.add_argument('--values', action='store_true',
                        help='also read every field and flag any that is entirely NaN (slow for large runs)')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='print only problems and the summary')
    args = parser.parse_args()
    return 1 if check_run(args.run_dir, values=args.values, quiet=args.quiet) else 0


if __name__ == '__main__':
    sys.exit(main())
