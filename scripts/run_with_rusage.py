#!/usr/bin/env python3
"""
Run a command and report its peak resident set size.

Exists because DPS does not give us the number.  getJobMetrics() is documented
to return max_mem_usage and job_duration_seconds, and it DOES for a job that
failed -- but for the first tile that succeeded (fdc4d767, 2026-09-08) it came
back an empty dict, and getJob().retrieve_attributes() populated only `status`.
Sizing the production queue needs peak memory per tile, so the job measures
itself rather than waiting on the platform.

RUSAGE_CHILDREN is the whole subtree, so this reports the high-water mark of
the solve including anything it forks.  ru_maxrss is kilobytes on Linux.

Usage:  run_with_rusage.py <label> <command> [args...]
Exits with the command's own exit status, so it is transparent in a pipeline.
"""
import os
import resource
import subprocess
import sys
import time


def main(argv):
    if len(argv) < 3:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    label, cmd = argv[1], argv[2:]

    before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    t0 = time.time()
    status = subprocess.call(cmd)
    elapsed = time.time() - t0
    after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

    # ru_maxrss is a high-water mark over ALL children ever reaped by this
    # process, so it never decreases; `after` is the peak including this child,
    # and `before` is what earlier children had already reached.  Report both,
    # since only the pair distinguishes "this step was the expensive one" from
    # "an earlier step still holds the record".
    print(f'=== rusage [{label}]: elapsed {elapsed:.1f} s, '
          f'peak RSS {after/1048576:.2f} GiB '
          f'(before this step: {before/1048576:.2f} GiB), exit {status}',
          flush=True)
    return status


if __name__ == '__main__':
    sys.exit(main(sys.argv))
