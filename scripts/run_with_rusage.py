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

AND WHERE THE CPU WENT (docs/plan_dps_speed.sh D1).  A second line per step:
  cpu_time   the subtree's user+system CPU seconds, and that over the wall
             time = the cores the step actually got on average.  4 threads
             that get 1.5 cores say the machine is short of CPU, whatever
             nproc claims.
  load1      /proc/loadavg sampled every 10 s while the step runs: min,
             median, max.  Host-wide in a container, so it counts other
             jobs' threads on the same node.
  steal      /proc/stat steal time over the step, summed over all CPUs:
             the hypervisor giving our vCPUs to someone else.
  throttled  cgroup cpu.stat throttled time and count over the step: a CPU
             quota below what the threads ask for.
Any of them is "na" where the file does not exist; none can fail the step.

Usage:  run_with_rusage.py <label> <command> [args...]
Exits with the command's own exit status, so it is transparent in a pipeline.
"""
import os
import resource
import statistics
import subprocess
import sys
import threading
import time

LOAD_SAMPLE_S = 10


def read(path):
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return None


def steal_s():
    """Host steal time so far, CPU-seconds summed over CPUs, or None."""
    text = read('/proc/stat')
    if not text or not text.startswith('cpu '):
        return None
    fields = text.split('\n', 1)[0].split()
    if len(fields) < 9:
        return None
    return int(fields[8]) / os.sysconf('SC_CLK_TCK')


def throttled():
    """(throttled seconds, times throttled) for this cgroup so far, or None."""
    text = read('/sys/fs/cgroup/cpu.stat')          # cgroup v2
    if text:
        stat = dict(line.split() for line in text.splitlines() if len(line.split()) == 2)
        if 'throttled_usec' in stat:
            return int(stat['throttled_usec']) / 1e6, int(stat.get('nr_throttled', 0))
    text = read('/sys/fs/cgroup/cpu/cpu.stat')      # cgroup v1
    if text:
        stat = dict(line.split() for line in text.splitlines() if len(line.split()) == 2)
        if 'throttled_time' in stat:
            return int(stat['throttled_time']) / 1e9, int(stat.get('nr_throttled', 0))
    return None


def load1():
    text = read('/proc/loadavg')
    return float(text.split()[0]) if text else None


class LoadSampler(threading.Thread):
    """Samples the 1-minute load average until stopped."""
    def __init__(self):
        super().__init__(daemon=True)
        self.samples, self.done = [], threading.Event()

    def run(self):
        while True:
            value = load1()
            if value is not None:
                self.samples.append(value)
            if self.done.wait(LOAD_SAMPLE_S):
                return


def cpu_line(label, elapsed, cpu, loads, steal, thr):
    cores = f'{cpu / elapsed:.2f}' if elapsed > 0 else 'na'
    if loads:
        load = (f'{min(loads):.1f}/{statistics.median(loads):.1f}/{max(loads):.1f}'
                f' (n={len(loads)})')
    else:
        load = 'na'
    steal = f'{steal:.1f} s' if steal is not None else 'na'
    thr = f'{thr[0]:.1f} s (n={thr[1]})' if thr is not None else 'na'
    return (f'=== cpu [{label}]: cpu_time {cpu:.1f} s = {cores} cores avg, '
            f'load1 min/med/max {load}, steal {steal}, throttled {thr}')


def main(argv):
    if len(argv) < 3:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    label, cmd = argv[1], argv[2:]

    r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    before = r0.ru_maxrss
    steal0, thr0 = steal_s(), throttled()
    sampler = LoadSampler()
    sampler.start()
    t0 = time.time()
    status = subprocess.call(cmd)
    elapsed = time.time() - t0
    sampler.done.set()
    sampler.join()
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    after = r1.ru_maxrss
    cpu = (r1.ru_utime + r1.ru_stime) - (r0.ru_utime + r0.ru_stime)
    steal1, thr1 = steal_s(), throttled()
    steal = steal1 - steal0 if None not in (steal0, steal1) else None
    thr = ((thr1[0] - thr0[0], thr1[1] - thr0[1])
           if None not in (thr0, thr1) else None)

    # ru_maxrss is a high-water mark over ALL children ever reaped by this
    # process, so it never decreases; `after` is the peak including this child,
    # and `before` is what earlier children had already reached.  Report both,
    # since only the pair distinguishes "this step was the expensive one" from
    # "an earlier step still holds the record".
    print(f'=== rusage [{label}]: elapsed {elapsed:.1f} s, '
          f'peak RSS {after/1048576:.2f} GiB '
          f'(before this step: {before/1048576:.2f} GiB), exit {status}',
          flush=True)
    print(cpu_line(label, elapsed, cpu, sampler.samples, steal, thr), flush=True)
    return status


if __name__ == '__main__':
    sys.exit(main(sys.argv))
