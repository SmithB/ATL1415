#!/usr/bin/env python3
"""
Print one WORKER: line describing the machine this job runs on.

Exists because DPS does not say: get_job_metrics() returns machine_type,
architecture and machine_memory_size as null (read 2026-09-24), and the SPQR
solves run 2-5x slower on DPS than on the ADE for reasons the logs cannot
show (docs/plan_dps_speed.sh).  Every field is a fact about THIS machine at
job start:

  instance_type, instance_id, az,  EC2 instance metadata (IMDSv2, then v1).
  lifecycle                        spot or on-demand (instance-life-cycle):
                                   Ben expects DPS workers to be spot, which
                                   would explain why a queue's instance type
                                   varies from job to job.
                                   instance_id says which jobs shared a node.
                                   "unreachable" when the container cannot
                                   reach it (an IMDSv2 hop limit of 1 blocks
                                   a container behind a bridge network).
  cpu                              the CPU model, spaces as underscores
  vcpus                            logical CPUs the kernel shows
  affinity                         CPUs this process may run on (cpuset)
  threads_per_core                 2 = hyperthreads: vcpus/2 real cores
  cpu_quota                        cgroup cpu.max as cores ("max" = none)
  mem_gib, mem_limit_gib           MemTotal; cgroup memory.max
  blas                             library-version-architecture-threading,
                                   num_threads as OpenBLAS would start
  load1                            /proc/loadavg; host-wide in a container,
                                   so it counts other jobs' threads too

Every read is guarded: a missing fact prints "na" and never fails the job.
Usage:  worker_facts.py        (prints one line, exits 0)
"""
import os
import re
import urllib.request

IMDS = 'http://169.254.169.254/latest'


def imds(keys, timeout=2):
    """{key: value} from the instance metadata service, or {} if unreachable."""
    headers = {}
    try:
        req = urllib.request.Request(f'{IMDS}/api/token', method='PUT',
                                     headers={'X-aws-ec2-metadata-token-ttl-seconds': '60'})
        headers['X-aws-ec2-metadata-token'] = urllib.request.urlopen(req, timeout=timeout).read().decode()
    except Exception:
        pass    # IMDSv1 may still answer without a token
    out = {}
    for key in keys:
        try:
            req = urllib.request.Request(f'{IMDS}/meta-data/{key}', headers=headers)
            out[key] = urllib.request.urlopen(req, timeout=timeout).read().decode().strip()
        except Exception:
            return {}
    return out


def read(path):
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return None


def cpu_model():
    text = read('/proc/cpuinfo') or ''
    m = re.search(r'^model name\s*:\s*(.+)$', text, re.M)
    return m.group(1).strip() if m else 'na'


def threads_per_core():
    text = read('/sys/devices/system/cpu/cpu0/topology/thread_siblings_list')
    if not text:
        return 'na'
    n = 0
    for part in text.strip().split(','):
        lo, _, hi = part.partition('-')
        n += int(hi or lo) - int(lo) + 1
    return str(n)


def cpu_quota():
    text = read('/sys/fs/cgroup/cpu.max')           # cgroup v2: "<quota> <period>"
    if text:
        quota, period = text.split()[:2]
        return 'max' if quota == 'max' else f'{int(quota) / int(period):.2f}'
    quota = read('/sys/fs/cgroup/cpu/cpu.cfs_quota_us')   # cgroup v1
    period = read('/sys/fs/cgroup/cpu/cpu.cfs_period_us')
    if quota and period:
        return 'max' if int(quota) < 0 else f'{int(quota) / int(period):.2f}'
    return 'na'


def mem_gib():
    m = re.search(r'^MemTotal:\s*(\d+) kB', read('/proc/meminfo') or '', re.M)
    return f'{int(m.group(1)) / 1048576:.1f}' if m else 'na'


def mem_limit_gib():
    text = (read('/sys/fs/cgroup/memory.max')
            or read('/sys/fs/cgroup/memory/memory.limit_in_bytes') or '').strip()
    if not text:
        return 'na'
    if text == 'max' or int(text) > 2**60:
        return 'max'
    return f'{int(text) / 2**30:.1f}'


def blas():
    try:
        import numpy  # noqa: F401  -- loads the BLAS threadpoolctl inspects
        from threadpoolctl import threadpool_info
        for info in threadpool_info():
            if info.get('user_api') == 'blas':
                return (f"{info.get('internal_api')}-{info.get('version')}-"
                        f"{info.get('architecture')}-{info.get('threading_layer')}"
                        f"-n{info.get('num_threads')}")
        return 'none'
    except Exception as exc:
        return f'na({type(exc).__name__})'


def main():
    meta = imds(['instance-type', 'instance-id', 'placement/availability-zone',
                 'instance-life-cycle'])
    try:
        affinity = str(len(os.sched_getaffinity(0)))
    except Exception:
        affinity = 'na'
    load = (read('/proc/loadavg') or 'na').split()[0]
    fields = {
        'instance_type': meta.get('instance-type', 'unreachable'),
        'instance_id': meta.get('instance-id', 'unreachable'),
        'az': meta.get('placement/availability-zone', 'unreachable'),
        'lifecycle': meta.get('instance-life-cycle', 'unreachable'),
        'cpu': cpu_model().replace(' ', '_'),
        'vcpus': str(os.cpu_count()),
        'affinity': affinity,
        'threads_per_core': threads_per_core(),
        'cpu_quota': cpu_quota(),
        'mem_gib': mem_gib(),
        'mem_limit_gib': mem_limit_gib(),
        'blas': blas(),
        'load1': load,
    }
    print('WORKER: ' + ' '.join(f'{k}={v}' for k, v in fields.items()), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
