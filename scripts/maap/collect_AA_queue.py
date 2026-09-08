#!/usr/bin/env python3
"""Status + time + memory + input size for every job in a transect ledger."""
import csv, re, sys, subprocess
from maap.maap import MAAP

LEDGER = sys.argv[1] if len(sys.argv) > 1 else 'AA_transect_jobs.csv'
m = MAAP(maap_host='api.maap-project.org')

def stdout_for(job_id):
    """Fetch a job's _stdout.txt from whichever prefix DPS put it in."""
    try:
        res = m.getJobResult(job_id)
    except Exception:
        return ''
    for item in (res if isinstance(res, list) else [res]):
        if isinstance(item, str) and item.startswith('s3://'):
            base = item.rstrip('/')
            for name in ('_stdout.txt',):
                p = subprocess.run(['aws', 's3', 'cp', f'{base}/{name}', '-'],
                                   capture_output=True, text=True, timeout=120)
                if p.returncode == 0:
                    return p.stdout
    return ''

N_RE   = re.compile(r'decimate_data: N_target:[^,]+, N=(\d+)')
RUSAGE = re.compile(r'=== rusage \[(\w+)\]: elapsed ([\d.]+) s, peak RSS ([\d.]+) GiB')
FIT_RE = re.compile(r'initial: (\d+):')
ITER_RE= re.compile(r'starting qr solve for iteration (\d+)')

print(f"{'tile':>22} {'status':>10} {'secs':>7} {'max_mem_GiB':>12} "
      f"{'N_ATL11':>9} {'N_fit':>8} {'iters':>5} {'machine':>14}")
for row in csv.DictReader(open(LEDGER)):
    jid = row['job_id']
    if not jid or jid.startswith('<'):
        print(f"{row['identifier'][-22:]:>22} {'NOT SUBMITTED':>10}"); continue
    # getJob returns a plain dict for a FAILED job but a DPSJob for a succeeded
    # one, and the DPSJob's fields are empty until retrieve_attributes() runs.
    # getJobMetrics has been observed to return {} for a succeeded job, which
    # is why run.sh reports peak RSS itself -- see scripts/run_with_rusage.py.
    d = {}
    try:
        info = m.getJob(jid)
        if isinstance(info, dict):
            d = dict(info)
        else:
            info.id = jid
            try:
                info.retrieve_attributes()
            except Exception:
                pass
            d = {k: getattr(info, k, None) for k in
                 ('status', 'machine_type', 'machine_memory_size',
                  'job_duration_seconds', 'max_mem_usage', 'outputs')}
    except Exception as exc:
        print(f"{row['identifier'][-22:]:>22} <getJob failed: {exc}>"); continue
    try:
        mt = m.getJobMetrics(jid)
        if isinstance(mt, dict):
            d.update({k: v for k, v in mt.items() if v is not None})
    except Exception:
        pass
    out = stdout_for(jid)
    n   = N_RE.search(out); fit = FIT_RE.search(out)
    iters = ITER_RE.findall(out)
    # Prefer what the job measured about itself over what DPS reports.
    steps = {k: (float(t), float(g)) for k, t, g in RUSAGE.findall(out)}
    mem = d.get('max_mem_usage')
    if steps:
        mem = max(g for _, g in steps.values()) * 2**30
    try:    mem_s = f'{float(mem)/2**30:12.2f}'
    except (TypeError, ValueError): mem_s = f'{str(mem):>12}'
    dur = d.get('job_duration_seconds')
    print(f"{row['identifier'][-22:]:>22} {str(d.get('status')):>10} "
          f"{str(dur):>7} {mem_s} {(n.group(1) if n else '-'):>9} "
          f"{(fit.group(1) if fit else '-'):>8} "
          f"{(max(map(int,iters))+1 if iters else '-'):>5} "
          f"{str(d.get('machine_type')):>14}")
    for label, (secs, gib) in steps.items():
        print(f"{'':22} {'step ' + label:>10} {secs:7.0f} {gib:12.2f}")
