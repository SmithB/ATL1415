#!/usr/bin/env python3
"""Submit the AA transect, one job per tile center, and write a ledger."""
import csv, sys, time, datetime
from maap.maap import MAAP

XY     = sys.argv[1] if len(sys.argv) > 1 else 'AA_transect_xy.txt'
ARGS   = sys.argv[2] if len(sys.argv) > 2 else (
    's3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/south/AA/input_args_AA.txt')
QUEUE  = sys.argv[3] if len(sys.argv) > 3 else 'maap-dps-worker-32gb'
LEDGER = sys.argv[4] if len(sys.argv) > 4 else 'AA_transect_jobs.csv'

m = MAAP(maap_host='api.maap-project.org')
centers = [tuple(int(float(v)) for v in ln.split())
           for ln in open(XY) if ln.strip()]

with open(LEDGER, 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['identifier', 'x0', 'y0', 'queue', 'job_id', 'submitted_utc'])
    for x0, y0 in centers:
        ident = f'AA_transect_E{x0//1000}_N{y0//1000}'
        try:
            job = m.submitJob(identifier=ident,
                              algo_id='ATL1415_tile_solve', version='on_s3',
                              queue=QUEUE, queue_name=QUEUE,
                              x0=x0, y0=y0, step='prelim', args_file=ARGS)
            jid = getattr(job, 'id', None) or getattr(job, 'job_id', None)
        except Exception as exc:
            jid = f'<submit failed: {type(exc).__name__}: {exc}>'
        stamp = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')
        w.writerow([ident, x0, y0, QUEUE, jid, stamp]); fh.flush()
        print(f'{ident:28} {jid}')
        time.sleep(2)          # be gentle with the public queue's rate limit
print(f'\nledger: {LEDGER}')
