#!/usr/bin/env python3
"""
Submit the Antarctic cost-characterisation queue, one job per tile, and write a
ledger.

ROUTES EACH TILE TO THE RIGHT HALF.  Antarctica is solved as two halves with
DIFFERENT TILE WIDTHS -- 60 km north of the 400 km line, 44 km south of it --
selected by make_ATL1415_queue.py with --min_xy 360000 and --max_xy 440000.
Those limits deliberately overlap (Ben, 2026-09-08), so a tile whose max|xy|
falls in 360-440 km belongs to BOTH halves and is submitted TWICE, once per
geometry.  For a cost experiment that is a feature: the same center solved at
two widths is a direct read on how width alone drives time and memory.

Getting this wrong is silent -- the solve would simply run at the wrong width
and produce a plausible tile -- which is why the routing lives here rather than
in the caller's head.

Usage:
  submit_AA_queue.py [xy_file] [args_60km_url] [args_44km_url] [queue] [ledger]
"""
import csv
import datetime
import sys
import time

from maap.maap import MAAP

S3_RUN = ('s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/south/AA')

XY        = sys.argv[1] if len(sys.argv) > 1 else 'scripts/maap/AA_queue_xy.txt'
ARGS_60   = sys.argv[2] if len(sys.argv) > 2 else f'{S3_RUN}/input_args_AA.txt'
ARGS_44   = sys.argv[3] if len(sys.argv) > 3 else f'{S3_RUN}/input_args_AA_44km.txt'
QUEUE     = sys.argv[4] if len(sys.argv) > 4 else 'maap-dps-worker-32gb'
LEDGER    = sys.argv[5] if len(sys.argv) > 5 else 'AA_queue_jobs.csv'

MIN_XY, MAX_XY = 360000, 440000        # the two halves' limits, as in the howto


def halves_for(x0, y0):
    """Which half (or halves) claims this tile center."""
    extent = max(abs(x0), abs(y0))
    out = []
    if extent <= MAX_XY:
        out.append(('44km', ARGS_44))
    if extent >= MIN_XY:
        out.append(('60km', ARGS_60))
    return out


def main():
    maap = MAAP(maap_host='api.maap-project.org')
    centers = [tuple(int(float(v)) for v in line.split())
               for line in open(XY) if line.strip()]

    with open(LEDGER, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['identifier', 'x0', 'y0', 'half', 'queue',
                         'args_file', 'job_id', 'submitted_utc'])
        n = 0
        for x0, y0 in centers:
            for half, args_url in halves_for(x0, y0):
                ident = f'AA_cost_{half}_E{x0//1000}_N{y0//1000}'
                try:
                    job = maap.submitJob(
                        identifier=ident,
                        algo_id='ATL1415_tile_solve', version='on_s3',
                        queue=QUEUE, queue_name=QUEUE,
                        x0=x0, y0=y0, step='prelim', args_file=args_url)
                    job_id = getattr(job, 'id', None) or getattr(job, 'job_id', None)
                except Exception as exc:
                    job_id = f'<submit failed: {type(exc).__name__}: {exc}>'
                stamp = datetime.datetime.now(
                    datetime.timezone.utc).isoformat(timespec='seconds')
                writer.writerow([ident, x0, y0, half, QUEUE, args_url,
                                 job_id, stamp])
                fh.flush()
                print(f'{ident:34} {job_id}')
                n += 1
                time.sleep(2)     # gentle with the public queue's rate limit
    print(f'\n{n} jobs -> {LEDGER}')


if __name__ == '__main__':
    main()
