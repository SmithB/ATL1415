#!/usr/bin/env python3
"""
Submit N bench jobs (run.sh step=bench) and write a ledger collect_jobs.py
reads -- docs/plan_dps_speed.sh D4.

Each job times the same saved least-squares system (scripts/maap/bench_solve.py)
at 1, 2 and 4 threads and reports its worker (WORKER: line) and the cores
each step really got (run_with_rusage.py's cpu line).  N at once, on one
queue, is the node-sharing test: collect_jobs lists which EC2 instance each
job ran on.

Usage:
  submit_bench.py <n> [--queue <q>] [--source <s3 dir>] [--ledger <f>] [--dry-run]
    n         jobs to submit, all at once
    --queue   default maap-dps-worker-16gb, the queue the IS prelim ran on
    --source  the system's directory; default run.sh's bench_default
    --ledger  default ~/ATL14_processing/maap_ledgers/bench_<queue>_<n>x_<time>_jobs.csv
Then:  collect_jobs.py <ledger>
"""
import argparse
import csv
import datetime
import os
import sys
import time

from maap.maap import MAAP

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ogc_jobs import find_process, job_id_from, load_config  # noqa: E402
from submit_MAAP_jobs import LEDGER_COLUMNS  # noqa: E402

DEFAULT_SOURCE = 's3://maap-ops-workspace/ben_smith/ATL1415/bench/E1340_N-2420_it0'
LEDGER_DIR = os.path.expanduser('~/ATL14_processing/maap_ledgers')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    parser.add_argument('n', type=int)
    parser.add_argument('--queue', default='maap-dps-worker-16gb')
    parser.add_argument('--source', default=DEFAULT_SOURCE)
    parser.add_argument('--ledger')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if args.n < 1:
        parser.error('n must be at least 1')

    stamp = time.strftime('%Y%m%dT%H%M%S', time.gmtime())
    ledger = args.ledger or os.path.join(
        LEDGER_DIR, f'bench_{args.queue.rsplit("-", 1)[-1]}_{args.n}x_{stamp}_jobs.csv')
    if os.path.exists(ledger):
        print(f'{ledger} exists; it is the only record of its jobs.  Nothing submitted.',
              file=sys.stderr)
        sys.exit(2)

    config = load_config()
    name, version = config['algorithm_name'], config['algorithm_version']
    maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
    pid = find_process(maap, name, version).get('processID')
    inputs = {'x0': '0', 'y0': '0', 'step': 'bench', 'args_file': args.source}
    print(f'process {name}:{version} processID={pid}  queue={args.queue}  n={args.n}')
    print(f'system  {args.source}\nledger  {ledger}\n')
    if args.dry_run:
        print(f'--dry-run: would submit {args.n} x submit_job({pid}, {inputs}, {args.queue!r})')
        return

    with open(ledger, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(LEDGER_COLUMNS)
        for i in range(1, args.n + 1):
            ident = f'bench_{stamp}_{i:02d}'
            # dedup=False: every job has identical inputs, and each must run
            r = maap.submit_job(pid, inputs, args.queue, dedup=False, tag=ident)
            job_id = job_id_from(r) if 200 <= r.status_code < 300 else None
            job_id = job_id or f'<submit failed: HTTP {r.status_code}: {r.text[:200]}>'
            writer.writerow([ident, 0, 0, 'bench', args.queue, args.source, job_id,
                             datetime.datetime.now(datetime.timezone.utc)
                             .isoformat(timespec='seconds'), '-'])
            fh.flush()
            print(f'  {i:3}/{args.n}  {ident}  {job_id}')
    print(f'\nthen: collect_jobs.py {ledger}')


if __name__ == '__main__':
    main()
