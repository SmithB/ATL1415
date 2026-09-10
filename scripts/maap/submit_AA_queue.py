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

OGC SYSTEM, since 2026-09-10 (docs/howto_MAAP_ogc.sh O7).  maap-py 5.x
dropped submitJob; this finds the deployed process by name and version and
submits with submit_job(process_id, inputs, queue, dedup=False, tag=...).
The inputs are the CWL's -- strings, bound onto run.sh as --x0/--y0/--step/
--args_file -- and the ledger is unchanged, so collect_AA_queue.py reads old
and new ledgers alike.

dedup=False, EXPLICITLY, as in check_build_id.py: a tile resubmitted after a
rebuild has identical inputs, and a deduplicated job would hand back the OLD
image's result -- exactly what a post-rebuild rerun is for escaping.

Usage:
  submit_AA_queue.py [xy_file] [args_60km_url] [args_44km_url] [queue] [ledger]
                     [--dry-run]

--dry-run finds the process and prints every job it would submit, then stops.
"""
import csv
import datetime
import os
import sys
import time

from maap.maap import MAAP

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ogc_jobs import S3_RUN, find_process, job_id_from, load_config  # noqa: E402

# The name AND version come from algorithm_config.yml, never hard-coded: the
# version is the container tag and the branch the build clones, and a
# submitter naming an old one would silently run an old image.
CONFIG = load_config()
NAME, VERSION = CONFIG['algorithm_name'], CONFIG['algorithm_version']

ARGV = [a for a in sys.argv[1:] if a != '--dry-run']
DRY_RUN = '--dry-run' in sys.argv[1:]
XY        = ARGV[0] if len(ARGV) > 0 else 'scripts/maap/AA_queue_xy.txt'
ARGS_60   = ARGV[1] if len(ARGV) > 1 else f'{S3_RUN}/input_args_AA.txt'
ARGS_44   = ARGV[2] if len(ARGV) > 2 else f'{S3_RUN}/input_args_AA_44km.txt'
QUEUE     = ARGV[3] if len(ARGV) > 3 else 'maap-dps-worker-32gb'
LEDGER    = ARGV[4] if len(ARGV) > 4 else 'AA_queue_jobs.csv'

MIN_XY, MAX_XY = 360000, 440000        # the two halves' limits, as in the howto


def check_args():
    """
    Refuse arguments that are plainly in the wrong slot, before anything runs.

    Five positionals are easy to transpose, and it has happened: the howto's
    own commands (howto_MAAP_AA 3b and 3b-i, until 2026-09-10) passed four,
    which put the QUEUE NAME in the 44 km args-file slot and the LEDGER NAME
    in the queue slot.  Nothing would have noticed until MAAP rejected the
    queue -- or worse, accepted a job whose args_file was a queue name.
    """
    problems = []
    for label, value in (('args_60km', ARGS_60), ('args_44km', ARGS_44)):
        if not (value.startswith('s3://') or os.path.isfile(value)):
            problems.append(f'{label}={value!r} is neither an s3:// URI nor a file')
    if not QUEUE.startswith('maap-dps-'):
        problems.append(f'queue={QUEUE!r} does not look like a DPS queue (maap-dps-*)')
    if not LEDGER.endswith('.csv'):
        problems.append(f'ledger={LEDGER!r} is not a .csv file name')
    if problems:
        print('ARGUMENTS OUT OF PLACE -- usage: submit_AA_queue.py [xy_file]'
              ' [args_60km_url] [args_44km_url] [queue] [ledger]', file=sys.stderr)
        for problem in problems:
            print(f'  - {problem}', file=sys.stderr)
        sys.exit(2)


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
    check_args()
    maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
    process = find_process(maap, NAME, VERSION)
    pid = process.get('processID')
    print(f'process {NAME}:{VERSION}  processID={pid}'
          f"  (modified {process.get('lastModifiedTime')})  queue={QUEUE}\n")
    centers = [tuple(int(float(v)) for v in line.split())
               for line in open(XY) if line.strip()]

    if DRY_RUN:
        for x0, y0 in centers:
            for half, args_url in halves_for(x0, y0):
                print(f'would submit AA_cost_{half}_E{x0//1000}_N{y0//1000:<6}'
                      f' x0={x0} y0={y0} args={args_url.rsplit("/", 1)[-1]}')
        print('\n--dry-run: nothing submitted.')
        return

    with open(LEDGER, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['identifier', 'x0', 'y0', 'half', 'queue',
                         'args_file', 'job_id', 'submitted_utc'])
        n = 0
        for x0, y0 in centers:
            for half, args_url in halves_for(x0, y0):
                ident = f'AA_cost_{half}_E{x0//1000}_N{y0//1000}'
                # Strings: every input is `type: string` in the CWL.
                inputs = {'x0': str(x0), 'y0': str(y0), 'step': 'prelim',
                          'args_file': args_url}
                try:
                    r = maap.submit_job(pid, inputs, QUEUE, dedup=False, tag=ident)
                    job_id = job_id_from(r) if 200 <= r.status_code < 300 else None
                    if not job_id:
                        job_id = f'<submit failed: HTTP {r.status_code}: {r.text[:200]}>'
                except Exception as exc:
                    job_id = f'<submit failed: {type(exc).__name__}: {exc}>'
                stamp = datetime.datetime.now(
                    datetime.timezone.utc).isoformat(timespec='seconds')
                writer.writerow([ident, x0, y0, half, QUEUE, args_url,
                                 job_id, stamp])
                fh.flush()
                print(f'{ident:34} {job_id}')
                n += 1
                time.sleep(2)     # gentle with the queue's rate limit
    print(f'\n{n} jobs -> {LEDGER}')


if __name__ == '__main__':
    main()
