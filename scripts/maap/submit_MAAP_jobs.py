#!/usr/bin/env python3
"""
Submit one DPS job per tile center for ONE region, and write the ledger.

THE GENERAL SUBMITTER.  scripts/maap/submit_AA_queue.py stays as it is: it
carries Antarctica's two-width routing (60 km north of the 400 km line, 44 km
south of it, deliberately overlapping so a tile in the band is submitted
twice), which is a real feature of AA and of nothing else.  Every other region
has one geometry and one args file, and this is for those.
Decided 2026-09-12, docs/plan_IS_run.sh QI2.

WHAT IT SHARES WITH THE AA SUBMITTER, deliberately and through ogc_jobs.py:
the process lookup by name+version (never a hard-coded processID -- it changes
on every redeploy), submit_job(pid, inputs, queue, dedup=False, tag=...), and
the ledger columns -- plus ONE MORE here, tile_prefix, appended last.
collect_jobs.py and fetch_tiles.py read columns by name (csv.DictReader), so
they read the ledger this writes, and the AA submitter's, and ledgers written
before the column existed, without knowing which is which.

THE tile_prefix COLUMN (added 2026-09-15): QI5a made tile_prefix a job input
so that where a tile went is written down per job rather than implied by a
convention -- and the ledger is where it is written down.  "-" means the job
was sent none, which is the registered default and run.sh's "none".

dedup=False, EXPLICITLY, as everywhere else: a tile resubmitted after a
rebuild has byte-identical inputs, and a deduplicated job would hand back the
OLD image's result -- exactly what a post-rebuild rerun exists to escape.

THE LEDGER IS THE ONLY RECORD OF WHAT WAS SUBMITTED.  A job id that is not
written down is a worker-hour that cannot be collected, so:
  - it is written incrementally and flushed per row, and an interrupted run
    still leaves a readable ledger;
  - an existing ledger is NEVER overwritten without --replace;
  - a submit that fails is recorded, with the error in the job_id field, and
    the run continues (Q11: record-and-continue, not retry).

Usage:
  submit_MAAP_jobs.py (--tile_list <f> | --xy_file <f>)
                      --step prelim|matched --args_url <uri|path>
                      [--ledger <f>] [--queue <q>] [--tag <prefix>]
                      [--tile_prefix <s3://...>]
                      [--limit N] [--rate S] [--max_in_flight N]
                      [--replace] [--dry-run]

  --tile_list      THE FAN-OUT INPUT (docs/plan_tile_lists.sh TL2, Ben's AM8):
                   one tile file name per line, E<x km>_N<y km>.h5, e.g.
                   ATL1415/resources/IS/40km_tile_list.txt.  Pruned of
                   no-data centers by scripts/maap/prune_tile_list.py.
  --xy_file        one "<x0> <y0>" per line, in meters -- for the one-center
                   smoke and retry files.  The region_files/*_prelim_xy.txt
                   center lists are retired for fan-outs.
  --args_url       the composed input_args_<REGION>.txt, on the bucket
  --tag            identifier prefix; default "<REGION>_<step>", with REGION
                   read out of the args file's name
  --tile_prefix    where the job writes its tile / reads its neighbours.
                   Deployed since the 6978a8a registration; main() still
                   refuses it if the config or deployed CWL lacks it.
  --limit N        submit only the first N centers.  USE IT: one job on a
                   queue nobody has used before costs one job to find out,
                   and 29 to find out the expensive way.
  --rate S         seconds between submissions (default 2)
  --max_in_flight  hold at N un-finished jobs, polling until one finishes

MATCHED PRE-FLIGHT.  With --step matched, every center's own prelim tile must
already exist at <tile_prefix>/prelim/.  If any is missing, nothing is
submitted and the missing names are printed: each is either a no-data center
not yet pruned from the list, or a failed prelim under investigation, and a
matched job for it could only fail on DPS.

Written for the IS run (docs/plan_IS_run.sh I2) and intended for GL next.
"""
import argparse
import csv
import datetime
import os
import re
import subprocess
import sys
import time

import requests
from maap.maap import MAAP

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ogc_jobs import (DEFAULT_QUEUE, DONE, POLL_S, find_process,  # noqa: E402
                      job_id_from, load_config)

LEDGER_COLUMNS = ['identifier', 'x0', 'y0', 'step', 'queue', 'args_file',
                  'job_id', 'submitted_utc', 'tile_prefix']


def region_of(args_url):
    """IS, from .../input_args_IS.txt -- else None, and --tag is required."""
    m = re.search(r'input_args_(.+?)\.txt$', os.path.basename(args_url))
    return m.group(1) if m else None


def read_centers(path):
    """The xy file: "<x0> <y0>" per line, meters.

    NO COMMENT SYNTAX.  Every non-blank line must parse, and a line that does
    not is an error rather than a skip: silently dropping one would submit a
    region short by a tile and nothing downstream would notice.
    """
    centers = []
    for n, line in enumerate(open(path), 1):
        if not line.strip():
            continue
        try:
            x0, y0 = (int(float(v)) for v in line.split())
        except ValueError:
            print(f'{path}:{n}: cannot read a tile center from {line.strip()!r}',
                  file=sys.stderr)
            sys.exit(2)
        centers.append((x0, y0))
    return centers


TILE_NAME = re.compile(r'^E(-?\d+)_N(-?\d+)\.h5$')


def tile_name(x0, y0):
    """'E%d_N%d.h5', truncated toward zero -- ATL11_to_ATL15's own name."""
    return 'E%d_N%d.h5' % (int(x0 / 1000), int(y0 / 1000))


def read_tile_list(path):
    """The resource tile list: "E<x km>_N<y km>.h5" per line -> meters.

    The same rule as read_centers: every non-blank line must be a tile name,
    and one that is not is an error rather than a skip -- a list that has
    picked up something else (a directory name from an `ls`, say) is not a
    list to submit from.
    """
    centers = []
    for n, line in enumerate(open(path), 1):
        text = line.strip()
        if not text:
            continue
        m = TILE_NAME.match(text)
        if not m:
            print(f'{path}:{n}: not a tile name (E<x km>_N<y km>.h5): {text!r}',
                  file=sys.stderr)
            sys.exit(2)
        centers.append((int(m.group(1)) * 1000, int(m.group(2)) * 1000))
    return centers


def s3_names(prefix):
    """The .h5 names directly under an s3:// prefix, as a set.

    ONE listing for the whole region rather than one call per tile.  A
    listing that fails is an error, not an empty set: an empty set would
    read as "every prelim tile is missing".
    """
    p = subprocess.run(['aws', 's3', 'ls', prefix.rstrip('/') + '/'],
                       capture_output=True, text=True, timeout=600)
    if p.returncode not in (0, 1) or (p.returncode == 1 and p.stderr.strip()):
        print(f'could not list {prefix}: {p.stderr.strip()[:300]}', file=sys.stderr)
        sys.exit(2)
    return {line.split()[-1] for line in p.stdout.splitlines()
            if line.strip().endswith('.h5') and not line.lstrip().startswith('PRE')}


def missing_prelim_tiles(centers, tile_prefix, lister=s3_names):
    """Centers whose OWN prelim tile is not at <tile_prefix>/prelim/."""
    present = lister(f'{tile_prefix.rstrip("/")}/prelim')
    return [tile_name(x0, y0) for x0, y0 in centers
            if tile_name(x0, y0) not in present]


def config_declares(config, name):
    """Does algorithm_config.yml declare an input called `name`?

    THE CHEAP, OFFLINE, AUTHORITATIVE-ABOUT-INTENT CHECK, and the one that
    runs first: it needs no network, and an input absent from the config was
    never registered by anybody.
    """
    return any(i.get('name') == name for i in config.get('inputs') or [])


def deployed_declares(process, name):
    """Does the DEPLOYED CWL declare it?  True / False / None if unreadable.

    The config says what we meant to register; the CWL says what is actually
    deployed, and the two differ for exactly as long as a registration takes
    to build.  None is NOT treated as yes: repo.maap-project.org timed out
    while this was being written, and an earlier version of this function
    returned True on that failure -- turning a hard guard into a pass because
    a host was slow, at a cost of one failed job per tile.
    """
    link = process.get('cwlLink')
    if not link:
        return None
    try:
        text = requests.get(link, timeout=30).text
    except requests.RequestException as exc:
        print(f'NOTE: could not read the deployed CWL ({exc.__class__.__name__});'
              ' falling back to algorithm_config.yml alone.')
        return None
    return re.search(rf'^\s*{re.escape(name)}:', text, re.M) is not None


def wait_for_slot(maap, live, limit):
    """Block until fewer than `limit` of the jobs in `live` are unfinished."""
    while len(live) >= limit:
        still = []
        for jid in live:
            try:
                status = maap.get_job_status(jid).json().get('status', '?')
            except Exception:
                status = '?'          # a job not yet visible counts as live
            if str(status) not in DONE:
                still.append(jid)
        live[:] = still
        if len(live) >= limit:
            print(f'    {len(live)} in flight; waiting {POLL_S}s')
            time.sleep(POLL_S)


def main():
    parser = argparse.ArgumentParser(
        description='Submit one DPS job per tile center for one region.')
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--tile_list')
    source.add_argument('--xy_file')
    parser.add_argument('--step', required=True, choices=['prelim', 'matched'])
    parser.add_argument('--args_url', required=True)
    parser.add_argument('--ledger')
    parser.add_argument('--queue', default=DEFAULT_QUEUE)
    parser.add_argument('--tag')
    parser.add_argument('--tile_prefix')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--rate', type=float, default=2.0)
    parser.add_argument('--max_in_flight', type=int)
    parser.add_argument('--replace', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    # "-" is tile_prefix's registered default and means "none", exactly as in
    # run.sh -- so `--tile_prefix -` cannot slip past the matched check below.
    if args.tile_prefix in ('', '-'):
        args.tile_prefix = None

    region = region_of(args.args_url)
    tag = args.tag or (f'{region}_{args.step}' if region else None)
    if not tag:
        print(f'cannot read a region from {args.args_url!r} -- pass --tag',
              file=sys.stderr)
        sys.exit(2)
    ledger = args.ledger or f'{tag}_jobs.csv'

    if not (args.args_url.startswith('s3://') or os.path.isfile(args.args_url)):
        print(f'--args_url {args.args_url!r} is neither an s3:// URI nor a file',
              file=sys.stderr)
        sys.exit(2)
    if not args.queue.startswith('maap-dps-'):
        print(f'--queue {args.queue!r} does not look like a DPS queue',
              file=sys.stderr)
        sys.exit(2)
    if os.path.exists(ledger) and not (args.replace or args.dry_run):
        print(f'{ledger} exists.  It is the only record of the jobs it names;\n'
              '  overwriting it loses their ids.  Move it aside, pass a\n'
              '  different --ledger, or --replace if you mean to discard it.',
              file=sys.stderr)
        sys.exit(2)

    source = args.tile_list or args.xy_file
    centers = read_tile_list(args.tile_list) if args.tile_list \
        else read_centers(args.xy_file)
    if args.limit:
        centers = centers[:args.limit]

    maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
    config = load_config()
    name, version = config['algorithm_name'], config['algorithm_version']
    process = find_process(maap, name, version)
    pid = process.get('processID')

    inputs_common = {'step': args.step, 'args_file': args.args_url}
    if args.tile_prefix:
        # THE I7 INPUT, AND IT IS NOT DEPLOYED YET.  Refuse rather than submit
        # an input the process does not declare: the CWL runner's complaint
        # about an unexpected input is not obviously that, and it would be one
        # failed job per tile before anyone read a log.
        if not config_declares(config, 'tile_prefix'):
            print('--tile_prefix given, but algorithm_config.yml declares no'
                  ' such input.\n  It is docs/plan_IS_run.sh I7 (QI4/QI5):'
                  ' run.sh and algorithm_config.yml\n  need the change, then a'
                  ' rebuild and a registration by Ben.', file=sys.stderr)
            sys.exit(2)
        if deployed_declares(process, 'tile_prefix') is False:
            print(f'--tile_prefix is in algorithm_config.yml but {name}:{version}'
                  ' as DEPLOYED\n  does not have it: the registration has not'
                  ' built and deployed yet.\n  Check with'
                  ' scripts/maap/check_build_id.py.', file=sys.stderr)
            sys.exit(2)
        inputs_common['tile_prefix'] = args.tile_prefix
    if args.step == 'matched' and not args.tile_prefix:
        print('--step matched needs --tile_prefix: a matched job reads its 8'
              ' neighbours\'\n  prelim tiles and has no other way to find'
              ' them (plan_IS_run.sh I7).', file=sys.stderr)
        sys.exit(2)
    if args.step == 'matched':
        missing = missing_prelim_tiles(centers, args.tile_prefix)
        if missing:
            print(f'{len(missing)} of {len(centers)} centers have no prelim tile'
                  f' at {args.tile_prefix}/prelim/:\n    ' + '\n    '.join(missing) +
                  '\n  Nothing submitted.  Each is a no-data center not yet pruned'
                  '\n  (scripts/maap/prune_tile_list.py) or a failed prelim still'
                  '\n  to be investigated; its matched job could only fail.',
                  file=sys.stderr)
            sys.exit(2)

    print(f'process {name}:{version}  processID={pid}'
          f"  (modified {process.get('lastModifiedTime')})")
    print(f'{len(centers)} centers from {source}  step={args.step}'
          f'  queue={args.queue}')
    print(f'args  {args.args_url}')
    print(f'ledger {ledger}\n')

    if args.dry_run:
        for x0, y0 in centers:
            print(f'would submit {tag}_E{int(x0/1000)}_N{int(y0/1000)}'
                  f'  x0={x0} y0={y0}')
        print(f'\n--dry-run: nothing submitted ({len(centers)} jobs).')
        return

    live, failed = [], 0
    with open(ledger, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(LEDGER_COLUMNS)
        for n, (x0, y0) in enumerate(centers, 1):
            if args.max_in_flight:
                wait_for_slot(maap, live, args.max_in_flight)
            ident = f'{tag}_E{int(x0 / 1000)}_N{int(y0 / 1000)}'
            inputs = dict(inputs_common, x0=str(x0), y0=str(y0))
            try:
                r = maap.submit_job(pid, inputs, args.queue,
                                    dedup=False, tag=ident)
                job_id = job_id_from(r) if 200 <= r.status_code < 300 else None
                if not job_id:
                    job_id = f'<submit failed: HTTP {r.status_code}: {r.text[:200]}>'
            except Exception as exc:
                job_id = f'<submit failed: {type(exc).__name__}: {exc}>'
            if job_id.startswith('<'):
                failed += 1
            else:
                live.append(job_id)
            writer.writerow([ident, x0, y0, args.step, args.queue,
                             args.args_url, job_id,
                             datetime.datetime.now(datetime.timezone.utc)
                             .isoformat(timespec='seconds'),
                             args.tile_prefix or '-'])
            fh.flush()
            print(f'  {n:4}/{len(centers)}  {ident:34} {job_id}')
            if n < len(centers):
                time.sleep(args.rate)

    print(f'\n{len(centers) - failed}/{len(centers)} submitted -> {ledger}')
    if failed:
        print(f'{failed} SUBMIT FAILURES, recorded in the ledger; re-run for'
              ' those centers with a different --ledger.')
    print(f'Watch:  scripts/maap/collect_jobs.py {ledger}')
    print(f'Fetch:  scripts/maap/fetch_tiles.py {ledger} <region_dir>'
          f' --step {args.step}')


if __name__ == '__main__':
    main()
