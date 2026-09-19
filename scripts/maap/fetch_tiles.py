#!/usr/bin/env python3
"""
Fetch the SOLVED TILES a ledger's DPS jobs produced into the local region tree.

THE COLLECTOR'S SIBLING, AND NOT THE SAME THING.  collect_jobs.py reads the
jobs' LOGS and reports what each tile cost; this moves the ~240 MB .h5 files
themselves, so the matched and mosaic steps have something to read.  The two
were both called "the collector" until 2026-09-12, which is why they are named
apart now (docs/plan_IS_run.sh QI3).

WHY IT EXISTS.  The howtos all say

    aws s3 sync $s3_out/prelim/ $region_dir/prelim/

and that cannot work yet.  Q9's deterministic output prefix is unimplemented,
so a job's products land wherever DPS put them --

    s3://maap-ops-workspace/ben_smith/dps_output/<algo>_<n>/<version>/
        <yyyy>/<mm>/<dd>/<HH>/<MM>/<SS>/<usec>/

-- a different, timestamped prefix for every job, discoverable only through
the job id the submitter wrote into the ledger.  So: walk the ledger, ask
get_job_result where each job's prefix is, copy the tiles down.  No rebuild,
no new CWL input, nothing re-registered.

THE TWO STEPS UPLOAD DIFFERENTLY, and the difference is in run.sh, not here:

  prelim   run.sh passes --base_directory $PWD/output and ATL11_to_ATL15
           appends '/prelim', so the products are
               <prefix>/prelim/E<x>_N<y>.h5
               <prefix>/prelim/field_sizes/E<x>_N<y>_report.json
           VERIFIED 2026-09-11 on AA job 55c01ec3: run.sh's output/prelim/
           subdirectory survives the upload rather than being flattened.

  matched  run.sh passes --out_name $PWD/output/E<x>_N<y>.h5, so the tile is
           at the TOP of the prefix, beside _stdout.txt:
               <prefix>/E<x>_N<y>.h5
           READ FROM run.sh:419 AND UNVERIFIED -- no matched job has ever run.

Both places are tried for every row, so a fetch does not depend on which of
those is right, and the summary says which layout each tile actually came
from.  If matched tiles turn up under <prefix>/matched/ instead, nothing here
needs changing -- that is the third place tried.

A SUCCESSFUL PRELIM JOB WITH NO TILE MEANS NO DATA, not an error: the fit
(docs/plan_tile_lists.sh TL1) or the uncertainty step (plan_IS_run.sh I7a)
found no data, exited 0 and left no tile.  Those rows get the verdict
'no tile', counted apart from failures -- AND THEIR NAMES ARE SAVED to
    <region_dir>/prelim/no_data_tiles.txt
one E<x>_N<y>.h5 per line, the format of ATL1415/resources/<region>/
40km_tile_list.txt, merged with what the file already holds (sorted, unique),
so a retry ledger or a second fetch never drops a name an earlier one found.
Nothing prunes the tile list automatically: after the run, remove those names
from it and commit (Ben, 2026-09-19; plan_tile_lists.sh TL8).  --dry-run
writes nothing and says what it would add.
A MATCHED job with no tile is NOT normal -- since TL1 a matched fit with no
data exits 1 -- so there 'no tile' is listed under NOT FETCHED.

THE TILE NAME IS DERIVED FROM THE LEDGER, not from whatever is on the bucket,
and then required to match: 'E%d_N%d.h5' % (x0/1e3, y0/1e3), truncated toward
zero, which is what ATL11_to_ATL15 builds and what run.sh's awk reproduces.
A job whose prefix holds some other tile is a mis-submission, and silently
syncing it into the region tree would put a wrong tile into the mosaic.

Usage:
  fetch_tiles.py <ledger.csv> <region_dir> [--step prelim|matched]
                 [--replace] [--dry-run] [--any-status]

  <region_dir>   the local region tree, e.g.
                 /home/jovyan/ATL14_processing/rel006/north/IS
                 Tiles land in <region_dir>/<step>/ .
  --step         default prelim.  Also chooses the destination subdirectory.
  --replace      re-fetch a tile that is already present at the same size.
  --any-status   try jobs that are not 'successful' too (they rarely have
                 products; a failed job's prefix is under triaged_job/).
"""
NO_DATA_FILE = 'no_data_tiles.txt'
import argparse
import csv
import os
import subprocess
import sys

from maap.maap import MAAP

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ogc_jobs import s3_prefixes  # noqa: E402


def tile_name(x0, y0):
    """'E%d_N%d.h5', truncated toward zero -- ATL11_to_ATL15's own name."""
    return 'E%d_N%d.h5' % (int(float(x0) / 1000), int(float(y0) / 1000))


def json_or_empty(response):
    try:
        body = response.json()
    except ValueError:
        return {}
    return body if isinstance(body, dict) else {}


def s3_size(uri):
    """The size of exactly this key, or None if it does not exist.

    `aws s3 ls` matches by PREFIX, so E10_N-20.h5 would also match
    E10_N-200.h5: the basename is compared explicitly.
    """
    p = subprocess.run(['aws', 's3', 'ls', uri],
                       capture_output=True, text=True, timeout=120)
    if p.returncode != 0:
        return None
    want = uri.rsplit('/', 1)[-1]
    for line in p.stdout.splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[-1] == want:
            return int(fields[-2])
    return None


def s3_cp(src, dest, dry_run):
    if dry_run:
        return True
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    p = subprocess.run(['aws', 's3', 'cp', src, dest],
                       capture_output=True, text=True, timeout=3600)
    if p.returncode != 0:
        print(f'    cp FAILED: {p.stderr.strip()[:200]}')
        return False
    return True


def candidates(prefix, step, name):
    """Where this tile could be under one job prefix, most likely first."""
    return [(f'{prefix}/{step}/{name}', f'{step}/'),     # prelim, verified
            (f'{prefix}/{name}', 'top level'),           # matched, from run.sh
            (f'{prefix}/output/{step}/{name}', f'output/{step}/')]


def record_no_data(path, names, dry_run=False):
    """
    Merge `names` into the no-data list at `path`; return the names it adds.

    Sorted and unique, so the file is the union of every fetch into this
    region directory and re-running a fetch changes nothing.
    """
    have = set()
    if os.path.isfile(path):
        with open(path) as fh:
            have = {line.strip() for line in fh if line.strip()}
    added = sorted(set(names) - have)
    if added and not dry_run:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as fh:
            fh.writelines(f'{name}\n' for name in sorted(have | set(names)))
    return added


def fetch_row(maap, row, args):
    """One ledger row -> (verdict, bytes fetched)."""
    ident = row.get('identifier', '?')
    jid = row.get('job_id', '')
    if not jid or jid.startswith('<'):
        return 'NOT SUBMITTED', 0

    status = str(json_or_empty(maap.get_job_status(jid)).get('status', '?'))
    if status != 'successful' and not args.any_status:
        return status.upper(), 0

    name = tile_name(row['x0'], row['y0'])
    dest = os.path.join(args.region_dir, args.step, name)
    have = os.path.getsize(dest) if os.path.isfile(dest) else None

    result = json_or_empty(maap.get_job_result(jid))
    for prefix in dict.fromkeys(s3_prefixes(result)):
        for uri, where in candidates(prefix, args.step, name):
            size = s3_size(uri)
            if size is None:
                continue
            if have == size and not args.replace:
                return f'have it ({where})', 0
            print(f'  {ident}: {where}  {size/2**20:.0f} MiB')
            if not s3_cp(uri, dest, args.dry_run):
                return 'CP FAILED', 0
            # The field-size report the solve writes beside a prelim tile.
            stem = name[:-3]
            report = f'{uri.rsplit("/", 1)[0]}/field_sizes/{stem}_report.json'
            if s3_size(report) is not None:
                s3_cp(report, os.path.join(args.region_dir, args.step,
                                           'field_sizes', f'{stem}_report.json'),
                      args.dry_run)
            return ('would fetch' if args.dry_run else 'fetched'), size
    return 'no tile', 0


def main():
    parser = argparse.ArgumentParser(
        description='Copy a ledger\'s solved tiles from DPS output into the '
                    'local region tree.')
    parser.add_argument('ledger')
    parser.add_argument('region_dir')
    parser.add_argument('--step', default='prelim',
                        choices=['prelim', 'matched'])
    parser.add_argument('--replace', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--any-status', action='store_true')
    args = parser.parse_args()

    if not os.path.isdir(args.region_dir):
        print(f'no such region directory: {args.region_dir}', file=sys.stderr)
        sys.exit(2)

    maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
    verdicts, total, no_data = {}, 0, []
    rows = list(csv.DictReader(open(args.ledger)))
    print(f'{len(rows)} rows in {args.ledger} -> '
          f'{os.path.join(args.region_dir, args.step)}'
          f'{"  (--dry-run)" if args.dry_run else ""}\n')
    for row in rows:
        try:
            verdict, size = fetch_row(maap, row, args)
        except Exception as exc:
            verdict, size = f'ERROR {type(exc).__name__}: {exc}'[:60], 0
        verdicts.setdefault(verdict, []).append(row.get('identifier', '?'))
        if verdict == 'no tile':
            no_data.append(tile_name(row['x0'], row['y0']))
        total += size
        print(f'  {row.get("identifier", "?"):34} {verdict}')

    print(f'\n{total / 2**30:.2f} GiB')
    for verdict, idents in sorted(verdicts.items()):
        print(f'  {len(idents):4}  {verdict}')
    # 'no tile' is a normal outcome for PRELIM (no data) and is saved for the
    # cleanup; for matched it is unexpected.  Anything else odd is named.
    normal = ('fetched', 'would fetch') + (('no tile',) if args.step == 'prelim' else ())
    odd = {v: i for v, i in verdicts.items()
           if v not in normal and not v.startswith('have it')}
    if odd:
        print('\nNOT FETCHED:')
        for verdict, idents in sorted(odd.items()):
            print(f'  {verdict}: ' + ', '.join(idents))

    if args.step == 'prelim' and no_data:
        path = os.path.join(args.region_dir, 'prelim', NO_DATA_FILE)
        added = record_no_data(path, no_data, args.dry_run)
        print(f'\nNO DATA -- {len(no_data)} center(s) wrote no tile; nothing to'
              ' run there next time:')
        for name in sorted(no_data):
            print(f'  {name}{"" if name in added else "  (already listed)"}')
        print(f'  {"would add" if args.dry_run else "added"} {len(added)} to {path}')
        print('  AFTER THE RUN: remove these from the region\'s tile list and'
              ' commit, e.g.\n'
              f'    grep -vxFf {path} <tile_list> > t && mv t <tile_list>')


if __name__ == '__main__':
    main()
