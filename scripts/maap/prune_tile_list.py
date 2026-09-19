#!/usr/bin/env python3
"""
Remove no-data centers from a region's tile list, after a prelim fan-out.

docs/plan_tile_lists.sh TL3.  Ben, 2026-09-18 (plan_monthly_on_maap.sh AM5):
"If any tiles fail for lack of data on the prelim step, they should be deleted
from this list."  The lists are ATL1415/resources/<region>/40km_tile_list.txt,
and submit_MAAP_jobs.py --tile_list submits from them.

THE RULE, for each PRELIM row of the ledger:
  successful, and no tile at <tile_prefix>/prelim/<name>   -> PRUNE
  successful, tile present                                  -> keep
  failed / dismissed / never submitted                      -> INVESTIGATE, kept
  anything still running                                    -> refuse: too early

WHY THAT IS ENOUGH, AND NO LOG IS READ.  Both no-data exits -- the prelim fit's
(TL1) and the uncertainty step's (I7a, plan_IS_run.sh) -- end in a SUCCESSFUL
job that leaves no tile, while anything else that goes wrong still fails the
job.  So "successful, no tile" is exactly "no data", and a failed job is never
pruned: it is a real fault until someone has read its log.
CONSEQUENCE: a no-data job from a build older than TL1 reports as FAILED and is
listed to investigate, not pruned (the monthly E1020_N-2580 job 9266c3d7 is
one).  That is the old build's verdict, not this rule's.

Usage:
  prune_tile_list.py <prelim_ledger.csv> <tile_list.txt> [--write]

Without --write it only prints what it would remove.  --write rewrites the
list in place, keeping every other line and its order; commit the change so
it can be reviewed -- the list is not edited by hand.

Exit: 0 all prelim jobs accounted for; 1 some jobs need investigating (the
no-data centers were still pruned, with --write); 2 could not decide.
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ogc_jobs import DONE  # noqa: E402
from submit_MAAP_jobs import TILE_NAME, s3_names, tile_name  # noqa: E402


def json_or_empty(response):
    try:
        body = response.json()
    except ValueError:
        return {}
    return body if isinstance(body, dict) else {}


def classify(rows, status_of, present_under):
    """
    Sort prelim ledger rows by the rule above.

    rows           ledger rows (dicts), prelim only
    status_of      job_id -> status string
    present_under  tile_prefix -> the set of tile names under <prefix>/prelim/

    Returns {'prune': [...], 'keep': [...], 'investigate': [(name, why)],
             'running': [(name, status)], 'unknown_prefix': [name]}.
    """
    out = {'prune': [], 'keep': [], 'investigate': [], 'running': [],
           'unknown_prefix': []}
    for row in rows:
        name = tile_name(float(row['x0']), float(row['y0']))
        jid = row.get('job_id') or ''
        if not jid or jid.startswith('<'):
            out['investigate'].append((name, 'never submitted: ' + jid[:80]))
            continue
        status = status_of(jid)
        if status not in DONE:
            out['running'].append((name, status))
            continue
        if status != 'successful':
            out['investigate'].append((name, f'job {jid} {status}'))
            continue
        prefix = row.get('tile_prefix') or '-'
        if prefix == '-':
            # no prefix: nowhere to look for the tile, so no verdict on it
            out['unknown_prefix'].append(name)
            continue
        if name in present_under(prefix):
            out['keep'].append(name)
        else:
            out['prune'].append(name)
    return out


def rewrite(path, drop):
    """The list without the names in `drop`: every other line kept, in order."""
    with open(path) as fh:
        lines = fh.readlines()
    kept = [line for line in lines if line.strip() not in drop]
    with open(path, 'w') as fh:
        fh.writelines(kept)
    return len(lines) - len(kept)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Prune no-data centers from a tile list after a prelim fan-out.')
    parser.add_argument('ledger')
    parser.add_argument('tile_list')
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args(argv)

    with open(args.ledger, newline='') as fh:
        rows = list(csv.DictReader(fh))
    steps = {row.get('step') for row in rows}
    if steps != {'prelim'}:
        print(f'{args.ledger}: steps {sorted(map(str, steps))} -- this takes a PRELIM'
              ' ledger only.\n  Pruning is decided by the prelim step (AM5).',
              file=sys.stderr)
        return 2

    listed = [line.strip() for line in open(args.tile_list) if line.strip()]
    bad = [name for name in listed if not TILE_NAME.match(name)]
    if bad:
        print(f'{args.tile_list}: not tile names: {bad[:5]} -- fix the list first.',
              file=sys.stderr)
        return 2

    from maap.maap import MAAP
    maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
    listings = {}

    def present_under(prefix):
        if prefix not in listings:
            listings[prefix] = s3_names(f'{prefix.rstrip("/")}/prelim')
        return listings[prefix]

    def status_of(jid):
        return str(json_or_empty(maap.get_job_status(jid)).get('status', '?'))

    out = classify(rows, status_of, present_under)

    if out['running']:
        print(f"{len(out['running'])} jobs are not finished -- whether they write a"
              ' tile is not known yet.  Run again when they are:', file=sys.stderr)
        for name, status in out['running']:
            print(f'    {name}  {status}', file=sys.stderr)
        return 2
    if out['unknown_prefix']:
        print(f"{len(out['unknown_prefix'])} rows have no tile_prefix, so there is"
              ' nowhere to look for their tiles:\n    '
              + '\n    '.join(out['unknown_prefix']), file=sys.stderr)
        return 2

    in_list = set(listed)
    to_drop = [name for name in out['prune'] if name in in_list]
    already = [name for name in out['prune'] if name not in in_list]
    print(f'{args.ledger}: {len(rows)} prelim jobs')
    print(f"  kept, tile written       : {len(out['keep'])}")
    print(f'  no data, PRUNE           : {len(to_drop)}'
          + (''.join(f'\n      {n}' for n in to_drop)))
    if already:
        print(f'  no data, already absent  : {len(already)}'
              + ''.join(f'\n      {n}' for n in already))
    if out['investigate']:
        print(f"  INVESTIGATE, kept        : {len(out['investigate'])}"
              + ''.join(f'\n      {n}  ({why})' for n, why in out['investigate']))

    if to_drop and args.write:
        n = rewrite(args.tile_list, set(to_drop))
        print(f'\nremoved {n} from {args.tile_list} -- commit it for review.')
    elif to_drop:
        print(f'\n--write not given: {args.tile_list} unchanged.')
    return 1 if out['investigate'] else 0


if __name__ == '__main__':
    sys.exit(main())
