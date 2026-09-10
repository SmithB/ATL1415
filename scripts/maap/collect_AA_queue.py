#!/usr/bin/env python3
"""
Status, time, memory, input size and CROSSOVER COUNT for every job in a ledger.

OGC SYSTEM, since 2026-09-10 (docs/howto_MAAP_ogc.sh O7).  maap-py 5.x
dropped getJob/getJobResult/getJobMetrics; this uses get_job_status,
get_job_result and get_job_metrics, and reads the job's logs with
ogc_jobs.read_logs() -- BOTH _stdout.txt and _stderr.txt, because under the
CWL runner the solve's own lines (decimate_data, rusage) are in _stderr.txt,
and on the legacy system they were in _stdout.txt.  The OGC job endpoints
serve legacy jobs too, so old ledgers read the same way.

N_XO is a column now.  howto_MAAP_AA 3b-i / howto_MAAP_ogc O8 is decided by
it -- the crossover read works if N_XO > 0 -- and the collector used to leave
it to a hand-run grep.  It comes from the solve's
    Decimate_data: N_AT=<n>, N_XO=<n>
line, the FIT step's (the first): E220_N20 printed N_AT=935506, N_XO=0 before
the crossover fix.  N_ATL11 is still decimate_data's N=, which counts both.

WHICH BUILD RAN EACH TILE is a column too.  Since 2026-09-10 run.sh prints its
one-line BUILD_ID summary at the top of every tile job, because MAAP's runner
reuses a worker's cached image for a tag without re-pulling it -- so one
build_id job cannot vouch for every worker.  After the table the collector
WARNS if the ledger's tiles ran more than one build (a stale cached image, or
a rebuild landing mid-run) or any ran with MAAP_PGT unset (no NSIDC).  Jobs
from before per-tile stamping show '-'.

Time and peak memory are what the job reports about ITSELF
(scripts/run_with_rusage.py, one line per fit / error / matched step);
get_job_metrics() is used only for the wall clock, because its machine and
memory fields come back null.

Usage:
  collect_AA_queue.py [ledger.csv]        (default AA_queue_jobs.csv, the
                                            submitter's default)
"""
import csv
import os
import re
import sys

from maap.maap import MAAP

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ogc_jobs import read_logs  # noqa: E402

LEDGER = sys.argv[1] if len(sys.argv) > 1 else 'AA_queue_jobs.csv'

N_RE   = re.compile(r'decimate_data: N_target:[^,]+, N=(\d+)')
XO_RE  = re.compile(r'Decimate_data: N_AT=(\d+), N_XO=(\d+)')
RUSAGE = re.compile(r'=== rusage \[(\w+)\]: elapsed ([\d.]+) s, peak RSS ([\d.]+) GiB')
FIT_RE = re.compile(r'initial: (\d+):')
ITER_RE = re.compile(r'starting qr solve for iteration (\d+)')
BUILD_RE = re.compile(r'^BUILD_ID: (.*)$', re.M)


def json_or_empty(response):
    try:
        body = response.json()
    except ValueError:
        return {}
    return body if isinstance(body, dict) else {}


def collect(maap, row):
    """One ledger row -> a dict of the table's fields (strings, '-' if absent)."""
    jid = row['job_id']
    out = {'tile': row['identifier'], 'queue': row.get('queue', '-')}
    if not jid or jid.startswith('<'):
        out['status'] = 'NOT SUBMITTED'
        return out, {}
    out['status'] = str(json_or_empty(maap.get_job_status(jid)).get('status', '?'))
    secs = json_or_empty(maap.get_job_metrics(jid)).get('job_duration_seconds')
    out['secs'] = f'{float(secs):.0f}' if secs is not None else '-'

    text, _, _ = read_logs(json_or_empty(maap.get_job_result(jid)))
    n, xo, fit = N_RE.search(text), XO_RE.search(text), FIT_RE.search(text)
    iters = ITER_RE.findall(text)
    # Prefer what the job measured about itself over anything DPS reports.
    steps = {k: (float(t), float(g)) for k, t, g in RUSAGE.findall(text)}
    out['max_mem_GiB'] = (f'{max(g for _, g in steps.values()):.2f}'
                          if steps else '-')
    out['N_ATL11'] = n.group(1) if n else '-'
    out['N_AT'] = xo.group(1) if xo else '-'
    out['N_XO'] = xo.group(2) if xo else '-'
    out['N_fit'] = fit.group(1) if fit else '-'
    out['iters'] = str(max(map(int, iters)) + 1) if iters else '-'
    build = BUILD_RE.search(text)
    fields = dict(f.split('=', 1) for f in build.group(1).split() if '=' in f) if build else {}
    out['commit'] = fields.get('commit', '-')[:7]
    out['maap_pgt'] = fields.get('maap_pgt', '-')
    return out, steps


COLUMNS = (('tile', 26), ('status', 11), ('secs', 7), ('max_mem_GiB', 11),
           ('N_ATL11', 9), ('N_AT', 8), ('N_XO', 7), ('N_fit', 8), ('iters', 5),
           ('commit', 7), ('queue', 22))


def main():
    maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
    print(' '.join(f'{name:>{width}}' for name, width in COLUMNS))
    builds, no_pgt = {}, []
    for row in csv.DictReader(open(LEDGER)):
        try:
            fields, steps = collect(maap, row)
        except Exception as exc:
            print(f"{row['identifier'][-26:]:>26} <collect failed: {type(exc).__name__}: {exc}>")
            continue
        if fields.get('commit', '-') != '-':
            builds.setdefault(fields['commit'], []).append(fields['tile'])
        if fields.get('maap_pgt') == 'unset':
            no_pgt.append(fields['tile'])
        fields['tile'] = fields['tile'][-26:]
        print(' '.join(f'{fields.get(name, "-"):>{width}}' for name, width in COLUMNS))
        for label, (secs, gib) in steps.items():
            print(f"{'':26} {'step ' + label:>11} {secs:7.0f} {gib:11.2f}")

    if len(builds) > 1:
        print('\nWARNING: THESE TILES RAN DIFFERENT BUILDS -- a worker used a cached'
              ' older image, or a rebuild landed mid-run.  Results from different'
              ' builds must not be mixed:')
        for commit, tiles in builds.items():
            print(f'  {commit}: {", ".join(tiles)}')
    if no_pgt:
        print('\nWARNING: MAAP_PGT WAS UNSET on the workers for: ' + ', '.join(no_pgt) +
              ' -- they could not get NSIDC credentials.')


if __name__ == '__main__':
    main()
