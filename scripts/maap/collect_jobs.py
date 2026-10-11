#!/usr/bin/env python3
"""
Status, time, memory, input size and CROSSOVER COUNT for every job in a ledger.

REGION-GENERIC, and named so since 2026-09-12.  It was collect_jobs.py --
the name was the only AA-specific thing about it; it has always read whatever
ledger it was handed, and the AA records in howto_MAAP_ogc O7/O8 and
howto_MAAP_AA were written under the old name and still describe this file.

IT READS ABOUT TILES; IT DOES NOT MOVE THEM.  Copying the solved .h5 files out
of DPS output and into the region tree is fetch_tiles.py -- its sibling, and
the other half of what the howtos call "collect" (docs/plan_IS_run.sh QI3).

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

WHICH BUILD RAN EACH TILE is a column too, from the one-line BUILD_ID summary
run.sh prints at the top of every tile job -- the record that lets a run be
audited after the fact (howto_MAAP_ogc O11).  After the table the collector
lists the builds the ledger's tiles ran, with a tile count for each; a list,
not a warning, since a run that reruns patched tiles mixes builds on purpose
(what changed between them is the run's final step, O12b, done by hand).  It
WARNS only if a tile ran with MAAP_PGT unset (no NSIDC).  Jobs from before
per-tile stamping show '-'.

WHAT MACHINE RAN EACH TILE, and how much CPU each step really got
(docs/plan_dps_speed.sh D1): the instance type is a column, from run.sh's
WORKER: line; each step line adds the cores it averaged (CPU time / wall),
the median node load, and the hypervisor steal and cgroup throttling over
the step, from run_with_rusage.py's "=== cpu" line.  After the table, the
workers: one line per EC2 instance, with how many of the ledger's jobs ran
on it -- more than one means jobs shared a node.  A bench job's BENCH: line
is printed under it.  Jobs from before D1 show '-'.

Time and peak memory are what the job reports about ITSELF
(scripts/run_with_rusage.py, one line per fit / error / matched step);
get_job_metrics() is used only for the wall clock, because its machine and
memory fields come back null.

ONE ROW PER TILE, MACHINE-READABLE (--csv, plan_pack_tiles K5a): the same
numbers as the table, for one-tile and packed jobs alike -- the prediction
table pack_tiles.py reads, and the per-tile outcome a driver resubmits from.
Columns: tile (E<x>_N<y>, no .h5), step, status, secs, solve_secs,
max_mem_GiB, instance, lane, identifier, job_id.
  status  ok / nodata / failed for a tile that ran.  A one-tile job's tile is
          ok when the job is successful, nodata when its log says it wrote no
          tile; a tile that never reported (its job failed first, is still
          running, or was not submitted) has the job's status -- except that a
          successful packed job's unreported tile is failed, as run.sh says.
  secs    a one-tile job: the job's wall time (start-up and stage-in
          included); a packed tile: TILE_START to TILE_END (waiting for a
          lane's earlier tiles not included).  solve_secs: the sum of the
          tile's rusage steps, the same thing measured the same way in both.

Usage:
  collect_jobs.py <ledger.csv> [--csv <out.csv>]
      ledger  the CSV the submitter wrote.  Required: defaulting it invited a
              run against the wrong region's ledger, which reads as a table of
              jobs that inexplicably will not advance.
      --csv   also write one row per tile to <out.csv>.
"""
import csv
import datetime
import os
import re
import sys

from maap.maap import MAAP

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ogc_jobs import read_logs, parse_worker, parse_cpu, parse_bench  # noqa: E402



def parse_argv(argv):
    """(ledger, csv path or None); exits 2 with the usage on anything else."""
    if len(argv) == 1 and not argv[0].startswith('-'):
        return argv[0], None
    if len(argv) == 3 and argv[1] == '--csv' and not argv[0].startswith('-'):
        return argv[0], argv[2]
    print(__doc__.strip().rsplit('Usage:', 1)[-1], file=sys.stderr)
    sys.exit(2)


N_RE   = re.compile(r'decimate_data: N_target:[^,]+, N=(\d+)')
XO_RE  = re.compile(r'Decimate_data: N_AT=(\d+), N_XO=(\d+)')
RUSAGE = re.compile(r'=== rusage \[(\w+)\]: elapsed ([\d.]+) s, peak RSS ([\d.]+) GiB')
FIT_RE = re.compile(r'initial: (\d+):')
ITER_RE = re.compile(r'starting qr solve for iteration (\d+)')
BUILD_RE = re.compile(r'^BUILD_ID: (.*)$', re.M)
# a one-tile prelim job that wrote no tile (run.sh: both no-data exits)
NODATA_RE = re.compile(r'^no (fit written|tile to upload) for ', re.M)
TILE_ID_RE = re.compile(r'E-?\d+_N-?\d+$')
# a packed job (run.sh --tiles, plan_pack_tiles K4): each tile's log is printed
# after its banner, and the summary repeats one TILE_STATUS line per tile
TILE_LOG_RE = re.compile(r'^=+ TILE (\S+) log =+$', re.M)
TILE_STATUS_RE = re.compile(r'^TILE_STATUS (\S+) (ok|nodata|failed)$', re.M)
TILE_TIME_RE = re.compile(r'^TILE_(START|END) (\S+) lane (\d+) .*?(\S+Z)$', re.M)
PACK_SUMMARY = 'packed job summary'


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
    steps = solve_fields(text, out)
    out['_tiles'] = packed_tiles(text)
    out['_nodata'] = bool(NODATA_RE.search(text))
    build = BUILD_RE.search(text)
    fields = dict(f.split('=', 1) for f in build.group(1).split() if '=' in f) if build else {}
    out['commit'] = fields.get('commit', '-')[:7]
    out['commit_full'] = fields.get('commit', '-')
    out['version'] = fields.get('algorithm_version', '-')
    out['built'] = fields.get('built', '-')
    out['maap_pgt'] = fields.get('maap_pgt', '-')
    worker = parse_worker(text)
    out['instance'] = worker.get('instance_type', '-')
    out['_worker'], out['_cpu'], out['_bench'] = worker, parse_cpu(text), parse_bench(text)
    return out, steps


def solve_fields(text, out):
    """The solve's own lines in <text> -> fields in <out>; returns the steps
    {label: (secs, peak GiB)}.  Prefer what the job measured about itself
    over anything DPS reports."""
    n, xo, fit = N_RE.search(text), XO_RE.search(text), FIT_RE.search(text)
    iters = ITER_RE.findall(text)
    steps = {k: (float(t), float(g)) for k, t, g in RUSAGE.findall(text)}
    out['max_mem_GiB'] = (f'{max(g for _, g in steps.values()):.2f}'
                          if steps else '-')
    out['N_ATL11'] = n.group(1) if n else '-'
    out['N_AT'] = xo.group(1) if xo else '-'
    out['N_XO'] = xo.group(2) if xo else '-'
    out['N_fit'] = fit.group(1) if fit else '-'
    out['iters'] = str(max(map(int, iters)) + 1) if iters else '-'
    return steps


def packed_tiles(text):
    """A packed job's tiles, in its order: [(fields, steps)], [] for a
    one-tile job.  Each tile's fields come from its own log section; its
    status from the summary (a tile with no status line: 'failed', as run.sh
    reports it); secs from its TILE_START/TILE_END, i.e. the tile's own wall
    time, waiting included; lane is the lane it ran in."""
    if PACK_SUMMARY not in text:
        return []
    marks = list(TILE_LOG_RE.finditer(text))
    end_all = text.find(PACK_SUMMARY)
    status = dict(TILE_STATUS_RE.findall(text[end_all:]))
    times = {}
    for kind, name, lane, stamp in TILE_TIME_RE.findall(text):
        times.setdefault(name, {})[kind] = stamp
        times[name]['lane'] = lane
    tiles = []
    for i, m in enumerate(marks):
        name = m.group(1)
        stop = marks[i + 1].start() if i + 1 < len(marks) else end_all
        fields = {'tile': '  ' + name, 'status': status.get(name, 'failed')}
        steps = solve_fields(text[m.end():stop], fields)
        t = times.get(name, {})
        fields['secs'] = '-'
        if 'START' in t and 'END' in t:
            dt = (datetime.datetime.fromisoformat(t['END'].replace('Z', '+00:00')) -
                  datetime.datetime.fromisoformat(t['START'].replace('Z', '+00:00')))
            fields['secs'] = f'{dt.total_seconds():.0f}'
        fields['queue'] = f"lane {t['lane']}" if 'lane' in t else '-'
        tiles.append((fields, steps))
    return tiles


CSV_COLUMNS = ('tile', 'step', 'status', 'secs', 'solve_secs', 'max_mem_GiB',
               'instance', 'lane', 'identifier', 'job_id')


def ledger_tiles(row):
    """A ledger row's tiles as E<x>_N<y> names, in its order: the `tiles`
    column of a packed job, else the identifier's own tile."""
    packed = row.get('tiles', '-') or '-'
    if packed != '-':
        return ['E%d_N%d' % (int(float(x) / 1000), int(float(y) / 1000))
                for lane in packed.split('|') for pair in lane.split(';')
                for x, y in [pair.split(',')]]
    m = TILE_ID_RE.search(row['identifier'])
    return [m.group(0)] if m else [row['identifier']]


def tile_rows(row, fields, steps):
    """One ledger row and what collect() found -> the --csv rows, one per tile."""
    base = {'step': row.get('step', '-'), 'identifier': row['identifier'],
            'job_id': row.get('job_id', '-'), 'instance': fields.get('instance', '-')}
    job_status = fields.get('status', '?')
    names = ledger_tiles(row)
    if (row.get('tiles', '-') or '-') == '-':
        status = job_status
        if job_status == 'successful':
            status = 'nodata' if fields.get('_nodata') else 'ok'
        solve = f'{sum(t for t, _ in steps.values()):.0f}' if steps else '-'
        return [dict(base, tile=names[0], status=status, secs=fields.get('secs', '-'),
                     solve_secs=solve, max_mem_GiB=fields.get('max_mem_GiB', '-'), lane='-')]
    ran = {f['tile'].strip(): (f, st) for f, st in fields.get('_tiles', [])}
    rows = []
    for name in names:
        if name in ran:
            f, st = ran[name]
            solve = f'{sum(t for t, _ in st.values()):.0f}' if st else '-'
            rows.append(dict(base, tile=name, status=f['status'], secs=f.get('secs', '-'),
                             solve_secs=solve, max_mem_GiB=f.get('max_mem_GiB', '-'),
                             lane=f.get('queue', '-').replace('lane ', '')))
        else:
            rows.append(dict(base, tile=name,
                             status='failed' if job_status == 'successful' else job_status,
                             secs='-', solve_secs='-', max_mem_GiB='-', lane='-'))
    return rows


COLUMNS = (('tile', 26), ('status', 11), ('secs', 7), ('max_mem_GiB', 11),
           ('N_ATL11', 9), ('N_AT', 8), ('N_XO', 7), ('N_fit', 8), ('iters', 5),
           ('commit', 7), ('queue', 22), ('instance', 12))


def main():
    ledger, csv_path = parse_argv(sys.argv[1:])
    maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
    print(' '.join(f'{name:>{width}}' for name, width in COLUMNS))
    builds, no_pgt, workers = {}, [], {}
    out_csv = None
    if csv_path:
        csv_file = open(csv_path, 'w', newline='')
        out_csv = csv.DictWriter(csv_file, fieldnames=CSV_COLUMNS)
        out_csv.writeheader()
    for row in csv.DictReader(open(ledger)):
        try:
            fields, steps = collect(maap, row)
        except Exception as exc:
            print(f"{row['identifier'][-26:]:>26} <collect failed: {type(exc).__name__}: {exc}>")
            if out_csv:
                out_csv.writerows(tile_rows(row, {'status': 'collect failed'}, {}))
            continue
        if out_csv:
            out_csv.writerows(tile_rows(row, fields, steps))
        if fields.get('commit_full', '-') != '-':
            key = (fields['commit_full'], fields['version'], fields['built'])
            builds[key] = builds.get(key, 0) + 1
        if fields.get('maap_pgt') == 'unset':
            no_pgt.append(fields['tile'])
        worker = fields.get('_worker', {})
        if worker.get('instance_id', 'unreachable') != 'unreachable':
            workers.setdefault(worker['instance_id'], [worker, []])[1].append(fields['tile'])
        fields['tile'] = fields['tile'][-26:]
        print(' '.join(f'{fields.get(name, "-"):>{width}}' for name, width in COLUMNS))
        cpu = fields.get('_cpu', {})
        for label, (secs, gib) in steps.items():
            line = f"{'':26} {'step ' + label:>11} {secs:7.0f} {gib:11.2f}"
            if label in cpu:
                c = cpu[label]
                line += (f"   cores {c['cores']}, load {c['load']}, steal {c['steal_s']} s,"
                         f" throttled {c['throttled_s']} s")
            print(line)
        for tile_fields, tile_steps in fields.get('_tiles', []):
            print(' '.join(f'{tile_fields.get(name, "-"):>{width}}' for name, width in COLUMNS))
            for label, (secs, gib) in tile_steps.items():
                print(f"{'':26} {'step ' + label:>11} {secs:7.0f} {gib:11.2f}")
        if fields.get('_bench'):
            print(f"{'':26} {'bench':>11} " +
                  ' '.join(f'{k}={v}' for k, v in fields['_bench'].items()))

    if builds:
        print('\nBuilds that ran these tiles:')
        for (commit, version, built), n in builds.items():
            print(f'  {commit[:7]}  {version}  built {built}  {n} tile(s)')
    if workers:
        print('\nWorkers (EC2 instance: jobs from this ledger that ran on it):')
        for iid, (w, tiles) in workers.items():
            print(f"  {iid}  {w.get('instance_type', '-')}  {w.get('cpu', '-')}"
                  f"  vcpus={w.get('vcpus', '-')} threads_per_core={w.get('threads_per_core', '-')}"
                  f"  cpu_quota={w.get('cpu_quota', '-')}  blas={w.get('blas', '-')}"
                  f"  {len(tiles)} job(s)")
    if out_csv:
        csv_file.close()
    if no_pgt:
        print('\nWARNING: MAAP_PGT WAS UNSET on the workers for: ' + ', '.join(no_pgt) +
              ' -- they could not get NSIDC credentials.')


if __name__ == '__main__':
    main()
