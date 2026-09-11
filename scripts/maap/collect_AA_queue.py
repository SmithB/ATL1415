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

MIXED BUILDS CAN BE INTENDED (howto_MAAP_ogc O12b).  A production run that
patches a bug and reruns only the affected tiles is mixed ON PURPOSE.  Such a
change is declared, whenever convenient, in the run's run_notes.txt -- plain
text, beside the run's args files on the bucket, one line per change:
    # from          to             note (free text to the end of the line)
    on_s3-8935494   on_s3-3f2c1a7  tide-mask edge fix; grounding-line tiles rerun
A build is named by its algorithm_version or by a commit prefix of 7+
characters.  Builds linked by notes, directly or through a chain, are reported
as INTENDED with their notes; anything unlinked still warns, and the warning
prints the line that would declare it.  A name that matches more than one of
the ledger's builds is ignored as ambiguous -- which is what `on_s3` does until
every build has its own tag (O11) -- so a note cannot wave through a stale
worker by accident.

Time and peak memory are what the job reports about ITSELF
(scripts/run_with_rusage.py, one line per fit / error / matched step);
get_job_metrics() is used only for the wall clock, because its machine and
memory fields come back null.

Usage:
  collect_AA_queue.py [ledger.csv] [--notes <uri|path>]
      ledger  default AA_queue_jobs.csv, the submitter's default
      notes   default <directory of the ledger's args_file>/run_notes.txt
"""
import csv
import os
import re
import subprocess
import sys

from maap.maap import MAAP

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ogc_jobs import read_logs  # noqa: E402

_ARGV = list(sys.argv[1:])
NOTES = None
if '--notes' in _ARGV:
    _i = _ARGV.index('--notes')
    NOTES = _ARGV[_i + 1]
    del _ARGV[_i:_i + 2]
LEDGER = _ARGV[0] if _ARGV else 'AA_queue_jobs.csv'

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
    out['commit_full'] = fields.get('commit', '-')
    out['version'] = fields.get('algorithm_version', '-')
    out['maap_pgt'] = fields.get('maap_pgt', '-')
    return out, steps


def read_text(where):
    """A small text file from s3:// or disk, or None if it is not there."""
    if where.startswith('s3://'):
        p = subprocess.run(['aws', 's3', 'cp', where, '-'],
                           capture_output=True, text=True, timeout=120)
        return p.stdout if p.returncode == 0 else None
    return open(where).read() if os.path.isfile(where) else None


def parse_notes(text):
    """run_notes.txt -> [(from, to, note)]; '#' comments and blank lines skipped."""
    notes = []
    for line in (text or '').splitlines():
        if line.lstrip().startswith('#'):
            continue            # a '#' later in a line is note text, and kept
        parts = line.split(None, 2)
        if len(parts) >= 2:
            notes.append((parts[0], parts[1], parts[2] if len(parts) > 2 else ''))
    return notes


def resolve(name, builds):
    """The ledger builds a note's name refers to: by version, or commit prefix."""
    return {key for key, (commit, version) in builds.items()
            if name == version or (len(name) >= 7 and commit.startswith(name))}


def report_builds(builds, tiles_by_build, notes, notes_where):
    """
    Say whether the ledger's mix of builds is declared intended.

    builds: key -> (full commit, version); tiles_by_build: key -> [tiles].
    Builds are linked by notes (union-find), so a chain A->B->C covers A..C.
    """
    if len(builds) < 2:
        return
    parent = {key: key for key in builds}

    def root(key):
        while parent[key] != key:
            key = parent[key]
        return key

    used, ambiguous = [], []
    for a, b, note in notes:
        ra, rb = resolve(a, builds), resolve(b, builds)
        if len(ra) > 1 or len(rb) > 1:
            ambiguous.append((a, b))
            continue
        if ra and rb:
            parent[root(ra.pop())] = root(rb.pop())
            used.append((a, b, note))

    def label(key):
        commit, version = builds[key]
        return f'{commit[:7]} ({version})'

    groups = {}
    for key in builds:
        groups.setdefault(root(key), []).append(key)
    listing = ''.join(f'\n  {label(k)}: {len(tiles_by_build[k])} tile(s)'
                      for k in builds)
    if len(groups) == 1:
        print(f'\nMIXED BUILDS, ALL DECLARED INTENDED in {notes_where}:{listing}')
    else:
        print('\nWARNING: THESE TILES RAN DIFFERENT BUILDS -- a worker used a cached'
              ' older image, or a rebuild landed mid-run -- and run_notes.txt does'
              ' not declare the change intended.  Results from undeclared builds'
              ' must not be mixed:')
        for key in builds:
            print(f'  {label(key)}: {", ".join(tiles_by_build[key])}')
    for a, b, note in used:
        print(f'  INTENDED  {a} -> {b}: {note or "(no note text)"}')
    for a, b in ambiguous:
        print(f'  IGNORED   {a} -> {b}: a name there matches more than one build'
              ' in this ledger; use a commit prefix of 7+ characters.')
    if len(groups) > 1:
        # Offer the exact line for each unlinked group, against the largest.
        ordered = sorted(groups.values(), key=lambda g: -sum(len(tiles_by_build[k]) for k in g))
        main = ordered[0][0]
        where = notes_where or '<the run\'s run_notes.txt>'
        print(f'\n  If a change is intended, declare it in {where}, e.g.:')
        for group in ordered[1:]:
            other = group[0]
            names = [builds[k][1] for k in builds]
            name = (lambda k: builds[k][1] if names.count(builds[k][1]) == 1
                    else builds[k][0][:7])
            print(f'    {name(main)}  {name(other)}  <what changed, and which tiles were rerun>')


COLUMNS = (('tile', 26), ('status', 11), ('secs', 7), ('max_mem_GiB', 11),
           ('N_ATL11', 9), ('N_AT', 8), ('N_XO', 7), ('N_fit', 8), ('iters', 5),
           ('commit', 7), ('queue', 22))


def main():
    maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
    print(' '.join(f'{name:>{width}}' for name, width in COLUMNS))
    rows = list(csv.DictReader(open(LEDGER)))
    notes_where = NOTES
    if notes_where is None:
        args_files = [r.get('args_file', '') for r in rows if r.get('args_file')]
        if args_files:
            notes_where = args_files[0].rsplit('/', 1)[0] + '/run_notes.txt'
    builds, tiles_by_build, no_pgt = {}, {}, []
    for row in rows:
        try:
            fields, steps = collect(maap, row)
        except Exception as exc:
            print(f"{row['identifier'][-26:]:>26} <collect failed: {type(exc).__name__}: {exc}>")
            continue
        if fields.get('commit_full', '-') != '-':
            key = fields['commit_full']
            builds[key] = (fields['commit_full'], fields['version'])
            tiles_by_build.setdefault(key, []).append(fields['tile'])
        if fields.get('maap_pgt') == 'unset':
            no_pgt.append(fields['tile'])
        fields['tile'] = fields['tile'][-26:]
        print(' '.join(f'{fields.get(name, "-"):>{width}}' for name, width in COLUMNS))
        for label, (secs, gib) in steps.items():
            print(f"{'':26} {'step ' + label:>11} {secs:7.0f} {gib:11.2f}")

    if len(builds) > 1:
        text = read_text(notes_where) if notes_where else None
        if text is None:
            print(f'\n(no run notes found{" at " + notes_where if notes_where else ""})')
        report_builds(builds, tiles_by_build, parse_notes(text), notes_where)
    if no_pgt:
        print('\nWARNING: MAAP_PGT WAS UNSET on the workers for: ' + ', '.join(no_pgt) +
              ' -- they could not get NSIDC credentials.')


if __name__ == '__main__':
    main()
