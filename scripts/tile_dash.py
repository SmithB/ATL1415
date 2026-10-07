#!/usr/bin/env python3
"""
tile_dash.py -- terminal dashboard for a tile run, on MAAP (DPS ledgers) or on
discover (slurm task-queue run directories).  Adapted from parallel_boss's
pdash.py: the top box has one timeline bar per ledger / run directory (jobs in
submission order), and the rest of the screen is a map of the region, one
character per block of 40 km tile centers, colored by state.

  MAAP:      tile_dash.py --ledger 'maap_ledgers/GL_0332_maskv5_monthly_prelim_run*_jobs.csv'
             [--tile_list ATL1415/resources/GL/40km_tile_list.txt] [--interval 60]
  discover:  tile_dash.py --run_dir GL_prelim [--run_dir ...] [--tile_list ...] [--interval 10]

QUERIES (MAAP).  Never one call per job per refresh.  Each refresh lists the
user's 'accepted' and 'running' jobs (one paged call each) and checks only the
ledger jobs that have left those lists since the last refresh -- they finished
-- with one get_job_status each.  A job seen successful or failed is never
asked about again.  The first refresh has to learn the state of jobs that are
already finished: at most --max_status_calls of them per refresh, the rest
shown as unknown until later refreshes reach them.
discover reads only the run directory (file names; each task file once).

A tile's state is that of its LATEST job across the ledgers / run dirs given,
so a retry that succeeds turns a failed cell done.

Map: '.' not started (in --tile_list, no job), ':' queued (MAAP accepted /
discover queue), digit = number of tiles running in the block ('+' for 10 or
more), '#' done, 'X' failed (any failed tile in the block, unless something
there is running), '?' unknown.
"""
import argparse
import csv
import datetime
import glob
import os
import re
import shutil
import sys
import time

DONE_STATES = {'successful'}
FAIL_STATES = {'failed', 'dismissed', 'deleted', 'offline'}
TILE_RE = re.compile(r'E(-?\d+)_N(-?\d+)')
XY0_RE = re.compile(r'--xy0[ =]+(-?[\d.]+)[ ,]+(-?[\d.]+)')

C = {'reset': '\033[0m', 'dim': '\033[2m', 'yellow': '\033[33m', 'green': '\033[32;1m',
     'blue': '\033[34;1m', 'red': '\033[31;1m', 'magenta': '\033[35m'}


def color(text, name, use):
    return f'{C[name]}{text}{C["reset"]}' if use else text


# ------------------------------------------------------------------ sources

class Job:
    __slots__ = ('key', 'xy', 'order', 'source', 'state')

    def __init__(self, key, xy, order, source, state='unknown'):
        self.key, self.xy, self.order, self.source, self.state = key, xy, order, source, state


class MaapSource:
    """Jobs from submit_MAAP_jobs.py ledgers; states from a few list calls."""

    def __init__(self, ledger_globs, max_status_calls=50, page_size=500):
        from maap.maap import MAAP
        self.maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
        self.globs, self.max_status_calls, self.page_size = ledger_globs, max_status_calls, page_size
        self.jobs = {}          # job_id -> Job
        self.active_last = set()
        self.calls = 0
        self.note = ''

    def _load_ledgers(self):
        files = sorted({f for g in self.globs for f in glob.glob(g)})
        for f in files:
            name = os.path.basename(f).replace('_jobs.csv', '')
            with open(f) as fh:
                for n, r in enumerate(csv.DictReader(fh)):
                    jid = r['job_id']
                    if jid in self.jobs:
                        continue
                    # km, truncated toward zero: the tile name's own rule
                    xy = (int(float(r['x0']) / 1000), int(float(r['y0']) / 1000))
                    if xy == (0, 0) and r.get('task', '-') not in ('-', ''):
                        xy = None         # a mosaic / nc task: timeline only
                    state = 'submit_failed' if jid.startswith('<') else 'unknown'
                    self.jobs[jid] = Job(jid, xy, (r.get('submitted_utc', ''), n), name, state)
        return files

    def _list(self, status):
        ids, offset = set(), 0
        while True:
            r = self.maap.list_jobs(status=status, page_size=self.page_size, offset=offset,
                                    get_job_details=False)
            self.calls += 1
            r.raise_for_status()
            page = r.json().get('jobs', [])
            new = {j['jobID'] for j in page if 'jobID' in j} - ids
            # the API caps a page (250 seen 2026-10-07) whatever page_size asks
            # for: page on until a page brings nothing new
            if not new:
                return ids
            ids |= new
            offset += len(page)

    def refresh(self):
        self.calls, self.note = 0, ''
        files = self._load_ledgers()
        try:
            listed = {s: self._list(s) for s in ('accepted', 'running')}
        except Exception as e:
            self.note = f'list_jobs failed ({type(e).__name__}); states as of last refresh'
            return files
        active = listed['accepted'] | listed['running']
        for jid, job in self.jobs.items():
            if jid in listed['running']:
                job.state = 'running'
            elif jid in listed['accepted']:
                job.state = 'accepted'
        # jobs that left the active lists, then jobs never seen: ask, within budget
        left = [j for j in self.active_last - active
                if j in self.jobs and self.jobs[j].state in ('running', 'accepted')]
        unknown = [j for j, job in self.jobs.items()
                   if job.state == 'unknown' and j not in active]
        budget = max(self.max_status_calls, len(left))
        for jid in (left + unknown)[:budget]:
            try:
                st = self.maap.get_job_status(jid).json().get('status', 'unknown')
            except Exception:
                st = 'unknown'
            self.calls += 1
            self.jobs[jid].state = st if st in DONE_STATES | FAIL_STATES else 'unknown'
        self.active_last = active
        n_unknown = sum(job.state == 'unknown' for job in self.jobs.values())
        if n_unknown:
            self.note = f'{n_unknown} job states not yet known (learning {self.max_status_calls}/refresh)'
        return files


class DiscoverSource:
    """Jobs from slurm task-queue run directories (setup_slurm_run.py layout)."""

    def __init__(self, run_dirs):
        self.run_dirs = run_dirs
        self.xy_cache = {}      # (run_dir, task) -> list of (x, y) km
        self.jobs = {}
        self.calls = 0
        self.note = ''

    def _xy(self, run_dir, task, path):
        key = (run_dir, task)
        if key not in self.xy_cache:
            xys = []
            try:
                with open(path, errors='replace') as fh:
                    for line in fh:
                        m = XY0_RE.search(line)
                        if m:
                            xys.append((int(float(m.group(1)) / 1000), int(float(m.group(2)) / 1000)))
            except OSError:
                return []
            self.xy_cache[key] = xys
        return self.xy_cache[key]

    def refresh(self):
        self.jobs = {}
        for rd in self.run_dirs:
            name = os.path.basename(os.path.normpath(rd))
            errors = {re.search(r'task_(\d+)', f).group(1)
                      for f in glob.glob(os.path.join(rd, 'error_logs', 'task_*'))
                      if re.search(r'task_(\d+)', f)}
            for sub, state in (('queue', 'accepted'), ('running', 'running'), ('done', 'successful')):
                for f in glob.glob(os.path.join(rd, sub, 'task_*')):
                    m = re.match(r'task_(\d+)', os.path.basename(f))
                    if not m:
                        continue
                    task = m.group(1)
                    st = 'failed' if (task in errors and state == 'successful') else state
                    xys = self._xy(rd, task, f) or [None]
                    for k, xy in enumerate(xys):
                        self.jobs[(rd, task, k)] = Job((rd, task, k), xy, (int(task), k), name, st)
        return self.run_dirs


# ------------------------------------------------------------------ rendering

def timeline(jobs, width, use):
    """jobs in order -> bar: = done  # running  : queued  X failed  ? unknown"""
    if not jobs:
        return '-' * width
    sym = {'successful': ('=', 'blue'), 'running': ('#', 'green'), 'accepted': (':', 'yellow'),
           'unknown': ('?', 'magenta')}
    cells = []
    n = len(jobs)
    for c in range(width):
        chunk = jobs[c * n // width:max((c + 1) * n // width, c * n // width + 1)]
        states = [j.state for j in chunk]
        if any(s in FAIL_STATES or s == 'submit_failed' for s in states):
            cells.append(color('X', 'red', use))
        else:
            # the column's most common state
            s = max(set(states), key=states.count)
            ch, col = sym.get(s, ('?', 'magenta'))
            cells.append(color(ch, col, use))
    return ''.join(cells)


def render_map(tiles, all_xy, rows, cols, use):
    """tiles: {(x,y) km: state}; all_xy: every center to place (incl. not started)."""
    if not all_xy or rows < 3 or cols < 10:
        return []
    xs = sorted({x for x, _ in all_xy}); ys = sorted({y for _, y in all_xy})
    dx = min((b - a for a, b in zip(xs, xs[1:]) if b > a), default=40)
    dy = min((b - a for a, b in zip(ys, ys[1:]) if b > a), default=40)
    nx = (xs[-1] - xs[0]) // dx + 1
    ny = (ys[-1] - ys[0]) // dy + 1
    # tiles per character: a character is ~2x as tall as wide, so a square
    # block of tiles is k wide and k tall in tiles, drawn 1 char per k tiles
    # across and 1 char per 2k tiles down
    k = 1
    while (nx + k - 1) // k > cols or (ny + 2 * k - 1) // (2 * k) > rows:
        k += 1
    kx, ky = k, 2 * k
    if (nx + kx - 1) // kx * 2 <= cols:   # room to draw each block two chars wide
        wide = True
    else:
        wide = False
    blocks = {}
    for xy in all_xy:
        bx = ((xy[0] - xs[0]) // dx) // kx
        by = ((ys[-1] - xy[1]) // dy) // ky      # north at the top
        blocks.setdefault((bx, by), []).append(tiles.get(xy, 'not_started'))
    out = []
    for by in range(max(b[1] for b in blocks) + 1):
        line = []
        for bx in range(max(b[0] for b in blocks) + 1):
            st = blocks.get((bx, by))
            if st is None:
                ch, col = ' ', None
            else:
                n_run = sum(s == 'running' for s in st)
                if n_run:
                    ch, col = (str(n_run) if n_run < 10 else '+'), 'green'
                elif any(s in FAIL_STATES or s == 'submit_failed' for s in st):
                    ch, col = 'X', 'red'
                elif all(s in DONE_STATES for s in st):
                    ch, col = '#', 'blue'
                elif any(s == 'accepted' for s in st):
                    ch, col = ':', 'yellow'
                elif any(s == 'unknown' for s in st):
                    ch, col = '?', 'magenta'
                elif any(s in DONE_STATES for s in st):
                    ch, col = '#', 'blue'          # part done, rest not started
                else:
                    ch, col = '.', 'dim'
            cell = (ch * 2 if wide else ch)
            line.append(color(cell, col, use) if col else cell)
        out.append(''.join(line).rstrip())
    out.append(f'  1 char = {kx}x{ky} tiles ({kx * dx} x {ky * dy} km)'
               + ('' if not wide else ', drawn 2 wide'))
    return out


def render(src, sources, tile_list_xy, interval, use):
    term = shutil.get_terminal_size((100, 40))
    W = max(60, term.columns)
    inner = W - 4
    jobs = list(src.jobs.values())
    # latest job per tile
    tiles = {}
    for j in sorted(jobs, key=lambda j: (j.order, str(j.source))):
        if j.xy is not None:
            tiles[j.xy] = j.state
    counts = {}
    for s in tiles.values():
        counts[s] = counts.get(s, 0) + 1
    n_tiles = len(set(tiles) | set(tile_list_xy))
    n_done = sum(s in DONE_STATES for s in tiles.values())
    n_fail = sum(s in FAIL_STATES or s == 'submit_failed' for s in tiles.values())
    div = '+' + '-' * (W - 2) + '+'

    def row(text, raw_len=None):
        pad = inner - (raw_len if raw_len is not None else len(text))
        return f'| {text}{" " * max(0, pad)} |'

    now = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%SZ')
    lines = [div, row(f'tile_dash   {now}   refresh {interval:.0f}s   API calls last refresh: {src.calls}'), div]
    lines.append(row(f'tiles {n_tiles}:  done {n_done}   running {counts.get("running", 0)}   '
                     f'queued {counts.get("accepted", 0)}   failed {n_fail}   '
                     f'unknown {counts.get("unknown", 0)}   not started {n_tiles - len(tiles)}'))
    name_w = min(34, max((len(str(s)) for s in sources), default=8))
    bar_w = max(10, inner - name_w - 14)
    by_src = {}
    for j in jobs:
        by_src.setdefault(j.source, []).append(j)
    for name in sorted(by_src, key=lambda s: min(j.order for j in by_src[s])):
        js = sorted(by_src[name], key=lambda j: j.order)
        bar = timeline(js, bar_w, use)
        label = f'{str(name)[-name_w:]:<{name_w}}'
        lines.append(row(f'{label} [{bar}] {len(js):>6}', raw_len=name_w + bar_w + 10))
    if src.note:
        lines.append(row(src.note[:inner]))
    lines.append(div)
    lines.append('  = # : X ?  done running queued failed unknown   map: . not started, digit = running in block')
    rows_left = term.lines - len(lines) - 2
    all_xy = set(tiles) | set(tile_list_xy)
    lines += render_map(tiles, all_xy, rows_left, W - 2, use)
    return '\n'.join(lines)


def read_tile_list(path):
    out = set()
    if path:
        for line in open(path):
            m = TILE_RE.search(line)
            if m:
                out.add((int(m.group(1)), int(m.group(2))))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--ledger', action='append', help='ledger CSV or glob (MAAP); repeatable')
    g.add_argument('--run_dir', action='append', help='slurm task-queue run directory (discover); repeatable')
    p.add_argument('--tile_list', help='E<x>_N<y> list: every center of the region, to show not-started cells')
    p.add_argument('--interval', type=float, help='seconds between refreshes (default 60 MAAP, 10 discover)')
    p.add_argument('--max_status_calls', type=int, default=50,
                   help='MAAP: most get_job_status calls per refresh for jobs of unknown state (default 50)')
    p.add_argument('--once', action='store_true', help='print once and exit')
    p.add_argument('--no-color', action='store_true')
    args = p.parse_args(argv)
    use = not args.no_color and sys.stdout.isatty()
    if args.ledger:
        src = MaapSource(args.ledger, max_status_calls=args.max_status_calls)
        interval = args.interval or 60
    else:
        src = DiscoverSource(args.run_dir)
        interval = args.interval or 10
    tile_list_xy = read_tile_list(args.tile_list)
    try:
        while True:
            sources = src.refresh()
            if not src.jobs and not tile_list_xy:
                print('tile_dash: no jobs found in ' + ', '.join(args.ledger or args.run_dir), file=sys.stderr)
                return 1
            text = render(src, [os.path.basename(str(s)).replace('_jobs.csv', '') for s in sources],
                          tile_list_xy, interval, use)
            if not args.once:
                sys.stdout.write('\033[2J\033[H')
            print(text, flush=True)
            if args.once:
                return 0
            time.sleep(interval)
    except KeyboardInterrupt:
        return 0


if __name__ == '__main__':
    sys.exit(main())
