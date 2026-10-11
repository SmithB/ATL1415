#!/usr/bin/env python3
"""
Run every tile of one prelim or matched step as PACKED DPS jobs, in rounds,
resubmitting PER TILE (plan_pack_tiles K5c; QK9: one reusable driver instead
of a per-run copy in maap_ledgers).

Each round:
  1. packs its tiles (pack_tiles.py rules, the --predict tables) into
     <stem>_r<n>_pack.txt; round 0 is every listed tile (matched: those with
     their own prelim tile -- the others are logged and dropped, since a
     matched job cannot run without one);
  2. records its start time in <stem>_r<n>_start.txt and submits with
     submit_MAAP_jobs.py --pack_file, ledger <stem>_r<n>_jobs.csv;
  3. waits until every job is final (ogc_jobs.JobStates: list calls, not a
     status call per job per poll; job ids in <stem>_cancelled.txt count as
     final -- MAAP has left cancelled jobs 'running');
  4. collects (collect_jobs.py --csv -> <stem>_r<n>_tiles.csv) and decides
     PER TILE: DONE if its .h5 is under <tile_prefix>/<step>/ with a
     LastModified after the round's start, or (prelim) collect says nodata.
     Every other tile goes into the next round, packed again: a job that died
     part-way loses only its unfinished tiles.  A tile that was not done
     after sharing a node (a 2-lane job) runs ALONE next time (one lane, as
     if its peak were --node_gib): an out-of-memory kill from two tiles that
     were each predicted to fit must not repeat until the retries run out.
After --max_retries rounds the tiles still not done are written to
<stem>_still_failed_tile_list.txt.  The log ends with ALL_DONE, or STOPPED and
the reason.

RESTARTABLE: a round whose ledger exists is not submitted again, only waited
for; tiles of its pack file missing from its ledgers (the ADE died while it was
submitting) are submitted to <stem>_r<n>_resume<k>_jobs.csv first.  Run it
detached, with the python that has maap-py:
  cd ~/git_repos/ATL1415 && setsid nohup /srv/conda/envs/notebook/bin/python -u \\
      scripts/maap/drive_packed.py <options> >> <stem>_driver.log 2>&1 &

--dry-run packs round 0, prints the summary and stops: nothing is submitted.
"""
import argparse
import csv
import datetime
import glob
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import pack_tiles  # noqa: E402

FINAL = {'successful', 'failed', 'dismissed', 'deleted'}


def log(*a):
    print(time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), *a, flush=True)


def stop(why):
    log(f'STOPPED: {why}')
    sys.exit(1)


def ledger_tiles(row):
    """E<x>_N<y> names a ledger row submitted (its `tiles` column)."""
    packed = row.get('tiles', '-') or '-'
    if packed == '-':
        return []
    return ['E%d_N%d' % (int(float(x) / 1000), int(float(y) / 1000))
            for lane in packed.split('|') for pair in lane.split(';') for x, y in [pair.split(',')]]


def pack_tiles_of(path):
    """tile names in a pack file, in order, without .h5"""
    return [pack_tiles.tile_key(t) for line in open(path) for t in line.replace('|', ' ').split()]


def shared_node_tiles(ledger_rows):
    """tiles that ran in a job with more than one lane"""
    return {t for r in ledger_rows if '|' in (r.get('tiles', '-') or '-') for t in ledger_tiles(r)}


def outcome(tiles, new_h5, tile_rows, step):
    """(done, todo): a tile is done if its .h5 is new since the round
    started, or (prelim) collect says nodata."""
    nodata = {r['tile'] for r in tile_rows if r['status'] == 'nodata'} if step == 'prelim' else set()
    done = [t for t in tiles if f'{t}.h5' in new_h5 or t in nodata]
    return done, [t for t in tiles if t not in set(done)]


def s3_mtimes(prefix):
    """{name.h5: LastModified (aware, UTC)} directly under an s3:// prefix.
    A failed listing stops the driver: an empty answer would read as "no
    tile was written" and resubmit everything."""
    import s3fs
    fs = s3fs.S3FileSystem()
    path = prefix.replace('s3://', '').rstrip('/')
    try:
        entries = fs.ls(path, detail=True, refresh=True)
    except FileNotFoundError:
        return {}
    return {os.path.basename(e['name']): e['LastModified'] for e in entries
            if e['name'].endswith('.h5') and 'LastModified' in e}


class Round:
    def __init__(self, a, n):
        self.a, self.n = a, n
        self.stem = f'{a.stem}_r{n}'
        self.pack = f'{self.stem}_pack.txt'

    def ledgers(self):
        return ([f'{self.stem}_jobs.csv'] if os.path.exists(f'{self.stem}_jobs.csv') else []) + \
            sorted(glob.glob(f'{self.stem}_resume*_jobs.csv'))

    def rows(self):
        return [r for f in self.ledgers() for r in csv.DictReader(open(f))]

    def start(self):
        return datetime.datetime.fromisoformat(open(f'{self.stem}_start.txt').read().strip())

    def submit(self, pack_file, ledger_stem, tag):
        a = self.a
        cmd = [sys.executable, f'{HERE}/submit_MAAP_jobs.py', '--pack_file', pack_file,
               '--step', a.step, '--args_url', a.args_url, '--tile_prefix', a.tile_prefix,
               '--queue', a.queue, '--tag', tag, '--ledger', f'{ledger_stem}_jobs.csv',
               '--max_in_flight', str(a.max_in_flight)]
        t0 = time.time()
        p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
        open(f'{ledger_stem}_submit.log', 'a').write(p.stdout + p.stderr)
        log(f'round {self.n}: submit {os.path.basename(pack_file)} exit {p.returncode}'
            f' in {time.time() - t0:.0f} s; {(p.stdout.strip().splitlines() or ["-"])[-3:]}')
        if not os.path.exists(f'{ledger_stem}_jobs.csv'):
            stop(f'round {self.n} wrote no ledger; see {ledger_stem}_submit.log')

    def run_submission(self, tiles, pred):
        if not os.path.exists(f'{self.stem}_jobs.csv'):
            write_pack(self.pack, tiles, pred, self.a)
            open(f'{self.stem}_start.txt', 'w').write(
                datetime.datetime.now(datetime.timezone.utc).isoformat())
            self.submit(self.pack, self.stem, f'{self.a.tag}_r{self.n}')
            return
        log(f'round {self.n}: ledger exists, not submitting again')
        have = {t for r in self.rows() for t in ledger_tiles(r)}
        lines = [line for line in open(self.pack) if line.strip()]
        missing = [line for line in lines
                   if not {pack_tiles.tile_key(t) for t in line.replace('|', ' ').split()} & have]
        if missing:
            k = len(self.ledgers())
            resume = f'{self.stem}_resume{k}'
            open(f'{resume}_pack.txt', 'w').writelines(missing)
            log(f'round {self.n}: {len(missing)} of {len(lines)} jobs never reached a ledger;'
                f' submitting them as resume{k}')
            self.submit(f'{resume}_pack.txt', resume, f'{self.a.tag}_r{self.n}_resume{k}')

    def wait(self, maap):
        from ogc_jobs import JobStates
        ids = [r['job_id'] for r in self.rows() if not r['job_id'].startswith('<')]
        js, deadline = JobStates(maap), time.time() + self.a.round_timeout
        cancelled_file = f'{self.a.stem}_cancelled.txt'
        while True:
            cancelled = ({line.strip() for line in open(cancelled_file) if line.strip()}
                         if os.path.exists(cancelled_file) else set())
            by = js.poll(ids)
            by.update({j: 'dismissed' for j in ids if j in cancelled})
            counts = {}
            for s in by.values():
                counts[s] = counts.get(s, 0) + 1
            log(f'round {self.n}: {len(ids)} jobs {dict(sorted(counts.items()))}; {js.calls} API calls'
                + (f'; LISTING FAILED: {js.error[:100]}' if js.error else ''))
            if all(by[j] in FINAL for j in ids):
                return
            if time.time() > deadline:
                stop(f'round {self.n} not finished after {self.a.round_timeout} s')
            time.sleep(self.a.poll_secs)

    def collect(self):
        rows = []
        for f in self.ledgers():
            out = f.replace('_jobs.csv', '_tiles.csv')
            p = subprocess.run([sys.executable, f'{HERE}/collect_jobs.py', f, '--csv', out],
                               cwd=REPO, capture_output=True, text=True)
            open(f.replace('_jobs.csv', '_collect.txt'), 'w').write(p.stdout + p.stderr)
            if os.path.exists(out):
                rows += list(csv.DictReader(open(out)))
        return rows


def write_pack(path, tiles, pred, a):
    jobs, _, big = pack_tiles.pack(tiles, pred, a.target_secs, a.big_gib, a.big_secs, a.slowdown)
    with open(path, 'w') as f:
        for job in jobs:
            f.write(' | '.join(' '.join(f'{t}.h5' for t in lane) for lane in job) + '\n')
    walls = [pack_tiles.job_wall(j, pred, a.slowdown) for j in jobs]
    log(f'packed {len(tiles)} tiles -> {len(jobs)} jobs ({len(big)} big tiles in one lane);'
        f' predicted wall per job max {max(walls):.0f} s;'
        f' node peak max {max(pack_tiles.job_peak(j, pred) for j in jobs):.1f} GiB -> {path}')
    return jobs


def predictions(a, tiles):
    pred = pack_tiles.read_predictions(a.predict, a.step, a.slowdown, a.overhead_secs)
    missing = [t for t in tiles if t not in pred]
    if missing and a.default_gib is None:
        stop(f'{len(missing)} tiles have no `ok` {a.step} row in --predict: {missing[:20]};'
             ' give --default_gib and --default_secs')
    measured = [t for t in tiles if t in pred]
    if a.big_secs is None and len(measured) > 1:
        import statistics
        a.big_secs = statistics.quantiles([pred[t][1] for t in measured], n=10)[-1]
    for t in missing:
        pred[t] = (a.default_gib, a.default_secs)
    return pred


def parse(argv):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--step', required=True, choices=('prelim', 'matched'))
    p.add_argument('--tile_list', required=True)
    p.add_argument('--args_url', required=True)
    p.add_argument('--tile_prefix', required=True)
    p.add_argument('--queue', required=True)
    p.add_argument('--tag', required=True)
    p.add_argument('--stem', required=True, help='ledger path stem, e.g. maap_ledgers/GL_0332_prelim')
    p.add_argument('--predict', required=True, nargs='+')
    p.add_argument('--default_gib', type=float)
    p.add_argument('--default_secs', type=float)
    p.add_argument('--target_secs', type=float, default=3600.)
    p.add_argument('--node_gib', type=float, default=26.)
    p.add_argument('--big_gib', type=float)
    p.add_argument('--big_secs', type=float)
    p.add_argument('--slowdown', type=float, default=1.5)
    p.add_argument('--overhead_secs', type=float, default=0.)
    p.add_argument('--max_in_flight', type=int, default=100)
    p.add_argument('--max_retries', type=int, default=2)
    p.add_argument('--round_timeout', type=float, default=8 * 3600.)
    p.add_argument('--poll_secs', type=float, default=120.)
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args(argv)
    a.big_gib = a.node_gib / 2 if a.big_gib is None else a.big_gib
    if 2 * a.big_gib > a.node_gib:
        p.error(f'--big_gib {a.big_gib} would let two lanes need more than --node_gib {a.node_gib}')
    if (a.default_gib is None) != (a.default_secs is None):
        p.error('--default_gib and --default_secs go together')
    return a


def main(argv=None):
    a = parse(argv)
    tiles = pack_tiles.read_tile_names(a.tile_list)
    if os.path.exists(f'{a.stem}_r0_jobs.csv'):
        # a restart: round 0's pack file is the tile set (matched drops done)
        tiles = pack_tiles_of(f'{a.stem}_r0_pack.txt')
    elif a.step == 'matched':
        have = set(s3_mtimes(f'{a.tile_prefix.rstrip("/")}/prelim'))
        dropped = [t for t in tiles if f'{t}.h5' not in have]
        if dropped:
            log(f'{len(dropped)} tiles have no prelim tile, not submitted: {dropped[:20]}'
                f'{" ..." if len(dropped) > 20 else ""}')
            open(f'{a.stem}_no_prelim_tile_list.txt', 'w').write(
                ''.join(f'{t}.h5\n' for t in dropped))
        tiles = [t for t in tiles if t not in set(dropped)]
    if not tiles:
        log('ALL_DONE: no tiles to run')
        return
    pred = predictions(a, tiles)
    if a.dry_run:
        write_pack(f'{a.stem}_r0_pack.txt', tiles, pred, a)
        log('--dry-run: nothing submitted')
        return

    from maap.maap import MAAP
    maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
    for n in range(a.max_retries + 1):
        r = Round(a, n)
        if os.path.exists(r.pack):
            tiles = pack_tiles_of(r.pack)
        r.run_submission(tiles, pred)
        r.wait(maap)
        rows = r.collect()
        recent = {name for name, t in s3_mtimes(f'{a.tile_prefix.rstrip("/")}/{a.step}').items()
                  if t >= r.start()}
        done, todo = outcome(tiles, recent, rows, a.step)
        log(f'round {n}: {len(done)} of {len(tiles)} tiles done, {len(todo)} not: {todo[:20]}'
            f'{" ..." if len(todo) > 20 else ""}')
        if not todo:
            log('ALL_DONE: every tile is done')
            return
        alone = sorted(set(todo) & shared_node_tiles(r.rows()))
        for t in alone:
            pred[t] = (a.node_gib, pred[t][1])
        if alone:
            log(f'round {n}: {len(alone)} of them shared a node; they run alone next round')
        tiles = todo
    open(f'{a.stem}_still_failed_tile_list.txt', 'w').write(''.join(f'{t}.h5\n' for t in tiles))
    log(f'ALL_DONE with {len(tiles)} tiles not done after {a.max_retries} retries:'
        f' {a.stem}_still_failed_tile_list.txt')


if __name__ == '__main__':
    main()
