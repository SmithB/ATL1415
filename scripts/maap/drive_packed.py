#!/usr/bin/env python3
"""
Run every tile of one prelim or matched step as PACKED DPS jobs, resubmitting
PER TILE as soon as its job is final (plan_pack_tiles K5c; K9: rolling
resubmission instead of rounds; QK9: one reusable driver instead of a per-run
copy in maap_ledgers).

START: every listed tile is packed once (pack_tiles.py rules, the --predict
tables) into <stem>_plan_pack.txt, one-lane (big, longest) jobs first.
Matched: only tiles with their own prelim tile -- the others are logged to
<stem>_no_prelim_tile_list.txt and dropped, since a matched job cannot run
without one.

EACH POLL (every --poll_secs):
  1. job states (ogc_jobs.JobStates: list calls, not a status call per job per
     poll; job ids in <stem>_cancelled.txt count as final -- MAAP has left
     cancelled jobs 'running');
  2. every job that went final since the last poll is collected
     (collect_jobs rows, appended to <stem>_tiles.csv: a --predict table for
     later runs) and decided PER TILE in <stem>_decisions.csv: DONE if its
     .h5 is under <tile_prefix>/<step>/ with a LastModified after THAT JOB was
     submitted, NODATA (prelim) if collect says so, else RETRY -- or GAVE_UP
     once the tile has had 1 + --max_retries jobs that STARTED it (K10): a
     failed submit, or a job whose log shows no TILE_START for the tile (a
     start-up failure: the runner could not reach the MAAP API), is recorded
     as 'NOT SUBMITTED' / 'not started' and does not count -- up to
     --max_unstarted such jobs per tile, so an outage cannot resubmit a tile
     forever.  A job with no readable log counts.  A tile that was not done
     after sharing a node (a 2-lane job) runs ALONE from then on (one lane, as
     if its peak were --node_gib): an out-of-memory kill from two tiles that
     were each predicted to fit must not repeat;
  3. the tiles to retry are packed together and go FIRST, then the plan's
     jobs not yet submitted, as many as keep --max_in_flight jobs unfinished:
     one batch, <stem>_b<k>_{pack.txt,jobs.csv,submit.log}.
So a start-up failure is resubmitted within a poll of going final, not after
the slowest job of a round (K7c: near-pole tiles waited ~90 min for that).

SECOND QUEUE (Ben 2026-10-18: "make a second queue for jobs likely to need the
64gb queue"): with --big_queue <queue> --big_queue_gib <G>, a tile predicted
to peak above G GiB runs ALONE (one tile per job, its retries too) on <queue>;
every other tile is packed and submitted as before on --queue.  A batch that
holds both kinds is submitted as two ledgers, one per queue.  (AA 2026-10-17:
two 44 km tiles peaked at 30.5-30.7 GiB, over the ~30.4 GiB a 32gb node has.)
The log ends ALL_DONE (every tile done, or "with N tiles not done": their
names in <stem>_still_failed_tile_list.txt) or STOPPED and the reason.

RESTARTABLE: the state is rebuilt from the batch ledgers and the decisions
file at every poll.  A job in a ledger is never submitted again, only waited
for and decided; pack lines that never reached a ledger (the ADE died while
a batch was submitting) are simply still pending.  A stem the round driver
ran before K9 (<stem>_r0_jobs.csv) is refused.  Run it detached, with the
python that has maap-py:
  cd ~/git_repos/ATL1415 && setsid nohup /srv/conda/envs/notebook/bin/python -u \\
      scripts/maap/drive_packed.py <options> >> <stem>_driver.log 2>&1 &

STOPPING: setsid makes the driver a process-group leader; stop the whole group
with `kill -- -<pgid>` (ps -o pgid= -p <pid>).  Killing only the driver's pid
leaves a submit_MAAP_jobs.py child submitting.  Then dismiss what it
submitted.

--dry-run packs the plan into <stem>_dryrun_pack.txt, prints the summary and
stops: nothing is submitted.
"""
import argparse
import csv
import datetime
import glob
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import pack_tiles  # noqa: E402

FINAL = {'successful', 'failed', 'dismissed', 'deleted'}
DECISION_COLUMNS = ('job', 'job_id', 'tile', 'attempt', 'decision', 'tile_status',
                    'shared', 'decided_utc')
TILE_COLUMNS = ('tile', 'step', 'status', 'secs', 'solve_secs', 'max_mem_GiB',
                'instance', 'lane', 'identifier', 'job_id')     # = collect_jobs.CSV_COLUMNS
# a decision's tile_status when the job never ran the tile: not an attempt (K10)
NOT_STARTED = 'not started'
NOT_RUN = ('NOT SUBMITTED', NOT_STARTED)


def log(*a):
    print(time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), *a, flush=True)


def stop(why):
    log(f'STOPPED: {why}')
    sys.exit(1)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc)


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


def line_tiles(line):
    return [pack_tiles.tile_key(t) for t in line.replace('|', ' ').split()]


def submit_failed(row):
    return row['job_id'].startswith('<')


def decide_tile(tile, row, tile_status, mtimes, step):
    """'done' / 'nodata' / None (not done) for one tile of a final job: done if
    its .h5 is new since THIS job was submitted, or (prelim) collect says
    nodata.  A matched job that writes no tile failed: matched has no nodata."""
    if not submit_failed(row):
        submitted = datetime.datetime.fromisoformat(row['submitted_utc'])
        m = mtimes.get(f'{tile}.h5')
        if m is not None and m >= submitted:
            return 'done'
    if step == 'prelim' and tile_status == 'nodata':
        return 'nodata'
    return None


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


def collect_rows(maap, rows):
    """collect_jobs' --csv rows (one per tile) for each ledger row, each with
    `started`: False if the job's log was read and has no TILE_START for the
    tile, else True (an unreadable log counts as started)."""
    import collect_jobs
    out = []
    for r in rows:
        try:
            fields, steps = collect_jobs.collect(maap, r)
        except Exception as exc:
            log(f'collect failed for {r["identifier"]}: {type(exc).__name__}: {exc}')
            fields, steps = {'status': 'collect failed'}, {}
        started = fields.get('_started')
        out.append([dict(x, started=started is None or x['tile'] in started)
                    for x in collect_jobs.tile_rows(r, fields, steps)])
    return out


def append_csv(path, columns, rows):
    new = not os.path.exists(path)
    with open(path, 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        if new:
            w.writeheader()
        w.writerows(rows)


def batch_num(path):
    return int(re.search(r'_b(\d+)_jobs\.csv$', path).group(1))


def pack_lines(jobs):
    return [' | '.join(' '.join(f'{t}.h5' for t in lane) for lane in job) + '\n' for job in jobs]


def big_queue_tiles(a, pred):
    """the tiles that go to --big_queue: predicted peak over --big_queue_gib"""
    if a.big_queue is None:
        return set()
    return {t for t, (gib, _) in pred.items() if gib > a.big_queue_gib}


def pack_split(tiles, pred, a, bigq):
    """jobs for `tiles`: the --big_queue ones one tile per job, first; the
    rest by pack_tiles.pack.  Returns (jobs, the packer's big tiles)."""
    mine = sorted((t for t in tiles if t in bigq), key=lambda t: (-pred[t][1], t))
    rest = [t for t in tiles if t not in bigq]
    jobs, big = [], []
    if rest:
        jobs, _, big = pack_tiles.pack(rest, pred, a.target_secs, a.big_gib, a.big_secs, a.slowdown)
    return [[[t]] for t in mine] + jobs, big


class Run:
    """One step's state, rebuilt from the files at every poll."""

    def __init__(self, a, pred, maap=None):
        self.a, self.pred, self.maap = a, pred, maap
        self.plan = [line for line in open(f'{a.stem}_plan_pack.txt') if line.strip()]
        self.decisions_file = f'{a.stem}_decisions.csv'
        self.bigq = big_queue_tiles(a, pred)
        for line in self.plan:
            self.queue_of(line)          # a plan packed without --big_queue stops here

    def queue_of(self, line):
        """the queue a pack line goes to; a line that mixes --big_queue tiles
        with others cannot be placed, and stops the driver"""
        ts = line_tiles(line)
        n = sum(t in self.bigq for t in ts)
        if n == 0:
            return self.a.queue
        if n == len(ts) == 1:
            return self.a.big_queue
        stop(f'pack line {line.strip()!r} puts a tile predicted over --big_queue_gib'
             f' {self.a.big_queue_gib:g} in a job with others: the plan was packed without'
             ' this --big_queue setting.  Use a new --stem.')

    def jobs(self):
        """[(key, row)] for every ledger row, oldest batch first; key =
        '<ledger name>#<row>' (a submit failure has no job id)."""
        out = []
        for f in sorted(glob.glob(f'{self.a.stem}_b*_jobs.csv'), key=batch_num):
            out += [(f'{os.path.basename(f)}#{i}', r) for i, r in enumerate(csv.DictReader(open(f)))]
        return out

    def decided(self):
        """{job key: {tile: decision row}}"""
        out = {}
        if os.path.exists(self.decisions_file):
            for r in csv.DictReader(open(self.decisions_file)):
                out.setdefault(r['job'], {})[r['tile']] = r
        return out

    def tiles(self, jobs, decided):
        """{tile: (attempts, its latest job's decision row or None, unstarted)}:
        attempts = its jobs so far, except those decided as never having run
        it; unstarted = those (K10).  A job not decided yet is an attempt."""
        latest = {}
        for key, r in jobs:
            for t in ledger_tiles(r):
                n, _, free = latest.get(t, (0, None, 0))
                d = decided.get(key, {}).get(t)
                if d is not None and d['tile_status'] in NOT_RUN:
                    free += 1
                else:
                    n += 1
                latest[t] = (n, key, free)
        return {t: (n, decided.get(key, {}).get(t), free) for t, (n, key, free) in latest.items()}

    def alone(self):
        """tiles not done after sharing a node: one lane from then on"""
        if not os.path.exists(self.decisions_file):
            return set()
        return {r['tile'] for r in csv.DictReader(open(self.decisions_file))
                if r['decision'] in ('retry', 'gave_up') and r['shared'] == '1'
                and r['tile_status'] not in NOT_RUN}

    def decide(self, new, tiles):
        """collect the newly final jobs and decide each of their tiles"""
        rows = collect_rows(self.maap, [r for _, r in new])
        append_csv(f'{self.a.stem}_tiles.csv', TILE_COLUMNS, [x for rs in rows for x in rs])
        mtimes = s3_mtimes(f'{self.a.tile_prefix.rstrip("/")}/{self.a.step}')
        out = []
        for (key, r), trows in zip(new, rows):
            status = {x['tile']: x['status'] for x in trows}
            started = {x['tile']: x.get('started', True) for x in trows}
            shared = '|' in (r.get('tiles', '-') or '-')
            for t in ledger_tiles(r):
                d = decide_tile(t, r, status.get(t, '-'), mtimes, self.a.step)
                n, _, free = tiles[t]
                tile_status = status.get(t, '-')
                if d is None:
                    if submit_failed(r) or not started.get(t, True):
                        # never ran here: not an attempt (tiles() counted this
                        # job as one while it was undecided)
                        n, free = n - 1, free + 1
                        tile_status = 'NOT SUBMITTED' if submit_failed(r) else NOT_STARTED
                    d = ('gave_up' if n > self.a.max_retries or free > self.a.max_unstarted
                         else 'retry')
                out.append({'job': key, 'job_id': r['job_id'], 'tile': t, 'attempt': n,
                            'decision': d, 'tile_status': tile_status,
                            'shared': int(shared), 'decided_utc': utcnow().isoformat(timespec='seconds')})
        append_csv(self.decisions_file, DECISION_COLUMNS, out)
        for (key, r) in new:
            mine = [x for x in out if x['job'] == key]
            redo = [x['tile'] for x in mine if x['decision'] in ('retry', 'gave_up')]
            if redo:
                log(f'{r["identifier"]} ({r["job_id"][:36]}) final: {len(mine) - len(redo)} done,'
                    f' {len(redo)} not: ' + ', '.join(f'{x["tile"]} {x["tile_status"]} -> {x["decision"]}'
                                                       for x in mine if x['tile'] in redo))
        return out

    def pending(self, tiles):
        """pack lines to submit, retries first, then the plan's untouched jobs"""
        retry = sorted(t for t, (_, d, _) in tiles.items() if d is not None and d['decision'] == 'retry')
        lines = []
        if retry:
            pred = dict(self.pred)
            for t in self.alone() & set(retry):
                pred[t] = (self.a.node_gib, pred[t][1])
            jobs, _ = pack_split(retry, pred, self.a, self.bigq)
            lines = pack_lines(jobs)
        seen = set(tiles)
        return lines + [line for line in self.plan if not seen.intersection(line_tiles(line))], set(retry)

    def submit(self, lines, retried):
        """submit pack lines, one batch (ledger) per queue; `retried` = the
        tiles among them being retried, for the log"""
        a = self.a
        queues = [a.queue] + ([a.big_queue] if a.big_queue else [])
        for queue in queues:
            mine = [line for line in lines if self.queue_of(line) == queue]
            if not mine:
                continue
            n_retry = sum(t in retried for line in mine for t in line_tiles(line))
            k = len(glob.glob(f'{a.stem}_b*_jobs.csv'))
            stem = f'{a.stem}_b{k}'
            open(f'{stem}_pack.txt', 'w').writelines(mine)
            cmd = [sys.executable, f'{HERE}/submit_MAAP_jobs.py', '--pack_file', f'{stem}_pack.txt',
                   '--step', a.step, '--args_url', a.args_url, '--tile_prefix', a.tile_prefix,
                   '--queue', queue, '--tag', f'{a.tag}_b{k}', '--ledger', f'{stem}_jobs.csv']
            t0 = time.time()
            p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
            open(f'{stem}_submit.log', 'a').write(p.stdout + p.stderr)
            log(f'batch {k}: {len(mine)} jobs ({n_retry} retried tiles first)'
                f'{f" on {queue}" if queue != a.queue else ""} exit {p.returncode}'
                f' in {time.time() - t0:.0f} s; {(p.stdout.strip().splitlines() or ["-"])[-3:]}')
            if not os.path.exists(f'{stem}_jobs.csv'):
                stop(f'batch {k} wrote no ledger; see {stem}_submit.log')

    def tick(self, states):
        """One poll, given {job_id: status}.  Returns (finished, changed)."""
        jobs = self.jobs()
        decided = self.decided()
        # a job is decided when every one of its tiles has a row (a crash while
        # appending leaves it undecided, and it is decided again in full)
        whole = lambda: {k for k, r in jobs if set(ledger_tiles(r)) <= set(decided.get(k, {}))}
        done_jobs = whole()
        tiles = self.tiles(jobs, decided)
        new = [(k, r) for k, r in jobs if k not in done_jobs
               and (submit_failed(r) or states.get(r['job_id']) in FINAL)]
        if new:
            self.decide(new, tiles)
            decided = self.decided()
            done_jobs = whole()
            tiles = self.tiles(jobs, decided)
        in_flight = [r['job_id'] for k, r in jobs if k not in done_jobs]
        lines, retried = self.pending(tiles)
        counts = {}
        for _, d, _ in tiles.values():
            key = d['decision'] if d is not None else 'running'
            counts[key] = counts.get(key, 0) + 1
        seen = set(tiles)
        counts['unsubmitted'] = sum(len(line_tiles(line)) for line in self.plan
                                    if not seen.intersection(line_tiles(line)))
        log(f'{len(jobs)} jobs, {len(in_flight)} unfinished; tiles {dict(sorted(counts.items()))}')
        free = self.a.max_in_flight - len(in_flight)
        if lines and free > 0:
            self.submit(lines[:free], retried)
            return False, True
        if not lines and not in_flight:
            failed = sorted(t for t, (_, d, _) in tiles.items() if d['decision'] == 'gave_up')
            if failed:
                open(f'{self.a.stem}_still_failed_tile_list.txt', 'w').write(
                    ''.join(f'{t}.h5\n' for t in failed))
                log(f'ALL_DONE with {len(failed)} tiles not done after {self.a.max_retries} retries:'
                    f' {self.a.stem}_still_failed_tile_list.txt')
            else:
                log('ALL_DONE: every tile is done')
            return True, bool(new)
        return False, bool(new)


def write_pack(path, tiles, pred, a):
    bigq = big_queue_tiles(a, pred) & set(tiles)
    jobs, big = pack_split(tiles, pred, a, bigq)
    open(path, 'w').writelines(pack_lines(jobs))
    rest = jobs[len(bigq):]
    if rest:
        walls = [pack_tiles.job_wall(j, pred, a.slowdown) for j in rest]
        log(f'packed {len(tiles) - len(bigq)} tiles -> {len(rest)} jobs ({len(big)} big tiles in one lane);'
            f' predicted wall per job max {max(walls):.0f} s;'
            f' node peak max {max(pack_tiles.job_peak(j, pred) for j in rest):.1f} GiB -> {path}')
    if a.big_queue is not None:
        log(f'{len(bigq)} tiles predicted over {a.big_queue_gib:g} GiB -> one tile per job on'
            f' {a.big_queue}' + (f': {", ".join(f"{t} {pred[t][0]:g} GiB" for t in sorted(bigq))}'
                                 if bigq else ''))
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
    p.add_argument('--big_queue', help='second queue, for tiles predicted over --big_queue_gib'
                                       ' (one tile per job), e.g. maap-dps-worker-64gb')
    p.add_argument('--big_queue_gib', type=float,
                   help='predicted peak (GiB) above which a tile goes to --big_queue')
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
    p.add_argument('--max_retries', type=int, default=2,
                   help='jobs per tile after its first, counting only jobs that started it')
    p.add_argument('--max_unstarted', type=int, default=5,
                   help='jobs per tile that never started it (failed submit or start-up) before it'
                        ' gives up anyway')
    p.add_argument('--stall_secs', type=float, default=8 * 3600.,
                   help='STOP if jobs are unfinished and none has ended or been submitted for this long')
    p.add_argument('--poll_secs', type=float, default=120.)
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args(argv)
    a.big_gib = a.node_gib / 2 if a.big_gib is None else a.big_gib
    if 2 * a.big_gib > a.node_gib:
        p.error(f'--big_gib {a.big_gib} would let two lanes need more than --node_gib {a.node_gib}')
    if (a.default_gib is None) != (a.default_secs is None):
        p.error('--default_gib and --default_secs go together')
    if (a.big_queue is None) != (a.big_queue_gib is None):
        p.error('--big_queue and --big_queue_gib go together')
    if a.big_queue is not None and a.big_queue == a.queue:
        p.error(f'--big_queue is --queue ({a.queue})')
    return a


def main(argv=None):
    a = parse(argv)
    plan = f'{a.stem}_plan_pack.txt'
    if os.path.exists(f'{a.stem}_r0_jobs.csv'):
        stop(f'{a.stem} was run by the round driver (before K9): {a.stem}_r0_jobs.csv;'
             ' finish it with that version, or use a new --stem')
    if os.path.exists(plan) and not a.dry_run:
        tiles = pack_tiles_of(plan)          # a restart: the plan is the tile set
        log(f'restart: {len(tiles)} tiles in {plan}')
    else:
        tiles = pack_tiles.read_tile_names(a.tile_list)
        if a.step == 'matched':
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
        write_pack(f'{a.stem}_dryrun_pack.txt', tiles, pred, a)
        log('--dry-run: nothing submitted')
        return
    if not os.path.exists(plan):
        write_pack(plan, tiles, pred, a)

    from maap.maap import MAAP
    from ogc_jobs import JobStates
    maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
    run, js = Run(a, pred, maap), JobStates(maap)
    cancelled_file = f'{a.stem}_cancelled.txt'
    last_change = time.time()
    while True:
        ids = [r['job_id'] for _, r in run.jobs() if not submit_failed(r)]
        cancelled = ({line.strip() for line in open(cancelled_file) if line.strip()}
                     if os.path.exists(cancelled_file) else set())
        states = js.poll(ids)
        states.update({j: 'dismissed' for j in ids if j in cancelled})
        if js.error:
            log(f'LISTING FAILED: {js.error[:100]}; statuses are the last known')
        finished, changed = run.tick(states)
        if finished:
            return
        if changed:
            last_change = time.time()
        elif time.time() - last_change > a.stall_secs:
            stop(f'no job ended and none was submitted for {a.stall_secs:.0f} s')
        time.sleep(a.poll_secs)


if __name__ == '__main__':
    main()
