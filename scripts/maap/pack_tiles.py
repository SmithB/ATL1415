#!/usr/bin/env python3
"""
Pack a tile list into DPS jobs: write the pack file that
submit_MAAP_jobs.py --pack_file reads (plan_pack_tiles K5b).

PREDICTIONS come from collect_jobs.py --csv tables of earlier runs of the SAME
STEP (prelim numbers never stand in for matched).  Per tile: its last `ok` row
in the order the tables are given, so list the newest run last.  Time is
solve_secs (the tile's rusage steps, measured the same way in one-tile and
packed jobs); a row from a job that ran 2 lanes is divided by --slowdown, so
every prediction is a one-lane time; --overhead_secs is added to each, for
what a tile spends outside its solve steps (K6b matched: ~25 s fetching its
neighbours; prelim ~5 s).  The instance type is not corrected for:
the 32gb queue mixes r5.xlarge, m5.2xlarge and c5.4xlarge, ~1.5x apart.

RULES (Ben 2026-10-11, QK7: "2 lanes per job, but run the slowest, largest
jobs in a single lane"):
  BIG tiles -- predicted peak > --big_gib (default half of --node_gib, so any
  two other tiles fit on one node together) or predicted time > --big_secs
  (default: the 90th percentile of the listed tiles' predictions) -- go in
  ONE-LANE jobs, in series.  The rest go in 2-LANE jobs; a lane's time is its
  tiles' predictions x --slowdown (K6 G4: 1.52x per tile on r5.xlarge).  Both
  kinds fill longest tile first, each tile into the first job with room
  (the lighter lane of a 2-lane job), up to --target_secs of predicted wall
  time per job (QK1: about an hour).  A tile longer than that alone gets a job
  of its own.  One-lane jobs are written first: they are the longest.

A TILE WITH NO PREDICTION is refused, by name, unless --default_gib and
--default_secs say what to assume for it (QK10: a silent guess would decide
the packing unseen).

Usage:
  pack_tiles.py <tile_list> --step prelim|matched --predict <csv> [<csv> ...]
                --out <pack_file> [options]
"""
import argparse
import csv
import statistics
import sys

STEPS = ('prelim', 'matched')


def tile_key(name):
    """'E160_N-1640.h5' or 'E160_N-1640' -> 'E160_N-1640'."""
    return name[:-3] if name.endswith('.h5') else name


def read_tile_names(path):
    names = [line.strip() for line in open(path) if line.strip()]
    keys = [tile_key(n) for n in names]
    dup = sorted({k for k in keys if keys.count(k) > 1})
    if dup:
        sys.exit(f'{path}: tiles listed more than once: {dup}')
    return keys


def read_predictions(paths, step, slowdown, overhead_secs=0.):
    """{tile: (peak GiB, one-lane secs)} from collect_jobs --csv tables."""
    rows = [r for p in paths for r in csv.DictReader(open(p))]
    lanes = {}
    for r in rows:
        if r.get('lane', '-') != '-':
            lanes.setdefault(r['identifier'], set()).add(r['lane'])
    pred = {}
    for r in rows:
        if r['step'] != step or r['status'] != 'ok':
            continue
        try:
            gib, secs = float(r['max_mem_GiB']), float(r['solve_secs'])
        except ValueError:
            continue
        if len(lanes.get(r['identifier'], ())) > 1:
            secs /= slowdown
        pred[r['tile']] = (gib, secs + overhead_secs)
    return pred


def fill(tiles, pred, n_lanes, factor, target):
    """First fit, longest first: [[lane, ...], ...] of tile names; a tile goes
    to the lighter lane of the first job that has room for it."""
    jobs = []
    for t in sorted(tiles, key=lambda t: (-pred[t][1], t)):
        cost = pred[t][1] * factor
        for job in jobs:
            i = min(range(n_lanes), key=lambda i: job['load'][i])
            if job['load'][i] + cost <= target:
                break
        else:
            job = {'lanes': [[] for _ in range(n_lanes)], 'load': [0.] * n_lanes}
            jobs.append(job)
            i = 0
        job['lanes'][i].append(t)
        job['load'][i] += cost
    return [[lane for lane in job['lanes'] if lane] for job in jobs]


def pack(tiles, pred, target_secs=3600., big_gib=13., big_secs=None, slowdown=1.5):
    """tiles (all with a prediction) -> (jobs, big_secs, big tiles); a job is
    a list of lanes of tile names; the big tiles' one-lane jobs come first."""
    if big_secs is None:
        big_secs = statistics.quantiles([pred[t][1] for t in tiles], n=10)[-1] \
            if len(tiles) > 1 else float('inf')
    big = [t for t in tiles if pred[t][0] > big_gib or pred[t][1] > big_secs]
    rest = [t for t in tiles if t not in set(big)]
    jobs = fill(big, pred, 1, 1., target_secs) + fill(rest, pred, 2, slowdown, target_secs)
    return jobs, big_secs, big


def job_wall(job, pred, slowdown):
    factor = slowdown if len(job) > 1 else 1.
    return max(sum(pred[t][1] for t in lane) * factor for lane in job)


def job_peak(job, pred):
    """The lanes' largest tiles all at once: an upper bound on the node."""
    return sum(max(pred[t][0] for t in lane) for lane in job)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[1],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('tile_list')
    p.add_argument('--step', required=True, choices=STEPS)
    p.add_argument('--predict', required=True, nargs='+')
    p.add_argument('--out', required=True)
    p.add_argument('--target_secs', type=float, default=3600.)
    p.add_argument('--node_gib', type=float, default=26.,
                   help='memory two lanes may use together (default 26 of a 32 GiB node)')
    p.add_argument('--big_gib', type=float, default=None, help='default --node_gib / 2')
    p.add_argument('--big_secs', type=float, default=None, help='default: p90 of the predictions')
    p.add_argument('--slowdown', type=float, default=1.5)
    p.add_argument('--overhead_secs', type=float, default=0.)
    p.add_argument('--default_gib', type=float, default=None)
    p.add_argument('--default_secs', type=float, default=None)
    a = p.parse_args(argv)
    big_gib = a.node_gib / 2 if a.big_gib is None else a.big_gib
    if 2 * big_gib > a.node_gib:
        p.error(f'--big_gib {big_gib} would let two lanes need more than --node_gib {a.node_gib}')
    if (a.default_gib is None) != (a.default_secs is None):
        p.error('--default_gib and --default_secs go together')

    tiles = read_tile_names(a.tile_list)
    pred = read_predictions(a.predict, a.step, a.slowdown, a.overhead_secs)
    missing = [t for t in tiles if t not in pred]
    if missing and a.default_gib is None:
        sys.exit(f'{len(missing)} of {len(tiles)} tiles have no `ok` {a.step} row in --predict:'
                 f' {missing[:20]}{" ..." if len(missing) > 20 else ""}\n'
                 f'Give --default_gib and --default_secs to assume a size for them.')
    for t in missing:
        pred[t] = (a.default_gib, a.default_secs)
    # big_secs from the measured tiles only: defaults must not move the cut
    measured = [t for t in tiles if t not in missing]
    if a.big_secs is None and len(measured) > 1:
        a.big_secs = statistics.quantiles([pred[t][1] for t in measured], n=10)[-1]
    jobs, big_secs, big = pack(tiles, pred, a.target_secs, big_gib, a.big_secs, a.slowdown)

    with open(a.out, 'w') as f:
        for job in jobs:
            f.write(' | '.join(' '.join(f'{t}.h5' for t in lane) for lane in job) + '\n')
    one = [j for j in jobs if len(j) == 1]   # big tiles' jobs, and a leftover lone tile
    walls = [job_wall(j, pred, a.slowdown) for j in jobs]
    print(f'{len(tiles)} {a.step} tiles -> {len(jobs)} jobs in {a.out}'
          f' ({len(one)} one-lane, {len(jobs) - len(one)} two-lane);'
          f' {len(missing)} tiles at the default size')
    print(f'big: peak > {big_gib:g} GiB or time > {big_secs:.0f} s'
          f' -> {len(big)} tiles in one-lane jobs')
    print(f'tiles per job: median {statistics.median(sum(map(len, j)) for j in jobs):g},'
          f' max {max(sum(map(len, j)) for j in jobs)}')
    print(f'predicted wall per job: median {statistics.median(walls):.0f} s,'
          f' max {max(walls):.0f} s (target {a.target_secs:.0f})')
    print(f'largest predicted node peak: {max(job_peak(j, pred) for j in jobs):.1f} GiB'
          f' (budget {a.node_gib:g})')


if __name__ == '__main__':
    main()
