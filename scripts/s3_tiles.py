#!/usr/bin/env python3
"""
Move solved tiles between a DPS worker and the canonical tile tree on S3.

WORKER-SIDE, and deliberately thin: it runs inside the image, beside
run_with_rusage.py, and imports nothing but s3fs.  The maap-py scripts in
scripts/maap/ are the ADE's; a worker has no business talking to the job API.

s3fs RATHER THAN THE aws CLI, for the same reason run.sh fetches its args file
that way: build-env.sh proves s3fs importable, and nothing guarantees an aws
binary on maap_base.  It uses the worker's own credential chain, as every
other bucket read in the solve does.

THE CANONICAL TREE (docs/plan_IS_run.sh QI4, answered 2026-09-12) is
    <tile_prefix>/{prelim,matched}/E<x>_N<y>.h5
where tile_prefix is passed to run.sh per job and carries the release, the
hemisphere-with-period suffix and the region, e.g.
    s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north/IS
    s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north_monthly/IS
So this file never builds that path itself -- it appends one directory to what
it is given.  Nothing here knows what a region or a release is.

TWO SUBCOMMANDS:
  put    <src> <tile_prefix> <step>      one solved tile (and its field-size
                                         report, if it is there) up
  get    <tile_prefix> <step> <x0> <y0> <spacing> <dest>
                                         the 3x3 neighbourhood down

A MISSING NEIGHBOUR IS NOT AN ERROR (QI5b).  On a small coastal region most
tiles have fewer than 8 neighbours, and a tile with too little data writes
nothing at all, so `get` fetches what exists, NAMES each key it did not find,
and exits 0 regardless.  The tile's OWN prelim file is not special-cased here
either: run.sh already refuses to run a matched solve without it, and one
guard in one place is better than two that can disagree.
"""
import argparse
import os
import sys

import s3fs


def tile_name(x0, y0):
    """'E%d_N%d.h5', truncated toward zero -- ATL11_to_ATL15's own name."""
    return 'E%d_N%d.h5' % (int(x0 / 1000), int(y0 / 1000))


def put(fs, src, tile_prefix, step):
    if not os.path.isfile(src):
        print(f's3_tiles: nothing to upload at {src}', file=sys.stderr)
        return 1
    dest = f'{tile_prefix.rstrip("/")}/{step}/{os.path.basename(src)}'
    print(f's3_tiles: {src} -> {dest}')
    fs.put(src, dest)
    # The field-size report ATL11_to_ATL15 writes beside a prelim tile.  Best
    # effort: its absence is not a reason to fail a solve that succeeded.
    report = os.path.join(os.path.dirname(src), 'field_sizes',
                          os.path.basename(src)[:-3] + '_report.json')
    if os.path.isfile(report):
        r_dest = (f'{tile_prefix.rstrip("/")}/{step}/field_sizes/'
                  f'{os.path.basename(report)}')
        print(f's3_tiles: {report} -> {r_dest}')
        try:
            fs.put(report, r_dest)
        except Exception as exc:
            print(f's3_tiles: NOTE: field-size report not uploaded ({exc})')
    return 0


def get(fs, tile_prefix, step, x0, y0, spacing, dest):
    os.makedirs(dest, exist_ok=True)
    got, missing = [], []
    for dx in (-spacing, 0, spacing):
        for dy in (-spacing, 0, spacing):
            name = tile_name(x0 + dx, y0 + dy)
            src = f'{tile_prefix.rstrip("/")}/{step}/{name}'
            target = os.path.join(dest, name)
            if os.path.isfile(target):
                got.append(name + ' (already here)')
                continue
            try:
                if not fs.exists(src):
                    missing.append(name)
                    continue
                fs.get(src, target)
                got.append(name)
            except Exception as exc:
                # Report and continue: one unreadable neighbour must not cost
                # the eight that are readable.
                missing.append(f'{name} ({type(exc).__name__})')
    print(f's3_tiles: {len(got)}/9 of the neighbourhood of '
          f'{tile_name(x0, y0)} localized into {dest}')
    for name in got:
        print(f'  got     {name}')
    for name in missing:
        print(f'  MISSING {name}')
    if missing:
        print('s3_tiles: missing neighbours are expected at a region edge and '
              'for tiles with too little data to fit; the solve continues '
              'with the priors it has.')
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    sub = parser.add_subparsers(dest='cmd', required=True)

    p_put = sub.add_parser('put')
    p_put.add_argument('src')
    p_put.add_argument('tile_prefix')
    p_put.add_argument('step')

    p_get = sub.add_parser('get')
    p_get.add_argument('tile_prefix')
    p_get.add_argument('step')
    p_get.add_argument('x0', type=float)
    p_get.add_argument('y0', type=float)
    p_get.add_argument('spacing', type=float)
    p_get.add_argument('dest')

    args = parser.parse_args()
    fs = s3fs.S3FileSystem()
    if args.cmd == 'put':
        return put(fs, args.src, args.tile_prefix, args.step)
    return get(fs, args.tile_prefix, args.step, args.x0, args.y0,
               args.spacing, args.dest)


if __name__ == '__main__':
    sys.exit(main())
