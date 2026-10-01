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

FOUR SUBCOMMANDS:
  put    <src> <tile_prefix> <step>      one solved tile (and its field-size
                                         report, if it is there) up
  get    <tile_prefix> <step> <x0> <y0> <spacing> <dest>
                                         the 3x3 neighbourhood down
  put_tree <src_dir> <prefix>            every file under a local directory up,
                                         keeping relative paths (a mosaic
                                         step's products: 200 km tiles,
                                         mosaics, netCDFs)
  get_glob <prefix> <pattern> <dest>     the files matching <prefix>/<pattern>
                                         down (the mosaics an nc job reads);
                                         --require fails if there are none

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


def put_tree(fs, src_dir, prefix):
    """upload every file under src_dir to prefix, keeping relative paths"""
    files = sorted(os.path.join(root, name) for root, _, names in os.walk(src_dir) for name in names)
    if not files:
        print(f's3_tiles: nothing to upload under {src_dir}', file=sys.stderr)
        return 1
    for src in files:
        dest = f'{prefix.rstrip("/")}/{os.path.relpath(src, src_dir)}'
        print(f's3_tiles: {src} -> {dest}')
        fs.put(src, dest)
    return 0


def get_glob(fs, prefix, pattern, dest, require=False):
    """download the files matching prefix/pattern into dest"""
    os.makedirs(dest, exist_ok=True)
    found = sorted(fs.glob(f'{prefix.rstrip("/")}/{pattern}'))
    for src in found:
        target = os.path.join(dest, os.path.basename(src))
        print(f's3_tiles: {src} -> {target}')
        fs.get(src, target)
    if not found:
        print(f's3_tiles: no {prefix.rstrip("/")}/{pattern}',
              file=sys.stderr if require else sys.stdout)
        return 1 if require else 0
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

    p_put_tree = sub.add_parser('put_tree')
    p_put_tree.add_argument('src_dir')
    p_put_tree.add_argument('prefix')

    p_get_glob = sub.add_parser('get_glob')
    p_get_glob.add_argument('prefix')
    p_get_glob.add_argument('pattern')
    p_get_glob.add_argument('dest')
    p_get_glob.add_argument('--require', action='store_true')

    args = parser.parse_args()
    fs = s3fs.S3FileSystem()
    if args.cmd == 'put':
        return put(fs, args.src, args.tile_prefix, args.step)
    if args.cmd == 'put_tree':
        return put_tree(fs, args.src_dir, args.prefix)
    if args.cmd == 'get_glob':
        return get_glob(fs, args.prefix, args.pattern, args.dest, require=args.require)
    return get(fs, args.tile_prefix, args.step, args.x0, args.y0,
               args.spacing, args.dest)


if __name__ == '__main__':
    sys.exit(main())
