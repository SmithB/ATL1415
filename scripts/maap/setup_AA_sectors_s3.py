#!/usr/bin/env python3
"""
The four Antarctic sectors (A1-A4) for the DPS mosaic and netCDF steps, from
the tiles on S3 (plan_AA_mosaic_dps.sh M1).  The discover counterpart is
setup_AA_sectors.py, which symlinks both halves' tiles into each sector; S3
has no links, so here each sector gets small files instead, and the 200 km
jobs write straight into the sector they belong to:

  <south_dir>/A<n>/bounds.txt               the sector's crop (setup_AA_sectors' rule)
  <south_dir>/A<n>/input_args_A<n>.txt      AA's args, -b and --region -> A<n> (its rule)
  <south_dir>/A<n>/prelim_tiles.txt         the sector's prelim tiles, one URI per line:
                                            the nc writers' --tiles_dir
  <south_dir>/A<n>/mosaic200_tasks_60km.txt the 200 km cells this sector takes from
  <south_dir>/A<n>/mosaic200_tasks_44km.txt   each half, '<x>_<y>' (submit --task)

THE RULES, each from where it already lives:
  - the 200 km cells: ATL1415/resources/AA/200km_tile_list.txt (422);
  - which half makes a cell: make_200km_tiles.select_200km_tiles with the
    halves' own limits (60 km --min_xy 360000, 44 km --max_xy 440000), which
    split the cells at 400 km, setup_AA_sectors' near_pole_radius;
  - a cell's sector: setup_AA_sectors.get_quadrant of its center;
  - a tile's cell (Ben 2026-10-18, QA1: "200km-ownership rule"):
    mosaic_groups.centers_200km, the rule tile_meta records by.  So a sector's
    tiles are, for each of its cells, the tiles of that cell's half whose
    center the cell holds -- exactly what the 200 km jobs' tile_meta records
    hold, which the nc writers check.
GATES (all, or nothing is written): the halves' cells partition the 422; the
manifests list each center that has a prelim tile in either half exactly
once; every tile's cell is a canonical cell.

The files are written locally; upload them (aws s3 cp) to <tiles root>/A<n>/
(bounds.txt, prelim_tiles.txt) and run_args/.../A<n>/ (input_args_A<n>.txt).

usage: setup_AA_sectors_s3.py <south_dir> --tiles_60km s3://.../rel006/south/AA
                              --tiles_44km s3://.../rel006/south/AA_44km [--dry-run]
  <south_dir> holds AA/input_args_AA.txt (the 60 km half's args, == S3's).
"""
import argparse
import os
import sys

import numpy as np

from ATL1415.mosaic_groups import centers_200km, read_200km_centers
from ATL1415.scripts.make_200km_tiles import select_200km_tiles
from ATL1415.scripts.setup_AA_sectors import (N_SECTORS, get_quadrant, write_sector_bounds,
                                               write_sector_input_args)

HALVES = {'60km': dict(min_xy=360000.), '44km': dict(max_xy=440000.)}   # howto_MAAP_AA.sh step 3
BOUNDS_PAD = 4.e4      # setup_AA_sectors.py --bounds_pad default


def cell_key(c):
    return (float(c[0]), float(c[1]))


def task_name(c):
    return f'{int(c[0])}_{int(c[1])}'


def plan(names, prefixes):
    """
    names {half: [tile file names]}, prefixes {half: tile prefix} ->
    {sector: {'tiles': [uri, ...], 'tasks': {half: [cell, ...]}}}, raising on any gate.
    """
    centers = [cell_key(c) for c in read_200km_centers('AA')]
    cells = {h: {cell_key(c) for c in select_200km_tiles(centers, **lim)} for h, lim in HALVES.items()}
    if cells['60km'] & cells['44km'] or (cells['60km'] | cells['44km']) != set(centers):
        raise RuntimeError('the halves\' 200 km cells do not partition the canonical list:'
                           f' {len(cells["60km"])} + {len(cells["44km"])} of {len(centers)},'
                           f' {len(cells["60km"] & cells["44km"])} in both')
    half_of = {c: h for h, cs in cells.items() for c in cs}
    sector_of = {c: int(get_quadrant(np.array(c))[0]) for c in centers}
    out = {n: {'tiles': [], 'tasks': {h: [] for h in HALVES}} for n in range(1, N_SECTORS + 1)}
    for c in centers:
        out[sector_of[c]]['tasks'][half_of[c]].append(c)
    stray = []
    for h in HALVES:
        for name in names[h]:
            owner = centers_200km([name])
            c = cell_key(owner[0]) if owner else None
            if c not in half_of:
                stray.append(f'{h}/{name}')
            elif half_of[c] == h:
                out[sector_of[c]]['tiles'].append(f'{prefixes[h].rstrip("/")}/prelim/{name}')
    if stray:
        raise RuntimeError(f'{len(stray)} tiles are in no canonical 200 km cell: {stray[:10]}')
    listed = [os.path.basename(t) for s in out.values() for t in s['tiles']]
    expected = set(names['60km']) | set(names['44km'])
    if len(listed) != len(set(listed)) or set(listed) != expected:
        raise RuntimeError(f'the sectors list {len(listed)} tiles ({len(set(listed))} names) but the'
                           f' halves have {len(expected)} centers: missing'
                           f' {sorted(expected - set(listed))[:10]}')
    for s in out.values():
        s['tiles'].sort()
    return out


def full_bounds(names_60km):
    """setup_AA_sectors.compute_full_bounds over the 60 km half's prelim tiles (by name)"""
    xy = np.array([[int(v) for v in os.path.splitext(n)[0][1:].split('_N')] for n in names_60km])
    pad_km = BOUNDS_PAD / 1000
    return [1000 * np.array(jj) for jj in ([xy[:, 0].min() - pad_km, xy[:, 0].max() + pad_km],
                                           [xy[:, 1].min() - pad_km, xy[:, 1].max() + pad_km])]


def write(south_dir, sectors, names_60km):
    for n, s in sectors.items():
        d = os.path.join(south_dir, f'A{n}')
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, 'prelim_tiles.txt'), 'w') as fh:
            fh.writelines(t + '\n' for t in s['tiles'])
        for h, cs in s['tasks'].items():
            with open(os.path.join(d, f'mosaic200_tasks_{h}.txt'), 'w') as fh:
                fh.writelines(task_name(c) + '\n' for c in cs)
    write_sector_bounds(south_dir, full_bounds(names_60km))
    write_sector_input_args(south_dir, 'AA')


def list_prelim(prefix):
    import pointCollection as pc
    found = pc.io_utils.glob_remote(f'{prefix.rstrip("/")}/prelim/E*.h5', fs=pc.io_utils.get_s3fs(daac=None))
    return [os.path.basename(f) for f in found]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('south_dir')
    p.add_argument('--tiles_60km', required=True)
    p.add_argument('--tiles_44km', required=True)
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args(argv)
    if not os.path.isfile(os.path.join(a.south_dir, 'AA', 'input_args_AA.txt')):
        sys.exit(f'{a.south_dir}/AA/input_args_AA.txt not found: the sectors\' args are made from it')
    prefixes = {'60km': a.tiles_60km, '44km': a.tiles_44km}
    names = {h: list_prelim(pfx) for h, pfx in prefixes.items()}
    for h in HALVES:
        if not names[h]:
            sys.exit(f'no prelim tiles under {prefixes[h]}/prelim/')
    sectors = plan(names, prefixes)
    for n, s in sectors.items():
        print(f'A{n}: {len(s["tiles"])} tiles; 200 km cells: '
              + ', '.join(f'{len(cs)} from {h}' for h, cs in s['tasks'].items()))
    print(f'{sum(len(s["tiles"]) for s in sectors.values())} tiles in all'
          f' ({len(names["60km"])} 60 km + {len(names["44km"])} 44 km prelim tiles listed)')
    if a.dry_run:
        print('--dry-run: nothing written')
        return 0
    write(a.south_dir, sectors, names['60km'])
    print(f'written under {a.south_dir}/A1..A{N_SECTORS}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
