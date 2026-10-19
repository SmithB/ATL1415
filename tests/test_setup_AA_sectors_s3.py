"""
scripts/maap/setup_AA_sectors_s3.py (plan_AA_mosaic_dps.sh M1): the Antarctic
sectors' manifests, 200 km task lists, bounds and args, from the halves' tile
names.  No S3: the names are made from the canonical 40 km list.
"""
import importlib.util
import os
import re

import numpy as np
import pytest

from ATL1415.scripts import setup_AA_sectors

HERE = os.path.dirname(__file__)
REPO = os.path.join(HERE, '..')
_spec = importlib.util.spec_from_file_location(
    'setup_AA_sectors_s3', os.path.join(REPO, 'scripts', 'maap', 'setup_AA_sectors_s3.py'))
s3s = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s3s)

PREFIXES = {'60km': 's3://b/rel006/south/AA', '44km': 's3://b/rel006/south/AA_44km'}


def xy_km(name):
    return tuple(int(v) for v in re.match(r'E(-?\d+)_N(-?\d+)', name).groups())


@pytest.fixture(scope='module')
def names():
    """every AA center in the half (or halves) howto_MAAP_AA.sh step 3 puts it in"""
    cs = [line.strip() for line in open(os.path.join(REPO, 'ATL1415', 'resources', 'AA', '40km_tile_list.txt'))
          if line.strip()]
    ext = lambda n: max(abs(v) for v in xy_km(n))
    return {'60km': [n for n in cs if ext(n) >= 360], '44km': [n for n in cs if ext(n) <= 440]}


def test_the_sectors_partition_the_centers_by_200km_ownership(names):
    sectors = s3s.plan(names, PREFIXES)
    tiles = [t for s in sectors.values() for t in s['tiles']]
    assert len(tiles) == 8944 and len({os.path.basename(t) for t in tiles}) == 8944
    cells = {n: {h: len(c) for h, c in s['tasks'].items()} for n, s in sectors.items()}
    assert cells == {1: {'60km': 127, '44km': 4}, 2: {'60km': 90, '44km': 4},
                     3: {'60km': 61, '44km': 4}, 4: {'60km': 128, '44km': 4}}
    half = {os.path.basename(t): ('44km' if '/AA_44km/' in t else '60km') for t in tiles}
    # QA1: the 39 centers at x or y = -400 km inside the near-pole box are the 44 km half's
    # (their cell is [-400, -200)); +400 belongs to the 60 km half's cells
    edge = [n for n in half if -400 in xy_km(n) and all(-400 <= v < 400 for v in xy_km(n))]
    assert len(edge) == 39 and {half[n] for n in edge} == {'44km'}
    assert half['E400_N0.h5'] == '60km' and half['E-400_N0.h5'] == '44km'
    assert half['E360_N360.h5'] == '44km' and half['E0_N-440.h5'] == '60km'
    # each tile is in the sector of its cell
    for n, s in sectors.items():
        for t in s['tiles']:
            x, y = xy_km(os.path.basename(t))
            c = np.floor(np.array([x, y]) / 200.) * 200 + 100
            assert setup_AA_sectors.get_quadrant(c)[0] == n


def test_a_missing_owner_tile_is_an_error(names):
    # E-400_N0's cell takes the 44 km tile: without it, the 60 km one must not stand in
    broken = dict(names, **{'44km': [n for n in names['44km'] if n != 'E-400_N0.h5']})
    with pytest.raises(RuntimeError, match='missing'):
        s3s.plan(broken, PREFIXES)


def test_a_tile_in_no_canonical_cell_is_an_error(names):
    with pytest.raises(RuntimeError, match='no canonical 200 km cell'):
        s3s.plan(dict(names, **{'60km': names['60km'] + ['E9000_N9000.h5']}), PREFIXES)


def test_bounds_and_args_are_setup_AA_sectors_own(names, tmp_path):
    # the discover script's bounds, from a directory of (empty) 60 km prelim tiles
    top = tmp_path / 'discover'
    (top / 'AA' / 'prelim').mkdir(parents=True)
    for n in names['60km']:
        (top / 'AA' / 'prelim' / n).touch()
    (top / 'AA' / 'input_args_AA.txt').write_text('-b=/x/rel006/south/AA\n--region=AA\n-W=60000\n')
    for q in range(1, 5):
        (top / f'A{q}').mkdir()
    setup_AA_sectors.write_sector_bounds(str(top), setup_AA_sectors.compute_full_bounds(str(top / 'AA'), 4.e4))
    setup_AA_sectors.write_sector_input_args(str(top), 'AA')

    mine = tmp_path / 'maap'
    (mine / 'AA').mkdir(parents=True)
    (mine / 'AA' / 'input_args_AA.txt').write_text((top / 'AA' / 'input_args_AA.txt').read_text())
    s3s.write(str(mine), s3s.plan(names, PREFIXES), names['60km'])
    for q in range(1, 5):
        for f in ('bounds.txt', f'input_args_A{q}.txt'):
            assert (mine / f'A{q}' / f).read_text() == (top / f'A{q}' / f).read_text(), (q, f)
    assert (mine / 'A1' / 'input_args_A1.txt').read_text().startswith('-b=/x/rel006/south/A1\n--region=A1\n')
    tasks = (mine / 'A3' / 'mosaic200_tasks_44km.txt').read_text().split()
    assert sorted(tasks) == ['-100000_-100000', '-100000_-300000', '-300000_-100000', '-300000_-300000']
    listed = (mine / 'A2' / 'prelim_tiles.txt').read_text().split()
    assert listed == sorted(listed) and all(t.startswith('s3://b/rel006/south/AA') for t in listed)


def test_one_partition_takes_every_cell_and_tile(names, tmp_path):
    # the monthly run (plan_AA_monthly_dps.sh AM1): 60 km everywhere, as setup_AA_sectors.py
    # --near_pole_radius 0
    every = sorted(set(names['60km']) | set(names['44km']))
    one = {'60km': every}
    sectors = s3s.plan(one, {'60km': PREFIXES['60km']}, s3s.ONE_PARTITION)
    cells = {n: {h: len(c) for h, c in s['tasks'].items()} for n, s in sectors.items()}
    assert cells == {1: {'60km': 131}, 2: {'60km': 94}, 3: {'60km': 65}, 4: {'60km': 132}}
    tiles = [t for s in sectors.values() for t in s['tiles']]
    assert len(tiles) == 8944 and all(t.startswith(PREFIXES['60km'] + '/prelim/') for t in tiles)
    for n, s in sectors.items():
        for t in s['tiles']:
            c = np.floor(np.array(xy_km(os.path.basename(t))) / 200.) * 200 + 100
            assert setup_AA_sectors.get_quadrant(c)[0] == n
    (tmp_path / 'AA').mkdir()
    (tmp_path / 'AA' / 'input_args_AA.txt').write_text('-b=/x/rel006/south_monthly/AA\n--region=AA\n')
    s3s.write(str(tmp_path), sectors, every)
    assert sorted(p.name for p in (tmp_path / 'A1').iterdir()) == [
        'bounds.txt', 'input_args_A1.txt', 'mosaic200_tasks_60km.txt', 'prelim_tiles.txt']

