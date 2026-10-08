"""
make_200km_tiles.py: the canonical 200 km tile list, filtered by partition.

Ben, 2026-09-19: both Antarctic partitions use ATL1415/resources/AA/
200km_tile_list.txt as the canonical listing of possible 200 km tiles, and
each decides by its own xy limits which to include.  Those limits are the ones
its 40 km tiles use (60 km half --min_xy 360000, 44 km half --max_xy 440000),
and they must reproduce setup_AA_sectors.py's split, which draws a 200 km tile
from the 44 km half when |x| and |y| of its center are both < 400 km.
"""
import importlib.util
import os

import pytest

HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    'make_200km_tiles', os.path.join(HERE, '..', 'ATL1415', 'scripts', 'make_200km_tiles.py'))
m2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m2)

AA_LIST = os.path.join(HERE, '..', 'ATL1415', 'resources', 'AA', '200km_tile_list.txt')
NEAR_POLE_RADIUS = 4.e5          # setup_AA_sectors.py's default


def key(xyc):
    return {tuple(xy) for xy in xyc}


# --- the real Antarctic list ---------------------------------------------------

def test_the_real_list_reads_as_422_centers():
    # 413, plus the 9 cells its tiles reach into but hold no tile center of
    # (plan_200km_footprint.sh, 2026-10-08)
    assert len(m2.read_200km_tile_list(AA_LIST)) == 422


def test_the_real_list_covers_every_cell_its_tiles_reach():
    import re
    from ATL1415.mosaic_groups import footprint_centers_200km
    names = [ln.strip() for ln in open(os.path.join(HERE, '..', 'ATL1415', 'resources', 'AA',
                                                      '40km_tile_list.txt')) if ln.strip()]
    ext = lambda n: max(abs(int(v)) for v in re.search(r'E(-?\d+)_N(-?\d+)', n).groups()) * 1000
    need = key(footprint_centers_200km([n for n in names if ext(n) >= 360000], half_width=30e3)) \
        | key(footprint_centers_200km([n for n in names if ext(n) <= 440000], half_width=22e3))
    assert need <= key(m2.read_200km_tile_list(AA_LIST))


def test_each_half_takes_its_own_tiles_and_together_they_take_all():
    xyc = m2.read_200km_tile_list(AA_LIST)
    north = key(m2.select_200km_tiles(xyc, min_xy=360000))
    south = key(m2.select_200km_tiles(xyc, max_xy=440000))
    assert (len(north), len(south)) == (406, 16)
    assert not north & south                      # no tile built twice
    assert north | south == key(xyc)              # no tile missed


def test_the_halves_split_where_setup_AA_sectors_does():
    # a sector draws a 200 km tile from the 44 km half exactly when it is
    # near the pole by setup_AA_sectors' rule -- so each half builds exactly
    # the tiles it will be asked for
    xyc = m2.read_200km_tile_list(AA_LIST)
    near_pole = key(xy for xy in xyc
                    if abs(xy[0]) < NEAR_POLE_RADIUS and abs(xy[1]) < NEAR_POLE_RADIUS)
    assert key(m2.select_200km_tiles(xyc, max_xy=440000)) == near_pole
    assert key(m2.select_200km_tiles(xyc, min_xy=360000)) == key(xyc) - near_pole


def test_no_limits_keeps_everything():
    # the monthly run is one partition, with no limits
    xyc = m2.read_200km_tile_list(AA_LIST)
    assert m2.select_200km_tiles(xyc) == xyc


# --- the limits' meaning, as make_ATL1415_queue.py's -----------------------------

@pytest.mark.parametrize('xy, north, south', [
    ([100000., 300000.], False, True),
    ([300000., -300000.], False, True),
    ([500000., 100000.], True, False),     # one coordinate past the line is enough
    ([-100000., -2700000.], True, False),
])
def test_limits_follow_max_abs_xy(xy, north, south):
    assert bool(m2.select_200km_tiles([xy], min_xy=360000)) == north
    assert bool(m2.select_200km_tiles([xy], max_xy=440000)) == south


# --- the canonical list overrides the per-directory cache ----------------------

def test_a_given_list_is_used_and_the_region_cache_is_left_alone(tmp_path):
    region = tmp_path / 'AA'
    region.mkdir()
    cache = region / '200km_tile_list.txt'
    cache.write_text('900000.0 900000.0\n')        # a stale cache
    canon = tmp_path / 'canon.txt'
    canon.write_text('100000.0 -100000.0\n-2700000.0 1500000.0\n')
    xyc = m2.make_200km_tiles(str(region), tile_list_file=str(canon))
    assert xyc == [[100000., -100000.], [-2700000., 1500000.]]
    assert cache.read_text() == '900000.0 900000.0\n'


def test_a_given_list_writes_no_cache(tmp_path):
    region = tmp_path / 'AA'
    region.mkdir()
    canon = tmp_path / 'canon.txt'
    canon.write_text('100000.0 -100000.0\n')
    m2.make_200km_tiles(str(region), tile_list_file=str(canon))
    assert not (region / '200km_tile_list.txt').exists()


def test_a_bad_line_is_an_error_naming_it(tmp_path):
    canon = tmp_path / 'canon.txt'
    canon.write_text('100000.0 -100000.0\nfield_sizes\n')
    with pytest.raises(ValueError, match=':2:'):
        m2.read_200km_tile_list(str(canon))


# with no list: every cell the tiles' 60 km squares reach, not only those
# holding a tile center (E420_N20 spans x 390..450, y -10..50 km: four cells;
# E-620_N-1340 spans x -650..-590 km: two)
FOOTPRINT = {(300000., 100000.), (300000., -100000.), (500000., 100000.), (500000., -100000.),
             (-700000., -1300000.), (-500000., -1300000.)}


def test_without_a_list_the_cells_come_from_the_prelim_footprints(tmp_path):
    region = tmp_path / 'XX'                       # a region with no resource list
    (region / 'prelim').mkdir(parents=True)
    for name in ['E420_N20.h5', 'E460_N20.h5', 'E-620_N-1340.h5']:
        (region / 'prelim' / name).write_bytes(b'')
    xyc = key(m2.make_200km_tiles(str(region)))
    assert xyc == FOOTPRINT
    assert (region / '200km_tile_list.txt').exists()


@pytest.mark.parametrize('region', ['GL', 'AA'])
def test_the_region_resource_list_is_used_by_default(tmp_path, region):
    from ATL1415.mosaic_groups import read_200km_centers
    xyc = m2.make_200km_tiles(str(tmp_path), region=region)
    assert xyc == read_200km_centers(region)
    assert not (tmp_path / '200km_tile_list.txt').exists()


# --- DPS: tiles read in place, one 200 km tile per job (plan_dps_mosaic D3a-1) --

import sys

import fsspec
import pointCollection as pc

PRELIM = ['E420_N20.h5', 'E460_N20.h5', 'E-620_N-1340.h5']


@pytest.fixture
def remote_region(monkeypatch):
    """the prelim tile NAMES under a memory:// region prefix (listing only)"""
    fs = fsspec.filesystem('memory')
    monkeypatch.setattr(pc.io_utils, 'get_s3fs', lambda daac=None, **kw: fs)
    prefix = 'memory://bucket/rel006/south/AA'
    for name in PRELIM:
        fs.pipe(f'{prefix}/prelim/{name}', b'')
    yield prefix
    fs.rm('memory://bucket', recursive=True)


def test_centers_from_a_remote_tiles_base(tmp_path, remote_region):
    region = tmp_path / 'AA'
    region.mkdir()
    xyc = key(m2.make_200km_tiles(str(region), tiles_base=remote_region))
    assert xyc == FOOTPRINT
    # the cache is written locally, beside the outputs
    assert (region / '200km_tile_list.txt').exists()


def run_main(monkeypatch, tmp_path, argv):
    # the slurm file is not what these test, and its template lookup
    # (importlib.resources on the ATL1415 package) fails once
    # test_setup_region.py has put a stub ATL1415 in sys.modules
    monkeypatch.setattr(m2.ATL1415, 'make_slurm_file', lambda *args, **kwargs: None)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, 'argv', ['make_200km_tiles.py'] + argv)
    m2.main()
    return sorted((tmp_path / 'tile_run_AA' / 'queue').iterdir())


def test_one_center_one_task_reading_from_tiles_base(tmp_path, monkeypatch, remote_region):
    region = tmp_path / 'AA'
    region.mkdir()
    tasks = run_main(monkeypatch, tmp_path, [str(region), 'AA', '--dzdt_lags', '1,4',
                                             '--tiles_base', remote_region, '--center', '500000', '100000'])
    assert [t.name for t in tasks] == ['task_1']
    lines = [ln for ln in tasks[0].read_text().splitlines() if ln.startswith('make_mosaic.py')]
    assert lines and all(f'-d {remote_region} ' in ln for ln in lines)
    # outputs stay local, under region_dir/200km_tiles, named by the tile's bounds
    assert all(f'-O {region}/200km_tiles/' in ln for ln in lines)
    assert all(ln.split(' -O ')[1].split()[0].endswith('400_600_0_200.h5') for ln in lines)
    # the search window is the 200 km square plus 10 km
    assert all('-r 390000.0 610000.0 -10000.0 210000.0' in ln for ln in lines)


def test_a_center_not_in_the_region_is_refused(tmp_path, monkeypatch, remote_region):
    region = tmp_path / 'AA'
    region.mkdir()
    # the check is against the region's list (AA: the resource), not the tiles
    with pytest.raises(SystemExit, match='not one of this region'):
        run_main(monkeypatch, tmp_path, [str(region), 'AA', '--dzdt_lags', '1',
                                         '--tiles_base', remote_region, '--center', '9100000', '9100000'])
