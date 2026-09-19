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

def test_the_real_list_reads_as_413_centers():
    assert len(m2.read_200km_tile_list(AA_LIST)) == 413


def test_each_half_takes_its_own_tiles_and_together_they_take_all():
    xyc = m2.read_200km_tile_list(AA_LIST)
    north = key(m2.select_200km_tiles(xyc, min_xy=360000))
    south = key(m2.select_200km_tiles(xyc, max_xy=440000))
    assert (len(north), len(south)) == (397, 16)
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


def test_without_a_list_the_centers_are_still_derived_from_prelim(tmp_path):
    # discover's behaviour, unchanged: 200 km cells covering the prelim tiles
    region = tmp_path / 'AA'
    (region / 'prelim').mkdir(parents=True)
    for name in ['E420_N20.h5', 'E460_N20.h5', 'E-620_N-1340.h5']:
        (region / 'prelim' / name).write_bytes(b'')
    xyc = key(m2.make_200km_tiles(str(region)))
    assert xyc == {(500000., 100000.), (-700000., -1300000.)}
    assert (region / '200km_tile_list.txt').exists()
