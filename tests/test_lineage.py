"""
ATL14/15 lineage: read at solve time, stored in the tile, used by the netCDF step.

docs/plan_lineage_at_solve_time.sh.  The solve records each granule's
file-only attributes in meta/lineage/<granule>; the netCDF step reads them
from the tiles and never opens ATL11, marking as 'NOT_SET' whatever the tiles
do not carry.  No network: granules and tiles here are synthetic.
"""
import importlib.util
import os
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from ATL1415.ATL1415_attrs_meta import (FILE_ONLY_LINEAGE_ATTRS,
                                        as_lineage_text,
                                        attributes_for_ATL11_file, set_lineage)
from ATL1415.ATL11_to_ATL15 import write_lineage

# read_ATL11.py directly, as tests/test_read_ATL11.py does, to keep LSsurf out
_spec = importlib.util.spec_from_file_location(
    'read_ATL11', os.path.join(os.path.dirname(__file__), '..', 'ATL1415', 'read_ATL11.py'))
read_ATL11_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(read_ATL11_mod)
lineage_attributes = read_ATL11_mod.lineage_attributes
lineage_for_granules = read_ATL11_mod.lineage_for_granules

AT = 'ATL11_023003_0332_007_05.h5'
XO = 'ATL11XO_AR_E1200_N-2600_c01_007_03.h5'
# the values the real granules carry (probed 2026-09-17, plan Background)
AT_ATTRS = {'uuid': 'a805a367-98cb-4a88-b3d4-1ac9fd969f8f',
            'start_geoseg': 354552, 'end_geoseg': 446050,
            'start_orbit': 3205, 'end_orbit': 43428}
XO_ATTRS = {'uuid': '11a11d42-62df-4034-a81c-51c228b9cfac',
            'start_geoseg': 352080, 'end_geoseg': 649385,
            'start_rgt': 238, 'end_rgt': 1381}


def make_granule(path, attrs, uuid_as_bytes=True):
    """A file shaped like an ATL11 granule, holding only what lineage reads."""
    with h5py.File(path, 'w') as h5f:
        if 'uuid' in attrs:
            uuid = attrs['uuid']
            h5f.require_group('METADATA/DatasetIdentification').attrs['uuid'] = (
                uuid.encode('utf-8') if uuid_as_bytes else uuid)
        anc = h5f.require_group('ancillary_data')
        for key, value in attrs.items():
            if key != 'uuid':
                anc.create_dataset(key, data=np.array([value], dtype='int32'))
    return str(path)


# ---------------------------------------------------------------- L2, reading

def test_along_track_granule_attributes(tmp_path):
    path = make_granule(tmp_path / AT, AT_ATTRS)
    with h5py.File(path, 'r') as h5f:
        attrs = lineage_attributes(h5f)
    assert attrs == AT_ATTRS
    assert isinstance(attrs['uuid'], str)


def test_crossover_granule_has_rgts_and_no_orbit(tmp_path):
    path = make_granule(tmp_path / XO, XO_ATTRS)
    with h5py.File(path, 'r') as h5f:
        attrs = lineage_attributes(h5f)
    assert attrs == XO_ATTRS
    assert 'start_orbit' not in attrs and 'end_orbit' not in attrs


def test_uuid_as_str_and_missing_pieces_are_absent(tmp_path):
    path = make_granule(tmp_path / AT, {'uuid': 'u', 'start_geoseg': 1}, uuid_as_bytes=False)
    with h5py.File(path, 'r') as h5f:
        assert lineage_attributes(h5f) == {'uuid': 'u', 'start_geoseg': 1}
    # a granule with neither METADATA nor ancillary_data gives {}, not an error
    with h5py.File(tmp_path / 'empty.h5', 'w') as h5f:
        pass
    with h5py.File(tmp_path / 'empty.h5', 'r') as h5f:
        assert lineage_attributes(h5f) == {}


def test_pair_suffixes_and_duplicates_collapse(tmp_path):
    path = make_granule(tmp_path / AT, AT_ATTRS)
    lineage = lineage_for_granules([path + ':pair1', path + ':pair2', path, None])
    assert list(lineage) == [AT]
    assert lineage[AT] == AT_ATTRS


def test_an_unreadable_granule_is_reported_and_skipped(tmp_path, capsys):
    good = make_granule(tmp_path / AT, AT_ATTRS)
    lineage = lineage_for_granules([str(tmp_path / 'gone.h5'), good])
    assert list(lineage) == [AT]          # the solve goes on with what it has
    assert 'could not read lineage attributes' in capsys.readouterr().out


# ----------------------------------------------------------------- L3, writing

def test_write_lineage_puts_one_group_per_granule(tmp_path):
    # the round trip the solve makes: granule -> lineage_attributes -> tile
    granules = tmp_path / 'granules'
    granules.mkdir()
    make_granule(granules / AT, AT_ATTRS)
    make_granule(granules / XO, XO_ATTRS)
    lineage = lineage_for_granules([str(granules / AT), str(granules / XO)])
    with h5py.File(tmp_path / 'tile.h5', 'w') as h5f:
        write_lineage(h5f, lineage)
    with h5py.File(tmp_path / 'tile.h5', 'r') as h5f:
        assert sorted(h5f['meta/lineage']) == sorted([AT, XO])
        stored = dict(h5f[f'meta/lineage/{XO}'].attrs)
        assert stored['end_rgt'] == 1381 and stored['uuid'] == XO_ATTRS['uuid']
        assert 'start_orbit' not in stored
        # the tile keeps the granule's own types; only the netCDF forces strings
        assert stored['end_rgt'].dtype == np.int32
        assert as_lineage_text(stored['end_rgt']) == '1381'


@pytest.mark.parametrize('lineage', [None, {}])
def test_write_lineage_writes_nothing_without_lineage(tmp_path, lineage):
    with h5py.File(tmp_path / 'tile.h5', 'w') as h5f:
        write_lineage(h5f, lineage)
        assert 'meta' not in h5f


# ------------------------------------------------------- L4, the netCDF reader

def test_stored_attributes_fill_the_file_only_ones():
    fa, this_format = attributes_for_ATL11_file(
        AT, stored={k: as_lineage_text(v) for k, v in AT_ATTRS.items()})
    assert this_format == 'along-track'
    assert (fa['uuid'], fa['start_geoseg'], fa['end_orbit']) == \
        (AT_ATTRS['uuid'], '354552', '43428')
    # the name stays the authority for the rest
    assert (fa['shortName'], fa['start_rgt'], fa['end_rgt'], fa['start_region'],
            fa['start_cycle'], fa['end_cycle'], fa['release'], fa['version']) == \
        ('ATL11', '0230', '0230', '03', '03', '32', '007', '05')


def test_a_crossover_keeps_both_of_its_rgts():
    # the bug this fixes: end_rgt used to be overwritten with start_rgt
    fa, this_format = attributes_for_ATL11_file(
        XO, stored={k: as_lineage_text(v) for k, v in XO_ATTRS.items()})
    assert this_format == 'xo'
    assert (fa['start_rgt'], fa['end_rgt']) == ('238', '1381')


def test_without_stored_attributes_everything_file_only_is_invalid():
    for name, fmt in [(AT, 'along-track'), (XO, 'xo')]:
        fa, this_format = attributes_for_ATL11_file(name)
        assert this_format == fmt
        for attr in FILE_ONLY_LINEAGE_ATTRS[fmt]:
            assert fa[attr] == 'NOT_SET'


def test_unrecognized_name_raises():
    with pytest.raises(ValueError, match='not_an_atl11.h5'):
        attributes_for_ATL11_file('not_an_atl11.h5')


class FakeGroup:
    def __init__(self):
        self.attrs = {}

    def setncattr(self, name, value):
        self.attrs[name] = value


def write_tile(path, input_files, lineage=None):
    with h5py.File(path, 'w') as h5f:
        h5f.require_group('meta').attrs['input_files'] = input_files.encode('ascii')
        write_lineage(h5f, lineage)


def run_set_lineage(tiles_dir):
    group = FakeGroup()
    set_lineage({'METADATA/Lineage/ATL11': group}, {},
                SimpleNamespace(tiles_dir=str(tiles_dir)))
    return group.attrs


def test_set_lineage_uses_what_the_tiles_recorded(tmp_path, capsys):
    write_tile(tmp_path / 'E1.h5', ','.join([AT, AT, XO]), {AT: AT_ATTRS, XO: XO_ATTRS})
    write_tile(tmp_path / 'E2.h5', ','.join([XO, AT]), {AT: AT_ATTRS, XO: XO_ATTRS})
    write_tile(tmp_path / 'E3.h5', '')          # a matched tile
    attrs = run_set_lineage(tmp_path)
    assert attrs['fileName'] == sorted([AT, XO])
    order = sorted([AT, XO]).index(AT)
    assert attrs['uuid'][order] == AT_ATTRS['uuid']
    assert attrs['end_rgt'][sorted([AT, XO]).index(XO)] == '1381'
    assert 'NOT_SET' not in attrs['uuid']
    assert 'INVALID' not in capsys.readouterr().out


def test_every_attribute_is_a_string_valid_or_not(tmp_path, capsys):
    write_tile(tmp_path / 'E1.h5', AT, {AT: AT_ATTRS})
    all_valid = run_set_lineage(tmp_path)
    write_tile(tmp_path / 'E2.h5', XO)          # no stored attributes: invalid
    mixed = run_set_lineage(tmp_path)
    for attrs in (all_valid, mixed):
        for field, values in attrs.items():
            assert all(isinstance(value, str) for value in values), field


def test_a_tile_without_stored_lineage_still_warns(tmp_path, capsys):
    write_tile(tmp_path / 'E1.h5', ','.join([AT, XO]))
    attrs = run_set_lineage(tmp_path)
    assert set(attrs['uuid']) == {'NOT_SET'}
    out = capsys.readouterr().out
    assert 'INVALID for 1 along-track files' in out and 'INVALID for 1 xo files' in out


def test_partly_recorded_lineage_warns_only_about_what_is_missing(tmp_path, capsys):
    write_tile(tmp_path / 'E1.h5', AT, {AT: {'uuid': 'u', 'start_geoseg': 1}})
    attrs = run_set_lineage(tmp_path)
    assert attrs['uuid'] == ['u'] and attrs['end_geoseg'] == ['NOT_SET']
    out = capsys.readouterr().out
    assert 'INVALID for 1 along-track files' in out
    assert 'end_geoseg, end_orbit, start_orbit are NOT_SET' in out


def test_the_same_granule_must_agree_across_tiles(tmp_path):
    write_tile(tmp_path / 'E1.h5', AT, {AT: AT_ATTRS})
    write_tile(tmp_path / 'E2.h5', AT, {AT: dict(AT_ATTRS, uuid='another-uuid')})
    with pytest.raises(ValueError, match=r'E[12]\.h5.*E[12]\.h5'):
        run_set_lineage(tmp_path)


def test_set_lineage_names_the_tile_for_a_bad_name(tmp_path):
    write_tile(tmp_path / 'E9.h5', ','.join([AT, 'junk.h5']))
    with pytest.raises(ValueError, match=r'junk\.h5.*E9\.h5'):
        run_set_lineage(tmp_path)


def test_an_unreadable_tile_is_skipped(tmp_path, capsys):
    write_tile(tmp_path / 'E1.h5', AT, {AT: AT_ATTRS})
    (tmp_path / 'E2.h5').write_bytes(b'not hdf5')
    with h5py.File(tmp_path / 'E3.h5', 'w') as h5f:    # no meta/input_files
        h5f.require_group('meta')
    attrs = run_set_lineage(tmp_path)
    assert attrs['fileName'] == [AT]
    assert 'failed to open tile file' in capsys.readouterr().out
