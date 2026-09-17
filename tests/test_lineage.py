"""
ATL14/15 lineage from the tiles' meta/input_files, with no granule opened.

TEMPORARY (docs/plan_IS_run.sh I9g2): attributes that need the granule are
'NOT_SET' -- invalid -- until the prelim step records them in the tiles.
"""
from types import SimpleNamespace

import h5py
import pytest

from ATL1415.ATL1415_attrs_meta import (FILE_ONLY_LINEAGE_ATTRS,
                                        attributes_for_ATL11_file, set_lineage)

AT = 'ATL11_023003_0331_007_04.h5'
XO = 'ATL11XO_AR_E1200_N-2600_c01_007_03.h5'


def test_along_track_name(tmp_path, monkeypatch):
    # a relative name that exists nowhere: nothing may be opened
    monkeypatch.chdir(tmp_path)
    fa, this_format = attributes_for_ATL11_file(AT)
    assert this_format == 'along-track'
    assert (fa['fileName'], fa['shortName'], fa['start_rgt'], fa['end_rgt'],
            fa['start_region'], fa['end_region'], fa['start_cycle'],
            fa['end_cycle'], fa['release'], fa['version']) == \
        (AT, 'ATL11', '0230', '0230', '03', '03', '03', '31', '007', '04')
    for attr in FILE_ONLY_LINEAGE_ATTRS['along-track']:
        assert fa[attr] == 'NOT_SET'


def test_xover_name(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fa, this_format = attributes_for_ATL11_file(XO)
    assert this_format == 'xo'
    assert (fa['fileName'], fa['shortName'], fa['start_cycle'], fa['end_cycle'],
            fa['release'], fa['version']) == (XO, 'ATL11XO', '01', '01', '007', '03')
    for attr in FILE_ONLY_LINEAGE_ATTRS['xo'] + ['start_region', 'end_region',
                                                 'start_orbit', 'end_orbit']:
        assert fa[attr] == 'NOT_SET'


def test_unrecognized_name_raises():
    with pytest.raises(ValueError, match='not_an_atl11.h5'):
        attributes_for_ATL11_file('not_an_atl11.h5')


class FakeGroup:
    def __init__(self):
        self.attrs = {}

    def setncattr(self, name, value):
        self.attrs[name] = value


def write_tile(path, input_files):
    with h5py.File(path, 'w') as h5f:
        h5f.create_group('meta').attrs['input_files'] = input_files.encode('ascii')


def run_set_lineage(tiles_dir):
    group = FakeGroup()
    set_lineage({'METADATA/Lineage/ATL11': group}, {},
                SimpleNamespace(tiles_dir=str(tiles_dir)))
    return group.attrs


def test_set_lineage_dedupes_skips_empty_and_warns(tmp_path, capsys):
    write_tile(tmp_path / 'E1.h5', ','.join([AT, AT, AT, XO]))
    write_tile(tmp_path / 'E2.h5', ','.join([XO, AT]))
    write_tile(tmp_path / 'E3.h5', '')          # a matched tile
    attrs = run_set_lineage(tmp_path)
    assert attrs['fileName'] == sorted([AT, XO])
    assert attrs['uuid'] == ['NOT_SET', 'NOT_SET']
    out = capsys.readouterr().out
    assert 'INVALID for 1 along-track files' in out
    assert 'INVALID for 1 xo files' in out


def test_set_lineage_names_the_tile_for_a_bad_name(tmp_path):
    write_tile(tmp_path / 'E9.h5', ','.join([AT, 'junk.h5']))
    with pytest.raises(ValueError, match=r'junk\.h5.*E9\.h5'):
        run_set_lineage(tmp_path)
