"""
ATL1415.tile_meta: the 200 km jobs save each prelim tile's tile_stats and
lineage inputs; the netCDF writers read those instead of every tile
(plan_tile_meta.sh).  The products must not change: everything here compares
the meta path against reading the tiles.  No network.
"""
import json
import os
import time
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from ATL1415 import tile_meta
from test_lineage import AT, AT_ATTRS, XO, XO_ATTRS, FakeGroup
from test_tiles_on_s3 import bucket, tile_stats, write_stats_tile  # noqa: F401 (fixture)
from ATL1415.ATL11_to_ATL15 import write_lineage
from ATL1415.ATL1415_attrs_meta import set_lineage

NAMES = ['E0_N-1000.h5', 'E40_N-1000.h5', 'E0_N-960.h5', 'E200_N-1000.h5', 'E240_N-760.h5']


def make_tiles(d, rng):
    """prelim tiles carrying both what tile_stats and what lineage read"""
    d.mkdir(exist_ok=True)
    for n, name in enumerate(NAMES):
        write_stats_tile(d / name, rng)
        with h5py.File(d / name, 'a') as h5f:
            files = [AT, XO] if n % 2 else [AT, AT]
            h5f.require_group('meta').attrs['input_files'] = ','.join(files).encode('ascii')
            write_lineage(h5f, {AT: AT_ATTRS, XO: XO_ATTRS} if n != 3 else None)
    return d


def lineage(args):
    group = FakeGroup()
    set_lineage({'METADATA/Lineage/ATL11': group}, {}, args)
    return group.attrs


def write_by_center(tiles, meta):
    from ATL1415.mosaic_groups import centers_200km
    n = 0
    for c in centers_200km(os.listdir(tiles)):
        n += tile_meta.write_meta(str(tiles), str(meta / f'{c[0]:.0f}_{c[1]:.0f}.json'), center=c)
    return n


def test_centers_partition_the_tiles(tmp_path):
    tiles = make_tiles(tmp_path / 'prelim', np.random.default_rng(0))
    meta = tmp_path / 'meta'
    assert write_by_center(tiles, meta) == len(NAMES)
    names = [r['name'] for f in sorted(meta.iterdir())
             for r in json.loads(f.read_text())['tiles']]
    assert sorted(names) == sorted(NAMES)
    assert len(list(meta.iterdir())) == 3     # (0..200) x (-1000..-800), (200..400) x2


def test_records_from_meta_equal_records_from_tiles(tmp_path):
    tiles = make_tiles(tmp_path / 'prelim', np.random.default_rng(1))
    write_by_center(tiles, tmp_path / 'meta')
    direct = tile_meta.tile_records(SimpleNamespace(tiles_dir=str(tiles)))
    saved = tile_meta.tile_records(SimpleNamespace(tiles_dir=str(tiles), tile_meta_dir=str(tmp_path / 'meta')))
    assert [p for p, _ in direct] == [p for p, _ in saved]
    for (_, a), (_, b) in zip(direct, saved):
        b = {k: v for k, v in b.items() if k != 'stamp'}
        assert a == b


def test_tile_stats_group_is_the_same_from_meta(tmp_path):
    tiles = make_tiles(tmp_path / 'prelim', np.random.default_rng(2))
    write_by_center(tiles, tmp_path / 'meta')
    (m1, d1) = tile_stats(tiles, str(tmp_path / 'a.nc'))
    tile_meta._CACHE.clear()
    from netCDF4 import Dataset
    from ATL1415.make_tile_stats_group import make_tile_stats_group
    with Dataset(tmp_path / 'b.nc', 'w') as nc:
        make_tile_stats_group(nc, SimpleNamespace(tiles_dir=str(tiles), region='GL',
                                                  tile_meta_dir=str(tmp_path / 'meta')))
    with Dataset(tmp_path / 'b.nc') as nc:
        grp = nc['tile_stats']
        for name in grp.variables:
            assert np.array_equal(np.ma.getmaskarray(grp[name][:]), m1[name]), name
            assert np.array_equal(np.ma.getdata(grp[name][:]), d1[name]), name


def test_lineage_is_the_same_from_meta(tmp_path, capsys):
    tiles = make_tiles(tmp_path / 'prelim', np.random.default_rng(3))
    write_by_center(tiles, tmp_path / 'meta')
    direct = lineage(SimpleNamespace(tiles_dir=str(tiles)))
    saved = lineage(SimpleNamespace(tiles_dir=str(tiles), tile_meta_dir=str(tmp_path / 'meta')))
    assert direct == saved and direct['fileName'] == sorted([AT, XO])


def test_unreadable_and_odd_tiles_behave_as_before(tmp_path, capsys):
    tiles = make_tiles(tmp_path / 'prelim', np.random.default_rng(4))
    (tiles / 'notes.h5').write_bytes(b'not hdf5')     # no E/N name: skipped by both readers
    tile_meta.write_meta(str(tiles), str(tmp_path / 'meta' / 'all.json'))
    recs = dict((os.path.basename(p), r) for p, r in tile_meta.tile_records(
        SimpleNamespace(tiles_dir=str(tiles), tile_meta_dir=str(tmp_path / 'meta'))))
    assert recs['notes.h5']['stats'] is None and recs['notes.h5']['lineage'] is None
    lineage(SimpleNamespace(tiles_dir=str(tiles), tile_meta_dir=str(tmp_path / 'meta')))
    assert 'failed to open tile file' in capsys.readouterr().out


@pytest.mark.parametrize('change', ['missing', 'extra', 'stale', 'duplicate', 'empty'])
def test_meta_that_does_not_match_the_tiles_stops(tmp_path, change):
    tiles = make_tiles(tmp_path / 'prelim', np.random.default_rng(5))
    meta = tmp_path / 'meta'
    write_by_center(tiles, meta)
    if change == 'missing':
        write_stats_tile(tiles / 'E600_N-1000.h5', np.random.default_rng(9))
    elif change == 'extra':
        os.remove(tiles / NAMES[0])
    elif change == 'stale':
        time.sleep(0.01)
        write_stats_tile(tiles / NAMES[1], np.random.default_rng(9))
    elif change == 'duplicate':
        tile_meta.write_meta(str(tiles), str(meta / 'zz_all.json'))
    elif change == 'empty':
        for f in meta.iterdir():
            f.unlink()
    with pytest.raises(RuntimeError, match='tile_meta'):
        tile_meta.tile_records(SimpleNamespace(tiles_dir=str(tiles), tile_meta_dir=str(meta)))


def test_tiles_are_read_once_per_run(tmp_path, monkeypatch):
    tiles = make_tiles(tmp_path / 'prelim', np.random.default_rng(6))
    opened = []
    real = tile_meta.tile_record
    monkeypatch.setattr(tile_meta, 'tile_record', lambda p: opened.append(p) or real(p))
    args = SimpleNamespace(tiles_dir=str(tiles))
    for _ in range(4):           # ATL15: four files, each tile_stats + lineage
        tile_meta.tile_records(args)
        lineage(args)
    assert len(opened) == len(NAMES)


def test_remote_tiles_and_remote_meta(tmp_path, bucket):  # noqa: F811
    tiles = make_tiles(tmp_path / 'prelim', np.random.default_rng(7))
    prefix = bucket(tiles)
    meta = prefix.rsplit('/prelim', 1)[0] + '/tile_meta'
    # --workers 1: spawned processes would not see the patched in-memory fs
    for out, cx, cy in [('0_-900000', '100000', '-900000'), ('rest', '300000', '-900000'),
                        ('r2', '300000', '-700000')]:
        assert tile_meta.main(['write', prefix, f'{meta}/{out}.json', '--center', cx, cy,
                               '--workers', '1']) == 0
    direct = lineage(SimpleNamespace(tiles_dir=prefix))
    saved = lineage(SimpleNamespace(tiles_dir=prefix, tile_meta_dir=meta))
    assert direct == saved


def test_a_cell_with_no_tile_center_writes_an_empty_file(tmp_path):
    # the footprint rule makes 200 km tiles for cells tiles only reach into
    tiles = make_tiles(tmp_path / 'prelim', np.random.default_rng(8))
    assert tile_meta.main(['write', str(tiles), str(tmp_path / 'x.json'),
                           '--center', '-100000', '-900000', '--workers', '1']) == 0
    assert json.loads((tmp_path / 'x.json').read_text())['tiles'] == []


def test_write_with_no_tiles_at_all_fails(tmp_path):
    (tmp_path / 'prelim').mkdir()
    assert tile_meta.main(['write', str(tmp_path / 'prelim'), str(tmp_path / 'x.json')]) == 1


def test_process_pool_gives_the_same_records(tmp_path):
    tiles = make_tiles(tmp_path / 'prelim', np.random.default_rng(10))
    tile_meta.write_meta(str(tiles), str(tmp_path / 'a' / 'all.json'), workers=1)
    tile_meta.write_meta(str(tiles), str(tmp_path / 'b' / 'all.json'), workers=3)
    a = json.loads((tmp_path / 'a' / 'all.json').read_text())['tiles']
    b = json.loads((tmp_path / 'b' / 'all.json').read_text())['tiles']
    assert a == b and len(a) == len(NAMES)


def test_a_manifest_gathers_a_sector_from_two_prefixes(tmp_path, bucket):  # noqa: F811
    """plan_AA_mosaic_dps.sh M0: an Antarctic sector's tiles are under the 60 km and the
    44 km prefixes; its manifest lists them, and the records of the 200 km jobs (one
    meta dir per sector, written from each half) must match it exactly"""
    tiles = make_tiles(tmp_path / 'prelim', np.random.default_rng(11))
    # E0_N-1000, E40_N-1000, E0_N-960 from one half; E200_N-1000, E240_N-760 from the other
    half_a, half_b = tmp_path / 'a', tmp_path / 'b'
    for d, names in ((half_a, NAMES[:3]), (half_b, NAMES[3:])):
        d.mkdir()
        for n in names:
            os.link(tiles / n, d / n)
    pa = bucket(half_a, 'memory://bucket/rel006/south/AA/prelim')
    pb = bucket(half_b, 'memory://bucket/rel006/south/AA_44km/prelim')
    meta = 'memory://bucket/rel006_0332_testing/south/A4/tile_meta'
    assert tile_meta.main(['write', pa, f'{meta}/100000_-900000.json', '--center', '100000', '-900000',
                           '--workers', '1']) == 0
    for c in ('300000 -900000', '300000 -700000'):
        cx, cy = c.split()
        assert tile_meta.main(['write', pb, f'{meta}/{cx}_{cy}.json', '--center', cx, cy,
                               '--workers', '1']) == 0
    import pointCollection as pc
    fs = pc.io_utils.get_s3fs()
    manifest = 'memory://bucket/rel006/south/A4/prelim_tiles.txt'
    fs.pipe(manifest, ''.join(f'{p}/{n}\n' for p, ns in ((pb, NAMES[3:]), (pa, NAMES[:3]))
                              for n in ns).encode())
    from ATL1415.paths import list_tiles
    assert [os.path.basename(t) for t in list_tiles(manifest)] == \
        sorted(NAMES[:3]) + sorted(NAMES[3:])            # sorted by full path: AA/ before AA_44km/
    direct = lineage(SimpleNamespace(tiles_dir=str(tiles)))
    assert lineage(SimpleNamespace(tiles_dir=manifest, tile_meta_dir=meta)) == direct
    # a tile the manifest lists but no record holds stops the writer
    tile_meta._CACHE.clear()
    fs.pipe(manifest, fs.cat(manifest) + f'{pa}/E999_N999.h5\n'.encode())
    fs.pipe(f'{pa}/E999_N999.h5', (tiles / NAMES[0]).read_bytes())
    with pytest.raises(RuntimeError, match='1 with no record: E999_N999.h5'):
        tile_meta.tile_records(SimpleNamespace(tiles_dir=manifest, tile_meta_dir=meta))


def test_a_manifest_naming_a_tile_twice_is_refused(tmp_path):
    from ATL1415.paths import list_tiles
    m = tmp_path / 'prelim_tiles.txt'
    m.write_text('/x/AA/prelim/E0_N0.h5\n/x/AA_44km/prelim/E0_N0.h5\n')
    with pytest.raises(ValueError, match='more than once'):
        list_tiles(str(m))
