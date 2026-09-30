"""
The netCDF writers read the prelim tiles in place on the bucket
(docs/plan_dps_mosaic.sh D2): --tiles_dir may be an s3:// prefix.
ATL1415.paths.list_tiles/open_tile do the listing and opening, and
set_lineage and make_tile_stats_group use them.

No network: remote tiles live in fsspec's in-memory filesystem, returned by a
patched get_s3fs.
"""
import os
from types import SimpleNamespace

import fsspec
import h5py
import numpy as np
import pytest
import pointCollection as pc

from ATL1415.paths import list_tiles, open_tile
from test_lineage import AT, AT_ATTRS, XO, XO_ATTRS, run_set_lineage, write_tile


@pytest.fixture
def bucket(monkeypatch):
    """a memory filesystem standing in for s3; yields a function that copies a
    local directory's .h5 files to a prefix and returns the prefix"""
    fs = fsspec.filesystem('memory')
    monkeypatch.setattr(pc.io_utils, 'get_s3fs', lambda daac=None, **kw: fs)

    def upload(local_dir, prefix='memory://bucket/rel006/north/GL/prelim'):
        for name in os.listdir(local_dir):
            if name.endswith('.h5'):
                with open(os.path.join(local_dir, name), 'rb') as fh:
                    fs.pipe(f'{prefix}/{name}', fh.read())
        return prefix
    yield upload
    if fs.exists('memory://bucket'):
        fs.rm('memory://bucket', recursive=True)


def test_list_tiles_is_sorted_local_and_remote(tmp_path, bucket):
    for name in ['E40_N0.h5', 'E0_N0.h5', 'E0_N40.h5', 'notes.txt']:
        (tmp_path / name).write_bytes(b'x')
    local = list_tiles(str(tmp_path))
    assert [os.path.basename(f) for f in local] == ['E0_N0.h5', 'E0_N40.h5', 'E40_N0.h5']
    remote = list_tiles(bucket(tmp_path))
    assert [os.path.basename(f) for f in remote] == ['E0_N0.h5', 'E0_N40.h5', 'E40_N0.h5']
    assert all(f.startswith('memory://') for f in remote)


def test_open_tile_reads_the_same_either_way(tmp_path, bucket):
    with h5py.File(tmp_path / 'E0_N0.h5', 'w') as h5f:
        h5f['RMS/data'] = 0.25
        h5f['data/three_sigma_edit'] = np.array([1, 0, 1])
    prefix = bucket(tmp_path)
    for path in [str(tmp_path / 'E0_N0.h5'), prefix + '/E0_N0.h5']:
        with open_tile(path) as h5f:
            assert h5f['RMS/data'][()] == 0.25
            assert h5f['data/three_sigma_edit'][:].sum() == 2


def test_set_lineage_from_a_remote_tiles_dir(tmp_path, bucket, capsys):
    write_tile(tmp_path / 'E1.h5', ','.join([AT, XO]), {AT: AT_ATTRS, XO: XO_ATTRS})
    write_tile(tmp_path / 'E2.h5', XO, {XO: XO_ATTRS})
    write_tile(tmp_path / 'E3.h5', '')          # a matched tile
    assert run_set_lineage(bucket(tmp_path)) == run_set_lineage(tmp_path)
    assert 'failed to open' not in capsys.readouterr().out


def write_stats_tile(path, rng):
    """the datasets make_tile_stats_group reads from a prelim tile"""
    with h5py.File(path, 'w') as h5f:
        h5f['data/three_sigma_edit'] = rng.random(20) > 0.3
        for key in ['data', 'grad2_z0', 'd2z_dt2', 'grad2_dzdt']:
            h5f[f'RMS/{key}'] = rng.random()
        h5f['bias/val'] = rng.normal(size=4)
        h5f['bias/expected'] = 1 + rng.random(4)
        for key in ['d2z0_dx2', 'd2z_dt2', 'd3z_dx2dt']:
            h5f[f'E_RMS/{key}'] = rng.random()


def tile_stats(tiles_dir, out):
    from netCDF4 import Dataset
    from ATL1415.make_tile_stats_group import make_tile_stats_group
    with Dataset(out, 'w') as nc:
        make_tile_stats_group(nc, SimpleNamespace(tiles_dir=str(tiles_dir), region='GL'))
    with Dataset(out) as nc:
        grp = nc['tile_stats']
        return {name: np.ma.getmaskarray(grp[name][:]) for name in grp.variables}, \
               {name: np.ma.getdata(grp[name][:]) for name in grp.variables}


def test_tile_stats_from_a_remote_tiles_dir(tmp_path, bucket):
    rng = np.random.default_rng(1)
    tiles = tmp_path / 'prelim'
    tiles.mkdir()
    for name in ['E0_N-1000.h5', 'E40_N-1000.h5', 'E0_N-960.h5']:
        write_stats_tile(tiles / name, rng)
    local = tile_stats(tiles, str(tmp_path / 'local.nc'))
    remote = tile_stats(bucket(tiles), str(tmp_path / 'remote.nc'))
    (local_mask, local_data), (remote_mask, remote_data) = local, remote
    assert local_data.keys() == remote_data.keys() and 'N_data' in local_data
    for name in local_data:
        assert np.array_equal(local_mask[name], remote_mask[name]), name
        assert np.array_equal(local_data[name], remote_data[name]), name
