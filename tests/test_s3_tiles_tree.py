"""
scripts/s3_tiles.py put_tree / get_glob: how a DPS mosaic or nc job moves its
products (plan_dps_mosaic D3-1).  An in-memory fsspec filesystem stands in
for s3fs.
"""
import importlib.util
import os

import fsspec
import pytest

HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location('s3_tiles', os.path.join(HERE, '..', 'scripts', 's3_tiles.py'))
s3_tiles = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s3_tiles)


@pytest.fixture
def fs():
    fs = fsspec.filesystem('memory')
    yield fs
    if fs.exists('/bucket'):
        fs.rm('/bucket', recursive=True)


def test_put_tree_keeps_relative_paths(tmp_path, fs):
    (tmp_path / '200km_tiles' / 'dz').mkdir(parents=True)
    (tmp_path / '200km_tiles' / 'dz' / 'dz0_200_-1000_-800.h5').write_bytes(b'a')
    (tmp_path / '200km_tiles' / 'z0').mkdir()
    (tmp_path / '200km_tiles' / 'z0' / 'z00_200_-1000_-800.h5').write_bytes(b'b')
    assert s3_tiles.put_tree(fs, str(tmp_path / '200km_tiles'), '/bucket/GL/200km_tiles/') == 0
    assert fs.cat('/bucket/GL/200km_tiles/dz/dz0_200_-1000_-800.h5') == b'a'
    assert fs.cat('/bucket/GL/200km_tiles/z0/z00_200_-1000_-800.h5') == b'b'


def test_put_tree_of_nothing_fails(tmp_path, fs):
    assert s3_tiles.put_tree(fs, str(tmp_path), '/bucket/GL') == 1


def test_get_glob_takes_only_matches(tmp_path, fs):
    for name in ['z0.h5', 'dz.h5', 'ATL14_GL.nc']:
        fs.pipe(f'/bucket/GL/{name}', name.encode())
    fs.pipe('/bucket/GL/200km_tiles/dz/x.h5', b'not this')
    assert s3_tiles.get_glob(fs, '/bucket/GL', '*.h5', str(tmp_path)) == 0
    assert sorted(os.listdir(tmp_path)) == ['dz.h5', 'z0.h5']


@pytest.mark.parametrize('require, status', [(False, 0), (True, 1)])
def test_get_glob_with_no_match(tmp_path, fs, require, status):
    assert s3_tiles.get_glob(fs, '/bucket/GL', 'bounds.txt', str(tmp_path), require=require) == status
