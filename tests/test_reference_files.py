"""
--ATL14_reference_file: a single file, a local glob, or (Ben 2026-10-18, AA monthly
on S3) a URI glob over the Antarctic sectors' ATL14s, which ATL11_to_ATL15 mosaics
with pc.grid.mosaic().from_list.  No network: the "bucket" is fsspec's memory
filesystem.
"""
import fsspec
import numpy as np
import pointCollection as pc
import pytest

from ATL1415.ATL11_to_ATL15 import _expand_reference_files

PATTERN = 'memory://bucket/rel006_0332_testing/south/A?/ATL14_A?_0332_100m_006_02.nc'


def sector(x, y, h0):
    """a small ATL14-like grid: h = h0 + x/1000, h_sigma = 1"""
    xx, yy = np.meshgrid(x, y)
    return pc.grid.data().from_dict({'x': x, 'y': y, 'h': h0 + xx / 1000.,
                                     'h_sigma': np.ones_like(xx)})


@pytest.fixture
def bucket(monkeypatch, tmp_path):
    """four sectors (quadrants) at 1 km, written to a memory bucket"""
    fs = fsspec.filesystem('memory')
    monkeypatch.setattr(pc.io_utils, 'get_s3fs', lambda daac=None, **kw: fs)
    neg, pos = np.arange(-20000., 0., 1000.), np.arange(0., 20001., 1000.)
    for n, (x, y) in enumerate([(pos, pos), (neg, pos), (neg, neg), (pos, neg)], start=1):
        local = str(tmp_path / f'A{n}.nc')
        sector(x, y, 100. * n).to_nc(local, replace=True)
        with open(local, 'rb') as fh:
            fs.pipe(PATTERN.replace('?', str(n)), fh.read())
    yield fs
    fs.rm('memory://bucket', recursive=True)


def test_a_uri_pattern_lists_the_sectors_sorted(bucket):
    files = _expand_reference_files(PATTERN)
    assert files == sorted(files) and len(files) == 4
    assert [f.split('/')[-2] for f in files] == ['A1', 'A2', 'A3', 'A4']
    assert all(f.startswith('memory://') for f in files)


def test_a_uri_pattern_that_matches_nothing_is_an_error(bucket):
    with pytest.raises(RuntimeError, match='matched no files'):
        _expand_reference_files(PATTERN.replace('0332', '0331'))


def test_a_single_uri_passes_through_unlisted():
    uri = 's3://bucket/rel006_0332_testing/north/GL/ATL14_GL_0332_100m_006_02.nc'
    assert _expand_reference_files(uri) == [uri]


def test_a_local_pattern_is_globbed_sorted(tmp_path):
    for n in (2, 1):
        (tmp_path / f'ATL14_A{n}.nc').write_text('')
    assert _expand_reference_files(str(tmp_path / 'ATL14_A?.nc')) == \
        [str(tmp_path / 'ATL14_A1.nc'), str(tmp_path / 'ATL14_A2.nc')]
    with pytest.raises(RuntimeError):
        _expand_reference_files(str(tmp_path / 'nothing*.nc'))


def test_the_reference_read_mosaics_the_sectors_a_tile_overlaps(bucket):
    """the call ATL11_to_ATL15 makes, for a tile across all four sectors"""
    bounds = [[-5000., 5000.], [-3000., 3000.]]
    ref = pc.grid.mosaic().from_list(_expand_reference_files(PATTERN), group='',
                                     bounds=bounds, fields=['h', 'h_sigma'])
    assert np.array_equal(ref.x, np.arange(-5000., 5001., 1000.))
    assert np.array_equal(ref.y, np.arange(-3000., 3001., 1000.))
    xx, yy = np.meshgrid(ref.x, ref.y)
    quadrant = np.where(xx >= 0, np.where(yy >= 0, 1, 4), np.where(yy >= 0, 2, 3))
    assert np.array_equal(ref.h, 100. * quadrant + xx / 1000.)
    assert np.all(ref.h_sigma == 1)
