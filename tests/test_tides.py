"""
Tests for ATL1415.tides, which sends a local --tide_directory to
pyTMD.compute.tide_elevations unchanged and an s3:// one to the zarr stores
pyTMD publishes.

The default tests need no network and no pyTMD install: pyTMD is stubbed, so
what is under test is the dispatch and the arguments handed on, not pyTMD.
Set ATL1415_TIDE_NETWORK_TESTS=1 to additionally run the live read against
s3://pytmd (needs pyTMD, zarr, s3fs and timescale).

tides.py is imported directly via importlib, bypassing ATL1415/__init__.py.
"""
import os
import sys
import types
import importlib.util
import numpy as np
import pytest

_HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    'atl1415_tides', os.path.join(_HERE, '..', 'ATL1415', 'tides.py'))
tides = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tides)


# ---------------------------------------------------------------------------
# which path a tide_directory takes
# ---------------------------------------------------------------------------

def test_is_remote_tide_directory():
    assert tides.is_remote_tide_directory('s3://pytmd')
    assert tides.is_remote_tide_directory('s3://pytmd/stores/')
    assert not tides.is_remote_tide_directory('/home/jovyan/tide_models')
    assert not tides.is_remote_tide_directory(None)


def test_store_path_matches_pytmd_naming():
    """pyTMD publishes its stores as '<bucket>/<model>.zarr'"""
    assert tides._store_path('s3://pytmd', 'CATS2008-v2023') == 'pytmd/CATS2008-v2023.zarr'
    assert tides._store_path('s3://pytmd/', 'Gr1km-v2') == 'pytmd/Gr1km-v2.zarr'


# ---------------------------------------------------------------------------
# a local directory must still go to pyTMD unchanged
# ---------------------------------------------------------------------------

@pytest.fixture
def fake_pytmd(monkeypatch):
    """stand in for pyTMD and record how compute.tide_elevations was called"""
    calls = []

    def tide_elevations(x, y, delta_time, **kwargs):
        calls.append(kwargs)
        return np.zeros_like(np.asarray(x, dtype=float))

    mod = types.ModuleType('pyTMD')
    mod.compute = types.SimpleNamespace(tide_elevations=tide_elevations)
    monkeypatch.setitem(sys.modules, 'pyTMD', mod)
    return calls


def test_local_directory_delegates_to_pytmd(fake_pytmd):
    x = np.array([0., 1.]); y = np.array([0., 1.]); dt = np.array([0., 1.])
    out = tides.tide_elevations(x, y, dt, tide_directory='/local/tide_models',
                                tide_model='CATS2008', crs=3031,
                                extrapolate=True, cutoff=10.)
    assert out.shape == x.shape
    assert len(fake_pytmd) == 1
    kw = fake_pytmd[0]
    # the arguments the production call has always passed must survive intact
    assert kw['directory'] == '/local/tide_models'
    assert kw['model'] == 'CATS2008'
    assert kw['crs'] == 3031
    assert kw['type'] == 'drift'
    assert kw['standard'] == 'GPS'
    assert kw['epoch'] == (2018, 1, 1, 0, 0, 0)
    assert kw['extrapolate'] is True and kw['cutoff'] == 10.


def test_remote_directory_does_not_call_compute(fake_pytmd, monkeypatch):
    """an s3:// directory must not fall through to the local reader"""
    def boom(*a, **k):
        raise AssertionError('the remote path must not open a local model')
    monkeypatch.setattr(tides, 'open_tide_dataset', boom)
    with pytest.raises(AssertionError):
        tides.tide_elevations(np.array([0.]), np.array([0.]), np.array([0.]),
                              tide_directory='s3://pytmd', tide_model='Gr1km-v2')
    assert fake_pytmd == []


# ---------------------------------------------------------------------------
# anonymous access
# ---------------------------------------------------------------------------

def test_defaults_to_anonymous_access(monkeypatch):
    """
    pyTMD's bucket is public but REJECTS signed requests from other accounts,
    so a DPS worker's own credentials would get AccessDenied.  The default has
    to be anonymous, and it has to be overridable.
    """
    seen = {}

    class FakeS3FS:
        def __init__(self, **kw): seen.update(kw)

    fake_s3fs = types.ModuleType('s3fs'); fake_s3fs.S3FileSystem = FakeS3FS
    fake_zarr = types.ModuleType('zarr')
    fake_zarr.storage = types.SimpleNamespace(FsspecStore=lambda fs, path=None: ('store', path))
    fake_xr = types.ModuleType('xarray')
    fake_xr.open_zarr = lambda store, group=None, zarr_format=None: 'ds'
    fake_pytmd = types.ModuleType('pyTMD')
    fake_pytmd.io = types.SimpleNamespace(
        model=lambda **kw: types.SimpleNamespace(from_database=lambda n: 'model'))
    for name, mod in [('s3fs', fake_s3fs), ('zarr', fake_zarr),
                      ('xarray', fake_xr), ('pyTMD', fake_pytmd)]:
        monkeypatch.setitem(sys.modules, name, mod)

    monkeypatch.delenv('ATL1415_TIDE_ANON', raising=False)
    tides.open_tide_dataset('s3://pytmd', 'Gr1km-v2')
    assert seen == {'anon': True}

    seen.clear()
    monkeypatch.setenv('ATL1415_TIDE_ANON', '0')
    tides.open_tide_dataset('s3://pytmd', 'Gr1km-v2')
    assert seen == {'anon': False}


# ---------------------------------------------------------------------------
# live read (opt in)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(os.environ.get('ATL1415_TIDE_NETWORK_TESTS') != '1',
                    reason='set ATL1415_TIDE_NETWORK_TESTS=1 to read s3://pytmd')
@pytest.mark.parametrize('model,epsg,xc,yc,expected', [
    ('CATS2008-v2023', 3031, -4.0e5, -9.6e5, +0.00923),
    ('Gr1km-v2',       3413, -2.0e5, -2.2e6, -0.00677),
])
def test_live_zarr_read(model, epsg, xc, yc, expected):
    rng = np.random.default_rng(0)
    x = xc + rng.uniform(-3.e4, 3.e4, 5000)
    y = yc + rng.uniform(-3.e4, 3.e4, 5000)
    dt = rng.uniform(0., 7*86400., 5000)
    h = tides.tide_elevations(x, y, dt, tide_directory='s3://pytmd',
                              tide_model=model, crs=epsg)
    assert np.isfinite(h).all()
    assert np.abs(h).max() < 10.
