"""
Tests for scripts/setup_ATL1415_region.py, which composes the several
default_args files into the single input_args_<REGION>.txt that every later
stage reads.

No network access is required.  The module is loaded directly via importlib
behind a stub ATL1415 package, so ATL1415/__init__.py (LSsurf/sparseqr, not
installed in every environment this runs in) is never executed.
"""
import os
import sys
import types
import importlib.util
import pytest

_HERE = os.path.dirname(__file__)
_ATL1415_DIR = os.path.abspath(os.path.join(_HERE, '..', 'ATL1415'))


def _load(name, relpath):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(_ATL1415_DIR, relpath))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_pkg = types.ModuleType('ATL1415')
_pkg.__path__ = [_ATL1415_DIR]
sys.modules['ATL1415'] = _pkg
_pkg.infer_dzdt_lags = _load('ATL1415.lags', 'lags.py').infer_dzdt_lags
_pkg.paths = _load('ATL1415.paths', 'paths.py')

setup_region = _load('setup_ATL1415_region',
                     os.path.join('scripts', 'setup_ATL1415_region.py'))


def _compose(tmp_path, lines, extra_files=()):
    """
    run main() over a defaults file built from `lines`, and return the
    composed input_args_<REGION>.txt as a list of stripped lines
    """
    root = tmp_path / 'ATL14_processing'
    root.mkdir(exist_ok=True)
    defaults = tmp_path / 'defaults.txt'
    defaults.write_text('\n'.join(lines) + '\n')

    argv = ['setup_ATL1415_region.py', str(defaults), *map(str, extra_files)]
    old_argv = sys.argv
    sys.argv = argv
    try:
        setup_region.main()
    finally:
        sys.argv = old_argv

    out = root / 'rel006' / 'south' / 'AA' / 'input_args_AA.txt'
    return [ln.strip() for ln in out.read_text().splitlines() if ln.strip()]


def _base_lines(tmp_path, index_path):
    return [
        f'--ATL14_root={tmp_path}/ATL14_processing',
        '--region=AA',
        '--Release=006',
        '--Hemisphere=-1',
        '--mask_file=Antarctic/mask.tif',
        f'--ATL11_index={index_path}',
    ]


# ---------------------------------------------------------------------------
# bare store_true flags
# ---------------------------------------------------------------------------

def test_bare_flags_survive_composition(tmp_path):
    """
    Regression: the key=value regex does not match a bare '--flag' line, so
    a store_true option asked for in a defaults file used to be dropped from
    the composed file without any warning -- the run then silently proceeded
    with the option off.
    """
    index = tmp_path / 'GeoIndex.h5'
    index.write_text('')
    lines = _base_lines(tmp_path, index) + ['--tide_adjustment']
    composed = _compose(tmp_path, lines)
    assert '--tide_adjustment' in composed


def test_value_containing_equals_round_trips(tmp_path):
    """
    the key is captured non-greedily, so only the FIRST '=' separates key
    from value -- otherwise a proj4 string would be split at its last '='
    """
    index = tmp_path / 'GeoIndex.h5'
    index.write_text('')
    lines = _base_lines(tmp_path, index) + ['--SRS_proj4=+proj=stere +lat_0=-90']
    composed = _compose(tmp_path, lines)
    assert '--SRS_proj4=+proj=stere +lat_0=-90' in composed


def test_surrounding_whitespace_is_stripped(tmp_path):
    index = tmp_path / 'GeoIndex.h5'
    index.write_text('')
    lines = _base_lines(tmp_path, index) + ['-W=60000   ']
    composed = _compose(tmp_path, lines)
    assert '-W=60000' in composed


# ---------------------------------------------------------------------------
# --mask_dir joining
# ---------------------------------------------------------------------------

def test_mask_dir_joins_onto_s3_prefix(tmp_path):
    """os.path.join would have produced '<mask_dir>/s3://...' for a URI value"""
    index = tmp_path / 'GeoIndex.h5'
    index.write_text('')
    lines = _base_lines(tmp_path, index) + [
        '--mask_dir=s3://maap-ops-workspace/ben_smith/ATL1415/masks',
        '--tide_mask_file=Antarctic/tide_mask.tif']
    composed = _compose(tmp_path, lines)
    assert ('--mask_file=s3://maap-ops-workspace/ben_smith/ATL1415/masks/'
            'Antarctic/mask.tif') in composed
    assert ('--tide_mask_file=s3://maap-ops-workspace/ben_smith/ATL1415/masks/'
            'Antarctic/tide_mask.tif') in composed
    # --mask_dir is consumed here and must not reach the composed file
    assert not any(ln.startswith('--mask_dir') for ln in composed)


def test_mask_dir_leaves_an_absolute_uri_value_alone(tmp_path):
    index = tmp_path / 'GeoIndex.h5'
    index.write_text('')
    lines = _base_lines(tmp_path, index) + [
        '--mask_dir=s3://bucket/masks',
        '--geoid_file=s3://other-bucket/geoid/EGM2008_geoid_h.nc']
    composed = _compose(tmp_path, lines)
    assert '--geoid_file=s3://other-bucket/geoid/EGM2008_geoid_h.nc' in composed


# ---------------------------------------------------------------------------
# cloud ATL11 index
# ---------------------------------------------------------------------------

def test_cloud_index_root_may_be_a_directory(tmp_path):
    """
    Regression: with --ATL11_earthaccess the index is the ROOT of a
    per-granule geoIndex tree, i.e. a directory.  The old os.path.isfile()
    check rejected it outright, so no cloud args file could be composed.
    """
    index_root = tmp_path / 'ATL11_index'
    (index_root / 'ATL11_index_0331_007_04').mkdir(parents=True)
    lines = _base_lines(tmp_path, index_root) + ['--ATL11_earthaccess']
    composed = _compose(tmp_path, lines)
    assert f'--ATL11_index={index_root}' in composed
    assert '--ATL11_earthaccess' in composed


def test_cloud_index_is_not_synthesized_from_release(tmp_path):
    """
    the local layout '<root>/ATL11_<release>/<hemi>/index/GeoIndex.h5' does not
    exist for a cloud run, so --ATL11_index must be given rather than derived
    """
    lines = [ln for ln in _base_lines(tmp_path, '') if '--ATL11_index' not in ln]
    lines += ['--ATL11_release=007_cycle_03_31_v04', '--ATL11_earthaccess']
    with pytest.raises(OSError, match='must be given explicitly'):
        _compose(tmp_path, lines)


def test_local_index_must_still_exist(tmp_path):
    lines = _base_lines(tmp_path, tmp_path / 'no_such_index.h5')
    with pytest.raises(OSError, match='does not exist'):
        _compose(tmp_path, lines)
