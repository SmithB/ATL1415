"""
Tests for ATL1415.paths, which keeps s3:// URIs intact where the usual
local-path normalization would corrupt them.

No network access is required.  paths.py is imported directly via importlib,
bypassing ATL1415/__init__.py (which pulls in LSsurf/sparseqr, not installed
in every environment this runs in, and unrelated to this file).
"""
import os
import importlib.util

_HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    'atl1415_paths', os.path.join(_HERE, '..', 'ATL1415', 'paths.py'))
paths = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(paths)


# ---------------------------------------------------------------------------
# path_or_uri
# ---------------------------------------------------------------------------

def test_path_or_uri_leaves_uri_intact():
    """
    abspath() would rewrite this as '<cwd>/s3:/bucket/key.tif' -- it collapses
    the '//' and prefixes the working directory, so the URI must not reach it.
    """
    uri = 's3://maap-ops-workspace/ben_smith/masks/Antarctic/mask.tif'
    assert paths.path_or_uri(uri) == uri


def test_path_or_uri_normalizes_local_paths():
    assert paths.path_or_uri('~/x.tif') == os.path.join(os.path.expanduser('~'), 'x.tif')
    assert os.path.isabs(paths.path_or_uri('relative/x.tif'))


def test_path_or_uri_passes_none_through():
    assert paths.path_or_uri(None) is None


# ---------------------------------------------------------------------------
# join_path_or_uri
# ---------------------------------------------------------------------------

def test_join_relative_onto_uri_root():
    assert paths.join_path_or_uri('s3://bucket/masks', 'Antarctic/m.tif') \
        == 's3://bucket/masks/Antarctic/m.tif'
    # a trailing separator on the root must not double up
    assert paths.join_path_or_uri('s3://bucket/masks/', 'Antarctic/m.tif') \
        == 's3://bucket/masks/Antarctic/m.tif'


def test_join_leaves_absolute_values_alone():
    """a URI or an absolute path is already resolved and must win outright"""
    assert paths.join_path_or_uri('/local/masks', 's3://bucket/m.tif') == 's3://bucket/m.tif'
    assert paths.join_path_or_uri('s3://bucket/masks', '/abs/m.tif') == '/abs/m.tif'


def test_join_relative_onto_local_root_is_unchanged_behaviour():
    assert paths.join_path_or_uri('/local/masks', 'Antarctic/m.tif') \
        == os.path.join('/local/masks', 'Antarctic/m.tif')


# ---------------------------------------------------------------------------
# exists
# ---------------------------------------------------------------------------

def test_exists_local(tmp_path):
    f = tmp_path / 'x.tif'
    f.write_text('')
    assert paths.exists(str(f))
    assert paths.exists(str(tmp_path))          # a directory counts
    assert not paths.exists(str(tmp_path / 'no_such_file'))
