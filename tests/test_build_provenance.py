"""
Tiles record the build that wrote them (docs/howto_MAAP_ogc.sh O12a).

run.sh exports ATL1415_BUILD_COMMIT / _VERSION / _COMPLETED from the image's
build stamp; write_build_provenance() turns them into /meta attributes.
"""
import h5py
import pytest

from ATL1415.ATL11_to_ATL15 import BUILD_PROVENANCE_ENV, write_build_provenance

VALUES = {'ATL1415_BUILD_COMMIT': '8935494' + 'a' * 33,
          'ATL1415_BUILD_VERSION': 'on_s3',
          'ATL1415_BUILD_COMPLETED': '2026-09-10T20:47:45Z'}


def text(value):
    """
    An attribute as text.  Written as ascii bytes (like meta/input_files), but
    h5py 3.16 hands fixed-length ascii back as str -- and the existing reader,
    ATL1415_attrs_meta.py, already does str(...) on it.
    """
    return value.decode('ascii') if isinstance(value, bytes) else value


@pytest.fixture
def meta(tmp_path):
    with h5py.File(tmp_path / 'tile.h5', 'w') as h5f:
        yield h5f.require_group('meta')


def test_writes_every_field_as_ascii(meta, monkeypatch):
    for env, value in VALUES.items():
        monkeypatch.setenv(env, value)
    write_build_provenance(meta)
    for attr, env in BUILD_PROVENANCE_ENV.items():
        assert text(meta.attrs[attr]) == VALUES[env]


def test_prefix_keeps_the_error_step_separate(meta, monkeypatch):
    for env, value in VALUES.items():
        monkeypatch.setenv(env, value)
    write_build_provenance(meta, prefix='errors_')
    assert text(meta.attrs['errors_build_version']) == 'on_s3'
    assert 'build_version' not in meta.attrs


def test_unset_means_absent_not_empty(meta, monkeypatch):
    # discover, or any run outside MAAP: no stamp, so no attributes at all
    for env in BUILD_PROVENANCE_ENV.values():
        monkeypatch.delenv(env, raising=False)
    write_build_provenance(meta)
    assert not any(a.startswith('build_') for a in meta.attrs)


def test_empty_value_is_not_written(meta, monkeypatch):
    # run.sh exports an empty string when the image has no stamp
    for env in BUILD_PROVENANCE_ENV.values():
        monkeypatch.setenv(env, '')
    write_build_provenance(meta)
    assert not any(a.startswith('build_') for a in meta.attrs)
