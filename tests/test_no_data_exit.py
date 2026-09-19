"""
ATL11_to_ATL15.main(): which no-data outcomes exit 0, and which stay failures.

docs/plan_tile_lists.sh TL1.  smooth_fit has two no-data exits
(smooth_fit.py:485, data None; :513, data empty after masking) and returns
normally from both, with an empty 'm' and 'E'.  main() turns that into:
  - prelim FIT, no data       -> 0, no tile   (TL1, new)
  - error step, no data       -> 0, tile gone (I7a, plan_IS_run.sh)
  - MATCHED fit, no data      -> 1            (unexpected: its prelim had data)
  - data present, nothing fit -> 1            (a real fault, never silenced)
The fit itself is stubbed, so no data, LSsurf solve or network is involved.
"""
import importlib
import os
from types import SimpleNamespace

import pytest

mod = importlib.import_module('ATL1415.ATL11_to_ATL15')

EMPTY = SimpleNamespace(size=0)
SOME = SimpleNamespace(size=327)


def run(monkeypatch, tmp_path, S, prelim=True, calc_error=False, stale=False):
    """main() with the fit returning S; returns (status, tile, report)."""
    step_dir = tmp_path / 'prelim'
    step_dir.mkdir()
    tile = step_dir / 'E1020_N-2580.h5'
    report = step_dir / 'field_sizes' / 'E1020_N-2580_report.json'
    if stale or calc_error:
        # the error step always starts from the tile the fit wrote
        report.parent.mkdir()
        tile.write_bytes(b'tile')
        report.write_text('{}')
    args = SimpleNamespace(
        xy0=[1020000, -2580000],
        base_directory=str(tmp_path), out_name=str(tile),
        write_data_only=False, prelim=prelim, matched=not prelim,
        calc_error_file=str(tile) if calc_error else None,
        error_res_scale=[5, 2], dzdt_lags=[1], reference_epoch=0)
    saved = []
    monkeypatch.setattr(mod, 'parse_args', lambda: args)
    monkeypatch.setattr(mod, 'resolve_run_config',
                        lambda a: {'dest_dir': str(step_dir)})
    monkeypatch.setattr(mod, 'build_fit_kwargs', lambda a, c: {})
    monkeypatch.setattr(mod, 'ATL11_to_ATL15', lambda xy0, **kw: S)
    monkeypatch.setattr(mod, 'save_fit_to_file',
                        lambda S, f, **kw: (saved.append('fit'), open(f, 'wb').close()))
    monkeypatch.setattr(mod, 'save_errors_to_file', lambda S, f, **kw: saved.append('errors'))
    monkeypatch.setattr(mod, 'save_field_size_report', lambda f: None)
    monkeypatch.setattr(mod, 'interp_ds', lambda ds, scale: ds)
    return mod.main(), tile, report, saved


def fit(data, m=None, E=None):
    return {'m': m or {}, 'E': E or {}, 'data': data}


# --- TL1: the new branch ----------------------------------------------------

@pytest.mark.parametrize('data', [None, EMPTY], ids=['data_None_485', 'data_empty_513'])
def test_no_data_prelim_fit_exits_0_and_leaves_no_tile(monkeypatch, tmp_path, data):
    status, tile, report, saved = run(monkeypatch, tmp_path, fit(data))
    assert status == 0
    assert not tile.exists() and not report.exists()
    assert saved == []


def test_no_data_prelim_fit_removes_a_stale_tile_and_report(monkeypatch, tmp_path):
    # the verdict is "no tile": nothing left over may contradict it
    status, tile, report, _ = run(monkeypatch, tmp_path, fit(None), stale=True)
    assert status == 0
    assert not tile.exists() and not report.exists()


# --- the boundaries that must NOT move ---------------------------------------

@pytest.mark.parametrize('data', [None, EMPTY], ids=['data_None', 'data_empty'])
def test_no_data_matched_fit_still_fails(monkeypatch, tmp_path, data):
    status, *_ = run(monkeypatch, tmp_path, fit(data), prelim=False)
    assert status == 1


def test_data_present_but_nothing_fit_still_fails(monkeypatch, tmp_path):
    # data survived to the solve and still no model came back: a real fault
    status, *_ = run(monkeypatch, tmp_path, fit(SOME))
    assert status == 1


# --- the paths that already worked, unchanged --------------------------------

def test_normal_prelim_fit_saves_and_exits_0(monkeypatch, tmp_path):
    status, tile, _, saved = run(monkeypatch, tmp_path, fit(SOME, m={'z0': 1}))
    assert status == 0
    assert saved == ['fit'] and tile.exists()


@pytest.mark.parametrize('data', [None, EMPTY], ids=['data_None', 'data_empty'])
def test_no_data_error_step_still_removes_the_tile_I7a(monkeypatch, tmp_path, data):
    status, tile, report, _ = run(monkeypatch, tmp_path, fit(data), calc_error=True)
    assert status == 0
    assert not tile.exists() and not report.exists()


def test_error_step_with_errors_saves_them(monkeypatch, tmp_path):
    E = {'sigma_z0': 1, 'sigma_dz': 1, 'sigma_dzdt_lag1': 1}
    status, tile, _, saved = run(monkeypatch, tmp_path, fit(SOME, E=E), calc_error=True)
    assert status == 0
    assert saved == ['errors'] and tile.exists()
