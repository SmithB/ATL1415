"""
fetch_tiles.py saves the no-data centers for the cleanup after a run.

docs/plan_tile_lists.sh TL8 (Ben, 2026-09-19: flag, don't prune).  A
SUCCESSFUL prelim job that left no tile found no data; its name goes to
<region_dir>/prelim/no_data_tiles.txt, merged with what is there.  For matched,
"no tile" is unexpected and is reported, not saved.  fetch_row and MAAP are
stubbed: no network.
"""
import importlib.util
import os
import sys
from types import SimpleNamespace

import pytest

HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    'fetch_tiles', os.path.join(HERE, '..', 'scripts', 'maap', 'fetch_tiles.py'))
fetch = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fetch)

LEDGER_HEAD = 'identifier,x0,y0,step,queue,args_file,job_id,submitted_utc,tile_prefix\n'


# --- record_no_data ------------------------------------------------------------

def test_a_new_file_gets_the_names_sorted(tmp_path):
    path = str(tmp_path / 'prelim' / 'no_data_tiles.txt')
    added = fetch.record_no_data(path, ['E1180_N-2380.h5', 'E1020_N-2580.h5'])
    assert added == ['E1020_N-2580.h5', 'E1180_N-2380.h5']
    assert open(path).read() == 'E1020_N-2580.h5\nE1180_N-2380.h5\n'


def test_names_merge_and_nothing_is_dropped(tmp_path):
    # a retry ledger, or a second fetch, must not lose what the first found
    path = tmp_path / 'no_data_tiles.txt'
    path.write_text('E1020_N-2580.h5\n')
    added = fetch.record_no_data(str(path), ['E1180_N-2380.h5', 'E1020_N-2580.h5'])
    assert added == ['E1180_N-2380.h5']
    assert path.read_text() == 'E1020_N-2580.h5\nE1180_N-2380.h5\n'


def test_refetching_changes_nothing(tmp_path):
    path = tmp_path / 'no_data_tiles.txt'
    path.write_text('E1020_N-2580.h5\n')
    before = os.stat(path).st_mtime_ns
    assert fetch.record_no_data(str(path), ['E1020_N-2580.h5']) == []
    assert os.stat(path).st_mtime_ns == before          # not even rewritten


def test_dry_run_writes_nothing(tmp_path):
    path = str(tmp_path / 'prelim' / 'no_data_tiles.txt')
    assert fetch.record_no_data(path, ['E1020_N-2580.h5'], dry_run=True) == ['E1020_N-2580.h5']
    assert not os.path.exists(path)


# --- main: prelim saves, matched reports -------------------------------------------

def run_main(tmp_path, monkeypatch, step, verdicts, dry_run=False):
    region = tmp_path / 'IS'
    region.mkdir()
    ledger = tmp_path / 'ledger.csv'
    rows = [(1340000, -2460000), (1020000, -2580000), (1220000, -2460000)]
    ledger.write_text(LEDGER_HEAD + ''.join(
        f'IS_{step}_E{x // 1000}_N{y // 1000},{x},{y},{step},q,a,j{i},t,s3://b/IS\n'
        for i, (x, y) in enumerate(rows)))
    it = iter(verdicts)
    monkeypatch.setattr(fetch, 'fetch_row', lambda maap, row, args: next(it))
    monkeypatch.setattr(fetch, 'MAAP', lambda **k: None)
    monkeypatch.setattr('sys.argv', ['fetch_tiles.py', str(ledger), str(region),
                                     '--step', step] + (['--dry-run'] if dry_run else []))
    fetch.main()
    return region


def test_prelim_saves_only_the_no_tile_rows(tmp_path, monkeypatch, capsys):
    region = run_main(tmp_path, monkeypatch, 'prelim',
                      [('fetched', 10), ('no tile', 0), ('FAILED', 0)])
    out = capsys.readouterr().out
    assert open(region / 'prelim' / 'no_data_tiles.txt').read() == 'E1020_N-2580.h5\n'
    assert 'NO DATA' in out and 'E1020_N-2580.h5' in out
    assert 'FAILED: IS_prelim_E1220_N-2460' in out      # a failure stays a failure
    assert 'no tile:' not in out                         # not under NOT FETCHED


def test_prelim_dry_run_reports_but_saves_nothing(tmp_path, monkeypatch, capsys):
    region = run_main(tmp_path, monkeypatch, 'prelim',
                      [('would fetch', 10), ('no tile', 0), ('would fetch', 10)], dry_run=True)
    assert 'would add 1' in capsys.readouterr().out
    assert not (region / 'prelim' / 'no_data_tiles.txt').exists()


def test_matched_no_tile_is_reported_not_saved(tmp_path, monkeypatch, capsys):
    region = run_main(tmp_path, monkeypatch, 'matched',
                      [('fetched', 10), ('no tile', 0), ('fetched', 10)])
    out = capsys.readouterr().out
    assert 'NOT FETCHED' in out and 'no tile: IS_matched_E1020_N-2580' in out
    assert not (region / 'prelim' / 'no_data_tiles.txt').exists()
    assert 'NO DATA' not in out


def test_a_clean_prelim_fetch_writes_no_file(tmp_path, monkeypatch, capsys):
    region = run_main(tmp_path, monkeypatch, 'prelim',
                      [('fetched', 10), ('fetched', 10), ('fetched', 10)])
    assert not (region / 'prelim' / 'no_data_tiles.txt').exists()
    assert 'NO DATA' not in capsys.readouterr().out
