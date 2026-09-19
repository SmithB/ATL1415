"""
prune_tile_list.py: which prelim outcomes remove a center from the tile list.

docs/plan_tile_lists.sh TL3.  The rule: successful and no tile -> prune;
failed -> keep and investigate; unfinished -> refuse.  MAAP and S3 are stubbed.
"""
import csv
import importlib.util
import os
import sys
from types import SimpleNamespace

import pytest

HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    'prune_tile_list', os.path.join(HERE, '..', 'scripts', 'maap', 'prune_tile_list.py'))
prune = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(prune)

PREFIX = 's3://b/rel006/north/IS'
COLS = ['identifier', 'x0', 'y0', 'step', 'queue', 'args_file', 'job_id',
        'submitted_utc', 'tile_prefix']


def row(x0, y0, job_id, step='prelim', prefix=PREFIX):
    return {'identifier': f'IS_prelim_E{x0 // 1000}_N{y0 // 1000}', 'x0': str(x0),
            'y0': str(y0), 'step': step, 'queue': 'q', 'args_file': 'a',
            'job_id': job_id, 'submitted_utc': 't', 'tile_prefix': prefix}


# --- the rule -----------------------------------------------------------------

def test_the_rule_over_every_outcome():
    rows = [row(1340000, -2460000, 'ok-tile'),       # successful, tile there
            row(1020000, -2580000, 'ok-notile'),     # successful, no tile: NO DATA
            row(1220000, -2460000, 'bad'),           # failed
            row(1260000, -2460000, 'gone'),          # dismissed
            row(1300000, -2460000, '<submit failed: HTTP 500>'),
            row(1380000, -2460000, 'busy')]          # still running
    status = {'ok-tile': 'successful', 'ok-notile': 'successful', 'bad': 'failed',
              'gone': 'dismissed', 'busy': 'running'}
    out = prune.classify(rows, status.get, lambda p: {'E1340_N-2460.h5'})
    assert out['keep'] == ['E1340_N-2460.h5']
    assert out['prune'] == ['E1020_N-2580.h5']
    assert [n for n, _ in out['investigate']] == ['E1220_N-2460.h5', 'E1260_N-2460.h5',
                                                   'E1300_N-2460.h5']
    assert out['running'] == [('E1380_N-2460.h5', 'running')]


def test_a_failed_no_data_job_from_an_old_build_is_not_pruned():
    # 9266c3d7: the monthly E1020 no-data fit, before TL1, exited 1
    out = prune.classify([row(1020000, -2580000, '9266c3d7')],
                         lambda j: 'failed', lambda p: set())
    assert out['prune'] == []
    assert out['investigate'][0][0] == 'E1020_N-2580.h5'


def test_a_row_without_a_tile_prefix_gets_no_verdict():
    out = prune.classify([row(1020000, -2580000, 'j', prefix='-')],
                         lambda j: 'successful',
                         lambda p: pytest.fail('looked for tiles with no prefix'))
    assert out['unknown_prefix'] == ['E1020_N-2580.h5']
    assert out['prune'] == []


# --- rewriting the list --------------------------------------------------------

def test_rewrite_drops_only_the_named_and_keeps_order(tmp_path):
    path = tmp_path / 'list.txt'
    path.write_text('E1020_N-2420.h5\nE1020_N-2580.h5\nE1140_N-2500.h5\n')
    assert prune.rewrite(str(path), {'E1020_N-2580.h5'}) == 1
    assert path.read_text() == 'E1020_N-2420.h5\nE1140_N-2500.h5\n'


# --- main, end to end ------------------------------------------------------------

def setup(tmp_path, monkeypatch, rows, status, present):
    ledger = tmp_path / 'ledger.csv'
    with open(ledger, 'w', newline='') as fh:
        w = csv.DictWriter(fh, COLS)
        w.writeheader()
        w.writerows(rows)
    listing = tmp_path / '40km_tile_list.txt'
    listing.write_text('E1020_N-2580.h5\nE1340_N-2460.h5\nE1220_N-2460.h5\n')

    class FakeMAAP:
        def __init__(self, **kw):
            pass

        def get_job_status(self, jid):
            return SimpleNamespace(json=lambda: {'status': status[jid]})
    import maap.maap
    monkeypatch.setattr(maap.maap, 'MAAP', FakeMAAP)
    monkeypatch.setattr(prune, 's3_names', lambda prefix: present)
    return str(ledger), str(listing)


def test_write_prunes_the_no_data_center(tmp_path, monkeypatch, capsys):
    ledger, listing = setup(tmp_path, monkeypatch,
                            [row(1020000, -2580000, 'a'), row(1340000, -2460000, 'b')],
                            {'a': 'successful', 'b': 'successful'}, {'E1340_N-2460.h5'})
    assert prune.main([ledger, listing, '--write']) == 0
    assert open(listing).read() == 'E1340_N-2460.h5\nE1220_N-2460.h5\n'


def test_without_write_the_list_is_untouched(tmp_path, monkeypatch):
    ledger, listing = setup(tmp_path, monkeypatch, [row(1020000, -2580000, 'a')],
                            {'a': 'successful'}, set())
    before = open(listing).read()
    assert prune.main([ledger, listing]) == 0
    assert open(listing).read() == before


def test_a_failure_to_investigate_exits_1_but_still_prunes(tmp_path, monkeypatch, capsys):
    ledger, listing = setup(tmp_path, monkeypatch,
                            [row(1020000, -2580000, 'a'), row(1220000, -2460000, 'b')],
                            {'a': 'successful', 'b': 'failed'}, set())
    assert prune.main([ledger, listing, '--write']) == 1
    assert 'E1020_N-2580.h5' not in open(listing).read()
    assert 'E1220_N-2460.h5' in open(listing).read()          # kept
    assert 'INVESTIGATE' in capsys.readouterr().out


def test_unfinished_jobs_refuse_and_change_nothing(tmp_path, monkeypatch):
    ledger, listing = setup(tmp_path, monkeypatch,
                            [row(1020000, -2580000, 'a'), row(1340000, -2460000, 'b')],
                            {'a': 'successful', 'b': 'running'}, set())
    before = open(listing).read()
    assert prune.main([ledger, listing, '--write']) == 2
    assert open(listing).read() == before


def test_a_matched_ledger_is_refused(tmp_path, monkeypatch):
    ledger, listing = setup(tmp_path, monkeypatch,
                            [row(1020000, -2580000, 'a', step='matched')],
                            {'a': 'successful'}, set())
    assert prune.main([ledger, listing, '--write']) == 2


def test_a_list_that_is_not_tile_names_is_refused(tmp_path, monkeypatch):
    ledger, listing = setup(tmp_path, monkeypatch, [row(1020000, -2580000, 'a')],
                            {'a': 'successful'}, set())
    with open(listing, 'a') as fh:
        fh.write('field_sizes\n')                   # the real AA list's line 8945
    assert prune.main([ledger, listing, '--write']) == 2
