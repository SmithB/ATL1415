"""collect_jobs.py --csv: one row per tile (plan_pack_tiles K5a)."""
import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def cj():
    spec = importlib.util.spec_from_file_location(
        'collect_jobs', REPO / 'scripts' / 'maap' / 'collect_jobs.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def one_tile_row(**kw):
    return dict({'identifier': 'GL_x_E160_N-1640', 'step': 'prelim', 'job_id': 'j1',
                 'tiles': '-'}, **kw)


PACKED = {'identifier': 'GL_x_E0_N0_n4', 'step': 'prelim', 'job_id': 'j2',
          'tiles': '0,0;40000,0|0,40000;-40000,-40000'}


def test_one_tile_ok(cj):
    fields = {'status': 'successful', 'secs': '590', 'max_mem_GiB': '10.47',
              'instance': 'r5.xlarge', '_nodata': False}
    steps = {'fit': (294.0, 10.47), 'error': (278.0, 6.12)}
    [row] = cj.tile_rows(one_tile_row(), fields, steps)
    assert row == {'tile': 'E160_N-1640', 'step': 'prelim', 'status': 'ok', 'secs': '590',
                   'solve_secs': '572', 'max_mem_GiB': '10.47', 'instance': 'r5.xlarge',
                   'lane': '-', 'identifier': 'GL_x_E160_N-1640', 'job_id': 'j1'}


def test_one_tile_nodata_and_failed(cj):
    [row] = cj.tile_rows(one_tile_row(), {'status': 'successful', '_nodata': True}, {})
    assert row['status'] == 'nodata' and row['solve_secs'] == '-'
    [row] = cj.tile_rows(one_tile_row(), {'status': 'failed'}, {})
    assert row['status'] == 'failed'
    # an old ledger has no `tiles` column
    old = one_tile_row()
    del old['tiles']
    [row] = cj.tile_rows(old, {'status': 'running'}, {})
    assert (row['tile'], row['status']) == ('E160_N-1640', 'running')


def packed_fields(job_status):
    ran = [({'tile': '  E0_N0', 'status': 'ok', 'secs': '982', 'max_mem_GiB': '10.48',
             'queue': 'lane 1'}, {'fit': (400.0, 10.48), 'error': (577.0, 6.1)}),
           ({'tile': '  E0_N40', 'status': 'failed', 'secs': '50', 'max_mem_GiB': '0.70',
             'queue': 'lane 2'}, {'fit': (49.0, 0.7)})]
    return {'status': job_status, 'instance': 'c5.4xlarge', '_tiles': ran}


def test_packed_tiles_in_ledger_order(cj):
    rows = cj.tile_rows(PACKED, packed_fields('failed'), {})
    got = [(r['tile'], r['status'], r['secs'], r['solve_secs'], r['lane']) for r in rows]
    # tiles that never reported take the job's status
    assert got == [('E0_N0', 'ok', '982', '977', '1'),
                   ('E40_N0', 'failed', '-', '-', '-'),
                   ('E0_N40', 'failed', '50', '49', '2'),
                   ('E-40_N-40', 'failed', '-', '-', '-')]
    assert {r['instance'] for r in rows} == {'c5.4xlarge'}
    assert {r['identifier'] for r in rows} == {'GL_x_E0_N0_n4'}


def test_packed_unreported_tile(cj):
    # a successful job whose tile never reported: failed, as run.sh says
    rows = cj.tile_rows(PACKED, packed_fields('successful'), {})
    assert rows[1]['status'] == 'failed'
    # a job still running: its unreported tiles say so
    rows = cj.tile_rows(PACKED, {'status': 'running'}, {})
    assert [r['status'] for r in rows] == ['running'] * 4


def test_collect_failed_row(cj):
    rows = cj.tile_rows(PACKED, {'status': 'collect failed'}, {})
    assert [r['status'] for r in rows] == ['collect failed'] * 4


def test_nodata_lines(cj):
    assert cj.NODATA_RE.search('x\nno fit written for E0_N0.h5 (insufficient data); skipping\n')
    assert cj.NODATA_RE.search('no tile to upload for E0_N0.h5 (removed: no data ...)')
    assert not cj.NODATA_RE.search('uploaded E0_N0.h5')


@pytest.mark.parametrize('argv,want', [
    (['l.csv'], ('l.csv', None)),
    (['l.csv', '--csv', 'o.csv'], ('l.csv', 'o.csv')),
])
def test_parse_argv(cj, argv, want):
    assert cj.parse_argv(argv) == want


@pytest.mark.parametrize('argv', [[], ['--csv', 'o.csv'], ['l.csv', '--cs', 'o.csv'],
                                  ['l.csv', '--csv'], ['-h']])
def test_parse_argv_refuses(cj, argv):
    with pytest.raises(SystemExit) as e:
        cj.parse_argv(argv)
    assert e.value.code == 2
