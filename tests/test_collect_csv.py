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


class Answer:
    def __init__(self, body):
        self.body = body

    def json(self):
        return self.body


class Maap:
    def get_job_status(self, jid):
        return Answer({'status': 'failed'})

    def get_job_metrics(self, jid):
        return Answer({})

    def get_job_result(self, jid):
        return Answer({})


@pytest.mark.parametrize('text,names,want', [
    # a start-up failure: the container's log was read and no tile started
    ('Traceback ...\nConnectTimeout: api.maap-project.org\n', ['_stdout.txt', '_stderr.txt'], set()),
    ('TILE_START E0_N0 lane 1 2026-10-16T19:00:00Z\nTILE_START E40_N0 lane 2 2026-10-16T19:00:01Z\n'
     'TILE_END E0_N0 lane 1 rc 1 2026-10-16T19:05:00Z\n', ['_stderr.txt'], {'E0_N0', 'E40_N0'}),
    # only the runner's chatter, or nothing at all: nobody knows what ran
    ('runner chatter\n', ['_stdout.txt'], None),
    ('', [], None)])
def test_started_tiles(cj, monkeypatch, text, names, want):
    """plan_pack_tiles K10: the tiles a job started, from its TILE_START lines."""
    monkeypatch.setattr(cj, 'read_logs', lambda result: (text, 's3://x', names))
    fields, _ = cj.collect(Maap(), PACKED)
    assert fields['_started'] == want


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
