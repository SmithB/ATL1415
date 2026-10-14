"""drive_packed.py: per-tile decisions, rolling resubmission, restart
(plan_pack_tiles K5c, K9).  No MAAP, no S3: submission, collect and the S3
listing are stubbed."""
import csv
import datetime
import importlib.util
import os
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, '..', 'scripts', 'maap'))
_spec = importlib.util.spec_from_file_location(
    'drive_packed', os.path.join(HERE, '..', 'scripts', 'maap', 'drive_packed.py'))
dp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dp)

NOW = datetime.datetime(2026, 10, 14, 12, 0, tzinfo=datetime.timezone.utc)


def test_ledger_tiles():
    assert dp.ledger_tiles({'tiles': '160000,-1640000;160000,-1600000|200000,-1640000'}) == \
        ['E160_N-1640', 'E160_N-1600', 'E200_N-1640']
    assert dp.ledger_tiles({'tiles': '-'}) == [] and dp.ledger_tiles({}) == []


def test_pack_tiles_of(tmp_path):
    p = tmp_path / 'p.txt'
    p.write_text('E0_N0.h5 E0_N40.h5 | E40_N0.h5\nE80_N0.h5\n')
    assert dp.pack_tiles_of(str(p)) == ['E0_N0', 'E0_N40', 'E40_N0', 'E80_N0']


def test_decide_tile():
    row = {'job_id': 'j', 'submitted_utc': NOW.isoformat()}
    new, old = NOW + datetime.timedelta(seconds=5), NOW - datetime.timedelta(hours=1)
    assert dp.decide_tile('E0_N0', row, 'ok', {'E0_N0.h5': new}, 'prelim') == 'done'
    # a tile written before THIS job was submitted is an earlier attempt's
    assert dp.decide_tile('E0_N0', row, 'ok', {'E0_N0.h5': old}, 'prelim') is None
    assert dp.decide_tile('E0_N0', row, 'nodata', {}, 'prelim') == 'nodata'
    # matched has no nodata: a matched job that writes no tile failed
    assert dp.decide_tile('E0_N0', row, 'nodata', {}, 'matched') is None
    failed = {'job_id': '<submit failed: HTTP 502>', 'submitted_utc': NOW.isoformat()}
    assert dp.decide_tile('E0_N0', failed, '-', {'E0_N0.h5': new}, 'prelim') is None


def args(tmp_path, *extra):
    return dp.parse(['--step', 'prelim', '--tile_list', 'x', '--args_url', 's3://a/args.txt',
                     '--tile_prefix', 's3://a/t', '--queue', 'q', '--tag', 'T',
                     '--stem', str(tmp_path / 'run'), '--predict', 'p.csv', *extra])


class Fake:
    """DPS stand-in: submit writes a ledger (ids j1, j2, ...) at self.now;
    collect reports self.status per job id; S3 holds self.h5 {name: mtime}."""

    def __init__(self, monkeypatch):
        self.n, self.now, self.status, self.h5, self.batches = 0, NOW, {}, {}, []
        monkeypatch.setattr(dp.subprocess, 'run', self.submit)
        monkeypatch.setattr(dp, 'collect_rows', self.collect)
        monkeypatch.setattr(dp, 's3_mtimes', lambda prefix: dict(self.h5))

    def submit(self, cmd, **kw):
        pack, ledger = cmd[cmd.index('--pack_file') + 1], cmd[cmd.index('--ledger') + 1]
        lines = [line for line in open(pack) if line.strip()]
        self.batches.append(lines)
        xy = lambda t: '%d,%d' % (int(t.split('_N')[0][1:]) * 1000, int(t.split('_N')[1]) * 1000)
        with open(ledger, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=('identifier', 'job_id', 'submitted_utc', 'tiles'))
            w.writeheader()
            for line in lines:
                self.n += 1
                tiles = '|'.join(';'.join(xy(t) for t in lane.split())
                                 for lane in line.replace('.h5', '').split('|'))
                w.writerow({'identifier': f'T_j{self.n}', 'job_id': f'j{self.n}',
                            'submitted_utc': self.now.isoformat(), 'tiles': tiles})

        class P:
            returncode, stdout, stderr = 0, 'ok\n', ''
        return P()

    def collect(self, maap, rows):
        return [[{'tile': t, 'status': self.status.get(r['job_id'], 'ok'), 'step': 'prelim',
                  'identifier': r['identifier'], 'job_id': r['job_id']}
                 for t in dp.ledger_tiles(r)] for r in rows]

    def write(self, tile, ago=-5):
        self.h5[f'{tile}.h5'] = self.now - datetime.timedelta(seconds=ago)


def plan(tmp_path, text):
    (tmp_path / 'run_plan_pack.txt').write_text(text)


def decisions(tmp_path):
    return [(r['tile'], r['attempt'], r['decision'], r['shared'])
            for r in csv.DictReader(open(tmp_path / 'run_decisions.csv'))]


def test_failures_are_resubmitted_while_other_jobs_run(tmp_path, monkeypatch):
    a = args(tmp_path, '--max_in_flight', '2', '--max_retries', '1')
    plan(tmp_path, 'E0_N0.h5 | E40_N0.h5\nE80_N0.h5\nE120_N0.h5\n')
    pred = {t: (10., 500.) for t in ('E0_N0', 'E40_N0', 'E80_N0', 'E120_N0')}
    f = Fake(monkeypatch)
    run = dp.Run(a, pred)

    assert run.tick({}) == (False, True)                       # fills the 2 slots
    assert f.batches[0] == ['E0_N0.h5 | E40_N0.h5\n', 'E80_N0.h5\n']

    # j1 (2 lanes) fails while j2 runs: its tiles go next, ALONE (one lane),
    # ahead of the plan's E120_N0, in the one free slot
    f.now += datetime.timedelta(minutes=10)
    assert run.tick({'j1': 'failed', 'j2': 'running'}) == (False, True)
    assert f.batches[1] == ['E0_N0.h5 E40_N0.h5\n']
    assert decisions(tmp_path) == [('E0_N0', '1', 'retry', '1'), ('E40_N0', '1', 'retry', '1')]

    # a restart rebuilds the state from the files: nothing is sent twice
    assert dp.Run(a, pred).tick({'j1': 'failed', 'j2': 'running', 'j3': 'running'}) == (False, False)
    assert len(f.batches) == 2

    f.now += datetime.timedelta(minutes=10)
    f.write('E80_N0')
    assert run.tick({'j1': 'failed', 'j2': 'successful', 'j3': 'running'}) == (False, True)
    assert f.batches[2] == ['E120_N0.h5\n']

    # j3: E40_N0 written, E0_N0 not -> out of retries.  j4: E120_N0's file is
    # older than j4's submission (an earlier run's) -> not done, retried
    f.now += datetime.timedelta(minutes=10)
    f.write('E40_N0')
    f.write('E120_N0', ago=3600)
    assert run.tick({'j1': 'failed', 'j2': 'successful', 'j3': 'failed',
                     'j4': 'successful'}) == (False, True)
    assert f.batches[3] == ['E120_N0.h5\n']
    assert decisions(tmp_path)[-3:] == [('E0_N0', '2', 'gave_up', '0'), ('E40_N0', '2', 'done', '0'),
                                        ('E120_N0', '1', 'retry', '0')]

    f.now += datetime.timedelta(minutes=10)
    f.write('E120_N0')
    states = {'j1': 'failed', 'j2': 'successful', 'j3': 'failed', 'j4': 'successful',
              'j5': 'successful'}
    assert run.tick(states) == (True, True)
    assert open(tmp_path / 'run_still_failed_tile_list.txt').read() == 'E0_N0.h5\n'
    # the collect rows are kept as a --predict table for later runs
    assert len(list(csv.DictReader(open(tmp_path / 'run_tiles.csv')))) == 7


def test_nodata_and_submit_failures(tmp_path, monkeypatch):
    a = args(tmp_path, '--max_in_flight', '5')
    plan(tmp_path, 'E0_N0.h5\nE40_N0.h5\n')
    pred = {'E0_N0': (10., 500.), 'E40_N0': (10., 500.)}
    f = Fake(monkeypatch)
    run = dp.Run(a, pred)
    run.tick({})
    # E0_N0 has no data; E40_N0's submission failed (no job id): resubmitted at once
    ledger = tmp_path / 'run_b0_jobs.csv'
    rows = list(csv.DictReader(open(ledger)))
    rows[1]['job_id'] = '<submit failed: HTTP 502>'
    with open(ledger, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    f.status['j1'] = 'nodata'
    assert run.tick({'j1': 'successful'}) == (False, True)
    assert f.batches[1] == ['E40_N0.h5\n']
    assert decisions(tmp_path) == [('E0_N0', '1', 'nodata', '0'), ('E40_N0', '1', 'retry', '0')]


def test_a_half_written_decision_is_decided_again(tmp_path, monkeypatch):
    a = args(tmp_path)
    plan(tmp_path, 'E0_N0.h5 | E40_N0.h5\n')
    pred = {'E0_N0': (10., 500.), 'E40_N0': (10., 500.)}
    f = Fake(monkeypatch)
    run = dp.Run(a, pred)
    run.tick({})
    f.write('E0_N0')
    f.write('E40_N0')
    dp.append_csv(run.decisions_file, dp.DECISION_COLUMNS,
                  [{'job': 'run_b0_jobs.csv#0', 'tile': 'E0_N0', 'decision': 'done', 'shared': 1}])
    assert run.tick({'j1': 'successful'}) == (True, True)
    assert [d[:3] for d in decisions(tmp_path)][-2:] == [('E0_N0', '1', 'done'), ('E40_N0', '1', 'done')]


def test_a_round_driver_stem_is_refused(tmp_path):
    (tmp_path / 'run_r0_jobs.csv').write_text('identifier\n')
    with pytest.raises(SystemExit):
        dp.main(['--step', 'prelim', '--tile_list', 'x', '--args_url', 's3://a/x',
                 '--tile_prefix', 's3://a/t', '--queue', 'q', '--tag', 'T',
                 '--stem', str(tmp_path / 'run'), '--predict', 'p.csv'])


@pytest.mark.parametrize('extra', [['--big_gib', '14'], ['--default_secs', '600']])
def test_bad_options(tmp_path, extra):
    with pytest.raises(SystemExit) as e:
        args(tmp_path, *extra)
    assert e.value.code == 2


def test_dry_run_packs_and_submits_nothing(tmp_path, capsys):
    lst = tmp_path / 'list.txt'
    lst.write_text('E0_N0.h5\nE40_N0.h5\nE80_N0.h5\n')
    pc = tmp_path / 'p.csv'
    with open(pc, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=('tile', 'step', 'status', 'solve_secs',
                                          'max_mem_GiB', 'lane', 'identifier'))
        w.writeheader()
        for t in ('E0_N0', 'E40_N0'):
            w.writerow({'tile': t, 'step': 'prelim', 'status': 'ok', 'solve_secs': '500',
                        'max_mem_GiB': '10', 'lane': '-', 'identifier': t})
    stem = tmp_path / 'run'
    base = ['--step', 'prelim', '--tile_list', str(lst), '--args_url', 's3://a/x',
            '--tile_prefix', 's3://a/t', '--queue', 'q', '--tag', 'T', '--stem', str(stem),
            '--predict', str(pc), '--dry-run']
    with pytest.raises(SystemExit):          # E80_N0 has no prediction
        dp.main(base)
    assert 'E80_N0' in capsys.readouterr().out
    dp.main(base + ['--default_gib', '10', '--default_secs', '500'])
    assert sorted(dp.pack_tiles_of(f'{stem}_dryrun_pack.txt')) == ['E0_N0', 'E40_N0', 'E80_N0']
    assert not os.path.exists(f'{stem}_plan_pack.txt')
    assert not os.path.exists(f'{stem}_b0_jobs.csv')


def test_tiles_that_shared_a_node_run_alone_next(tmp_path):
    # as the driver sets them: a peak of node_gib makes the packer put them in one lane
    import pack_tiles
    pred = {'E0_N0': (26., 500.), 'E40_N0': (10., 500.), 'E120_N0': (10., 500.)}
    jobs, _, big = pack_tiles.pack(sorted(pred), pred, big_gib=13., big_secs=1e9)
    assert big == ['E0_N0'] and jobs[0] == [['E0_N0']]
