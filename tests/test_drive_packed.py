"""drive_packed.py: per-tile outcome, restart and resume (plan_pack_tiles K5c).
No MAAP, no S3: submission is stubbed."""
import csv
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

LEDGER_COLS = ('identifier', 'step', 'job_id', 'tiles')


def test_ledger_tiles():
    assert dp.ledger_tiles({'tiles': '160000,-1640000;160000,-1600000|200000,-1640000'}) == \
        ['E160_N-1640', 'E160_N-1600', 'E200_N-1640']
    assert dp.ledger_tiles({'tiles': '-'}) == [] and dp.ledger_tiles({}) == []


def test_pack_tiles_of(tmp_path):
    p = tmp_path / 'p.txt'
    p.write_text('E0_N0.h5 E0_N40.h5 | E40_N0.h5\nE80_N0.h5\n')
    assert dp.pack_tiles_of(str(p)) == ['E0_N0', 'E0_N40', 'E40_N0', 'E80_N0']


def test_outcome_prelim_and_matched():
    tiles = ['E0_N0', 'E0_N40', 'E40_N0', 'E80_N0']
    rows = [{'tile': 'E0_N40', 'status': 'nodata'}, {'tile': 'E40_N0', 'status': 'failed'},
            {'tile': 'E80_N0', 'status': 'ok'}]
    # E80_N0 says ok but no new .h5 (an old file, or the upload failed): not done
    assert dp.outcome(tiles, {'E0_N0.h5'}, rows, 'prelim') == \
        (['E0_N0', 'E0_N40'], ['E40_N0', 'E80_N0'])
    # matched has no nodata: a matched job that writes no tile failed
    assert dp.outcome(tiles, {'E0_N0.h5'}, rows, 'matched') == \
        (['E0_N0'], ['E0_N40', 'E40_N0', 'E80_N0'])


def args(tmp_path, *extra):
    return dp.parse(['--step', 'prelim', '--tile_list', 'x', '--args_url', 's3://a/args.txt',
                     '--tile_prefix', 's3://a/t', '--queue', 'q', '--tag', 'T',
                     '--stem', str(tmp_path / 'run'), '--predict', 'p.csv', *extra])


def write_ledger(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=LEDGER_COLS)
        w.writeheader()
        w.writerows(rows)


def test_a_round_is_packed_and_submitted_once(tmp_path, monkeypatch):
    a = args(tmp_path)
    r = dp.Round(a, 0)
    sent = []

    def fake_submit(self, pack_file, ledger_stem, tag):
        sent.append((os.path.basename(pack_file), os.path.basename(ledger_stem), tag))
        write_ledger(f'{ledger_stem}_jobs.csv', [{'identifier': 'j', 'step': 'prelim',
                                                  'job_id': 'id1', 'tiles': '0,0|40000,0'}])
    monkeypatch.setattr(dp.Round, 'submit', fake_submit)
    pred = {'E0_N0': (10., 500.), 'E40_N0': (10., 500.)}
    r.run_submission(['E0_N0', 'E40_N0'], pred)
    assert sent == [('run_r0_pack.txt', 'run_r0', 'T_r0')]
    assert os.path.exists(f'{r.stem}_start.txt')
    # a restart: the ledger covers the pack file, nothing more is sent
    r.run_submission(['E0_N0', 'E40_N0'], pred)
    assert len(sent) == 1


def test_a_partly_submitted_round_resumes(tmp_path, monkeypatch):
    a = args(tmp_path)
    r = dp.Round(a, 1)
    open(r.pack, 'w').write('E0_N0.h5 | E40_N0.h5\nE80_N0.h5 | E120_N0.h5\n')
    write_ledger(f'{r.stem}_jobs.csv', [{'identifier': 'j', 'step': 'prelim', 'job_id': 'id1',
                                         'tiles': '0,0|40000,0'}])
    sent = []

    def fake_submit(self, pack_file, ledger_stem, tag):
        sent.append((open(pack_file).read(), os.path.basename(ledger_stem), tag))
        write_ledger(f'{ledger_stem}_jobs.csv', [])
    monkeypatch.setattr(dp.Round, 'submit', fake_submit)
    r.run_submission(None, None)
    assert sent == [('E80_N0.h5 | E120_N0.h5\n', 'run_r1_resume1', 'T_r1_resume1')]
    assert [os.path.basename(f) for f in r.ledgers()] == ['run_r1_jobs.csv',
                                                          'run_r1_resume1_jobs.csv']


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
    assert sorted(dp.pack_tiles_of(f'{stem}_r0_pack.txt')) == ['E0_N0', 'E40_N0', 'E80_N0']
    assert not os.path.exists(f'{stem}_r0_jobs.csv')
