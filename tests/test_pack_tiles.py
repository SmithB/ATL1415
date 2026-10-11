"""pack_tiles.py: tile list + collect_jobs --csv predictions -> pack file
(plan_pack_tiles K5b).  No MAAP, no S3."""
import csv
import importlib.util
import os

import pytest

HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    'pack_tiles', os.path.join(HERE, '..', 'scripts', 'maap', 'pack_tiles.py'))
pk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pk)
_spec = importlib.util.spec_from_file_location(
    'submit_MAAP_jobs', os.path.join(HERE, '..', 'scripts', 'maap', 'submit_MAAP_jobs.py'))
sub = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sub)

COLUMNS = ('tile', 'step', 'status', 'secs', 'solve_secs', 'max_mem_GiB',
           'instance', 'lane', 'identifier', 'job_id')


def write_csv(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, restval='-')
        w.writeheader()
        w.writerows(rows)
    return str(path)


def row(tile, gib, secs, step='prelim', status='ok', lane='-', ident=None):
    return {'tile': tile, 'step': step, 'status': status, 'solve_secs': str(secs),
            'max_mem_GiB': str(gib), 'lane': lane, 'identifier': ident or f'one_{tile}'}


def test_predictions_last_ok_row_wins_and_steps_kept_apart(tmp_path):
    a = write_csv(tmp_path / 'a.csv', [row('E0_N0', 10, 600), row('E0_N40', 9, 500),
                                       row('E0_N80', 9, 500, status='failed')])
    b = write_csv(tmp_path / 'b.csv', [row('E0_N0', 11, 700),
                                       row('E0_N40', 5, 100, step='matched'),
                                       row('E0_N40', 8, 400, status='nodata')])
    pred = pk.read_predictions([a, b], 'prelim', 1.5)
    assert pred == {'E0_N0': (11., 700.), 'E0_N40': (9., 500.)}
    assert pk.read_predictions([a, b], 'matched', 1.5) == {'E0_N40': (5., 100.)}


def test_two_lane_rows_are_divided_by_the_slowdown(tmp_path):
    a = write_csv(tmp_path / 'a.csv', [
        row('E0_N0', 10, 900, lane='1', ident='J2'), row('E0_N40', 10, 600, lane='2', ident='J2'),
        row('E0_N80', 10, 900, lane='1', ident='J1'), row('E0_N120', 10, 600, lane='1', ident='J1')])
    pred = pk.read_predictions([a], 'prelim', 1.5)
    assert pred['E0_N0'] == (10., 600.) and pred['E0_N40'] == (10., 400.)
    # a one-lane packed job: its times stand as they are
    assert pred['E0_N80'] == (10., 900.) and pred['E0_N120'] == (10., 600.)


def names(jobs):
    return sorted(t for job in jobs for lane in job for t in lane)


def test_big_tiles_get_one_lane_jobs_first():
    pred = {f'E{i}_N0': (10., 500.) for i in range(8)}
    pred['E100_N0'] = (14., 500.)     # largest
    pred['E200_N0'] = (10., 2000.)    # slowest
    jobs, big_secs, big = pk.pack(sorted(pred), pred, target_secs=3600., big_gib=13.,
                                  big_secs=1000., slowdown=1.5)
    assert sorted(big) == ['E100_N0', 'E200_N0']
    assert jobs[0] == [['E200_N0', 'E100_N0']]          # one lane, longest first
    assert all(len(job) == 2 for job in jobs[1:])
    assert names(jobs) == sorted(pred)


def test_lanes_respect_the_target_with_the_slowdown():
    pred = {f'E{i}_N0': (10., 600.) for i in range(12)}
    jobs, _, big = pk.pack(sorted(pred), pred, target_secs=3600., big_gib=13.,
                           big_secs=1e9, slowdown=1.5)
    assert big == []
    # 600 s x 1.5 = 900 s a tile: 4 a lane, 8 a job
    assert [sorted(map(len, job)) for job in jobs] == [[4, 4], [2, 2]]
    assert all(pk.job_wall(j, pred, 1.5) <= 3600. for j in jobs)
    assert all(pk.job_peak(j, pred) <= 26. for j in jobs)


def test_a_tile_longer_than_the_target_runs_alone():
    pred = {'E0_N0': (10., 5000.), 'E40_N0': (10., 100.)}
    jobs, _, big = pk.pack(sorted(pred), pred, target_secs=3600., big_gib=13., big_secs=4000.)
    assert jobs == [[['E0_N0']], [['E40_N0']]]


def test_big_secs_defaults_to_p90():
    pred = {f'E{i}_N0': (10., 100. * (i + 1)) for i in range(10)}
    _, big_secs, big = pk.pack(sorted(pred), pred)
    assert 900. < big_secs < 1000. and big == ['E9_N0']


def run_main(tmp_path, tiles, rows, *extra):
    lst = tmp_path / 'list.txt'
    lst.write_text(''.join(f'{t}.h5\n' for t in tiles))
    out = tmp_path / 'pack.txt'
    pk.main([str(lst), '--step', 'prelim', '--predict', write_csv(tmp_path / 'p.csv', rows),
             '--out', str(out), *extra])
    return out


def test_pack_file_reads_back(tmp_path, capsys):
    tiles = [f'E{40 * i}_N-1600' for i in range(10)]
    rows = [row(t, 10 + (i == 3) * 4, 500 + 10 * i) for i, t in enumerate(tiles)]
    out = run_main(tmp_path, tiles, rows)
    jobs = sub.read_pack_file(str(out))
    got = sorted(pk.tile_key(sub.tile_name(x, y)) for job in jobs for lane in job for x, y in lane)
    assert got == sorted(tiles)
    assert '10 prelim tiles ->' in capsys.readouterr().out


def test_unpredicted_tiles_are_refused_by_name(tmp_path):
    with pytest.raises(SystemExit) as e:
        run_main(tmp_path, ['E0_N0', 'E40_N0'], [row('E0_N0', 10, 500)])
    assert 'E40_N0' in str(e.value) and '--default_gib' in str(e.value)


def test_unpredicted_tiles_take_the_defaults(tmp_path):
    out = run_main(tmp_path, ['E0_N0', 'E40_N0'], [row('E0_N0', 10, 500)],
                   '--default_gib', '20', '--default_secs', '900')
    # 20 GiB > 13: a big tile, in its own one-lane job
    assert out.read_text().splitlines()[0] == 'E40_N0.h5'


@pytest.mark.parametrize('extra', [['--default_gib', '20'], ['--big_gib', '14']])
def test_bad_options(tmp_path, extra):
    with pytest.raises(SystemExit) as e:
        run_main(tmp_path, ['E0_N0'], [row('E0_N0', 10, 500)], *extra)
    assert e.value.code == 2


def test_duplicate_tiles_refused(tmp_path):
    with pytest.raises(SystemExit) as e:
        run_main(tmp_path, ['E0_N0', 'E0_N0'], [row('E0_N0', 10, 500)])
    assert 'more than once' in str(e.value)


def test_overhead_is_added_after_the_slowdown(tmp_path):
    a = write_csv(tmp_path / 'a.csv', [row('E0_N0', 10, 900, lane='1', ident='J'),
                                       row('E0_N40', 10, 600, lane='2', ident='J')])
    assert pk.read_predictions([a], 'prelim', 1.5, 25.)['E0_N0'] == (10., 625.)
