"""
Tests for ATL1415/scripts/check_mosaic_outputs.py.

Every case builds a small synthetic run directory under tmp_path, so nothing
reads or writes a real ATL14_processing tree.  The script is loaded from its
file rather than through the ATL1415 package, so only h5py and numpy are needed.
"""
import importlib.util
import io
import os

import h5py
import numpy as np
import pytest

_HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    'check_mosaic_outputs',
    os.path.join(_HERE, '..', 'ATL1415', 'scripts', 'check_mosaic_outputs.py'))
checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(checker)


def make_run(tmp_path, tasks, state='done'):
    run = tmp_path / 'run'
    for sub in ('queue', 'running', 'done', 'logs', 'error_logs'):
        (run / sub).mkdir(parents=True, exist_ok=True)
    for number, lines in tasks.items():
        (run / state / f'task_{number}').write_text(
            'source activate ATL14;\n' + '\n'.join(lines) + '\n')
    return run


def write_h5(path, fields, group='dz', value=1.0):
    with h5py.File(path, 'a') as h5:
        for field in fields:
            data = np.full((4, 5, 3), value)
            h5.create_dataset(f'{group}/{field}', data=data)


def check(run, values=False, quiet=False):
    out = io.StringIO()
    problems = checker.check_run(str(run), values=values, quiet=quiet, out=out)
    return problems, out.getvalue()


def test_good_run_passes(tmp_path):
    base = tmp_path / 'IS'
    base.mkdir()
    write_h5(base / 'dz.h5', ['dz', 'cell_area', 'sigma_dz'])
    run = make_run(tmp_path, {1: [
        f"make_mosaic.py  -R -w -d {base} -g 'matched/*.h5' -p 5000 -f 10000 "
        f"-O {base}/dz.h5 --in_group dz/ -F dz cell_area",
        f"make_mosaic.py  -w -d {base} -g 'prelim/*.h5' -p 5000 -f 10000 "
        f"-O {base}/dz.h5 --in_group dz/ -F sigma_dz",
    ]})
    problems, text = check(run)
    assert problems == 0, text
    assert '1 output files, 3 fields checked' in text


def test_missing_sigma_on_non_last_line_is_caught(tmp_path):
    """The failure exit codes hide: the sigma line fails, the task still exits 0."""
    base = tmp_path / 'IS'
    base.mkdir()
    write_h5(base / 'dz.h5', ['dz', 'cell_area'])
    run = make_run(tmp_path, {1: [
        f"make_mosaic.py -R -d {base} -O {base}/dz.h5 --in_group dz/ -F dz cell_area",
        f"make_mosaic.py -d {base} -O {base}/dz.h5 --in_group dz/ -F sigma_dz",
        f"make_mosaic.py -d {base} -O {base}/dz.h5 --in_group dz/ -F cell_area",
    ]})
    problems, text = check(run)
    assert problems == 1
    assert 'MISSING FIELD' in text and 'dz/sigma_dz' in text


def test_missing_file_counts_every_field(tmp_path):
    run = make_run(tmp_path, {1: [
        f"make_mosaic.py -R -d {tmp_path} -O {tmp_path}/z0.h5 --in_group z0/ -F z0 mask"]})
    problems, text = check(run)
    assert problems == 2
    assert 'MISSING FILE' in text


def test_all_nan_field_is_a_problem_with_values(tmp_path):
    write_h5(tmp_path / 'dz.h5', ['dz'], value=np.nan)
    run = make_run(tmp_path, {1: [
        f"make_mosaic.py -R -d {tmp_path} -O {tmp_path}/dz.h5 --in_group dz/ -F dz"]})
    problems, text = check(run, values=True)
    assert problems == 1
    assert 'ALL NAN' in text


def test_default_does_not_read_values(tmp_path, monkeypatch):
    """Metadata only by default: an all-NaN field passes, and no data is read."""
    write_h5(tmp_path / 'dz.h5', ['dz'], value=np.nan)
    run = make_run(tmp_path, {1: [
        f"make_mosaic.py -R -d {tmp_path} -O {tmp_path}/dz.h5 --in_group dz/ -F dz"]})

    def fail(dataset):
        raise AssertionError('finite_fraction called without --values')
    monkeypatch.setattr(checker, 'finite_fraction', fail)
    problems, text = check(run)
    assert problems == 0, text
    assert 'metadata only' in text


def test_relative_output_resolves_against_directory(tmp_path):
    """make_mosaic.py writes os.path.join(-d, -O); the checker must look there."""
    base = tmp_path / 'IS'
    base.mkdir()
    write_h5(base / 'dz.h5', ['dz'])
    run = make_run(tmp_path, {1: [
        f"make_mosaic.py -R -d {base} -O dz.h5 --in_group dz/ -F dz"]})
    problems, text = check(run)
    assert problems == 0, text


def test_out_group_overrides_in_group(tmp_path):
    write_h5(tmp_path / 'out.h5', ['dz'], group='renamed')
    run = make_run(tmp_path, {1: [
        f"make_mosaic.py -R -d {tmp_path} -O {tmp_path}/out.h5 "
        f"--in_group dz/ --out_group renamed/ -F dz"]})
    problems, text = check(run)
    assert problems == 0, text


def test_unfinished_tasks_and_error_logs_are_problems(tmp_path):
    write_h5(tmp_path / 'dz.h5', ['dz'])
    line = f"make_mosaic.py -R -d {tmp_path} -O {tmp_path}/dz.h5 --in_group dz/ -F dz"
    run = make_run(tmp_path, {1: [line]}, state='running')
    (run / 'queue' / 'task_2').write_text(line + '\n')
    (run / 'error_logs' / 'task_1.log').write_text('exit code=1\n')
    problems, text = check(run)
    assert problems == 3
    assert 'NOT DONE      running/task_1' in text
    assert 'NOT DONE      queue/task_2' in text
    assert 'ERROR LOG' in text


def test_empty_run_is_a_problem(tmp_path):
    run = make_run(tmp_path, {})
    problems, text = check(run)
    assert problems == 1
    assert 'NO TASK FILES' in text


def test_finite_fraction_reads_in_blocks(tmp_path, monkeypatch):
    monkeypatch.setattr(checker, 'BLOCK_ROWS', 3)
    data = np.ones((10, 4, 2))
    data[7:, :, :] = np.nan
    with h5py.File(tmp_path / 'f.h5', 'w') as h5:
        h5['x'] = data
        assert checker.finite_fraction(h5['x']) == pytest.approx(0.7)


def test_quiet_prints_only_problems(tmp_path):
    write_h5(tmp_path / 'dz.h5', ['dz'])
    run = make_run(tmp_path, {1: [
        f"make_mosaic.py -R -d {tmp_path} -O {tmp_path}/dz.h5 --in_group dz/ -F dz"]})
    problems, text = check(run, quiet=True)
    assert problems == 0
    assert not any(line.startswith('ok ') for line in text.splitlines())
