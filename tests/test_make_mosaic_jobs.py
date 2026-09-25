"""
Tests for ATL1415/scripts/make_mosaic_jobs.py.

The z0 tasks must pass -w.  Without it make_mosaic.py ignores -p/-f and
overlapping tiles resolve in glob (directory, i.e. fetch) order: on IS
2026-09-25 the same 28 tiles, fetched in a different order, gave z0 mosaics
585 m apart on ice (plan_rerun_timing QT6).
"""
import importlib.util
import os

_HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    'make_mosaic_jobs', os.path.join(_HERE, '..', 'ATL1415', 'scripts', 'make_mosaic_jobs.py'))
mmj = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mmj)


def test_z0_tasks_are_weighted(tmp_path):
    base = tmp_path / 'IS'
    base.mkdir()
    run, _ = mmj.make_mosaic_jobs(str(base), 'IS', [1], run_name=str(tmp_path / 'run'))
    lines = [line for line in open(os.path.join(run, 'queue', 'task_1'))
             if line.startswith('make_mosaic.py') and 'z0.h5' in line]
    # six matched fields and sigma_z0 from prelim
    assert len(lines) == 7
    for line in lines:
        assert ' -w ' in line, line
