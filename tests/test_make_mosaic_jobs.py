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


# --- one group, tiles read elsewhere: the direct path as a DPS job --------------
# (docs/plan_dps_mosaic.sh D3b-2)

def _task_files(run):
    return sorted(os.listdir(os.path.join(run, 'queue')))


def test_defaults_read_and_write_under_base(tmp_path):
    base = tmp_path / 'IS'
    base.mkdir()
    run, n = mmj.make_mosaic_jobs(str(base), 'IS', [1, 2], run_name=str(tmp_path / 'run'))
    # z0 + dz + 2 dzdt + 3 x (avg_dz + 2 avg_dzdt)
    assert n == 13 and len(_task_files(run)) == 13
    text = open(os.path.join(run, 'queue', 'task_1')).read()
    assert text.startswith('source activate IS2;')
    assert f'-d {base} ' in text and ' -j ' not in text


def test_one_group_is_task_1_with_tiles_read_from_tiles_base(tmp_path):
    base = tmp_path / 'IS'
    base.mkdir()
    run, n = mmj.make_mosaic_jobs(str(base), 'IS', [1, 2], run_name=str(tmp_path / 'run'),
                                  tiles_base='s3://b/IS', group='avg_dzdt_20000m_lag2',
                                  environment='', workers=4)
    assert n == 1 and _task_files(run) == ['task_1']
    lines = open(os.path.join(run, 'queue', 'task_1')).read().splitlines()
    # the matched line, then the prelim sigma line; no activate line
    assert len(lines) == 2 and all(line.startswith('make_mosaic.py') for line in lines)
    assert "-d s3://b/IS -g 'matched/*.h5'" in lines[0]
    assert "-d s3://b/IS -g 'prelim/*.h5'" in lines[1] and '-F sigma_avg_dzdt_20000m_lag2' in lines[1]
    for line in lines:
        assert f'-O {base}/dzdt_20km_lag2.h5' in line and line.endswith(' -j 4')


def test_the_groups_are_make_fields_groups(tmp_path):
    # the submitter lists a region's jobs from make_fields; each must be a
    # group make_mosaic_jobs knows
    from ATL1415.mosaic_groups import make_fields
    base = tmp_path / 'IS'
    base.mkdir()
    for group in make_fields([1, 2])[0]:
        _, n = mmj.make_mosaic_jobs(str(base), 'IS', [1, 2], run_name=str(tmp_path / group),
                                    group=group)
        assert n == 1 and _task_files(str(tmp_path / group)) == ['task_1']


def test_an_unknown_group_is_refused(tmp_path):
    import pytest
    base = tmp_path / 'IS'
    base.mkdir()
    with pytest.raises(SystemExit, match="no group 'dzdt_lag99'"):
        mmj.make_mosaic_jobs(str(base), 'IS', [1], run_name=str(tmp_path / 'run'), group='dzdt_lag99')
