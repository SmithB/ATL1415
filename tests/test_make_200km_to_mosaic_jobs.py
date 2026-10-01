"""
make_200km_to_mosaic_jobs.py: stage 2 of the 200 km path, joining a region's
200 km tiles into its mosaics.  On DPS one job runs one group (--group),
reading the 200 km tiles in place (--in_base s3://...; plan_dps_mosaic D3a-2).

Restructured 2026-10-01 around mosaic_commands(); checked then that it writes
the same 41 GL tasks as before, apart from the z0 task, which now activates
the environment like the others.
"""
import fsspec
import pytest
import pointCollection as pc

from ATL1415.scripts import make_200km_to_mosaic_jobs as m2m

LAGS = [1, 4]


@pytest.fixture
def region(tmp_path):
    base = tmp_path / 'GL'
    (base / '200km_tiles' / 'z0').mkdir(parents=True)
    return str(base)


def test_groups_in_task_order(region):
    groups = list(m2m.mosaic_commands(region, LAGS))
    assert groups == ['dz', 'dzdt_lag1', 'dzdt_lag4',
                      'avg_dz_40000m', 'avg_dzdt_40000m_lag1', 'avg_dzdt_40000m_lag4',
                      'avg_dz_20000m', 'avg_dzdt_20000m_lag1', 'avg_dzdt_20000m_lag4',
                      'avg_dz_10000m', 'avg_dzdt_10000m_lag1', 'avg_dzdt_10000m_lag4',
                      'z0']


def test_z0_only_where_its_200km_tiles_exist(tmp_path):
    base = tmp_path / 'GL'
    (base / '200km_tiles' / 'dz').mkdir(parents=True)
    assert 'z0' not in m2m.mosaic_commands(str(base), LAGS)


def test_first_command_of_a_group_creates_the_file(region):
    for commands in m2m.mosaic_commands(region, LAGS).values():
        assert ' -R ' in commands[0]
        assert not any(' -R ' in c for c in commands[1:])


@pytest.mark.parametrize('group, out_file, n_fields', [
    ('dz', 'dz.h5', 7), ('z0', 'z0.h5', 7), ('dzdt_lag4', 'dzdt_lag4.h5', 3),
    ('avg_dz_40000m', 'dz_40km.h5', 3), ('avg_dzdt_20000m_lag1', 'dzdt_20km_lag1.h5', 3)])
def test_output_names_match_make_mosaic_jobs(region, group, out_file, n_fields):
    commands = m2m.mosaic_commands(region, LAGS)[group]
    assert len(commands) == n_fields
    assert all(f'-O {region}/{out_file} ' in c for c in commands)
    assert all(f'-d {region}/200km_tiles/{group} ' in c for c in commands)


def test_remote_in_base_reads_in_place_and_writes_locally(region, monkeypatch):
    fs = fsspec.filesystem('memory')
    monkeypatch.setattr(pc.io_utils, 'get_s3fs', lambda daac=None, **kw: fs)
    prefix = 'memory://bucket/rel006_0332_testing/north/GL'
    fs.pipe(f'{prefix}/200km_tiles/z0/z00_200_-1000_-800.h5', b'')
    try:
        commands = m2m.mosaic_commands(region, LAGS, in_base=prefix, workers=4)
    finally:
        fs.rm('memory://bucket', recursive=True)
    assert 'z0' in commands          # found on the bucket, not locally
    for group, these in commands.items():
        assert all(f'-d {prefix}/200km_tiles/{group} ' in c for c in these)
        assert all(f'-O {region}/' in c and c.endswith(' -j 4') for c in these)


def test_one_group_one_task(region, tmp_path, monkeypatch):
    monkeypatch.setattr(m2m.ATL1415, 'make_slurm_file', lambda *args, **kwargs: None)
    monkeypatch.chdir(tmp_path)
    m2m.make_mosaic_jobs(region, 'GL', LAGS, group='avg_dzdt_40000m_lag4', environment='')
    tasks = sorted((tmp_path / 'mosaic_run_GL' / 'queue').iterdir())
    assert [t.name for t in tasks] == ['task_1']
    lines = tasks[0].read_text().splitlines()
    assert len(lines) == 3 and all(ln.startswith('make_mosaic.py') for ln in lines)


def test_an_unknown_group_is_refused(region, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit, match='no group'):
        m2m.make_mosaic_jobs(region, 'GL', LAGS, group='dzdt_lag3')
