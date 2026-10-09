"""
submit_MAAP_jobs.py --tile_list, and matched skipping centers with no prelim tile.

docs/plan_tile_lists.sh TL2 (Ben's AM8: the resource lists drive prelim AND
matched submissions) and QT3 (matched skips, by name, rather than refusing).
No MAAP and no S3: the S3 listing is either stubbed or fed canned `aws s3 ls`
output, and argparse refuses before any network call.
"""
import importlib.util
import os
import subprocess
from types import SimpleNamespace

import pytest

HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    'submit_MAAP_jobs', os.path.join(HERE, '..', 'scripts', 'maap', 'submit_MAAP_jobs.py'))
sub = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sub)

IS_LIST = os.path.join(HERE, '..', 'ATL1415', 'resources', 'IS', '40km_tile_list.txt')
IS_MATCHED_XY = os.path.join(HERE, '..', 'region_files', 'IS_0332_monthly_matched_xy.txt')


def write(tmp_path, text, name='list.txt'):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


# --- reading the list ---------------------------------------------------------

def test_tile_names_become_meters(tmp_path):
    path = write(tmp_path, 'E1020_N-2420.h5\nE-1000_N0.h5\nE0_N-1040.h5\n')
    assert sub.read_tile_list(path) == [(1020000, -2420000), (-1000000, 0), (0, -1040000)]


def test_blank_lines_are_skipped(tmp_path):
    path = write(tmp_path, '\nE1020_N-2420.h5\n\n')
    assert sub.read_tile_list(path) == [(1020000, -2420000)]


@pytest.mark.parametrize('bad', ['field_sizes', 'E1020_N-2420', 'E1020_N-2420.nc',
                                 '1020000 -2420000', 'E10.5_N-24.h5'])
def test_a_line_that_is_not_a_tile_name_is_refused(tmp_path, capsys, bad):
    # 'field_sizes' is the real case: the directory name an `ls` of prelim/ adds
    path = write(tmp_path, f'E1020_N-2420.h5\n{bad}\n')
    with pytest.raises(SystemExit) as caught:
        sub.read_tile_list(path)
    assert caught.value.code == 2
    assert ':2:' in capsys.readouterr().err       # names the line


def test_names_round_trip_through_tile_name(tmp_path):
    names = ['E1020_N-2420.h5', 'E-1000_N0.h5', 'E-2700_N-1040.h5']
    path = write(tmp_path, '\n'.join(names) + '\n')
    assert [sub.tile_name(x, y) for x, y in sub.read_tile_list(path)] == names


def test_the_real_IS_list_is_the_28_centers_the_monthly_run_used():
    # ATL1415/resources/IS (Ben, 738bbd2) against the matched list built from
    # the tiles that existed after M7 -- the same 28, E1020_N-2580 absent
    from_list = set(sub.read_tile_list(IS_LIST))
    from_xy = set(sub.read_centers(IS_MATCHED_XY))
    assert len(from_list) == 28
    assert from_list == from_xy
    assert (1020000, -2580000) not in from_list


# --- matched: only centers with a prelim tile (QT3) ---------------------------

CENTERS = [(1020000, -2420000), (1140000, -2500000), (1020000, -2580000)]


def test_matched_skips_exactly_the_centers_without_a_prelim_tile():
    listed = {'E1020_N-2420.h5', 'E1140_N-2500.h5'}
    seen = []
    lister = lambda prefix: (seen.append(prefix), listed)[1]
    have, missing = sub.split_by_prelim(CENTERS, 's3://b/IS/', lister)
    assert have == [(1020000, -2420000), (1140000, -2500000)]    # order kept
    assert missing == ['E1020_N-2580.h5']
    assert seen == ['s3://b/IS/prelim']            # one listing, of prelim/


def test_matched_keeps_every_center_when_every_tile_is_there():
    listed = {sub.tile_name(x, y) for x, y in CENTERS}
    assert sub.split_by_prelim(CENTERS, 's3://b/IS', lambda p: listed) == (CENTERS, [])


def test_matched_with_no_prelim_tiles_at_all_keeps_nothing():
    assert sub.split_by_prelim(CENTERS, 's3://b/IS', lambda p: set()) == \
        ([], [sub.tile_name(x, y) for x, y in CENTERS])


# --- listing S3 --------------------------------------------------------------

def fake_run(returncode, stdout='', stderr=''):
    return lambda *a, **k: SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


def test_s3_names_reads_only_h5_keys(monkeypatch):
    out = ('                           PRE field_sizes/\n'
           '2026-09-18 17:32:03   60557867 E1340_N-2460.h5\n'
           '2026-09-18 17:30:00    4273533 E1180_N-2380.h5\n'
           '2026-09-18 17:30:00        101 notes.txt\n')
    monkeypatch.setattr(subprocess, 'run', fake_run(0, out))
    assert sub.s3_names('s3://b/IS/prelim') == {'E1340_N-2460.h5', 'E1180_N-2380.h5'}


def test_an_empty_prefix_is_an_empty_set(monkeypatch):
    # `aws s3 ls` exits 1 with no output when nothing matches
    monkeypatch.setattr(subprocess, 'run', fake_run(1))
    assert sub.s3_names('s3://b/IS/prelim') == set()


@pytest.mark.parametrize('rc', [1, 255])
def test_a_listing_that_fails_is_an_error_not_an_empty_set(monkeypatch, capsys, rc):
    # an empty set would read as "every prelim tile is missing"
    monkeypatch.setattr(subprocess, 'run', fake_run(rc, stderr='AccessDenied'))
    with pytest.raises(SystemExit) as caught:
        sub.s3_names('s3://b/IS/prelim')
    assert caught.value.code == 2
    assert 'AccessDenied' in capsys.readouterr().err


# --- the command line ---------------------------------------------------------

@pytest.mark.parametrize('argv', [
    ['--tile_list', 'a.txt', '--xy_file', 'b.txt'],     # both
    [],                                                  # neither
], ids=['both', 'neither'])
def test_exactly_one_of_tile_list_and_xy_file(monkeypatch, argv):
    monkeypatch.setattr('sys.argv', ['submit_MAAP_jobs.py', *argv,
                                     '--step', 'prelim', '--args_url', 's3://b/a.txt'])
    monkeypatch.setattr(sub, 'MAAP', lambda *a, **k: pytest.fail('reached MAAP'))
    with pytest.raises(SystemExit) as caught:
        sub.main()
    assert caught.value.code == 2


# --- the mosaic steps (docs/plan_dps_mosaic.sh D3-2) ----------------------------

GL_ARGS = '--region=GL\n-g=100,1000,0.25\n-t=2018.75,2026.5\n-W=60000\n--tile_spacing=40000\n'
PRELIM = ['E420_N20.h5', 'E460_N20.h5', 'E-620_N-1340.h5', 'E0_N-920.h5']


def test_mosaic200_tasks_are_the_region_list():
    # the canonical list, whatever the tile listing says: 81 GL cells, among
    # them the southern one no tile center is in (plan_200km_footprint.sh)
    tasks = sub.mosaic_tasks('mosaic200', GL_ARGS, 's3://b/GL', 's3://b/test/GL',
                             lister=lambda prefix: set(PRELIM))
    assert len(tasks) == 81 and '-100000_-3300000' in tasks


def test_the_gl_list_is_the_footprints_of_the_gl_tiles():
    from ATL1415.mosaic_groups import footprint_centers_200km, read_200km_centers
    names = [ln.strip() for ln in open(os.path.join(HERE, '..', 'ATL1415', 'resources', 'GL',
                                                      '40km_tile_list.txt')) if ln.strip()]
    assert read_200km_centers('GL') == footprint_centers_200km(names, half_width=30e3)


def test_the_centers_are_make_200km_tiles_own(tmp_path):
    # the submitter's (numpy-only) rule and the worker's must agree: the
    # worker refuses a --center that is not in its own list
    spec = importlib.util.spec_from_file_location(
        'make_200km_tiles', os.path.join(HERE, '..', 'ATL1415', 'scripts', 'make_200km_tiles.py'))
    m2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m2)
    worker = m2.make_200km_tiles(str(tmp_path), region='GL')
    tasks = sub.mosaic_tasks('mosaic200', GL_ARGS, 's3://b/GL', 's3://b/test/GL')
    assert tasks == [f'{int(x)}_{int(y)}' for x, y in worker]


@pytest.mark.parametrize('z0_tiles', [True, False])
def test_mosaic_tasks_are_the_groups(z0_tiles):
    tasks = sub.mosaic_tasks('mosaic', GL_ARGS, 's3://b/GL', 's3://b/test/GL',
                             exists=lambda uri: z0_tiles)
    # 9 lags: dz + 9 dzdt + 3 x (avg_dz + 9 avg_dzdt) = 40, + z0 where made
    assert len(tasks) == (41 if z0_tiles else 40)
    assert ('z0' in tasks) is z0_tiles
    assert {'dz', 'dzdt_lag28', 'avg_dz_40000m', 'avg_dzdt_10000m_lag1'} <= set(tasks)


IS_ARGS = GL_ARGS.replace('--region=GL', '--region=IS')


def test_only_greenland_and_antarctica_take_the_200km_step():
    # Ben 2026-10-01 (plan AD8)
    from ATL1415.mosaic_groups import uses_200km_tiles
    assert [r for r in ('GL', 'AA', 'A1', 'A4', 'IS', 'SV', 'CN', 'CS', 'RA', 'AK')
            if uses_200km_tiles(r)] == ['GL', 'AA', 'A1', 'A4']


def test_mosaic200_is_refused_for_a_region_without_the_step():
    with pytest.raises(ValueError, match='region IS has no 200 km step'):
        sub.mosaic_tasks('mosaic200', IS_ARGS, 's3://b/IS', 's3://b/test/IS',
                         lister=lambda prefix: set(PRELIM))


def test_direct_mosaic_tasks_do_not_wait_for_200km_z0():
    # no 200 km tiles to look for: z0 is mosaicked from the solve tiles
    def no_lookup(uri):
        raise AssertionError('looked for 200 km tiles in a region without them')
    tasks = sub.mosaic_tasks('mosaic', IS_ARGS, 's3://b/IS', 's3://b/test/IS', exists=no_lookup)
    assert len(tasks) == 41 and 'z0' in tasks


def test_mosaic_tasks_need_the_region():
    with pytest.raises(ValueError, match='--region='):
        sub.mosaic_tasks('mosaic', GL_ARGS.replace('--region=GL\n', ''), 's3://b/GL', 's3://b/t')


def test_nc_tasks():
    assert sub.mosaic_tasks('nc', GL_ARGS, 's3://b/GL', 's3://b/test/GL') == ['ATL14', 'ATL15']


def test_select_tasks_keeps_the_named_ones_in_step_order():
    tasks = ['z0', 'dz', 'dzdt_lag1', 'avg_dz_40000m']
    assert sub.select_tasks(tasks, ['avg_dz_40000m', 'z0']) == ['z0', 'avg_dz_40000m']


def test_select_tasks_refuses_a_name_the_step_does_not_have():
    with pytest.raises(ValueError, match='--task dz_lag1: not a task of this step'):
        sub.select_tasks(['z0', 'dz'], ['z0', 'dz_lag1'])


@pytest.mark.parametrize('argv, message', [
    (['--step', 'prelim', '--tile_list', 'x.txt', '--task', 'z0'], '--task is for the mosaic steps'),
    (['--step', 'mosaic', '--tile_list', 'x.txt', '--tile_prefix', 's3://b/GL'], 'lists its own jobs'),
    (['--step', 'nc'], 'needs --tile_prefix'),
    (['--step', 'prelim'], 'needs --tile_list, --xy_file or --pack_file'),
    (['--step', 'prelim', '--tile_list', 'x.txt', '--out_prefix', 's3://b/t'], 'for the mosaic steps'),
])
def test_mosaic_argument_refusals(argv, message, monkeypatch, capsys):
    monkeypatch.setattr(sub.sys, 'argv', ['submit_MAAP_jobs.py', '--args_url', 's3://b/a.txt'] + argv)
    with pytest.raises(SystemExit) as exit_info:
        sub.main()
    assert exit_info.value.code == 2
    assert message in capsys.readouterr().err


# --- packed jobs (plan_pack_tiles K4/K5) ---------------------------------------

def test_pack_file_jobs_lanes_and_tiles_input(tmp_path):
    path = write(tmp_path, 'E160_N-1640.h5 E160_N-1600.h5 | E200_N-1640.h5\n\n'
                           'E-240_N-960.h5\n')
    packs = sub.read_pack_file(path)
    assert packs == [[[(160000, -1640000), (160000, -1600000)], [(200000, -1640000)]],
                     [[(-240000, -960000)]]]
    assert sub.tiles_input(packs[0]) == '160000,-1640000;160000,-1600000|200000,-1640000'


@pytest.mark.parametrize('text, message', [
    ('E0_N0.h5 | | E40_N0.h5\n', 'an empty lane'),
    ('E0_N0.h5 E40_N0\n', 'not a tile name'),
    ('E0_N0.h5\nE40_N0.h5 | E0_N0.h5\n', 'already in line 1'),
], ids=['empty_lane', 'bad_name', 'duplicate_across_jobs'])
def test_a_bad_pack_file_is_refused(tmp_path, capsys, text, message):
    with pytest.raises(SystemExit) as caught:
        sub.read_pack_file(write(tmp_path, text))
    assert caught.value.code == 2 and message in capsys.readouterr().err


def test_matched_packs_drop_tiles_lanes_and_jobs_without_prelim():
    packs = [[[(0, 0), (40000, 0)], [(80000, 0)]], [[(120000, 0)]]]
    present = {'E0_N0.h5', 'E40_N0.h5'}
    kept, missing = sub.split_packs_by_prelim(packs, 's3://b/run', lister=lambda p: present)
    assert kept == [[[(0, 0), (40000, 0)]]]
    assert missing == ['E80_N0.h5', 'E120_N0.h5']


def test_pack_file_excludes_the_other_lists(monkeypatch):
    monkeypatch.setattr('sys.argv', ['submit_MAAP_jobs.py', '--pack_file', 'a.txt',
                                     '--tile_list', 'b.txt', '--step', 'prelim',
                                     '--args_url', 's3://b/a.txt'])
    monkeypatch.setattr(sub, 'MAAP', lambda *a, **k: pytest.fail('reached MAAP'))
    with pytest.raises(SystemExit) as caught:
        sub.main()
    assert caught.value.code == 2
