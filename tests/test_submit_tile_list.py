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
