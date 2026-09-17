"""
Tests for scripts/check_field_sizes.py.

The failure cases build synthetic step directories under tmp_path.  The tile
files are empty placeholders: the checker pairs tiles with reports by name and
never opens a tile.  The acceptance test reads the real IS reports when they
are on disk, and only reads them -- nothing is written into ATL14_processing.
"""
import importlib.util
import json
import os

import pytest

_HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    'check_field_sizes', os.path.join(_HERE, '..', 'scripts', 'check_field_sizes.py'))
checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(checker)

SHAPE = [61, 61, 32]
IS_DIR = os.path.expanduser('~/ATL14_processing/rel006/north/IS')

# the three lines that set the shape, among lines the checker must ignore
ARGS = """--ATL11_earthaccess
--tile_spacing=40000
-W=60000
-t=2018.75,2026.5
--t_crop=2019,2026.25
-g=100,1000,0.25
--Hemisphere=1
-b=/somewhere/IS
"""


def make_args(tmp_path, text=ARGS):
    path = tmp_path / 'input_args.txt'
    path.write_text(text)
    return str(path)


def shape_from(argv):
    args = checker.parse_args(['step_dir'] + argv)
    return checker.expected_shape(args.Width, args.time_span, args.grid_spacing)


def make_step(tmp_path, step, tiles):
    """A step dir with a tile and a well-formed report for each name."""
    step_dir = tmp_path / step
    (step_dir / 'field_sizes').mkdir(parents=True)
    for name in tiles:
        (step_dir / f'{name}.h5').write_bytes(b'')
        write_report(step_dir, name, SHAPE, SHAPE if step == 'prelim' else None)
    return step_dir


def write_report(step_dir, name, dz, sigma, **extra):
    report = {'file': f'/XyZabc/output/{name}.h5', 'dz/dz': dz, 'dz/sigma_dz': sigma}
    report.update(extra)
    (step_dir / 'field_sizes' / f'{name}_report.json').write_text(json.dumps(report))


def run(capsys, *argv):
    status = checker.main([str(a) for a in argv])
    out = capsys.readouterr()
    return status, out.out, out.err


@pytest.mark.parametrize('step', ['prelim', 'matched'])
def test_good_step_passes_and_prints_counts(tmp_path, capsys, step):
    step_dir = make_step(tmp_path, step, ['E1_N1', 'E2_N2'])
    status, out, _ = run(capsys, step_dir, '@' + make_args(tmp_path))
    assert status == 0
    assert f'step {step}' in out
    assert 'expected dz/dz [61, 61, 32]' in out
    assert '2 reports, 2 tiles, 2 of 2 passed, 0 problems' in out


def test_shape_comes_from_t_not_t_crop(tmp_path):
    shape, derivation = shape_from(['@' + make_args(tmp_path)])
    assert shape == SHAPE           # t_crop would give 30 epochs
    assert '-t=2018.75,2026.5' in derivation


def test_shape_follows_the_args_file(tmp_path):
    # the AA south half: -W 44000 at the same spacing
    shape, _ = shape_from(['@' + make_args(tmp_path, ARGS.replace('-W=60000', '-W=44000'))])
    assert shape == [45, 45, 32]


def test_long_option_names_are_read(tmp_path):
    text = '--Width=60000\n--time_span=2018.75,2026.5\n--grid_spacing=100,1000,0.25\n'
    assert shape_from(['@' + make_args(tmp_path, text)])[0] == SHAPE


def test_flags_can_be_given_directly(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'matched', ['E1_N1'])
    status, out, _ = run(capsys, step_dir, '-W=60000', '-g=100,1000,0.25', '-t=2018.75,2026.5')
    assert status == 0
    assert 'expected dz/dz [61, 61, 32]' in out


def test_nested_at_includes_are_followed(tmp_path):
    inner = tmp_path / 'inner.txt'
    inner.write_text('-g=100,1000,0.25\n')
    outer = make_args(tmp_path, ARGS.replace('-g=100,1000,0.25\n', f'@{inner}\n'))
    assert shape_from(['@' + outer])[0] == SHAPE


def test_a_typed_unknown_option_is_an_error_not_ignored(tmp_path, capsys):
    # the file's other lines are ignored; a typo on the command line must not be
    step_dir = make_step(tmp_path, 'prelim', ['E1_N1'])
    with pytest.raises(SystemExit) as e:
        checker.main([str(step_dir), '@' + make_args(tmp_path), '--stpe', 'matched'])
    assert e.value.code == 2
    assert '--stpe' in capsys.readouterr().err


def test_prelim_without_sigma_fails(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'prelim', ['E1_N1', 'E2_N2'])
    write_report(step_dir, 'E1_N1', SHAPE, None)
    status, out, _ = run(capsys, step_dir, '@' + make_args(tmp_path))
    assert status == 1
    assert 'PROBLEM E1_N1: prelim dz/sigma_dz is None' in out
    assert '1 of 2 passed, 1 problems' in out


def test_prelim_sigma_of_another_shape_fails(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'prelim', ['E1_N1'])
    write_report(step_dir, 'E1_N1', SHAPE, [61, 61, 30])
    status, out, _ = run(capsys, step_dir, '@' + make_args(tmp_path))
    assert status == 1
    assert 'prelim dz/sigma_dz is [61, 61, 30]' in out


def test_matched_with_sigma_fails(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'matched', ['E1_N1'])
    write_report(step_dir, 'E1_N1', SHAPE, SHAPE)
    status, out, _ = run(capsys, step_dir, '@' + make_args(tmp_path))
    assert status == 1
    assert 'matched dz/sigma_dz is [61, 61, 32], expected null' in out


def test_wrong_dz_shape_fails(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'matched', ['E1_N1'])
    write_report(step_dir, 'E1_N1', [61, 61, 30], None)
    status, out, _ = run(capsys, step_dir, '@' + make_args(tmp_path))
    assert status == 1
    assert 'dz/dz is [61, 61, 30], expected [61, 61, 32]' in out


def test_missing_dz_fails(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'matched', ['E1_N1'])
    write_report(step_dir, 'E1_N1', None, None)
    status, out, _ = run(capsys, step_dir, '@' + make_args(tmp_path))
    assert status == 1
    assert 'dz/dz is None' in out


def test_report_whose_tile_was_deleted_fails(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'prelim', ['E1_N1', 'E2_N2'])
    (step_dir / 'E2_N2.h5').unlink()
    status, out, _ = run(capsys, step_dir, '@' + make_args(tmp_path))
    assert status == 1
    assert 'PROBLEM E2_N2: report has no tile' in out
    assert '2 reports, 1 tiles, 1 of 2 passed' in out


def test_tile_without_report_fails(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'prelim', ['E1_N1', 'E2_N2'])
    (step_dir / 'field_sizes' / 'E2_N2_report.json').unlink()
    status, out, _ = run(capsys, step_dir, '@' + make_args(tmp_path))
    assert status == 1
    assert 'PROBLEM E2_N2: tile has no report' in out


def test_malformed_and_unexpected_reports_fail(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'matched', ['E1_N1', 'E2_N2', 'E3_N3'])
    (step_dir / 'field_sizes' / 'E1_N1_report.json').write_text('{not json')
    (step_dir / 'field_sizes' / 'E2_N2_report.json').write_text(json.dumps({'dz/dz': SHAPE}))
    write_report(step_dir, 'E3_N3', SHAPE, None, spare=1)
    status, out, _ = run(capsys, step_dir, '@' + make_args(tmp_path))
    assert status == 1
    assert 'PROBLEM E1_N1: cannot read' in out
    assert "PROBLEM E2_N2: missing keys ['dz/sigma_dz', 'file']" in out
    assert "PROBLEM E3_N3: unexpected keys ['spare']" in out
    assert '0 of 3 passed' in out


def test_report_file_value_is_never_used(tmp_path, capsys):
    # the worker path points nowhere on this machine, and that is fine
    step_dir = make_step(tmp_path, 'matched', ['E1_N1'])
    status, _, _ = run(capsys, step_dir, '@' + make_args(tmp_path))
    assert status == 0


def test_empty_field_sizes_is_not_a_pass(tmp_path, capsys):
    step_dir = tmp_path / 'prelim'
    (step_dir / 'field_sizes').mkdir(parents=True)
    (step_dir / 'E1_N1.h5').write_bytes(b'')
    status, out, err = run(capsys, step_dir, '@' + make_args(tmp_path))
    assert status == 2
    assert 'CHECK NOT DONE: no reports' in err and '1 tiles' in err
    assert out == ''


@pytest.mark.parametrize('case', ['no_dir', 'args_without_W', 'uneven_W'])
def test_cannot_check_exits_2(tmp_path, capsys, case):
    step_dir = make_step(tmp_path, 'prelim', ['E1_N1'])
    args_file = make_args(tmp_path)
    if case == 'no_dir':
        step_dir = tmp_path / 'nowhere' / 'prelim'
    elif case == 'args_without_W':
        args_file = make_args(tmp_path, ARGS.replace('-W=60000\n', ''))
    elif case == 'uneven_W':
        args_file = make_args(tmp_path, ARGS.replace('-W=60000', '-W=60500'))
    status, _, err = run(capsys, step_dir, '@' + args_file)
    assert status == 2
    assert 'CHECK NOT DONE' in err


def test_missing_args_file_exits_2(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'prelim', ['E1_N1'])
    with pytest.raises(SystemExit) as e:
        checker.main([str(step_dir), '@' + str(tmp_path / 'missing.txt')])
    assert e.value.code == 2


def test_no_shape_flags_at_all_exits_2(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'prelim', ['E1_N1'])
    status, _, err = run(capsys, step_dir)
    assert status == 2
    assert '-W and -t are both required' in err


def test_step_is_required_when_the_directory_name_does_not_say(tmp_path, capsys):
    step_dir = make_step(tmp_path, 'prelim', ['E1_N1'])
    renamed = tmp_path / 'tiles_copy'
    step_dir.rename(renamed)
    status, _, err = run(capsys, renamed, '@' + make_args(tmp_path))
    assert status == 2
    assert '--step' in err
    status, out, _ = run(capsys, renamed, '@' + make_args(tmp_path), '--step', 'prelim')
    assert status == 0


def test_explicit_step_overrides_the_directory_name(tmp_path, capsys):
    # a prelim directory checked as matched: every sigma is now a fault
    step_dir = make_step(tmp_path, 'prelim', ['E1_N1'])
    status, out, _ = run(capsys, step_dir, '@' + make_args(tmp_path), '--step', 'matched')
    assert status == 1
    assert 'expected null' in out


@pytest.mark.skipif(not os.path.isdir(os.path.join(IS_DIR, 'matched', 'field_sizes')),
                    reason='the IS run is not on this machine')
@pytest.mark.parametrize('step', ['prelim', 'matched'])
def test_the_real_IS_run_passes(capsys, step):
    # known good (plan_IS_run.sh I4, I8): any complaint means the checker is wrong
    status, out, _ = run(capsys, os.path.join(IS_DIR, step),
                         '@' + os.path.join(IS_DIR, 'input_args_IS.txt'))
    assert status == 0, out
    assert '28 reports, 28 tiles, 28 of 28 passed, 0 problems' in out
