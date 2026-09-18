"""
check_build_id.py refuses to treat an option as the args_file URL.

This script SUBMITS a DPS job.  Its positional args_url used to swallow
anything, so `check_build_id.py --help` submitted a real build_id job with
args_file='--help' on the default 32gb queue (job 2cfaac9e, 2026-09-18).
The guard must stop BEFORE load_config() or MAAP(), so these tests need no
credentials and no network.
"""
import importlib.util
import os

import pytest

_spec = importlib.util.spec_from_file_location(
    'check_build_id',
    os.path.join(os.path.dirname(__file__), '..', 'scripts', 'maap', 'check_build_id.py'))
checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(checker)


def run(monkeypatch, *argv):
    monkeypatch.setattr('sys.argv', ['check_build_id.py'] + list(argv))
    # if the guard ever stops firing first, these make the test fail loudly
    # rather than reach the network
    monkeypatch.setattr(checker, 'load_config',
                        lambda *a, **k: pytest.fail('load_config reached: the guard did not fire'))
    monkeypatch.setattr(checker, 'MAAP',
                        lambda *a, **k: pytest.fail('MAAP reached: the guard did not fire'))
    with pytest.raises(SystemExit) as caught:
        checker.main()
    return caught.value.code


@pytest.mark.parametrize('option', ['--help', '-h', '--bogus', '-x'])
def test_an_option_is_refused_and_nothing_is_submitted(monkeypatch, capsys, option):
    assert run(monkeypatch, option) == 2
    err = capsys.readouterr().err
    assert 'Nothing was submitted' in err
    assert repr(option) in err


def test_the_refusal_prints_the_usage(monkeypatch, capsys):
    run(monkeypatch, '--help')
    assert 'Usage:' in capsys.readouterr().err


def test_an_option_after_a_good_url_is_still_refused(monkeypatch, capsys):
    # --expect and --timeout are consumed earlier; a leftover flag is not
    assert run(monkeypatch, 's3://bucket/input_args_IS.txt', '--nope') == 2


def test_known_flags_are_consumed_and_do_not_trip_the_guard(monkeypatch):
    # --expect/--timeout/--dry-run are removed before the guard, so the guard
    # must NOT fire on them: it gets past and reaches load_config, which the
    # fixture turns into a clear failure -- so a SystemExit(2) here would mean
    # the guard wrongly rejected a valid invocation
    monkeypatch.setattr('sys.argv', ['check_build_id.py', 's3://b/a.txt',
                                     '--expect', 'abc1234', '--dry-run'])
    reached = []
    monkeypatch.setattr(checker, 'load_config',
                        lambda *a, **k: (reached.append(True),
                                         {'algorithm_name': 'n', 'algorithm_version': 'v'})[1])
    monkeypatch.setattr(checker, 'MAAP', lambda *a, **k: None)
    monkeypatch.setattr(checker, 'find_process', lambda *a, **k: {'processID': 1})
    monkeypatch.setattr(checker, 'cwl_facts', lambda *a, **k: ('abc1234', 'img'))
    monkeypatch.setattr(checker, 'origin_tip', lambda *a, **k: 'abc1234')
    checker.main()          # --dry-run returns without submitting
    assert reached, 'the guard rejected a valid invocation'
