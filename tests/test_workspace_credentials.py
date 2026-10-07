"""
scripts/workspace_credentials.py: the MAAP-brokered workspace keys run.sh
exports (docs/plan_workspace_credentials.sh W1).  The broker is a stand-in;
no test calls MAAP.
"""
import datetime
import importlib.util
import os

import pytest

HERE = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    'workspace_credentials', os.path.join(HERE, '..', 'scripts', 'workspace_credentials.py'))
wc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(wc)

# the shape measured on the ADE 2026-10-01 (the values are made up)
RESPONSE = {
    'credentials': {'aws_access_key_id': 'ASIAEXAMPLEKEYID0000',
                    'aws_secret_access_key': 'not/a+real$secret',
                    'aws_session_token': 'not-a-real-token==',
                    'expires_at': '2026-10-02T10:05:12+0000'},
    'authorized_s3_paths': [
        {'uri': 's3://maap-ops-workspace/ben_smith', 'access': 'read_write'},
        {'uri': 's3://maap-ops-workspace/shared/ben_smith', 'access': 'read_write'},
        {'uri': 's3://maap-ops-workspace/shared', 'access': 'read_only'},
        {'uri': 's3://maap-ops-workspace/dataset/triaged_job', 'access': 'read_only'}]}
SECRETS = [RESPONSE['credentials'][k] for k in wc.FIELDS]


@pytest.fixture
def brokered(monkeypatch):
    monkeypatch.setattr(wc, 'broker', lambda: RESPONSE)
    monkeypatch.setattr(wc, 'fetch', lambda: wc.validate(wc.broker()))


def test_exports_go_to_stdout_and_no_key_to_stderr(brokered, capsys):
    assert wc.main([]) == 0
    out, err = capsys.readouterr()
    assert out.splitlines() == [
        'export AWS_ACCESS_KEY_ID=ASIAEXAMPLEKEYID0000',
        "export AWS_SECRET_ACCESS_KEY='not/a+real$secret'",
        'export AWS_SESSION_TOKEN=not-a-real-token==',
        'export ATL1415_WORKSPACE_CREDENTIALS_EXPIRE=2026-10-02T10:05:12+0000']
    assert 'expire 2026-10-02T10:05:12+0000' in err
    assert 's3://maap-ops-workspace/ben_smith (read_write)' in err
    assert not any(secret in err for secret in SECRETS)


@pytest.mark.parametrize('field', list(wc.FIELDS) + ['expires_at'])
def test_a_missing_field_is_named(field):
    creds = {k: v for k, v in RESPONSE['credentials'].items() if k != field}
    with pytest.raises(wc.CredentialError, match=field):
        wc.validate({'credentials': creds})


def test_the_docstring_key_names_are_not_accepted():
    # maap-py's docstring says accessKeyId/...; the endpoint returns aws_* names
    with pytest.raises(wc.CredentialError, match='aws_access_key_id'):
        wc.validate({'credentials': {'accessKeyId': 'a', 'secretAccessKey': 'b',
                                     'sessionToken': 'c', 'expiration': 'd'}})


@pytest.mark.parametrize('uri', ['s3://maap-ops-workspace/ben_smith',
                                 's3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north/IS/',
                                 's3://maap-ops-workspace/shared/ben_smith/x'])
def test_check_accepts_a_writable_path(uri):
    wc.check_uri(uri, RESPONSE)


@pytest.mark.parametrize('uri', ['s3://maap-ops-workspace/ben_smithers/x',      # not a path boundary
                                 's3://maap-ops-workspace/shared/someone_else',  # read_only
                                 's3://maap-ops-workspace/dataset/triaged_job/x',
                                 's3://another-bucket/ben_smith'])
def test_check_refuses_what_the_keys_cannot_write(uri):
    with pytest.raises(wc.CredentialError, match='not under a path these credentials may write'):
        wc.check_uri(uri, RESPONSE)


def test_a_refused_prefix_stops_the_job_and_exports_nothing(brokered, capsys):
    assert wc.main(['--check', 's3://maap-ops-workspace/someone_else/out']) == 1
    out, err = capsys.readouterr()
    assert out == ''
    assert 'ERROR: workspace credentials: s3://maap-ops-workspace/someone_else/out' in err


class Clock:
    """A stand-in for time.monotonic that the stand-in sleep advances."""
    def __init__(self):
        self.t, self.naps = 0.0, []

    def __call__(self):
        return self.t

    def sleep(self, s):
        self.naps.append(s)
        self.t += s


def no_jitter(lo, hi):
    return 1.0


def test_a_failed_call_is_retried_after_growing_pauses(capsys):
    calls, clock = [], Clock()

    def flaky():
        calls.append(1)
        if len(calls) < 4:
            raise ConnectionError('timed out')
        return RESPONSE
    assert wc.fetch(call=flaky, pauses=(10, 20, 40, 60), budget=1000, sleep=clock.sleep,
                    clock=clock, rand=no_jitter) is RESPONSE
    assert len(calls) == 4 and clock.naps == [10, 20, 40]
    assert 'attempt 2 of 5 failed (ConnectionError: timed out)' in capsys.readouterr().err


def test_pauses_are_jittered_within_the_bounds():
    clock, draws = Clock(), []

    def rand(lo, hi):
        draws.append((lo, hi))
        return hi

    def dead():
        raise ConnectionError('refused')
    with pytest.raises(wc.CredentialError):
        wc.fetch(call=dead, pauses=(10, 20), jitter=0.5, budget=1000, sleep=clock.sleep,
                 clock=clock, rand=rand)
    assert draws == [(0.5, 1.5), (0.5, 1.5)] and clock.naps == [15, 30]


def test_no_pause_runs_past_the_budget(capsys):
    calls, clock = [], Clock()

    def dead():
        calls.append(1)
        raise ConnectionError('refused')
    with pytest.raises(wc.CredentialError, match='failed 3 times'):
        wc.fetch(call=dead, pauses=(10, 20, 40, 60), budget=50, sleep=clock.sleep,
                 clock=clock, rand=no_jitter)
    assert len(calls) == 3 and clock.naps == [10, 20]     # 30 s used; a 40 s pause would end at 70
    assert 'would pass the 50 s budget' in capsys.readouterr().err


def test_the_defaults_wait_at_most_the_budget():
    clock = Clock()

    def dead():
        raise ConnectionError('refused')
    with pytest.raises(wc.CredentialError):
        wc.fetch(call=dead, sleep=clock.sleep, clock=clock, rand=lambda lo, hi: hi)
    assert sum(clock.naps) <= wc.BUDGET_S


class Unauthorized(Exception):
    """Carries .response.status_code, as requests.HTTPError does."""
    def __init__(self):
        super().__init__('401 Client Error: UNAUTHORIZED')
        self.response = type('R', (), {'status_code': 401})()


def test_a_401_stops_at_once_and_names_the_token(capsys):
    calls, clock = [], Clock()

    def rejected():
        calls.append(1)
        raise Unauthorized()
    with pytest.raises(wc.CredentialError, match="HTTP 401: MAAP_PGT was rejected.*get_maap_pgt_token"):
        wc.fetch(call=rejected, sleep=clock.sleep, clock=clock)
    assert len(calls) == 1 and clock.naps == []


def test_another_http_error_is_retried():
    calls, clock = [], Clock()

    class Unavailable(Exception):
        response = type('R', (), {'status_code': 503})()

    def flaky():
        calls.append(1)
        if len(calls) < 2:
            raise Unavailable('503')
        return RESPONSE
    assert wc.fetch(call=flaky, sleep=clock.sleep, clock=clock, rand=no_jitter) is RESPONSE
    assert len(calls) == 2


def test_one_maap_client_is_kept_across_tries(monkeypatch):
    # building MAAP() is itself an API call; once one is built, later tries reuse it
    import sys
    import types
    built, asked = [], []

    class FakeMAAP:
        def __init__(self, maap_host):
            built.append(maap_host)
            if len(built) == 1:
                raise ConnectionError('environment/config refused')
            self.aws = self

        def workspace_bucket_credentials(self):
            asked.append(1)
            if len(asked) < 2:
                raise ConnectionError('timed out')
            return RESPONSE
    fake = types.ModuleType('maap.maap')
    fake.MAAP = FakeMAAP
    monkeypatch.setitem(sys.modules, 'maap', types.ModuleType('maap'))
    monkeypatch.setitem(sys.modules, 'maap.maap', fake)
    monkeypatch.setattr(wc, '_client', [])
    clock = Clock()
    assert wc.fetch(sleep=clock.sleep, clock=clock, rand=no_jitter) is RESPONSE
    assert len(built) == 2 and len(asked) == 2     # built again only after the failed build


def test_after_the_last_attempt_the_error_names_the_call(capsys):
    def dead():
        raise ConnectionError('timed out')
    with pytest.raises(wc.CredentialError, match=r'workspace_bucket_credentials\(\) failed 3 times;'
                                                 ' last error: ConnectionError: timed out'):
        wc.fetch(call=dead, pauses=(0, 0), sleep=lambda s: None)


def test_a_call_that_never_answers_is_cut_off():
    import time

    def hangs():
        time.sleep(5)
    with pytest.raises(wc.CredentialError, match='TimeoutError: no answer'):
        wc.fetch(call=hangs, pauses=(), timeout=1)


def test_main_stops_with_exit_1_when_the_broker_is_down(monkeypatch, capsys):
    def down():
        raise wc.CredentialError('maap.aws.workspace_bucket_credentials() failed 5 times; last error: x')
    monkeypatch.setattr(wc, 'fetch', down)
    assert wc.main([]) == 1
    out, err = capsys.readouterr()
    assert out == '' and err.startswith('ERROR: workspace credentials: maap.aws')


def test_hours_left():
    now = datetime.datetime(2026, 10, 1, 22, 5, 12, tzinfo=datetime.timezone.utc)
    assert wc.hours_left('2026-10-02T10:05:12+0000', now) == 12.0


def test_report_never_fails_the_build_id_job(monkeypatch, capsys):
    monkeypatch.setattr(wc, 'fetch', lambda: RESPONSE)
    assert wc.main(['--report']) == 0
    assert capsys.readouterr().out.startswith('workspace_credentials=ok lifetime_h=')

    def down():
        raise wc.CredentialError('no answer')
    monkeypatch.setattr(wc, 'fetch', down)
    assert wc.main(['--report']) == 0
    assert capsys.readouterr().out.strip() == 'workspace_credentials=FAILED (no answer)'
