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


def test_a_failed_call_is_retried(capsys):
    calls, naps = [], []

    def flaky():
        calls.append(1)
        if len(calls) < 3:
            raise ConnectionError('timed out')
        return RESPONSE
    assert wc.fetch(call=flaky, attempts=5, pause=10, sleep=naps.append) is RESPONSE
    assert len(calls) == 3 and naps == [10, 10]
    assert 'attempt 2 of 5 failed (ConnectionError: timed out)' in capsys.readouterr().err


def test_after_the_last_attempt_the_error_names_the_call(capsys):
    def dead():
        raise ConnectionError('timed out')
    with pytest.raises(wc.CredentialError, match=r'workspace_bucket_credentials\(\) failed 3 times;'
                                                 ' last error: ConnectionError: timed out'):
        wc.fetch(call=dead, attempts=3, pause=0, sleep=lambda s: None)


def test_a_call_that_never_answers_is_cut_off():
    import time

    def hangs():
        time.sleep(5)
    with pytest.raises(wc.CredentialError, match='TimeoutError: no answer'):
        wc.fetch(call=hangs, attempts=1, timeout=1)


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
