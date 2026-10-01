#!/usr/bin/env python3
"""
workspace_credentials.py -- MAAP-brokered keys for the workspace bucket, for
run.sh to export (docs/plan_workspace_credentials.sh, design A).

WHY.  A job used to reach s3://maap-ops-workspace/<user>/ with the worker's
IAM role, through botocore's default credential chain.  MAAP is removing that
role (admin, 2026-10-01); the supported route is
maap.aws.workspace_bucket_credentials(), which returns temporary keys for the
user's own workspace.  run.sh asks ONCE, before its first bucket read, and
exports the keys as the standard AWS variables, which botocore, s3fs and GDAL
all read first -- so no other code changes and no process looks anything up
from the instance again.

  eval "$(workspace_credentials.py [--check <s3 uri>]...)"
      STDOUT: `export` lines for AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY,
      AWS_SESSION_TOKEN and ATL1415_WORKSPACE_CREDENTIALS_EXPIRE -- for the
      eval, never for a log.  STDERR: one line, the expiry and the paths the
      keys cover.  NO KEY IS EVER WRITTEN TO STDERR.
      --check <uri> (repeatable): exit 1 unless <uri> is under a read_write
      path, so a prefix the job cannot write fails here and not at the upload.
  workspace_credentials.py --report
      make the call and print `workspace_credentials=ok lifetime_h=N.N` or
      `workspace_credentials=FAILED (<why>)`; always exit 0 (run.sh --build-id).

A CALL THAT FAILS IS RETRIED, THEN THE JOB STOPS (Ben, QW2).  The MAAP API
times out when many jobs start together (plan_GL_north.sh NM5: 9 runner
timeouts and 61 NSIDC broker failures in 556 jobs), so the call is tried
ATTEMPTS times, PAUSE_S apart, each given TIMEOUT_S.  After that: exit 1 with
the call and the last error named.  There is no fallback to the worker role;
it would hide the failure until the day the role is removed.
"""
import datetime
import os
import shlex
import signal
import sys
import time

ATTEMPTS, PAUSE_S, TIMEOUT_S = 5, 10, 30
FIELDS = {'aws_access_key_id': 'AWS_ACCESS_KEY_ID',
          'aws_secret_access_key': 'AWS_SECRET_ACCESS_KEY',
          'aws_session_token': 'AWS_SESSION_TOKEN'}
EXPIRE_VAR = 'ATL1415_WORKSPACE_CREDENTIALS_EXPIRE'


class CredentialError(Exception):
    pass


def broker():
    """One call to MAAP's workspace-credentials endpoint."""
    from maap.maap import MAAP
    return MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org')
                ).aws.workspace_bucket_credentials()


def _timeout(signum, frame):
    raise TimeoutError(f'no answer in {TIMEOUT_S} s')


def fetch(call=broker, attempts=ATTEMPTS, pause=PAUSE_S, timeout=TIMEOUT_S, sleep=time.sleep):
    """The broker's response, after up to `attempts` tries; CredentialError if none works."""
    last = None
    for attempt in range(1, attempts + 1):
        # maap-py's requests.get has no timeout of its own, and a connect
        # that never answers took ~130 s in the failed jobs' logs
        signal.signal(signal.SIGALRM, _timeout)
        signal.alarm(timeout)
        try:
            return validate(call())
        except Exception as exc:        # any failure of the call is retried the same way
            last = f'{type(exc).__name__}: {exc}'
            print(f'workspace credentials: attempt {attempt} of {attempts} failed ({last})',
                  file=sys.stderr)
        finally:
            signal.alarm(0)
        if attempt < attempts:
            sleep(pause)
    raise CredentialError(f'maap.aws.workspace_bucket_credentials() failed {attempts} times;'
                          f' last error: {last}')


def validate(response):
    """The response, if it carries every field run.sh exports; otherwise CredentialError."""
    creds = response.get('credentials') if isinstance(response, dict) else None
    if not isinstance(creds, dict):
        raise CredentialError('the response has no "credentials" block')
    missing = [k for k in list(FIELDS) + ['expires_at'] if not creds.get(k)]
    if missing:
        raise CredentialError(f'the response is missing {", ".join(missing)}')
    parse_time(creds['expires_at'])
    return response


def parse_time(text):
    """expires_at ('2026-10-02T10:05:12+0000') as an aware datetime."""
    try:
        return datetime.datetime.strptime(text, '%Y-%m-%dT%H:%M:%S%z')
    except (TypeError, ValueError):
        raise CredentialError(f'cannot read the expiry time {text!r}')


def hours_left(expires_at, now=None):
    now = now or datetime.datetime.now(datetime.timezone.utc)
    return (parse_time(expires_at) - now).total_seconds() / 3600


def writable_paths(response):
    return [p['uri'].rstrip('/') for p in response.get('authorized_s3_paths') or []
            if p.get('access') == 'read_write' and p.get('uri')]


def check_uri(uri, response):
    """CredentialError unless `uri` is at or under a path the keys may write."""
    paths = writable_paths(response)
    u = uri.rstrip('/')
    if not any(u == p or u.startswith(p + '/') for p in paths):
        raise CredentialError(f'{uri} is not under a path these credentials may write:'
                              f' {", ".join(paths) or "(none)"}')


def export_lines(response):
    creds = response['credentials']
    lines = [f'export {var}={shlex.quote(creds[key])}' for key, var in FIELDS.items()]
    lines.append(f'export {EXPIRE_VAR}={shlex.quote(creds["expires_at"])}')
    return lines


def summary(response):
    """One line for the job log: when the keys expire and what they cover -- no key."""
    creds = response['credentials']
    paths = ', '.join(f'{p.get("uri")} ({p.get("access")})'
                      for p in response.get('authorized_s3_paths') or [])
    return (f'workspace credentials: from maap.aws.workspace_bucket_credentials(), expire'
            f' {creds["expires_at"]} ({hours_left(creds["expires_at"]):.1f} h); paths: {paths}')


def main(argv):
    if argv == ['--report']:
        try:
            r = fetch()
            print(f'workspace_credentials=ok lifetime_h={hours_left(r["credentials"]["expires_at"]):.1f}')
        except CredentialError as exc:
            print(f'workspace_credentials=FAILED ({exc})')
        return 0
    checks = []
    while argv:
        if argv[0] == '--check' and len(argv) > 1:
            checks.append(argv[1])
            argv = argv[2:]
        else:
            print(f'ERROR: workspace_credentials.py: unknown argument {argv[0]!r}', file=sys.stderr)
            return 2
    try:
        response = fetch()
        for uri in checks:
            check_uri(uri, response)
    except CredentialError as exc:
        print(f'ERROR: workspace credentials: {exc}', file=sys.stderr)
        return 1
    print(summary(response), file=sys.stderr)
    print('\n'.join(export_lines(response)))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
