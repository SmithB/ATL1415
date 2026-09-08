#!/usr/bin/env python3
"""
Register the ATL1415 DPS algorithm with MAAP and print the build URL.

This is docs/howto_MAAP_staging.sh S5 as a script.  It replaces the four-line
recipe that otherwise gets pasted into a Python REPL:

    from maap.maap import MAAP
    maap = MAAP(maap_host='api.maap-project.org')
    response = maap.register_algorithm_from_yaml_file('algorithm_config.yml')
    response.json()['message']['job_web_url']

Two things make that recipe worth wrapping.

FIRST, the response is opaque.  register_algorithm_from_yaml_file() returns a
raw requests.Response -- maap-py's registerAlgorithm() just hands back whatever
requests_utils.make_request() gave it -- and the build URL is not part of
maap-py at all.  It arrives inside the server's own JSON, under
['message']['job_web_url'], which is not documented anywhere and was originally
found by poking at the response in the REPL.

SECOND, and the reason for the git checks below: REGISTERING TRIGGERS A BUILD
THAT CLONES FROM GITHUB.  DPS clones repository_url at algorithm_version and
bakes the result into a container; the ADE working copy has nothing to do with
it.  Anything uncommitted, unpushed, or committed after the last build is
simply not on the worker, and no job log will ever say so -- a stale image
surfaces as whatever the missing fix was meant to prevent.  On 2026-09-04 that
cost a cycle: the build predated the as_gdal_path, pyTMD AWS_NO_SIGN_REQUEST and
Q27 W1/W3/W4/W5 fixes, and would have failed the Iceland smoke test looking
exactly like a credentials problem.  This script refuses to register in that
state; --force overrides.

Run it from the ADE, with the notebook env's interpreter (see the import guard
below -- maap-py is not installed in the ATL14 env):

    ./register_algorithm.py [--dry-run] [--force] [algorithm_config.yml]
"""

import argparse
import json
import os
import subprocess
import sys

# maap-py is installed ONLY in the ADE's default notebook env.  It is not in the
# ATL14 conda env that every howto starts by activating, and it is not a
# dependency of this package -- registration is a client-side operation, so
# there is no reason for the DPS image to carry it.
ADE_PYTHON = '/srv/conda/envs/notebook/bin/python'

# The repo the script lives in, so cwd never matters -- this is both where
# algorithm_config.yml is found and which checkout the git checks run against.
REPO_DIR = os.path.dirname(os.path.realpath(__file__))

DEFAULT_HOST = 'api.maap-project.org'


def import_maap():
    '''
    Import maap-py, or exit 1 naming the interpreter that has it.

    Deliberately does NOT re-exec itself under ADE_PYTHON: a script that
    silently switches interpreters hides which env is in play, and the whole
    point of this file is that what runs where is worth being explicit about.
    '''
    try:
        from maap.maap import MAAP
        from maap.utils import algorithm_utils
        return MAAP, algorithm_utils
    except ImportError as exc:
        print(f"ERROR: 'maap' is not importable under {sys.executable}", file=sys.stderr)
        print(f"       ({type(exc).__name__}: {exc})", file=sys.stderr)
        print(file=sys.stderr)
        print('maap-py is in the ADE notebook env, and -- since 2026-09-08 --'
              ' in any ATL14', file=sys.stderr)
        print('env rebuilt from environment.yml.  An ATL14 env built before'
              ' that predates it,', file=sys.stderr)
        print('which is the usual cause.  The notebook env always has it:',
              file=sys.stderr)
        print(file=sys.stderr)
        script = os.path.relpath(os.path.realpath(__file__))
        if os.path.exists(ADE_PYTHON):
            print(f'    {ADE_PYTHON} {script}', file=sys.stderr)
        else:
            # Do not print a path that is not there -- on a fresh account or a
            # different image the notebook env lives somewhere else.
            print(f'    <notebook-env>/bin/python {script}', file=sys.stderr)
            print(file=sys.stderr)
            print(f'{ADE_PYTHON} does not exist here; find the env that has'
                  ' maap-py with', file=sys.stderr)
            print('`conda env list`, or install it with `pip install maap-py`.',
                  file=sys.stderr)
        sys.exit(1)


def git(*args, check=True):
    '''Run a git command in REPO_DIR and return its stripped stdout.'''
    result = subprocess.run(('git', '-C', REPO_DIR) + args,
                            capture_output=True, text=True)
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, args,
                                            result.stdout, result.stderr)
    return result.stdout.rstrip()


def check_build_will_contain_local_work(version, fetch=True):
    '''
    Report every way the build could come out without the local work in it.

    inputs:
        version (str): algorithm_version from the yaml -- the ref DPS clones
        fetch (bool): update the remote ref first, so the comparison is honest
    output:
        list of problem strings; empty means the build will match the checkout
    '''
    problems = []

    try:
        git('rev-parse', '--git-dir')
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f'NOTE: {REPO_DIR} is not a git checkout; skipping the push checks.')
        return problems

    if fetch:
        # Read-only as far as the working tree goes, but it is what makes
        # origin/<version> mean "what GitHub has" rather than "what it had the
        # last time anything fetched".  Offline is not fatal.
        try:
            git('fetch', '--quiet', 'origin', version)
        except subprocess.CalledProcessError as exc:
            print(f'NOTE: could not fetch origin/{version}'
                  f' ({exc.stderr.strip()}); comparing against the local copy.')

    dirty = git('status', '--porcelain')
    if dirty:
        listed = '\n'.join('      ' + line for line in dirty.splitlines())
        problems.append('uncommitted changes in the working tree:\n' + listed)

    try:
        remote_ref = git('rev-parse', '--verify', f'origin/{version}')
    except subprocess.CalledProcessError:
        problems.append(
            f'origin/{version} does not exist.  algorithm_config.yml names'
            f' algorithm_version: {version},\n'
            f'      so the build clones that ref -- push the branch first.')
        return problems

    branch = git('rev-parse', '--abbrev-ref', 'HEAD')
    if branch != version:
        # Not fatal on its own: HEAD may already be merged into the ref DPS
        # clones.  The ancestry check below decides.
        where = 'a detached HEAD' if branch == 'HEAD' else f'branch {branch}'
        print(f'NOTE: registering from {where}, but the build clones'
              f' origin/{version}.')

    missing = git('rev-list', f'origin/{version}..HEAD')
    if missing:
        lines = git('log', '--oneline', f'origin/{version}..HEAD').splitlines()
        listed = '\n'.join('      ' + line for line in lines)
        problems.append(
            f'{len(missing.splitlines())} commit(s) on HEAD are not in'
            f' origin/{version} ({remote_ref[:7]}):\n' + listed)

    return problems


def build_url(response):
    '''
    Pull the build URL out of the registration response.

    The shape belongs to the API, not to maap-py, so nothing here can be taken
    on faith: a successful registration whose JSON is laid out differently
    should still print something useful rather than raise a KeyError.
    '''
    try:
        payload = response.json()
    except ValueError:
        return None, response.text

    for path in (('message', 'job_web_url'), ('job_web_url',)):
        node = payload
        for key in path:
            if not isinstance(node, dict) or key not in node:
                node = None
                break
            node = node[key]
        if isinstance(node, str) and node:
            return node, payload

    return None, payload


def main():
    parser = argparse.ArgumentParser(
        description='Register the ATL1415 DPS algorithm and print the build URL.')
    parser.add_argument('config', nargs='?',
                        default=os.path.join(REPO_DIR, 'algorithm_config.yml'),
                        help='algorithm config yaml (default: the one beside this script)')
    parser.add_argument('--host', default=os.environ.get('MAAP_API_HOST', DEFAULT_HOST),
                        help='MAAP API host (default: $MAAP_API_HOST, else '
                             + DEFAULT_HOST + ')')
    parser.add_argument('--dry-run', action='store_true',
                        help='read the config and run the push checks, then stop')
    parser.add_argument('--force', action='store_true',
                        help='register even if the build would not contain the local work')
    parser.add_argument('--no-fetch', action='store_true',
                        help='do not fetch origin before comparing (offline)')
    args = parser.parse_args()

    MAAP, algorithm_utils = import_maap()

    if not os.path.exists(args.config):
        print(f'ERROR: no such config file: {args.config}', file=sys.stderr)
        sys.exit(1)

    # maap-py's own reader, so what is reported here is what gets registered.
    config = algorithm_utils.read_yaml_file(args.config)
    name = config.get('algorithm_name', '<unnamed>')
    version = config.get('algorithm_version')
    if not version:
        print(f'ERROR: {args.config} has no algorithm_version; the build has no'
              ' ref to clone.', file=sys.stderr)
        sys.exit(1)

    print(f'config:     {args.config}')
    print(f'algorithm:  {name}:{version}')
    print(f'repository: {config.get("repository_url", "<unset>")}')
    print(f'host:       {args.host}')
    print()

    sys.stdout.flush()
    problems = check_build_will_contain_local_work(version, fetch=not args.no_fetch)
    sys.stdout.flush()
    if problems:
        head = ('WARNING (--force given, registering anyway)' if args.force
                else 'REFUSING TO REGISTER')
        print(f'{head} -- the build would not contain your work:', file=sys.stderr)
        for problem in problems:
            print(f'  - {problem}', file=sys.stderr)
        print(file=sys.stderr)
        print(f'The build clones origin/{version} from GitHub, not this working'
              ' copy.', file=sys.stderr)
        if not args.force:
            print('Push (or commit and push) first, or re-run with --force.',
                  file=sys.stderr)
            sys.exit(1)
        print(file=sys.stderr)
    else:
        print(f'push check: OK -- origin/{version} contains everything in this'
              ' checkout.')

    if args.dry_run:
        print()
        print('--dry-run: stopping before registration.')
        return

    print()
    print(f'registering {name}:{version} ...')
    response = MAAP(maap_host=args.host).register_algorithm_from_yaml_file(args.config)
    print(f'HTTP {response.status_code}')

    url, payload = build_url(response)
    if url is None:
        print('No job_web_url in the response; the full body was:', file=sys.stderr)
        if isinstance(payload, (dict, list)):
            print(json.dumps(payload, indent=2), file=sys.stderr)
        else:
            print(payload, file=sys.stderr)

    if not 200 <= response.status_code < 300:
        print('Registration FAILED.', file=sys.stderr)
        sys.exit(2)

    if url is None:
        # Registered, but there is nothing to hand back to the caller.
        sys.exit(2)

    print()
    print('Build log (browser only -- it is not on this filesystem):')
    print(url)


if __name__ == '__main__':
    main()
