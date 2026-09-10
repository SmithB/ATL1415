#!/usr/bin/env python3
"""
Register the ATL1415 DPS algorithm with MAAP's OGC build service, and follow
the build into its deployment.

This is docs/howto_MAAP_ogc.sh step O3 -- and O4 is running it.

WHAT CHANGED, 2026-09-10.  The ADE's maap-py became 5.1.0a2, which has no
register_algorithm_from_yaml_file and none of the legacy job calls, and Ben
chose to move to MAAP's OGC system rather than stay on the ATL14 env's 4.2.0.
Registration is now a JSON POST to https://api.maap-project.org/api/build --
the same call MAAP's own Algorithm Catalog "Register New Algorithm" form makes
(maap_algorithms_jupyter_extension 1.0.2) -- and MAAP builds the image from
the repo with build-env.sh on maap_base, generates the CWL, and deploys it as
an OGC process.  So the model is unchanged; only the transport is new.

This script uses maap-py ONLY for its auth header, which has the same
signature in 4.2.0 and 5.x, so it runs under either -- the notebook env or
ATL14.  Neither version's algorithm methods are touched.

WHAT IT KEEPS: REGISTERING TRIGGERS A BUILD THAT CLONES FROM GITHUB.  The
build clones code_repository at algorithm_version and bakes the result into
an image; the ADE working copy has nothing to do with it.  Anything
uncommitted, unpushed, or committed after the last build is simply not on the
worker, and no job log will ever say so -- a stale image surfaces as whatever
the missing fix was meant to prevent.  On 2026-09-04 that cost a cycle.  This
script refuses to register in that state; --force overrides.

WHAT IT PRINTS: the build id, and EVERY URL in every response, labelled with
its key path.  Ben learned on 2026-09-10 that the one URL the legacy script
printed (job_web_url) was only the build log; the OGC service hands back a
build pipeline, a deployment job and a deployment pipeline as separate links.
It walks the JSON rather than naming keys, because nothing documents them.

Run it from the ADE:

    ./register_algorithm.py [--dry-run] [--force] [algorithm_config.yml]
    ./register_algorithm.py --status <build_id> [--wait]
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time

# The default interpreter in an ADE terminal.  Since 2026-09-10 either env
# works (see the docstring); this is only what the error message suggests.
ADE_PYTHON = '/srv/conda/envs/notebook/bin/python'

# The repo the script lives in, so cwd never matters -- this is both where
# algorithm_config.yml is found and which checkout the git checks run against.
REPO_DIR = os.path.dirname(os.path.realpath(__file__))

DEFAULT_HOST = 'api.maap-project.org'


def import_maap():
    '''
    Import maap-py, or exit 1 naming the interpreter that has it.

    Deliberately does NOT re-exec itself under ADE_PYTHON: a script that
    silently switches interpreters hides which env is in play.
    '''
    try:
        from maap.maap import MAAP
        from maap.utils import algorithm_utils
        return MAAP, algorithm_utils
    except ImportError as exc:
        print(f"ERROR: 'maap' is not importable under {sys.executable}", file=sys.stderr)
        print(f"       ({type(exc).__name__}: {exc})", file=sys.stderr)
        print(file=sys.stderr)
        print('Any maap-py works (4.2.0 or 5.x).  The ADE notebook env has it,'
              ' and so does', file=sys.stderr)
        print('an ATL14 env built from environment.yml:', file=sys.stderr)
        print(file=sys.stderr)
        script = os.path.relpath(os.path.realpath(__file__))
        if os.path.exists(ADE_PYTHON):
            print(f'    {ADE_PYTHON} {script}', file=sys.stderr)
        else:
            # Do not print a path that is not there -- on a fresh account or a
            # different image the notebook env lives somewhere else.
            print(f'    <notebook-env>/bin/python {script}', file=sys.stderr)
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


# ---------------------------------------------------------------------------
# The build service's schema, AS ITS OWN CLIENT ENFORCES IT.  Everything below
# is read out of maap_algorithms_jupyter_extension 1.0.2 (the Algorithm
# Catalog plugin's static/156.*.js), not from documentation, because there is
# none.  Checking here means a bad config fails in the ADE with a reason,
# rather than as an opaque rejection -- or worse, an accepted build of the
# wrong thing.
# ---------------------------------------------------------------------------

# The only keys the form ever sends.  Anything else in the yaml is comment-
# level metadata as far as the service is concerned, and is not sent.
BUILD_FIELDS = (
    'algorithm_name', 'algorithm_version', 'algorithm_description',
    'code_repository', 'run_command', 'build_command', 'ram_min', 'cores_min',
    'algorithm_container_url', 'base_container_url', 'author', 'contributor',
    'license', 'release_notes', 'citation', 'keywords', 'inputs', 'outputs',
    'outdir_max',
)
INPUT_FIELDS = ('name', 'label', 'doc', 'type', 'default')
INPUT_TYPES = ('string', 'int', 'File', 'Directory', 'long', 'float',
               'boolean', 'double')
NAME_RE = re.compile(r'^[a-z0-9_-]+$')
VERSION_RE = re.compile(r'^[a-zA-Z0-9_][a-zA-Z0-9._-]{0,127}$')

# Terminal states, as the plugin's own polling loop treats them.
BUILD_DONE = {'successful', 'failed', 'canceled', 'cancelled', 'dismissed'}
DEPLOY_DONE = {'deployed', 'successful', 'failed', 'error', 'canceled',
               'cancelled', 'dismissed', 'not found'}


def validate_config(config):
    """
    Every reason the build service's own client would refuse this config.

    output:
        list of problem strings; empty means the form would have accepted it
    """
    problems = []
    name = str(config.get('algorithm_name') or '')
    version = str(config.get('algorithm_version') or '')
    if not 2 <= len(name) <= 255 or not NAME_RE.match(name):
        problems.append(
            f"algorithm_name {name!r} must be 2-255 characters of lowercase"
            " letters, digits, '-' and '_' (^[a-z0-9_-]+$).")
    if not VERSION_RE.match(version):
        problems.append(
            f"algorithm_version {version!r} must match"
            " ^[a-zA-Z0-9_][a-zA-Z0-9._-]{0,127}$ -- it is the branch or tag"
            " the build clones.")
    for key in ('code_repository', 'run_command'):
        if not config.get(key):
            problems.append(f'{key} is required.')
    if not config.get('algorithm_container_url') and not (
            config.get('base_container_url') and config.get('build_command')):
        problems.append(
            "need either algorithm_container_url (a prebuilt image) or both"
            " base_container_url and build_command (build it here).")
    for key, cap in (('ram_min', 128), ('cores_min', 32)):
        value = config.get(key)
        if value is None:
            problems.append(f'{key} is required.')
            continue
        try:
            ok = 0 < float(value) <= cap
        except (TypeError, ValueError):
            ok = False
        if not ok:
            problems.append(f'{key}={value!r} must be a number in (0, {cap}].')
    inputs = config.get('inputs') or []
    if not isinstance(inputs, list):
        # The legacy schema's {positional: [...], file: [...]} -- iterating it
        # would walk the KEYS and report each as a nameless input.
        problems.append(
            'inputs must be a LIST of {name, label, doc, type, default}; the'
            ' legacy {positional:, file:} form is not accepted.')
        inputs = []
    for i, item in enumerate(inputs):
        if not isinstance(item, dict) or not item.get('name'):
            problems.append(f'inputs[{i}] has no name.')
            continue
        # The docs' own example uses 'string?' for an optional input; the
        # form's dropdown has no '?' variants, so accept both.
        kind = str(item.get('type') or '').rstrip('?')
        if kind not in INPUT_TYPES:
            problems.append(f"inputs[{i}] ({item['name']}): type"
                            f" {item.get('type')!r} is not one of {INPUT_TYPES}.")
    return problems


def build_body(config):
    """
    The JSON the form would send for this config.

    Mirrors the form's own serializer: unknown keys dropped, empty values
    dropped (including inside each input), outputs defaulting to the single
    Directory 'out', outdir_max defaulting to 20.  A folded yaml description
    ('>') ends in a newline, which is stripped.
    """
    body = {}
    for key in BUILD_FIELDS:
        value = config.get(key)
        if value is None or value == '' or value == []:
            continue
        if isinstance(value, str):
            value = value.strip()
        body[key] = value
    if 'inputs' in body:
        body['inputs'] = [
            {k: item[k] for k in INPUT_FIELDS if item.get(k) not in (None, '')}
            for item in body['inputs']]
    body.setdefault('outputs', [{'name': 'out', 'type': 'Directory'}])
    body.setdefault('outdir_max', 20)
    return body


def api(maap, host, method, path, body=None):
    """
    One authenticated call to the MAAP API, returning (http_status, payload).

    payload is the parsed JSON, or the raw text when the body is not JSON --
    an error page is still worth printing.
    """
    import requests      # maap-py depends on it, so it is always present
    url = f'https://{host}/api/{path.lstrip("/")}'
    headers = maap._get_api_header(content_type='application/json')
    response = requests.request(method, url, headers=headers, json=body,
                                timeout=120)
    try:
        payload = response.json()
    except ValueError:
        payload = response.text
    return response.status_code, payload


# A URL ends at whitespace, or at a quote or closing bracket, which is what
# ends one embedded in a sentence like 'see <https://...>'.  Trailing
# punctuation is stripped afterwards rather than excluded here, so a comma or
# full stop INSIDE a path survives.
URL_RE = re.compile(r'https?://[^\s\'"<>)\]]+')


def urls_in(node, path=()):
    '''
    Yield (key_path, url) for every http(s) URL anywhere in a parsed response.

    Walks dicts and lists all the way down, and finds URLs embedded inside a
    longer string as well as ones that are the whole value -- a status message
    that says "build started, see https://..." counts.

    inputs:
        node: the parsed JSON (dict, list or scalar)
        path (tuple): key path to node, used for the label
    output:
        generator of (str, str): dotted key path, e.g. 'message.job_web_url',
        with list indices as [i]; and the URL
    '''
    if isinstance(node, dict):
        for key, value in node.items():
            yield from urls_in(value, path + (str(key),))
    elif isinstance(node, list):
        for i, value in enumerate(node):
            yield from urls_in(value, path + (f'[{i}]',))
    elif isinstance(node, str):
        label = '.'.join(path).replace('.[', '[') or '<top level>'
        for url in URL_RE.findall(node):
            yield label, url.rstrip('.,;:')


def print_urls(payload, heading):
    """Every URL in a response, labelled by key path.  Silent when none."""
    found = list(urls_in(payload)) if isinstance(payload, (dict, list)) else []
    if not found:
        return
    width = max(len(label) for label, _ in found)
    print(f'{heading} ({len(found)}; browser-only):')
    for label, link in found:
        print(f'  {label:<{width}}  {link}')


def dump(payload, stream=sys.stdout):
    if isinstance(payload, (dict, list)):
        print(json.dumps(payload, indent=2), file=stream)
    else:
        print(payload, file=stream)


def deployment_id(build):
    """
    The deployment job id, from the build's deploymentLink.

    The link arrives as {'href': '.../deploymentJobs/<id>', ...}; the plugin
    pulls the id out with this same pattern.
    """
    link = build.get('deploymentLink') if isinstance(build, dict) else None
    href = link.get('href') if isinstance(link, dict) else link
    m = re.search(r'deploymentJobs/(\w+)', href or '')
    return m.group(1) if m else None


def report_processes(maap, host, name, version):
    """
    The processID(s) a deployed build became -- what submit_job() needs.

    Matched on the process's string id and version; submit_job() POSTs to
    /api/ogc/processes/<process_id>/execution, and whether it wants the
    numeric processID or the string id is not yet known (howto_MAAP_ogc F9),
    so both are printed.
    """
    code, payload = api(maap, host, 'GET', 'ogc/processes')
    procs = payload.get('processes', []) if isinstance(payload, dict) else []
    mine = [p for p in procs if p.get('id') == name
            or str(p.get('title', '')).lower() == name]
    if not mine:
        print(f'  no deployed process named {name!r} yet (HTTP {code}).')
        return
    for p in mine:
        mark = '   <- this version' if str(p.get('version')) == str(version) else ''
        print(f"  processID={p.get('processID')}  id={p.get('id')}"
              f"  version={p.get('version')}{mark}")


def show_status(maap, host, build_id, config, wait=False, poll_s=30):
    """
    Follow one build through to its deployment, printing every link.

    Returns 0 when the build and its deployment both succeeded, 2 when either
    ended any other way, and 3 when not finished (without --wait).
    """
    started = time.time()
    while True:
        code, build = api(maap, host, 'GET', f'build/{build_id}')
        print(f'GET build/{build_id} -> HTTP {code}')
        if not isinstance(build, dict) or not 200 <= code < 300:
            dump(build, sys.stderr)
            return 2
        b_status = str(build.get('status', '?'))
        print(f"  build status:      {b_status}")
        if build.get('deploymentError'):
            print(f"  deploymentError:   {build['deploymentError']}")
        print_urls(build, '  URLs in the build record')

        d_status, dep_id = None, deployment_id(build)
        if dep_id:
            dcode, dep = api(maap, host, 'GET', f'ogc/deploymentJobs/{dep_id}')
            print(f'GET ogc/deploymentJobs/{dep_id} -> HTTP {dcode}')
            if isinstance(dep, dict):
                d_status = str(dep.get('status', '?'))
                print(f'  deployment status: {d_status}')
                if dep.get('error'):
                    print(f"  deployment error:  {dep['error']}")
                print_urls(dep, '  URLs in the deployment record')
            else:
                dump(dep, sys.stderr)

        build_over = b_status.lower() in BUILD_DONE
        deploy_over = d_status is not None and d_status.lower() in DEPLOY_DONE
        if build_over and (deploy_over or b_status.lower() != 'successful'):
            break
        if not wait:
            print('\n(not finished -- rerun, or add --wait to poll)')
            return 3
        elapsed = int(time.time() - started)
        print(f'\n... {elapsed} s, polling again in {poll_s} s\n')
        time.sleep(poll_s)

    ok = (b_status.lower() == 'successful'
          and d_status is not None and d_status.lower() in ('deployed', 'successful'))
    print()
    print(f"Deployed processes for {config.get('algorithm_name')}:")
    report_processes(maap, host, config.get('algorithm_name'),
                     config.get('algorithm_version'))
    print('\nRESULT:', 'DEPLOYED' if ok else f'NOT DEPLOYED (build {b_status},'
          f' deployment {d_status})')
    return 0 if ok else 2


def main():
    parser = argparse.ArgumentParser(
        description='Register the ATL1415 DPS algorithm with the MAAP OGC build'
                    ' service, and follow it into deployment.')
    parser.add_argument('config', nargs='?',
                        default=os.path.join(REPO_DIR, 'algorithm_config.yml'),
                        help='algorithm config yaml (default: the one beside this script)')
    parser.add_argument('--host', default=os.environ.get('MAAP_API_HOST', DEFAULT_HOST),
                        help='MAAP API host (default: $MAAP_API_HOST, else '
                             + DEFAULT_HOST + ')')
    parser.add_argument('--dry-run', action='store_true',
                        help='validate, run the push checks, print the JSON that'
                             ' would be sent, then stop')
    parser.add_argument('--force', action='store_true',
                        help='register even if the build would not contain the local work')
    parser.add_argument('--no-fetch', action='store_true',
                        help='do not fetch origin before comparing (offline)')
    parser.add_argument('--status', metavar='BUILD_ID',
                        help='report a build and its deployment instead of registering')
    parser.add_argument('--wait', action='store_true',
                        help='with --status: poll until the build and deployment finish')
    args = parser.parse_args()

    MAAP, algorithm_utils = import_maap()

    if not os.path.exists(args.config):
        print(f'ERROR: no such config file: {args.config}', file=sys.stderr)
        sys.exit(1)

    # maap-py's own reader -- present in 4.2.0 and 5.x alike.
    config = algorithm_utils.read_yaml_file(args.config)
    name = config.get('algorithm_name', '<unnamed>')
    version = config.get('algorithm_version')

    if args.status:
        sys.exit(show_status(MAAP(maap_host=args.host), args.host, args.status,
                             config, wait=args.wait))

    print(f'config:     {args.config}')
    print(f'algorithm:  {name}:{version}')
    print(f'repository: {config.get("code_repository", "<unset>")}')
    print(f'host:       {args.host}')
    print()

    invalid = validate_config(config)
    if invalid:
        print('INVALID CONFIG -- the build service\'s own form would refuse it:',
              file=sys.stderr)
        for problem in invalid:
            print(f'  - {problem}', file=sys.stderr)
        sys.exit(1)

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

    body = build_body(config)
    if args.dry_run:
        print()
        print('JSON that would be POSTed to /api/build:')
        dump(body)
        print()
        print('--dry-run: stopping before registration.')
        return

    print()
    print(f'registering {name}:{version} -> POST /api/build ...')
    maap = MAAP(maap_host=args.host)
    code, payload = api(maap, args.host, 'POST', 'build', body)
    print(f'HTTP {code}')

    accepted = (200 <= code < 300 and isinstance(payload, dict)
                and str(payload.get('status', '')).lower() == 'accepted')
    if not accepted:
        print('Registration NOT accepted; the full response was:', file=sys.stderr)
        dump(payload, sys.stderr)
        sys.exit(2)

    build_id = payload.get('build_id')
    print(f'build_id:   {build_id}')
    print_urls(payload, 'URLs in the registration response')
    print()
    print('Follow the build and its deployment -- every link, as each appears:')
    print(f'    {os.path.relpath(os.path.realpath(__file__))} --status {build_id} --wait')


if __name__ == '__main__':
    main()
