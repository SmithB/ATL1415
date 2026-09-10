#!/usr/bin/env python3
"""
Submit ONE build_id job and say whether the deployed image is the commit the
build service says it built -- and whether its workers can read NSIDC.

WHY THIS EXISTS.  The build clones the repo and bakes it into an image, so
the ADE working copy has nothing to do with what a worker runs.  On
2026-09-09 a re-registered version ran the OLD code, and nothing in any log
said so; a stale image looked exactly like a fix that did not work.  This is
the cheap test for that: one job, a few worker-seconds, before spending
worker hours on a run whose results would have to be thrown away.

OGC SYSTEM, since 2026-09-10 (docs/howto_MAAP_ogc.sh O5).  maap-py 5.x
dropped submitJob/getJob/getJobResult; this uses submit_job, get_job_status
and get_job_result against the deployed process, found by name and version.

WHAT IT COMPARES.  Three commits, from three independent sources:
  image   the stamp build-env.sh wrote INSIDE the image at build time,
          printed by `run.sh --build-id` on a worker
  cwl     s:commitHash in the CWL the build service generated -- what the
          service says it cloned (for the first build, d401699)
  origin  the tip of origin/<algorithm_version> on GitHub, right now
image == cwl is the test: the image is what the service says it built.
origin is reported, not required -- GitHub moving on after a build is normal
(push after registering and it will differ), and comparing against it alone
would call a correct image stale.

AND THE WORKER'S CREDENTIALS.  The job also reports maap_pgt=set|unset.
pointCollection asks MAAP for NSIDC's temporary S3 credentials only when
MAAP_PGT is set, and otherwise falls back silently to earthaccess, which has
no credentials on a worker -- so unset means no tile here can read ATL11.

VERDICTS (exit status):
  0  MATCH     image == cwl, and maap_pgt=set.  Proceed.
  1  MISMATCH  image != cwl: the image is not what the service built (the
               2026-09-09 failure).  Do NOT revive on_s3_v2 -- move to an
               immutable per-build tag (8aad07d).
  1  NO STAMP  the image predates the build stamp, so it is stale by
               definition: every build since 2026-09-10 writes one.
  1  NO NSIDC  maap_pgt=unset, whatever the commits say.
  2  the check itself could not run: not deployed, submit refused, timed
     out, or the job log could not be read (the full status and result are
     printed, for howto_MAAP_ogc QD).

Usage:
  check_build_id.py [args_file_url] [queue] [--expect <sha>] [--timeout <s>]
                    [--dry-run]

--dry-run finds the process, reads its CWL and GitHub, prints what it would
submit, and stops -- no job.  --expect replaces the cwl commit as the one the
image must equal.

The args_file is a REQUIRED input of the process, so one must be named even
though run.sh exits before reading it; the AA release args file is the
default because it is already published.  The default queue is -32gb, not
-8gb or -16gb, because the CWL asks for ramMin 16 and a 16 GB worker may
have less than that to allocate.
"""
import json
import os
import re
import subprocess
import sys
import time

import requests
import yaml
from maap.maap import MAAP

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
S3_RUN = 's3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/south/AA'
DEFAULT_QUEUE = 'maap-dps-worker-32gb'
POLL_S = 15
# OGC API - Processes job states that do not change again.
DONE = {'successful', 'failed', 'dismissed', 'deleted'}


def load_config():
    with open(os.path.join(REPO, 'algorithm_config.yml')) as fh:
        config = yaml.safe_load(fh)
    for key in ('algorithm_name', 'algorithm_version'):
        if not config.get(key):
            print(f'algorithm_config.yml has no {key}', file=sys.stderr)
            sys.exit(2)
    return config


def find_process(maap, name, version):
    """
    The deployed process for name:version, as list_algorithms() reports it.

    Looked up every time rather than hard-coded: its numeric processID (64 for
    the first deployment) is what submit_job() puts in the URL -- the
    deployment record's processLocation is /ogc/processes/64 -- and it can
    change with every redeploy.
    """
    response = maap.list_algorithms()
    response.raise_for_status()
    procs = response.json().get('processes', [])
    mine = [p for p in procs
            if p.get('id') == name and str(p.get('version')) == str(version)]
    if not mine:
        # Exit 2, not SystemExit(message): a message exits 1, which is the
        # status reserved for a verdict against the image.
        print(f'NO DEPLOYED PROCESS {name}:{version} (of {len(procs)} listed).\n'
              '  Register with register_algorithm.py and wait for the build\n'
              '  pipeline it prints to finish deploying.', file=sys.stderr)
        sys.exit(2)
    if len(mine) > 1:
        print(f'NOTE: {len(mine)} processes match {name}:{version}; using the'
              ' most recently modified.')
    return max(mine, key=lambda p: str(p.get('lastModifiedTime', '')))


def cwl_facts(process):
    """
    s:commitHash and the image, from the CWL the build service generated.

    The CWL is public (a raw file in MAAP's GitLab), so this needs no auth.
    Either value is None if the CWL cannot be read or does not carry it.
    """
    link = process.get('cwlLink')
    if not link:
        return None, None
    try:
        text = requests.get(link, timeout=60).text
    except requests.RequestException as exc:
        print(f'NOTE: could not read the CWL at {link} ({exc})')
        return None, None
    commit = re.search(r'^s:commitHash:\s*(\S+)', text, re.M)
    image = re.search(r'dockerPull:\s*(\S+)', text)
    return (commit.group(1) if commit else None,
            image.group(1) if image else None)


def origin_tip(version):
    """The tip of origin/<version> on GitHub now -- a branch, else a tag."""
    for ref in (f'refs/heads/{version}', f'refs/tags/{version}'):
        p = subprocess.run(['git', '-C', REPO, 'ls-remote', 'origin', ref],
                           capture_output=True, text=True)
        if p.returncode == 0 and p.stdout.strip():
            return p.stdout.split()[0]
    return None


def normalize_s3(uri):
    """
    Drop the endpoint element from the s3 URI the job result hands back.

    It arrives as s3://s3-us-west-2.amazonaws.com:80/maap-ops-workspace/... --
    host:port first, bucket second -- and `aws s3 cp` reads the first element
    after the scheme as the bucket.  Same fix as collect_AA_queue.py.
    """
    rest = uri[len('s3://'):]
    head, _, tail = rest.partition('/')
    if 'amazonaws.com' in head:
        rest = tail
    return 's3://' + rest.rstrip('/')


def s3_prefixes(node):
    """Every s3:// string anywhere in a parsed job result, normalized."""
    if isinstance(node, dict):
        for value in node.values():
            yield from s3_prefixes(value)
    elif isinstance(node, list):
        for value in node:
            yield from s3_prefixes(value)
    elif isinstance(node, str) and node.startswith('s3://'):
        yield normalize_s3(node)


def read_stdout(result):
    """
    The job's _stdout.txt, from the output prefix in get_job_result().

    get_job_result() returns the prefix three ways (website, endpoint-style
    s3, console) under {"<name>": {"links": [{"href": ...}, ...]}}, and
    _stdout.txt sits at that prefix -- checked on a legacy job through these
    same OGC endpoints, 2026-09-10.  Returns (text, prefix) or ('', None).
    """
    for prefix in dict.fromkeys(s3_prefixes(result)):
        p = subprocess.run(['aws', 's3', 'cp', f'{prefix}/_stdout.txt', '-'],
                           capture_output=True, text=True, timeout=120)
        if p.returncode == 0:
            return p.stdout, prefix
    return '', None


def job_id_from(response):
    """The job id from a submit_job() response: body first, then Location."""
    try:
        body = response.json()
    except ValueError:
        body = {}
    if isinstance(body, dict):
        for key in ('jobID', 'job_id', 'id'):
            if body.get(key):
                return str(body[key])
    m = re.search(r'/jobs/([^/?#]+)', response.headers.get('Location', ''))
    return m.group(1) if m else None


def container_of(maap, tag):
    """
    Best effort: the container the job actually ran, from its job record.

    The record carries context.container_specification (url, digest) --
    seen on a legacy job, 2026-09-10.  A digest that stays the same across a
    rebuild would mean the image was reused, whatever the tag says.
    """
    try:
        jobs = maap.list_jobs(tag=tag, page_size=5).json().get('jobs', [])
    except Exception:
        return None
    for job in jobs:
        spec = (job.get('context') or {}).get('container_specification') or {}
        if spec:
            return spec
    return None


def parse_build_id_line(text):
    """The key=value fields of run.sh's one-line summary, or None."""
    line = re.search(r'^BUILD_ID: (.*)$', text, re.M)
    if not line:
        return None
    return dict(field.split('=', 1) for field in line.group(1).split()
                if '=' in field)


def verdict(fields, want, origin, image_label='image'):
    """
    (exit_status, lines) for the parsed BUILD_ID fields.

    want is the commit the image must equal (cwl, or --expect); origin is
    only reported.
    """
    got = fields.get('commit', 'unknown')
    lines, status = [], 0
    if got == 'unknown':
        lines.append('VERDICT: NO STAMP -- the image was built before the build'
                     ' stamp existed, so it is not a fresh build of current code.')
        status = 1
    elif want and got != want:
        lines.append(f'VERDICT: MISMATCH\n  {image_label:6} {got}\n  built  {want}\n'
                     '  The image is not the commit the build service says it'
                     ' built -- an image reused under the tag.  Do NOT revive'
                     ' on_s3_v2; move to an immutable per-build tag (8aad07d).')
        status = 1
    elif want:
        lines.append(f'VERDICT: MATCH -- the image is {got}, the commit the build'
                     ' service recorded.')
    else:
        lines.append(f'VERDICT: UNCHECKED -- the image is {got}, but there is no'
                     ' recorded commit to compare it with (no s:commitHash, no'
                     ' --expect).')
        status = 1
    if origin and got not in ('unknown', origin):
        lines.append(f'  NOTE: origin is now at {origin[:12]}, past this build.'
                     '  Normal after a push; register again to pick it up.')
    pgt = fields.get('maap_pgt', 'unknown')
    lines.append(f"  maap_py={fields.get('maap_py', 'unknown')}  maap_pgt={pgt}")
    if pgt != 'set':
        lines.append('VERDICT: NO NSIDC -- MAAP_PGT is not set on this worker, so'
                     ' pointCollection cannot get NSIDC credentials and no tile'
                     ' here can read ATL11.  Raise it with MAAP support.')
        status = 1
    return status, lines


def main():
    argv = list(sys.argv[1:])
    expect, timeout, dry_run = None, 1800, False
    if '--dry-run' in argv:
        argv.remove('--dry-run')
        dry_run = True
    for flag in ('--expect', '--timeout'):
        if flag in argv:
            i = argv.index(flag)
            value = argv[i + 1]
            del argv[i:i + 2]
            if flag == '--expect':
                expect = value
            else:
                timeout = int(value)

    args_url = argv[0] if argv else f'{S3_RUN}/input_args_AA.txt'
    queue = argv[1] if len(argv) > 1 else DEFAULT_QUEUE
    config = load_config()
    name, version = config['algorithm_name'], config['algorithm_version']

    maap = MAAP(maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org'))
    process = find_process(maap, name, version)
    pid = process.get('processID')
    cwl_commit, image = cwl_facts(process)
    origin = origin_tip(version)
    want = expect or cwl_commit

    print(f'process          : {name}:{version}  processID={pid}'
          f"  (modified {process.get('lastModifiedTime')})")
    print(f'image            : {image or "<not in the CWL>"}')
    print(f'built (cwl)      : {cwl_commit or "<no s:commitHash in the CWL>"}')
    print(f'origin/{version:9}: {origin or "<could not resolve>"}')
    if expect:
        print(f'--expect         : {expect}  (replaces the cwl commit)')
    print(f'args file        : {args_url}')
    print(f'queue            : {queue}\n')

    # tag: findable in the Jobs UI and in list_jobs(tag=...).
    tag = f'atl1415_build_id_{int(time.time())}'
    inputs = {'x0': '0', 'y0': '0', 'step': 'build_id', 'args_file': args_url}
    if dry_run:
        print('--dry-run: would submit_job(%s, %s, %r, dedup=False, tag=%r)'
              % (pid, json.dumps(inputs), queue, tag))
        return

    # dedup=False, EXPLICITLY: this job's inputs are identical every time, and
    # a service that deduplicated it would hand back the PREVIOUS image's
    # answer after a rebuild -- the one thing this check must never do.
    response = maap.submit_job(pid, inputs, queue, dedup=False, tag=tag)
    print(f'submit_job -> HTTP {response.status_code}')
    job_id = job_id_from(response)
    if not 200 <= response.status_code < 300 or not job_id:
        print('submit refused or returned no job id; the response was:')
        print(response.text[:2000])
        sys.exit(2)
    print(f'job id           : {job_id}   tag: {tag}\n')

    deadline, status, record = time.time() + timeout, None, {}
    while time.time() < deadline:
        r = maap.get_job_status(job_id)
        try:
            record = r.json()
        except ValueError:
            record = {'raw': r.text[:500]}
        status = str(record.get('status', f'HTTP {r.status_code}'))
        print(f'  [{time.strftime("%H:%M:%S")}] {status}')
        if status.lower() in DONE:
            break
        time.sleep(POLL_S)
    else:
        print(f'timed out after {timeout}s; last status {status}')
        sys.exit(2)

    r = maap.get_job_result(job_id)
    try:
        result = r.json()
    except ValueError:
        result = r.text
    out, prefix = read_stdout(result)
    if not out:
        # This is howto_MAAP_ogc QD failing: print everything, so the answer
        # to "where do OGC job logs go" can be read off and recorded.
        print(f'\njob finished {status}, but no _stdout.txt could be read.')
        print('get_job_status:'); print(json.dumps(record, indent=2))
        print(f'get_job_result (HTTP {r.status_code}):')
        print(json.dumps(result, indent=2) if not isinstance(result, str) else result)
        sys.exit(2)
    print(f'log              : {prefix}/_stdout.txt')

    spec = container_of(maap, tag)
    if spec:
        print(f"container        : {spec.get('url') or spec.get('id')}\n"
              f"digest           : {spec.get('digest')}")

    # The whole block: the per-field detail is what makes a mismatch
    # diagnosable rather than merely visible.
    start = out.find('  ATL1415 build id')
    print('\n' + (out[start:] if start >= 0 else out).rstrip() + '\n')

    fields = parse_build_id_line(out)
    if fields is None:
        print('VERDICT: NO STAMP -- no BUILD_ID line at all: this image predates'
              ' run.sh --build-id, so it did not pick up current code.')
        sys.exit(1)
    code, lines = verdict(fields, want, origin)
    print('\n'.join(lines))
    sys.exit(code)


if __name__ == '__main__':
    main()
