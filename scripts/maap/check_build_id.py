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
  1  the report was printed but the JOB then failed: the image verdict is
     shown, and the runner's errors after it -- not green, because every
     tile would fail the same way.
  2  the check itself could not run: not deployed, submit refused, timed
     out, no log readable, or the job failed BEFORE run.sh printed a
     report -- which says nothing about the image, and is never reported
     as a stale one.

Usage:
  check_build_id.py [args_file_url] [queue] [--expect <sha>] [--timeout <s>]
                    [--dry-run] [--job <job_id>]

--dry-run finds the process, reads its CWL and GitHub, prints what it would
submit, and stops -- no job.  --job re-reads a job that already ran (after a
timeout here, or one submitted another way) instead of submitting.  --expect replaces the cwl commit as the one the
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


def read_logs(result):
    """
    Everything the job left that could hold run.sh's report, joined.

    UNDER OGC THE REPORT IS IN _stderr.txt, NOT _stdout.txt.  The job runs in
    MAAP's CWL runner (container-maap-cwltool-executor), which launches our
    image with cwltool; the runner's own chatter is _stdout.txt, and our
    container's output comes out on cwltool's stderr.  The first OGC build_id
    job (db93c7f3..., 2026-09-10) printed a complete, correct report into
    _stderr.txt while _stdout.txt held nine lines of runner log -- and the
    first version of this script, reading only _stdout.txt, called a correct
    image stale.  So read both, plus build_id.txt, which run.sh also writes
    as a product since that job (the name it lands under is not yet seen, so
    both likely places are tried).

    get_job_result() gives the prefix three ways (website, endpoint-style s3,
    console); a failed job's prefix is under dataset/triaged_job/.
    Returns (text, prefix, [files read]); ('', None, []) when nothing reads.
    """
    for prefix in dict.fromkeys(s3_prefixes(result)):
        parts, names = [], []
        for name in ('_stdout.txt', '_stderr.txt', 'build_id.txt',
                     'output/build_id.txt'):
            p = subprocess.run(['aws', 's3', 'cp', f'{prefix}/{name}', '-'],
                               capture_output=True, text=True, timeout=120)
            if p.returncode == 0 and p.stdout:
                parts.append(p.stdout)
                names.append(name)
        if names:
            return '\n'.join(parts), prefix, names
    return '', None, []


def image_digest(text):
    """
    The digest of OUR image, from the docker pull the runner logs.

    The job record's container_specification is the RUNNER's container
    (maap-cwltool-executor), not ours -- so this, not container_of(), is the
    digest that says whether a rebuild produced a new image (QB).
    """
    m = re.search(r'^Digest: (sha256:[0-9a-f]+)', text, re.M)
    return m.group(1) if m else None


# Error lines the runner prints on EVERY job, failed or not, which read as
# alarming and mean nothing: docker inspects the image before pulling it
# ("No such object", followed by the pull), and cleans up a container that
# has already exited ("cannot kill container ... No such container").
BENIGN = (r'^Error: No such object: ', r'cannot kill container: .*No such container')


def runner_errors(text, keep=8):
    """The runner's error lines, for a job that failed, minus BENIGN ones."""
    lines, hits = text.splitlines(), []
    for i, line in enumerate(lines):
        if any(re.search(b, line) for b in BENIGN):
            continue
        if re.search(r'ERROR|permanentFail|Error|Traceback', line):
            hits.append(line)
            # cwltool puts the REASON on the line after "ERROR ... Job error:"
            if 'ERROR' in line and i + 1 < len(lines) and lines[i + 1] not in hits:
                hits.append(lines[i + 1])
    return list(dict.fromkeys(hits))[-keep:]


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
    Best effort: the container named in the job's record.

    On a LEGACY job that was the algorithm image.  On an OGC job it is the
    RUNNER -- container-maap-cwltool-executor -- which launches our image
    inside it; image_digest() reads ours from the log.  Reported, labelled.
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
    expect, timeout, dry_run, job_id = None, 1800, False, None
    if '--dry-run' in argv:
        argv.remove('--dry-run')
        dry_run = True
    for flag in ('--expect', '--timeout', '--job'):
        if flag in argv:
            i = argv.index(flag)
            value = argv[i + 1]
            del argv[i:i + 2]
            if flag == '--expect':
                expect = value
            elif flag == '--job':
                job_id = value
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

    if job_id:
        print(f'job id           : {job_id}   (--job: re-reading, not submitting)\n')
        tag = None
    else:
        # dedup=False, EXPLICITLY: this job's inputs are identical every time,
        # and a service that deduplicated it would hand back the PREVIOUS
        # image's answer after a rebuild -- the one thing this must never do.
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
        # A job is not visible to get_job_status for a few seconds after
        # submit_job accepts it: the first poll of the first OGC job got a
        # 404 problem document, and the next one 'accepted'.
        shown = ('not visible yet (404)' if status == '404' else status)
        print(f'  [{time.strftime("%H:%M:%S")}] {shown}')
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
    out, prefix, names = read_logs(result)
    if not out:
        # This is howto_MAAP_ogc QD failing: print everything, so the answer
        # to "where do OGC job logs go" can be read off and recorded.
        print(f'\njob finished {status}, but no log could be read.')
        print('get_job_status:'); print(json.dumps(record, indent=2))
        print(f'get_job_result (HTTP {r.status_code}):')
        print(json.dumps(result, indent=2) if not isinstance(result, str) else result)
        sys.exit(2)
    print(f"logs             : {prefix}/ {{{', '.join(names)}}}")
    digest = image_digest(out)
    print(f'image digest     : {digest or "<no docker pull in the log>"}')
    spec = container_of(maap, tag) if tag else None
    if spec:
        print(f"runner container : {spec.get('url') or spec.get('id')}")

    ok = status.lower() == 'successful'
    fields = parse_build_id_line(out)
    if fields is None:
        if ok:
            print('\nVERDICT: NO STAMP -- the job succeeded but printed no'
                  ' BUILD_ID line: this image predates run.sh --build-id.')
            sys.exit(1)
        # A failed job with no report says nothing about the image.  Never
        # turn "the log is missing" into "the image is stale".
        print(f'\nJOB {status.upper()} BEFORE run.sh PRINTED A BUILD ID -- no verdict'
              ' on the image.  The runner said:')
        print('  ' + '\n  '.join(runner_errors(out) or ['<no error lines found>']))
        sys.exit(2)

    # The whole block: the per-field detail is what makes a mismatch
    # diagnosable rather than merely visible.
    start = out.find('  ATL1415 build id')
    end = out.find('\n', out.find('BUILD_ID:'))
    print('\n' + out[start if start >= 0 else 0:end if end > 0 else None].rstrip())
    print('=' * 58 + '\n')

    code, lines = verdict(fields, want, origin)
    print('\n'.join(lines))
    if not ok:
        # The report is valid -- run.sh printed it -- but a job that fails
        # after a correct report still fails every tile the same way.  Not
        # green.
        print(f'\nBUT THE JOB ENDED {status.upper()} after the report.  The runner said:')
        print('  ' + '\n  '.join(runner_errors(out) or ['<no error lines found>']))
        code = max(code, 1)
    sys.exit(code)


if __name__ == '__main__':
    main()
