"""
Shared helpers for the ATL1415 scripts that talk to MAAP's OGC job system:
check_build_id.py, submit_AA_queue.py and collect_jobs.py.

One copy, because each of these facts was learned the hard way and a second
copy would drift: where a job's log is (_stderr.txt, not _stdout.txt, under
the CWL runner), how the job result spells its output prefix (endpoint-style
s3://), that the process id changes on every redeploy, and which runner
error lines mean nothing.  See docs/howto_MAAP_ogc.sh, QD and O5-O7.

Scripts in this directory import it directly (`from ogc_jobs import ...`);
running one as a file puts this directory on sys.path.
"""
import os
import re
import subprocess
import sys

import requests
import yaml

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
    after the scheme as the bucket.  Same fix as collect_jobs.py.
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
