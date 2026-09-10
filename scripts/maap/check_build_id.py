#!/usr/bin/env python3
"""
Submit ONE build_id job and say whether the image DPS would run is the commit
you think it is.

WHY THIS EXISTS.  DPS clones repository_url at algorithm_version and bakes the
result into a container, so the ADE working copy has nothing to do with what a
worker runs, and until run.sh grew --build-id no job log said which commit it
carried.  A stale image was therefore indistinguishable from a fix that did not
work: on 2026-09-09 an algorithm_version that had already been built once ran
the OLD code -- MAAP support has since confirmed that is not expected, and the
platform environment has been updated.  This is the cheap test of whether it
stayed fixed.  Run it after every register_algorithm.py, before spending worker
hours on a run whose results you would have to throw away.

READING THE VERDICT:
  MATCH        the image carries the commit you expected.  Proceed.
  MISMATCH     the build cloned a different commit -- the tag-reuse failure.
               Re-register, and if it persists bump algorithm_version.
  NO STAMP     the image predates the build stamp (2026-09-10).  Since every
               build from that commit on writes one, an unstamped image is
               ITSELF evidence of a stale container rather than a fresh build.

Usage:
  check_build_id.py [args_file_url] [queue] [--expect <sha>] [--timeout <s>]

The args_file is a required DPS input on this algorithm, so one must be named
even though run.sh exits before it ever looks in input/; the AA release args
file is the default simply because it is already published.
"""
import os
import re
import subprocess
import sys
import time

from maap.maap import MAAP

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
S3_RUN = 's3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/south/AA'
POLL_S = 15


def algorithm_version():
    """The container tag AND the git ref DPS clones; never hard-code it."""
    with open(os.path.join(REPO, 'algorithm_config.yml')) as fh:
        for line in fh:
            m = re.match(r'^algorithm_version:\s*(\S+)', line)
            if m:
                return m.group(1)
    raise RuntimeError('algorithm_config.yml has no algorithm_version')


def expected_commit(version):
    """
    What the build SHOULD have cloned: the remote tip of algorithm_version.

    The remote, not the local HEAD -- DPS clones from GitHub, so a commit that
    is local-only was never a candidate.  register_algorithm.py makes the same
    distinction when it refuses to register unpushed work.
    """
    for ref in (f'refs/heads/{version}', f'refs/tags/{version}'):
        p = subprocess.run(['git', '-C', REPO, 'ls-remote', 'origin', ref],
                           capture_output=True, text=True)
        if p.returncode == 0 and p.stdout.strip():
            return p.stdout.split()[0]
    return None


def normalize_s3(uri):
    """
    Drop the endpoint element from the URI getJobResult hands back.

    It arrives as s3://s3-us-west-2.amazonaws.com:80/maap-ops-workspace/... --
    host:port first, bucket second -- and `aws s3 cp` reads the first element
    after the scheme as the bucket.  Same fix as collect_AA_queue.py.
    """
    rest = uri[len('s3://'):]
    head, _, tail = rest.partition('/')
    if 'amazonaws.com' in head:
        rest = tail
    return 's3://' + rest.rstrip('/')


def stdout_for(maap, job_id):
    try:
        res = maap.getJobResult(job_id)
    except Exception:
        return ''
    for item in (res if isinstance(res, list) else [res]):
        if isinstance(item, str) and item.startswith('s3://'):
            base = normalize_s3(item)
            p = subprocess.run(['aws', 's3', 'cp', f'{base}/_stdout.txt', '-'],
                               capture_output=True, text=True, timeout=120)
            if p.returncode == 0:
                return p.stdout
    return ''


def status_of(maap, job_id):
    """getJob returns a dict for a failed job and a DPSJob for a live one."""
    try:
        info = maap.getJob(job_id)
    except Exception:
        return None
    if isinstance(info, dict):
        return info.get('status')
    info.id = job_id
    try:
        info.retrieve_attributes()
    except Exception:
        pass
    return getattr(info, 'status', None)


def main():
    argv = list(sys.argv[1:])
    expect, timeout = None, 1800
    for flag, cast in (('--expect', str), ('--timeout', int)):
        if flag in argv:
            i = argv.index(flag)
            value = cast(argv[i + 1])
            argv[i:i + 2] = []
            if flag == '--expect':
                expect = value
            else:
                timeout = value

    args_url = argv[0] if argv else f'{S3_RUN}/input_args_AA.txt'
    queue = argv[1] if len(argv) > 1 else 'maap-dps-worker-32gb'
    version = algorithm_version()
    want = expect or expected_commit(version)

    print(f'algorithm_version : {version}')
    print(f'expected commit   : {want or "<could not resolve from origin>"}')
    print(f'args file         : {args_url}')
    print(f'queue             : {queue}\n')

    maap = MAAP(maap_host='api.maap-project.org')
    job = maap.submitJob(
        identifier=f'ATL1415_build_id_{int(time.time())}',
        algo_id='ATL1415_tile_solve', version=version,
        queue=queue, queue_name=queue,
        x0=0, y0=0, step='build_id', args_file=args_url)
    job_id = getattr(job, 'id', None) or getattr(job, 'job_id', None)
    print(f'job_id            : {job_id}\n')
    if not job_id:
        sys.exit(f'submit returned no job id: {job!r}')

    deadline = time.time() + timeout
    status = None
    while time.time() < deadline:
        status = status_of(maap, job_id)
        print(f'  [{time.strftime("%H:%M:%S")}] {status}')
        if status and status.lower() in ('succeeded', 'failed', 'deleted',
                                         'dismissed'):
            break
        time.sleep(POLL_S)
    else:
        sys.exit(f'timed out after {timeout}s; last status {status}')

    out = stdout_for(maap, job_id)
    if not out:
        sys.exit(f'job finished {status} but its _stdout.txt could not be read')

    # Print the whole block: the per-field detail is what makes a mismatch
    # diagnosable rather than merely visible.
    start = out.find('  ATL1415 build id')
    print('\n' + (out[start:] if start >= 0 else out).rstrip())

    line = re.search(r'^BUILD_ID: .*$', out, re.M)
    if not line:
        sys.exit('\nNO BUILD_ID LINE -- this image predates run.sh --build-id, '
                 'which itself means the build did not pick up current code.')
    got = re.search(r'commit=(\S+)', line.group(0))
    got = got.group(1) if got else 'unknown'

    print()
    if got == 'unknown':
        print('VERDICT: NO STAMP -- the image was built before the build stamp '
              'existed, so it is not a fresh build of current code.')
        sys.exit(1)
    if want and got != want:
        print(f'VERDICT: MISMATCH\n  image  {got}\n  wanted {want}\n'
              '  The build cloned a different commit than origin/'
              f'{version} points at.  Re-register; if it recurs, the tag is '
              'being reused and algorithm_version has to be bumped.')
        sys.exit(1)
    print(f'VERDICT: MATCH -- image is at {got}')


if __name__ == '__main__':
    main()
