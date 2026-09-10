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
import sys
import time

from maap.maap import MAAP

# The shared OGC helpers live beside this file.  Put this directory on the
# path explicitly, so an importlib load (as the tests do) works too, not only
# `python scripts/maap/check_build_id.py`.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ogc_jobs import (DEFAULT_QUEUE, DONE, POLL_S, S3_RUN,  # noqa: E402
                      container_of, cwl_facts, find_process, image_digest,
                      job_id_from, load_config, origin_tip, read_logs,
                      runner_errors)


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
