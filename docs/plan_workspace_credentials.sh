# plan_workspace_credentials.sh -- jobs reach the workspace bucket with
# MAAP-brokered keys, not the worker's IAM role
#
# ############################################################################
# ##  WRITTEN 2026-10-01.  TENTATIVE: nothing below has been coded or run. ##
# ##  QW1-QW4 ANSWERED by Ben 2026-10-01 (all as recommended): design A.  ##
# ##  Coding (W1) waits for the GL north monthly run to finish.            ##
# ############################################################################
#
# Status tags per step: TODO / DONE / BLOCKED.  DECIDED = Ben said so;
# RECOMMENDATION = mine; STATEMENT = checked fact, with where it came from;
# ASSUMPTION = taken as true until a run says otherwise; QUESTION = open.
#
# WHY.
#   STATEMENT (MAAP admin, relayed by Ben 2026-10-01): "I would recommend
#     using the maap-py workspace credentials method to get your temp keys
#     https://docs.maap-project.org/en/latest/system_reference_guide/accessing_bucket_data.html.
#     This allows you to control the expiration and refreshing of
#     credentials as you see fit for your job.  The worker IAM role will soon
#     be deprecated (to keep jobs properly isolated and prevent unauthorized
#     access from within jobs)."
#   STATEMENT (code, grep 2026-10-01): every workspace-bucket read and write
#     in a job uses botocore's DEFAULT credential chain, which on a worker
#     ends at that IAM role (instance metadata service):
#       run.sh:375                   args file fetch, s3fs.S3FileSystem()
#       scripts/s3_tiles.py:167      tile and product get/put
#       ATL1415/paths.py:81,97,121   pc.io_utils.get_s3fs(daac=None)
#       pointCollection io_utils     get_s3fs(daac=None): open_remote,
#                                    glob_remote; tilingSchema; geoIndex
#                                    (ATL11 index); grid/data.py /vsis3/
#                                    (GDAL, same chain)
#     NOT on that chain, unchanged by this plan: NSIDC granules (keys from
#     maap.aws.earthdata_s3_credentials, passed to s3fs explicitly) and the
#     tide model s3://pytmd (anonymous).
#   STATEMENT: when the role goes, every one of those calls fails, in every
#     step (prelim, matched, mosaic200, mosaic, nc).
#   STATEMENT (plan_dps_mosaic.sh D5-D7): the intermittent botocore
#     NoCredentialsError of 2026-10-01 is the same route failing ~0.8% of
#     the time; run.sh's AWS_METADATA_SERVICE_* settings (0b29127) retry it.
#
# WHAT THE BROKER RETURNS.  STATEMENT, measured on the ADE 2026-10-01T22:05Z
#   with maap-py 5.1.0a2, `MAAP().aws.workspace_bucket_credentials()`:
#     credentials: aws_access_key_id, aws_secret_access_key,
#       aws_session_token, expires_at
#       (NOT accessKeyId/secretAccessKey/sessionToken as maap-py's own
#       docstring says -- those are the earthdata_s3_credentials names)
#     expires_at was 2026-10-02T10:05:12Z: 12 h 00 min after the call.
#     authorized_s3_paths:
#       s3://maap-ops-workspace/ben_smith            read_write
#       s3://maap-ops-workspace/shared/ben_smith     read_write
#       s3://maap-ops-workspace/shared               read_only
#       s3://maap-ops-workspace/dataset/triaged_job  read_only
#   STATEMENT (default_args/*.txt): every bucket path the args name is
#     under s3://maap-ops-workspace/ben_smith, apart from s3://pytmd.  So the
#     keys cover everything a job reads and writes.
#   NOT KNOWN: the lifetime of keys issued to a WORKER (the 12 h is the
#     ADE's); whether the call works on a worker at all (it needs MAAP_PGT,
#     which workers have -- check_build_id reports maap_pgt=set).
#   ASSUMPTION (Ben 2026-10-01): the keys outlive a single job.  One fetch
#     per job, no refresh.  Revisit if W6 shows otherwise.  (Longest job so
#     far: GL north ATL15 netCDF, 3752 s.)
#
# ===========================================================================
# QUESTIONS FOR BEN (answer inline)
# ===========================================================================
#   QW1  One fetch in run.sh, exported as the standard AWS variables?
#        RECOMMENDATION: yes (design A below).  run.sh asks the broker once,
#        before its first bucket read, and exports AWS_ACCESS_KEY_ID,
#        AWS_SECRET_ACCESS_KEY, AWS_SESSION_TOKEN.  botocore, s3fs and GDAL
#        all read those FIRST in their chains, so every site in the list
#        above uses them with no code change, and no process does a
#        metadata lookup again -- which also removes the NoCredentialsError
#        class.  The alternative (design B: broker inside pointCollection's
#        get_s3fs(daac=None), as it does for NSIDC) needs changes in pC,
#        run.sh and s3_tiles.py, does not cover GDAL, and only pays off if
#        keys must be refreshed mid-job, which the ASSUMPTION rules out.
#        QW1 answer (Ben 2026-10-01): as recommended.  DECIDED.
#
#   QW2  What if the broker call fails on a worker?
#        STATEMENT (triaged logs, 2026-10-01, GL north monthly, 556 jobs
#        submitted at once): of the first 156 failures, 60 are the NSIDC
#        broker call failing (-> 'NoneType' ... 'get_s3_filesystem'), 9 are
#        MAAP's own runner timing out on api.maap-project.org, 86 have no
#        logs.  The workspace call goes to the same API, so it will see the
#        same timeouts when many jobs start together.
#        RECOMMENDATION: retry the call a few times with a pause (say 5
#        tries, 10 s apart, each with a 30 s timeout); if it still fails,
#        STOP the job with one line naming the call and the last error --
#        no fallback to the worker role, which would hide the failure until
#        the role is removed (memory: fail loudly).  The job is then
#        resubmitted like any other failure.
#        QW2 answer (Ben 2026-10-01): as recommended.  DECIDED.
#
#   QW3  Turn the worker role OFF inside our jobs now?
#        RECOMMENDATION: yes -- export AWS_EC2_METADATA_DISABLED=true in
#        run.sh once the keys are exported.  Then a job that would have
#        fallen back to the worker role fails today, in our test, instead of
#        on the day MAAP removes the role; and W6 really tests the new
#        route.  The AWS_METADATA_SERVICE_* lines of 0b29127 then do nothing
#        and are removed.
#        QW3 answer (Ben 2026-10-01): as recommended.  DECIDED.
#
#   QW4  The ADE side (not jobs): scripts/maap/*.py, the howtos' `aws s3`
#        lines and the ~/my-private-bucket mount use the ADE's own role
#        (AWS_ROLE_ARN + web identity; env, 2026-10-01).  The admin's note is
#        about the WORKER role.  RECOMMENDATION: leave the ADE side alone;
#        ask the admin whether the ADE role is also going away.  (For you
#        to ask, if you want it settled.)
#        QW4 answer (Ben 2026-10-01): as recommended.  DECIDED.
#
# ===========================================================================
# DESIGN A (DECIDED, Ben 2026-10-01, QW1-QW3)
# ===========================================================================
#   scripts/workspace_credentials.py   NEW.  Calls
#     MAAP().aws.workspace_bucket_credentials() (retry per QW2), checks the
#     four credential fields are there, and prints to STDOUT three
#     `export AWS_...=` lines plus
#     `export ATL1415_WORKSPACE_CREDENTIALS_EXPIRE=<expires_at>`; to STDERR
#     one line: expiry and the authorized paths -- never a key.
#     --check <s3 uri>: exit 1, naming the uri and the authorized paths, if
#     the uri is not under a read_write path (run.sh gives it tile_prefix
#     and out_prefix), so a wrong prefix fails at the start of the job and
#     not at the upload.
#   run.sh, before the args-file fetch:
#       eval "$(conda run --no-capture-output -n "$env_name" \
#               python "${repo_dir}/scripts/workspace_credentials.py" ...)"
#     STDOUT goes into the eval, not the job log.  run.sh has no `set -x`;
#     the header block prints named inputs only, never the environment.
#     Then AWS_EC2_METADATA_DISABLED=true (QW3) and
#     AWS_DEFAULT_REGION=us-west-2 if unset (STATEMENT: WORKER lines show
#     az=us-west-2c; the ADE has AWS_REGION=us-west-2.  NOT TESTED whether
#     s3fs or GDAL's /vsis3/ need it with the metadata service off; W3b
#     shows it, and the line goes if they do not).
#   run.sh --build-id: also makes the call and prints
#     workspace_credentials=ok|FAILED and the lifetime in hours;
#     check_build_id.py reads it like maap_pgt -- anything but ok is a
#     NO-GO verdict, since no step could read the bucket.
#   run.sh, at the end of every step: one line,
#     "workspace credentials: N.N h left at exit" -- the lifetime data the
#     ASSUMPTION needs.
#   Nothing changes in ATL1415/, pointCollection, s3_tiles.py or the args.
#   OFF THE WORKER: on the ADE and on discover MAAP_PGT decides.  With
#     MAAP_PGT set (ADE) run.sh brokers as on a worker.  With it unset
#     (discover, a laptop) there is no broker: run.sh says so in one line
#     and leaves the environment alone -- local runs there read local paths.
#
# ===========================================================================
# STEPS
# ===========================================================================
# W0. [Ben] DONE 2026-10-01.  QW1-QW4 answered: all as recommended.
#
# W1. [code] DONE 2026-10-01.  scripts/workspace_credentials.py + tests
#     (tests/test_workspace_credentials.py, MAAP mocked): the export lines;
#     no key on stderr; a missing field -> exit 1 naming it; --check inside
#     and outside the authorized paths; retry then stop (QW2).
#
#     RESULT: 20 tests.  As built: no --left mode (run.sh works the exit
#     line out with `date`); the retry is 5 tries, 10 s apart, 30 s each
#     (SIGALRM, since maap-py's requests.get has no timeout).
#
# W2. [code] DONE 2026-10-01.  run.sh wiring (fetch, QW3 lines, exit line, --build-id
#     field) and check_build_id.py's verdict.  Suite passes.
#     RESULT: suite 251 passed, 2 skipped.  As built: use_workspace_credentials
#     runs before the args fetch (and before bench); the BUILD_ID line gains
#     workspace_credentials=ok|FAILED|unchecked; the report gains a
#     `workspace_credentials=ok lifetime_h=N.N` line; check_build_id.py adds
#     VERDICT: NO WORKSPACE CREDENTIALS (exit 1) unless the field is ok --
#     so an image built before this change now fails the check too.  The
#     AWS_METADATA_SERVICE_* lines are gone (QW3).  NO AWS_DEFAULT_REGION
#     line: W3 shows s3fs and GDAL work without a region.
#
# W3. [ADE] DONE 2026-10-01 (result below).  Prove the keys alone are enough, off the worker: in a
#     shell with the ADE role variables (AWS_ROLE_ARN,
#     AWS_WEB_IDENTITY_TOKEN_FILE) UNSET and AWS_EC2_METADATA_DISABLED=true,
#       a. control, no keys: an s3fs list of the IS prefix FAILS
#          (NoCredentialsError) -- shows the shell really has no other route
#       b. with the script's exports: list + get + put + delete under
#          s3://maap-ops-workspace/ben_smith/scratch/ all work; a GDAL
#          /vsis3/ read of one mask works
#       c. a put outside ben_smith is refused (AccessDenied)
#     Then run.sh --step mosaic for one IS group from that shell, to a
#     scratch out_prefix; compare with the D6 file; delete the scratch.
#     RESULT (shell with AWS_ROLE_ARN and AWS_WEB_IDENTITY_TOKEN_FILE unset,
#     AWS_EC2_METADATA_DISABLED=true, no AWS config files):
#       a. no keys: s3fs list -> NoCredentialsError.  The shell has no other
#          route.
#       b. with the exports: list (7 entries), put, get, delete under
#          .../ben_smith/scratch/; pointCollection glob_remote finds the 28
#          IS matched tiles; GDAL /vsis3/ opens and reads GL_Ed2z0dx2.tif
#          (35 x 63).  The same with AWS_REGION and AWS_DEFAULT_REGION unset.
#       c. put to .../shared/ and .../dataset/ -> PermissionError; --check
#          on a path outside the workspace -> exit 1, naming the path.
#       run.sh --build-id: workspace_credentials=ok lifetime_h=12.0.
#       run.sh --step mosaic, IS avg_dz_40000m, to a scratch out_prefix:
#          exit 0, 16 s; dz_40km.h5 identical in every dataset to D6's; log
#          has the summary line and "12.0 h left at exit"; no key in the log
#          (grep).  Scratch prefix deleted.
#       broker unreachable (MAAP_API_HOST pointed at a closed port): 5
#          attempts in 43 s, then exit 1 with
#          "ERROR: no workspace credentials (above); nothing was read or
#          written"; nothing uploaded.
#     NOT TESTED HERE: a worker (W4, W5); a broker call that hangs for the
#     whole 30 s on a real network (the unit test covers the cut-off).
#
# W4. [Ben] TODO.  Push; register when nothing is in flight;
#     check_build_id.py -> MATCH and workspace_credentials=ok.  RECORD the
#     lifetime a worker is given.
#
# W5. [DPS] TODO.  Smoke, one job per kind of bucket use, IS (small):
#       prelim  one tile   (args fetch, ATL11 index, masks via GDAL, NSIDC
#                           keys alongside, tile put)
#       matched one tile   (prelim tiles get, tile put)
#       mosaic  z0         (list + in-place reads, product put)
#       nc      ATL14      (mosaics get, prelim reads, product put)
#     to .../ATL14_processing/rel006_0332_testing/creds/IS.  GATE: all
#     successful; outputs identical to the canonical IS tiles and the D6
#     products (weighted mosaics at rounding); every log has the
#     "h left at exit" line.
#
# W6. [DPS] TODO.  The lifetime under a real run: the next fan-out that is
#     due anyway (not a run made for this).  RECORD min "h left at exit"
#     over all jobs, and broker-call failures by cause.  GATE for the
#     ASSUMPTION: no job ends with less than 1 h left; otherwise reopen
#     QW1 (refresh, design B).
#
# W7. [docs] TODO.  howto_MAAP_ogc.sh (credentials section),
#     Transition_to_maap.md, run.sh header; ask the admin the removal date
#     of the worker role and note it here.
