# ===========================================================================
# howto_MAAP_ogc.sh -- MOVING ATL1415 TO MAAP'S OGC/CWL ALGORITHM SYSTEM
# ===========================================================================
#
#   ********************************************************************
#   **  TENTATIVE.  Written 2026-09-10, BEFORE any of the code it     **
#   **  describes exists and before any step has been run.  It is the **
#   **  plan and the acceptance criteria, not a record of a run.      **
#   **  Every step carries its own status tag; revise as each lands.  **
#   ********************************************************************
#
# WHY.  On 2026-09-10 the ADE notebook env came back with maap-py 5.1.0a2 and
# register_algorithm.py died: AttributeError: 'MAAP' object has no attribute
# 'register_algorithm_from_yaml_file'.  Offered three ways forward -- stay on
# the legacy path via the ATL14 env's maap-py 4.2.0, migrate, or both in
# sequence -- BEN CHOSE TO MIGRATE: "New path, CWL + prebuilt image".
#
# Steps are numbered O1..O12 so they can be cited ("OGC step 3").  Status
# tags: [OK] done and verified, [UNTESTED] ready but not run, [NEEDS CODE: x].
# Where each runs: [ADE] this workspace, [DPS] a worker, [BEN] needs you.
#
# This file REPLACES the registration mechanism in howto_MAAP_staging.sh S5
# and the job calls in S5b and howto_MAAP_AA step 3b-i; those now point here.
#
#
# ---------------------------------------------------------------------------
# FINDINGS -- verified statements, each with how it was verified
# ---------------------------------------------------------------------------
# F1. THE CLIENT DROPPED THE WHOLE LEGACY API, not one method.  maap-py
#     5.1.0a2 (Development Status: Alpha; installed from git tag v5.1.0a2)
#     has no register_algorithm_from_yaml_file, and submitJob, getJob,
#     getJobResult, getJobMetrics, getQueues and listAlgorithms all raise
#     "no longer supported".  Replacements: submit_job(process_id, inputs,
#     queue), get_job_status, get_job_result, get_job_metrics, get_queues.
#     -- read in site-packages/maap/maap.py.  CONSEQUENCE: register_algorithm,
#     submit_AA_queue, collect_AA_queue and check_build_id all break under it.
#
# F2. THE LEGACY PATH STILL WORKS, SERVER-SIDE.  /api/mas/algorithm is alive
#     and still lists ATL1415_tile_solve:on_s3 and :on_s3_v2.  -- GET,
#     2026-09-10.  CLIENT-SIDE it now needs assembling: ATL14 moved to maap-py
#     5.1.0 the same day (QE), and register_algorithm.py was rewritten for OGC
#     (O3), so the fallback is the legacy script and config from git under a
#     throwaway 4.2.0 env:
#       python3 -m venv ~/maap42 && ~/maap42/bin/pip install maap-py==4.2.0
#       git show f3d5049:register_algorithm.py > /tmp/register_legacy.py
#       git show f3d5049:algorithm_config.yml  > /tmp/algorithm_config_legacy.yml
#       ~/maap42/bin/python /tmp/register_legacy.py /tmp/algorithm_config_legacy.yml
#     FROM /tmp ITS PUSH CHECKS DO NOT RUN: it takes REPO_DIR from its own
#     path, finds no git checkout there, and says "NOTE: ... skipping the push
#     checks."  So do them by hand first -- both must print nothing:
#       git -C ~/git_repos/ATL1415 status --short
#       git -C ~/git_repos/ATL1415 log --oneline origin/on_s3..on_s3
#
# F3. OUR ALGORITHM WAS NOT MIGRATED.  /api/ogc/processes lists 45 processes
#     and none is ATL1415.  -- list_algorithms(), 2026-09-10.
#
# F4. THERE ARE TWO WAYS TO GET AN IMAGE INTO THE OGC SYSTEM, and they split
#     the 45 cleanly.  -- every cwlLink fetched, 2026-09-10:
#       25  image on ghcr.io, CWL on raw.githubusercontent.com: the owner's
#           own GitHub Action (MAAP-Project/ogc-app-pack-generator) builds
#           from a Containerfile.  No build_command -- it would mean replacing
#           build-env.sh with a Dockerfile.
#       16  image at mas.maap-project.org/root/ogc-application-packages/
#           <user>/<name>:<version>, CWL on repo.maap-project.org: MAAP builds
#           it from the repo.  MAAP's own dps_tutorial is built this way from
#           a repo with a build-env.sh and NO Containerfile.
#
# F5. THE MAAP-BUILT PATH TAKES OUR CURRENT MODEL UNCHANGED.  It is what the
#     Algorithm Catalog's "Register New Algorithm" form drives, and its fields
#     are ours: code_repository, algorithm_version (the branch), build_command
#     (build-env.sh), base_container_url (maap_base), run_command (run.sh).
#     A config needs EITHER algorithm_container_url (prebuilt) OR
#     base_container_url + build_command.  -- docs.maap-project.org/en/ogc
#     dps_tutorial_demo, and the plugin's own source
#     (maap_algorithms_jupyter_extension 1.0.2, static/156.*.js).
#
# F6. AND IT IS SCRIPTABLE -- the form is a thin client over a REST call:
#       POST https://api.maap-project.org/api/build   body: JSON config
#            -> {"status": "accepted", "build_id": ...}
#       GET  .../api/build/{build_id}
#            -> status, pipelineLink, deploymentLink{href}, deploymentError
#       GET  .../api/ogc/deploymentJobs/{deployment_id}
#            -> status, pipeline.processPipelineLink, error
#     The JSON the form sends (empty fields dropped; outdir_max forced to 20;
#     outputs forced to one Directory named "out"):
#       algorithm_name, algorithm_version, algorithm_description,
#       code_repository, run_command, build_command, base_container_url,
#       ram_min, cores_min, author, contributor, license, release_notes,
#       citation, keywords,
#       inputs:  [{name, label, doc, type, default}],
#       outputs: [{name: out, type: Directory}],  outdir_max: 20
#     -- plugin source.  GET /api/build with maap-py's own _get_api_header()
#     answers HTTP 200 {"builds": []} from a plain script.  -- 2026-09-10.
#     This is also where Ben's "other URLs" come from: pipelineLink,
#     deploymentLink and deploymentPipelineLink are three separate pages.
#
# F7. OUR NAME IS INVALID THERE.  The plugin rejects any algorithm_name not
#     matching ^[a-z0-9_-]+$ -- ATL1415_tile_solve has capitals.
#     algorithm_version must match ^[a-zA-Z0-9_][a-zA-Z0-9._-]{0,127}$, which
#     on_s3 does.  ram_min is capped at 128 GB and cores_min at 32.
#     -- plugin validators and tooltips.
#
# F8. THE CALLING CONVENTION CHANGES.  In every generated CWL, each input is
#     bound as a PREFIXED option, e.g. shah_dps_tutorial's
#       inputBinding: {position: 1, prefix: --input_file}
#     so run.sh would be invoked as
#       run.sh --x0 220000 --y0 20000 --step prelim --args_file <path>
#     where today it expects `run.sh 220000 20000 prelim` and FINDS the args
#     file in input/.  -- the fetched CWLs.  CONSEQUENCE: run.sh has to change.
#     CONFIRMED by our own generated CWL (O4): baseCommand /app/ATL1415/run.sh,
#     inputs bound --x0, --y0, --step, --args_file at positions 1-4.
#     Outputs are collected by `glob: ./output*`, which run.sh's output/
#     already satisfies.
#
# F9. JOBS ADDRESS A PROCESS BY ID, NOT BY name:version.  submit_job() POSTs
#     to /api/ogc/processes/<process_id>/execution, and each deployed process
#     has a numeric processID beside its string id.  -- maap.py and
#     list_algorithms().  THE NUMERIC ONE, by the deployment record: its
#     processLocation is /ogc/processes/64, and the string id is not unique
#     across versions.  O6 is the first real submit.
#
# F10. THE GENERATED CWL RECORDS WHAT IT BUILT: s:codeRepository,
#     s:commitHash and s:version.  -- shah_dps_tutorial's CWL.  CONSEQUENCE: a
#     SERVER-SIDE build id to set beside our own stamp, which costs no job.
#
# F11. WORKERS USE maap-py FOR ONE CALL, maap.aws.earthdata_s3_credentials
#     (the NSIDC ATL11 credentials), and 5.1.0a2 still has it.
#     -- environment.yml comment and AWS.py.
#
#
# ---------------------------------------------------------------------------
# OPEN QUESTIONS -- each says who can answer it
# ---------------------------------------------------------------------------
# QA. [TEST, O4] Can a File input be an s3:// URI that the runner stages?
#     The tutorials only show https://.  RECOMMEND: sidestep it -- declare
#     args_file as a `string` and have run.sh fetch an s3:// URI itself with
#     `aws s3 cp`, so the answer does not matter.
#
# QB. [ANSWERED 2026-09-10 by O6 run 2: YES, the rebuild of on_s3 ran the
#     new commit, and again at run 3.  O11 closes the rest.]  Does a SECOND build of the same algorithm_version produce a
#     fresh image?  That is the bug that started today (8aad07d); nothing
#     says the OGC path is immune to it.  The build stamp answers it with one
#     job, and s:commitHash (F10) answers it with none.
#
# QC. [BEN / MAAP] Is the legacy /api/mas path being retired, and when?  It
#     decides whether the ATL14 fallback (F2) is safe to keep.
#
# QD. [ANSWERED 2026-09-10: failed job by O6 run 1, successful by run 2 --
#     see O6 for both paths]
#     ON THIS SYSTEM A JOB RUNS INSIDE MAAP'S CWL RUNNER
#     (container-maap-cwltool-executor:v1.1.0 -- which is also what the job
#     record's container_specification names, NOT our image).  The runner's
#     own log is _stdout.txt; OUR CONTAINER'S OUTPUT IS IN _stderr.txt.  A
#     failed job's files are under
#       s3://maap-ops-workspace/dataset/triaged_job/v1.4.0/
#         triaged_job-<job_id>_task-<uuid>/
#     -- which is what get_job_result returned for it.  Consequence for O7:
#     the collector's Decimate_data / rusage lines will be in _stderr.txt.
#     (Superseded draft of this entry:)  What do get_job_status /
#     get_job_result return, and is _stdout.txt still on the bucket?
#     FINDING: the OGC job endpoints serve the SAME backend as the legacy
#     jobs -- list_jobs() returns the legacy AA transect jobs, and for one of
#     them (37a86437...):
#       get_job_status -> {"jobID", "processID", "type", "status": "successful"}
#       get_job_result -> {"<id>": {"links": [{href: website}, {href:
#                          s3://s3-us-west-2.amazonaws.com:80/maap-ops-
#                          workspace/ben_smith/dps_output/...}, {href: console}]}}
#       and _stdout.txt, _stderr.txt sit at that prefix.
#     What is NOT yet seen: a job submitted THROUGH submit_job to an OGC
#     process.  The first build_id job (O6) is that test.
#
# QE. CLOSED 2026-09-10 by O6 run 1: an OGC worker HAS MAAP_PGT --
#     the runner starts cwltool with --preserve-environment MAAP_PGT and
#     --preserve-environment MAAP_API_HOST, and the job printed
#     maap_pgt=set, maap_py=5.1.0.
#     SETTLED 2026-09-10 (Ben): maap-py 5.1.0 in environment.yml, for the
#     worker AND the ATL14 env.  The first draft of this entry recommended
#     staying on 4.2.0; checked against the code, its reasons did not hold:
#     - FINDING: the worker's whole maap-py path -- MAAP(), _get_api_header,
#       aws.earthdata_s3_credentials, requests_utils -- is CODE-IDENTICAL in
#       4.2.0, 5.1.0a2 and 5.1.0 (ASTs compared, docstrings stripped).  The
#       one difference: config_reader looks up three OGC endpoints, which the
#       server already serves.  Requires-Dist is identical.
#     - FINDING: 5.1.0 FINAL is on PyPI (`pip index versions --pre maap-py`),
#       so "5.x is an alpha" was only true of the ADE's 5.1.0a2.
#     - so the worker cannot tell the difference, and the ADE needs 5.x: the
#       ported job scripts (O5, O7) use submit_job & co., and the howtos run
#       them after `conda activate ATL14`.
#     COST: ATL14 no longer has the legacy calls -- the fallback is in F2.
#     THE RISK THAT REMAINS IS VERSION-INDEPENDENT: pointCollection asks MAAP
#     for NSIDC credentials ONLY when $MAAP_PGT is set, and silently falls
#     back to earthaccess (which has no credentials on a worker) when it is
#     not.  Legacy workers set it; OGC workers are unverified.  run.sh
#     --build-id now reports maap_pgt=set|unset and maap_py=<version>, so the
#     O6 job answers it before any tile reads ATL11.
#
# QF. [BEN] The new name.  RECOMMEND atl1415_tile_solve: the old name
#     lowercased, so every log, ledger and doc stays greppable for it.
#
#
# ===========================================================================
# O1. [OK, 2026-09-10]  algorithm_config.yml in the new schema.   [ADE]
# ===========================================================================
# Same file, same information, the api/build field names (F6):
#   algorithm_name      ATL1415_tile_solve      ->  atl1415_tile_solve  (QF)
#   repository_url                              ->  code_repository
#   docker_container_url (maap_base:v6.0.0)     ->  base_container_url
#   build_command, run_command                      unchanged
#   queue: maap-dps-worker-32gb, disk_space      ->  ram_min / cores_min, and
#                                                   the queue moves to submit
#   inputs.positional + inputs.file             ->  one `inputs` list of
#                                                   {name, label, doc, type,
#                                                    default}:
#       x0, y0      string   (not int -- negative, and passed straight through)
#       step        string, default prelim
#       args_file   string   (QA: a URI that run.sh fetches, not a File)
# RECOMMEND the file is REWRITTEN rather than duplicated: it is exactly what
# the Algorithm Catalog's "Load Algorithm Configuration" button reads, so
# the script and the UI stay interchangeable.  The legacy form stays
# recoverable from git (e231ba4..f3d5049) for the F2 fallback.
# ACCEPTANCE: the plugin's own validators would accept it -- name and
# version regexes (F7), and base_container_url + build_command present.
# DONE: rewritten as above, name atl1415_tile_solve (QF, as recommended --
# say if you want another), ram_min 16 / cores_min 4 as FLOORS (the queue
# still sizes the worker; untested assumption, noted in the file).
# register_algorithm.py now runs the plugin's validators itself: the new file
# passes, and the legacy one (f3d5049) is refused with six reasons, the
# capitalised name first.


# ===========================================================================
# O2. [OK, 2026-09-10 -- locally, and ON A WORKER in both O6 runs]  The CWL calling convention.   [ADE]
# ===========================================================================
#   run.sh --x0 V --y0 V --step V --args_file V      (the CWL binding, F8)
#   run.sh V V V                                     (kept: local runs and
#                                                     the F2 fallback)
# args_file may be a local path OR an s3:// URI; a URI is copied into input/
# with `aws s3 cp` first, so everything downstream of that is unchanged.
# The build-id pre-scan is already token-based, so `--step build_id` is
# caught before any parsing, like every other form.
# ACCEPTANCE: locally, all three forms of build_id still exit 0; a
# prefixed prelim call reaches ATL11_to_ATL15 with the same argv the
# positional form produces today.
# MET, with the solver stubbed to print its argv: positional, prefixed with a
# local path, and prefixed with the REAL s3:// URI of the published AA args
# file (fetched with s3fs, 41 lines) give BYTE-IDENTICAL argv.  build_id exits
# 0 in all five spellings tried, and eight malformed calls exit 2 with a
# message.  Testing it turned up two bugs, both fixed:
#   - --step=build_id was refused with "must be ... 'build_id', got
#     'build_id'": the pre-scan only knew the bare token.  The step check now
#     handles build_id too.
#   - legacy positionals with no input/ exited 1 with no message: under
#     pipefail, find failing on a missing dir killed the assignment.  Older
#     than today; invisible on legacy DPS, where input/ always existed.


# ===========================================================================
# O3. [OK, 2026-09-10 -- USED FOR THE REAL REGISTRATION IN O4]  POST api/build.
# ===========================================================================
# Keep what the script is for -- REFUSE to register work that is not on
# GitHub -- and swap the transport:
#   - read algorithm_config.yml, validate it with the rules the Algorithm
#     Catalog form enforces (F7), send it as JSON to POST /api/build using
#     maap-py's _get_api_header() (F6); print the build_id.
#   - print the FIRST URL in the response on a line of its own, as the one to
#     open -- the build pipeline -- then any others, labelled by key path.
# Works under maap-py 4.2.0 AND 5.x: it uses only the auth header, not
# either version's algorithm methods.  --dry-run prints the exact JSON.
# Exit codes: 0 accepted, 1 invalid config or unpushed work, 2 rejected, or
# accepted with no URL to open.
#
# REVISED 2026-09-10, after the first real registration: the draft also had
# `--status <build_id> [--wait]`, which followed the build into its
# deployment job and printed every link at each stage.  REMOVED, per Ben:
# "The first url returned by release 5 of maap_py is the correct one.
# There's no need for the --status --wait step."  The pipeline page is where
# the build and its deployment are followed.  What --status also did --
# look up the processID of the deployed process -- moves to O5, the first
# code that needs it.


# ===========================================================================
# O4. [OK, 2026-09-10 -- Ben registered; built and DEPLOYED]  Register.   [ADE]
# ===========================================================================
/srv/conda/envs/notebook/bin/python register_algorithm.py
# and open the first URL it prints.  That is the whole step.
#
# RECORD of the first registration, read back from GET /api/build
# (read-only) just after Ben ran it:
#   build_id      18533aaa-19a5-4720-a66c-26c454c750ab
#   status        running   (updated 2026-09-10T17:59:09)
#   repository    https://github.com/SmithB/ATL1415.git @ on_s3
#   pipelineLink  https://repo.maap-project.org/root/build-ogc-app-pack/-/pipelines/20198
#                 -- "Link to build pipeline", rel=monitor: THE URL TO OPEN.
#                 The builds run in MAAP's GitLab project
#                 root/build-ogc-app-pack.
#   links         {href: /build/<id>, rel: self} -- RELATIVE, so not
#                 something a browser can open, and never printed.
#   no deploymentLink yet while the build is running.
# ODDITY, not acted on: the record says created 2026-09-04T17:15:15, though
# GET /api/build returned NO builds earlier on 2026-09-10.  Whatever `created`
# measures, it is not when this registration happened -- do not time builds
# with it.
# DEPLOYED, read back the same way:
#   build status  successful; still running at 17:59:09
#   deployment    ogc/deploymentJobs/140, created 18:07:06, successful --
#                 pipeline https://repo.maap-project.org/root/deploy-ogc-hysds/-/pipelines/20199
#   process       processID 64, id atl1415_tile_solve, version on_s3,
#                 lastModified 18:08:54, deployedBy ben_smith
#   CWL           https://repo.maap-project.org/api/v4/projects/137/repository/files/ben_smith%2Fatl1415_tile_solve%2Fon_s3%2Fprocess.cwl/raw
#   s:commitHash  d401699 -- exactly origin/on_s3 when Ben registered
#   image         mas.maap-project.org/root/ogc-application-packages/
#                 ben_smith/atl1415_tile_solve:on_s3 -- the TAG IS THE
#                 BRANCH, so a rebuild of on_s3 reuses it (why QB is live)
# How long the build itself took is not known: Ben's registration time was
# not captured and `created` is unreliable (above).  Registration to
# deployed was under ~10 minutes.


# ===========================================================================
# O5. [OK, 2026-09-10 -- dry-run and log-reading REAL; submit path mocked]
#     Port check_build_id.py to the OGC job calls.   [ADE]
# ===========================================================================
# FIRST, FIND THE PROCESS: submit_job() needs the deployed process's id
# (F9).  Look it up by name and version from list_algorithms() --
# id == atl1415_tile_solve, version == on_s3 -- rather than hard-coding a
# number that changes with every redeploy.  (This was in register_algorithm's
# --status until O3 was revised.)
# submit_job(process_id, {x0, y0, step: build_id, args_file}, queue), then
# get_job_status / get_job_result (QD).  TWO comparisons instead of one:
# the stamp in the image against origin/<algorithm_version>, AND against the
# s:commitHash in the deployed CWL (F10).  All three agreeing is MATCH.
# DONE -- AS BUILT, one change to the paragraph above: image == cwl is the
# test, and origin is only REPORTED.  Requiring origin too would call a
# correct image stale whenever GitHub moves on after a build -- which is
# exactly why the 65c09bf push is being held until O6 has run.
#   - the process is found by name+version every run (processID 64 today).
#   - submit_job(64, {x0: '0', y0: '0', step: build_id, args_file}, queue,
#     dedup=False, tag=atl1415_build_id_<time>).  dedup=False EXPLICITLY: the
#     inputs are identical every time, and a deduplicated job would hand back
#     the PREVIOUS image's answer after a rebuild.  Default queue -32gb: the
#     CWL's ramMin is 16, which a 16 GB worker may not have to allocate.
#   - the log comes from get_job_result(): walk it for s3:// hrefs, normalize
#     the endpoint-style one, read <prefix>/_stdout.txt with `aws s3 cp`.
#   - best effort: the job record's container url + digest, via
#     list_jobs(tag=...).  A digest unchanged across a rebuild would mean a
#     reused image, whatever the tag says.
#   - verdicts MATCH / MISMATCH / NO STAMP / NO NSIDC (maap_pgt unset) exit
#     0/1/1/1; 2 when the check itself could not run.  --dry-run stops before
#     submitting.
# TESTED: --dry-run for real under maap-py 5.1.0a2 and 5.1.0 (finds 64, reads
# d401699 from the CWL and from GitHub).  Log reading for real, against the
# legacy AA_cost_44km_E220_N20 job through the OGC endpoints: finds its
# prefix, reads its 5070-char _stdout.txt, and the tag lookup returns its
# container digest.  Verdicts on real run.sh output (stamp present, stamp
# removed, MAAP_PGT removed).  submit -> poll -> log -> verdict, a refused
# submit, a failed job without a log, and an undeployed process, all with
# maap mocked -- the one path not yet run for real is submit_job itself,
# which is O6.


# ===========================================================================
# O6. [OK, 2026-09-10 -- run 2 clean: MATCH, exit 0, job successful]
#     One build_id job.   [ADE -> DPS]
# ===========================================================================
/srv/conda/envs/notebook/bin/python scripts/maap/check_build_id.py
# Answers QB (fresh image?) and QD (how logs come back), and QE's remaining
# risk: the BUILD_ID line must say maap_pgt=set, or no tile on this system
# can read NSIDC.  Nothing below runs until this says MATCH and maap_pgt=set.
#
# RUN 1, 2026-09-10 -- job db93c7f3-bc2b-45ad-8343-d5d3e60e6139 on -32gb:
# queued 18:21:44, running 18:23:45, FAILED 18:25:31.
#   WHAT IT ESTABLISHED, read from the job's own log:
#   - IMAGE == BUILT: the stamp inside the image says commit=d401699, clean
#     tree, built 17:59:17-18:03:14 -- the same commit as the CWL's
#     s:commitHash.  Live git in the image agrees.  The worker PULLED it
#     fresh ("Status: Downloaded newer image"), digest
#     sha256:9a91a554d305b88351c756e78b4b7f638a9517363309d0da542cac735d844962.
#   - maap_py=5.1.0, maap_pgt=set: QE's risk is closed (see QE).
#   - run.sh got --x0 0 --y0 0 --step build_id --args_file ..., exactly O2's
#     convention -- F8 and O2 confirmed on a worker.
#   WHY IT FAILED -- A BUG IN run.sh, MINE, FIXED: the build_id path exited
#   before `mkdir -p output`, and the CWL collects `glob: ./output*`, so
#   cwltool ended a job with a complete, correct report as permanentFail:
#   "Did not find output file with glob pattern: ['./output*']".  The
#   build_id path now makes output/ and tees its report into
#   output/build_id.txt, which also makes the answer an uploaded product.
#   AND A BUG IN check_build_id.py, MINE, FIXED: it read only _stdout.txt,
#   which on this system is the RUNNER's nine lines; our report is in
#   _stderr.txt (QD).  Finding no BUILD_ID line it announced "NO STAMP" --
#   calling a correct image stale.  It now reads both logs and build_id.txt,
#   never gives an image verdict without a report, and on a job that failed
#   after the report prints the verdict AND the runner's error.  Re-reading
#   this job with the fixed script (--job, no new submit): "VERDICT: MATCH"
#   then "BUT THE JOB ENDED FAILED", quoting the glob error; exit 1.
#   ALSO SEEN: for one poll after submit, get_job_status answered 404 --
#   the job was not yet visible.  Harmless; the poller now says so.
#   cwltool also warned that it SKIPS the container --memory and --cpus
#   limits despite ramMin/coresMin: on this system they are not enforced on
#   the container, which can use the whole worker.
# RUN 2, 2026-09-10, after Ben re-registered at 8935494 (build 2cedbffb,
# deploymentJobs/141; the redeploy KEPT processID 64) -- job f537a824 on
# -32gb: submitted 20:53:47, successful 20:57:49, exit 0.
#   - VERDICT: MATCH.  The stamp says commit=8935494, built 20:47:45, the
#     CWL's s:commitHash is 8935494, maap_py=5.1.0, maap_pgt=set.
#   - The run.sh fix works: build_id.txt was uploaded as a product.
#   - A SUCCESSFUL job's files are under
#       s3://maap-ops-workspace/ben_smith/dps_output/atl1415_tile_solve_1786/
#         on_s3/<yyyy>/<mm>/<dd>/<HH>/<MM>/<SS>/<usec>/
#     -- _stdout.txt, _stderr.txt, and the output/ contents at the top
#     level (build_id.txt beside them).  QD is closed.
#   - QB is answered: a SECOND build of on_s3 ran the NEW code.
#   - No image pull was logged (run 1's was).  Not a stale image -- see run 3.
# RUN 3, 2026-09-11, after Ben re-registered at ab84687 (built 15:01:35-
# 15:04:35, process modified 15:08:55, processID STILL 64) -- job c0232572
# on -32gb: submitted 15:40:44, running 15:42:45, successful 15:44:16.
#   - VERDICT: MATCH -- stamp, CWL s:commitHash and origin/on_s3 are all
#     ab84687; clean tree; maap_py=5.1.0, maap_pgt=set.
#   - Again no pull in either log, though no earlier job had run this build:
#     the image reaches a worker by a route the job logs do not show, so a
#     missing pull line says nothing about staleness.  Nothing to act on (O11).


# ===========================================================================
# O7. [OK, 2026-09-10 -- all REAL: dry-runs, the collector on 12 legacy
#     jobs, and the submitter's first submissions (O8)]  submit_AA_queue.py, collect_AA_queue.py   [ADE]
# ===========================================================================
# Same two scripts, new calls: submit_job with inputs as a dict and the
# queue as an argument (F1, F9); the collector reads the log wherever O6
# found it (QD).  The half-routing and the ledger format do not change.
# FROM O6: the solve's own lines (Decimate_data N_XO, the rusage lines)
# will be in _stderr.txt, not _stdout.txt; reuse check_build_id's
# read_logs() rather than a third copy of the log-reading code.
# DONE:
#   - scripts/maap/ogc_jobs.py: ONE copy of the process lookup, the log
#     reader (_stdout + _stderr + build_id.txt), the s3 prefix normalizer,
#     the job-id parser and the runner-error filter.  check_build_id.py now
#     imports them; re-reading job db93c7f3 and --dry-run behave identically.
#   - submit_AA_queue.py: submit_job(process_id, string inputs, queue,
#     dedup=False, tag=identifier), process found by name+version, ledger
#     format unchanged.  --dry-run, run for real: the full transect routes
#     16 centers to 17 jobs (E420 is in the 360-440 km overlap, both halves),
#     the O8 pair to two 44 km jobs.
#   - collect_AA_queue.py: get_job_status / get_job_metrics (wall clock only;
#     its machine fields are null) / get_job_result -> read_logs, and a new
#     N_AT / N_XO pair of columns from the Decimate_data line.
#   REAL-DATA CHECK: pointed at the twelve legacy transect jobs of 2026-09-09
#   04:52 (found with list_jobs -- the OGC endpoints serve them), it
#   reproduces scripts/maap/AA_cost_results.csv EXACTLY: n_atl11, n_fit, and
#   each step's seconds and GiB, 12/12.  And it reads N_XO=0 on all twelve --
#   the pre-fix baseline, for the whole transect, not just E220_N20.
#   FOUND ON THE WAY, older than today: howto_MAAP_AA's two submit commands
#   (3b, 3b-i) passed four arguments to a five-positional script, putting the
#   queue name in the 44 km args slot and the ledger in the queue slot --
#   dry-run showed queue=AA_xo_check_jobs.csv, args=maap-dps-worker-32gb.
#   Both commands fixed; the script now refuses that arrangement (exit 2).
#   RUN FOR REAL at O8: both submit_job calls accepted, job ids in the ledger.
#   AND (O11): a per-tile commit column, from the BUILD_ID line run.sh now
#   prints in every tile job, then a list of the builds the ledger's tiles
#   ran, and a warning when a worker lacked MAAP_PGT.  (Until 2026-09-11 it
#   also warned on mixed builds; see O11, O12b.)  Tested: mocked ledger (two
#   builds, MAAP_PGT unset, an unstamped job) lists both builds and warns
#   once; the real O8 ledger gives the same numbers, commit '-', no list.


# ===========================================================================
# O8. [OK, 2026-09-11 -- N_XO > 0 on both tiles]  The two pole-hole tiles.   [ADE -> DPS]
# ===========================================================================
# howto_MAAP_AA step 3b-i, unchanged in intent: E220_N20 and E300_N20, and
# N_XO > 0 on both.  E220_N20 read N_XO=0 before the crossover fix.
# SUBMITTED on build 8935494 (the O6 run-2 image) with the corrected howto
# command; queue -32gb:
#   AA_cost_44km_E220_N20  d2e9bcba-cdf5-4966-9f6e-5c7e62321645
#   AA_cost_44km_E300_N20  3273c127-961c-4b93-80f3-22c0708e5b56
# Ledger: ~/ATL14_processing/maap_ledgers/AA_xo_check_jobs.csv -- OUTSIDE
# the checkout on purpose: an untracked file there makes
# register_algorithm.py refuse.  Read it with
#   scripts/maap/collect_AA_queue.py ~/ATL14_processing/maap_ledgers/AA_xo_check_jobs.csv
# These two were built before per-tile stamping (O11), so their
# commit column reads '-'.  A cached image cannot fake this test: every
# image ever built on this system (d401699, 8935494) has the crossover fix.
# PRE-FIX BASELINE, from the collector: E220_N20 N_AT=935506 N_XO=0;
# E300_N20 N_AT=1256488 N_XO=0.
# E300_N20 DONE (successful, 4696 s): N_XO=52882 -- CROSSOVERS ARE READ ON
# A WORKER.  N_AT unchanged at 1256488; N_ATL11 5430970 -> 5483852, which is
# +52882 exactly, so every added point is a crossover; N_fit 1256487 ->
# 1309369.  Cost: fit 2672 -> 3205 s (+20%), peak 15.70 -> 16.32 GiB, error
# step 1425 -> 1481 s.
# E220_N20 DONE (successful, 11032 s = 3.1 h): N_XO=171488.  N_AT unchanged
# at 935506; N_ATL11 12168633 -> 12340121, +171488 exactly; N_fit 743118 ->
# 893044.  Cost: fit 4646 -> 6599 s (+42%), error 3302 -> 4423 s (+34%),
# whole job 2.2 -> 3.1 h; peak 20.93 -> 21.40 GiB on a 32 GiB worker.
# THE ANSWER: ATL11XO is read on a MAAP worker, and the pole-hole edge --
# where tracks converge -- gains the most (18% more input points there,
# 1% at E300_N20).  CORRECTED 2026-09-11: those two divide by different
# things.  18% is N_XO/N_AT; 1% was N_XO/N_ATL11.  On one basis, N_XO/N_AT,
# E300_N20 is 4.2% -- and its N_fit rose by exactly N_XO.
# CONSEQUENCE FOR SIZING (S6): every number in AA_cost_results.csv is
# pre-crossover and now too low, most of all near the pole.  The transect
# (howto_MAAP_AA 3b) wants rerunning before the queue request is written;
# the submitter and collector for it are ready (O7).
# RERUN 2026-09-11 on ab84687, all 17 successful -- howto_MAAP_AA 3b.


# ===========================================================================
# O9. [OK, 2026-09-11]  Point everything that described the old path here.
# ===========================================================================
# howto_MAAP_staging.sh S5/S5b, howto_MAAP_AA 3b-i, Transition_to_maap.md,
# and the four region howtos' step 0.  Replace, do not duplicate, so two
# copies of the procedure cannot drift.
# DONE.  The rule applied: a PROCEDURE is replaced by the OGC one or a pointer
# here; a dated RECORD stays as written, marked as the legacy path where it
# could be mistaken for current.
#   - staging S5: the two commands, the rebuild rule restated for OGC (and
#     that register_algorithm.py now enforces the push checks), and the
#     legacy registrations as a dated history.  S5b: the OGC check_build_id,
#     its verdicts, and where the report is (_stderr.txt, build_id.txt).
#   - staging S6: the queue is submit_job's argument; sizing waits on the
#     post-crossover transect.  S7 kept as the legacy record, with a banner
#     giving the OGC submit call (checked up to the submit: finds 64).
#   - step 0 of GL, AA and the Arctic: register, then check_build_id, and a
#     pointer to S5/S5b.  (Three region howtos have a step 0; the fourth
#     howto, staging, is S5.)
#   - AA: the queue_name input, the 3b-i "BLOCKED" banner, getJobMetrics.
#     GL step 5 points at the OGC calls already written (O7).
#   - Transition_to_maap.md: the smoke test, crossover, build-id and OGC items
#     brought up to date; the fromfile and build-command items and the
#     Registration log marked legacy; the rebuild note restated.
#   - and the O11 leftovers the sweep found: check_build_id's MISMATCH advice
#     and algorithm_config.yml's comment no longer ask for a per-build tag.


# ===========================================================================
# O10. [BEN]  Retire what is left of the legacy path.
# ===========================================================================
#   - delete the on_s3_v2 branch on GitHub (still at e231ba4, unused).
#   - the legacy registrations ATL1415_tile_solve:on_s3 and :on_s3_v2 still
#     exist on /api/mas/algorithm (F2).  Keep them until QC is answered: they
#     ARE the fallback.


# ===========================================================================
# O11. [RESOLVED 2026-09-11 (Ben)]  Stale images: record the build, audit later.
# ===========================================================================
# Ben, 2026-09-11: once the ADE moved to maap-py 5 (the OGC path), workers
# stopped running stale images, and "I would consider the problem resolved if
# each job run by each worker records the git tag for its build, so that
# problems can be audited after the fact."  So no prevention: no per-build
# image tag, no request to MAAP for --force-docker-pull.  (Both were drafted
# here 2026-09-10; git has them, at c645cc3.)
#
# THE RECORD, every item in the image since ab84687:
#   - every job's log starts with run.sh's BUILD_ID line -- commit, build
#     time, algorithm_version, maap_pgt;
#   - every tile's /meta carries build_commit / build_version /
#     build_completed, and errors_build_* for the error step (O12a) -- the
#     copy that outlives the logs;
#   - collect_AA_queue.py prints each tile's commit and lists the builds a
#     ledger's tiles ran.  It does not warn on a mix: a run that reruns
#     patched tiles mixes builds on purpose, and O12b explains the mix.
# A commit, not a tag, names the build: on_s3 is a branch every build reuses,
# and the build time separates two builds of one commit.
# CHECKED ON REAL TILE JOBS 2026-09-11: all 17 of howto_MAAP_AA 3b's jobs
# log BUILD_ID commit=ab84687, and all 16 tiles written carry it in /meta.


# ===========================================================================
# O12. [a: OK ON DPS 2026-09-11 -- all 16 transect tiles (howto_MAAP_AA 3b)
#      carry build_commit and errors_build_commit ab84687 in /meta;
#      b: A SUGGESTION, NO SOFTWARE -- revised 2026-09-11]  Tiles record
#      their build; a run ends with an annotated build history.
# ===========================================================================
# THE SCENARIO, Ben's: a full Antarctic/Greenland/Arctic RUN is under way, a
# bug turns up that affects a few tiles, the bug is patched, those tiles are
# rerun -- and the thousands it did not affect are kept.  From then on the
# run's tiles come from two builds ON PURPOSE, and someone auditing the run
# later needs to know why.
#
# O12a. [OK LOCALLY; IN THE IMAGE SINCE ab84687]  THE TILE KNOWS ITS BUILD.
#   run.sh exports the stamp's commit, algorithm_version and build_completed
#   as ATL1415_BUILD_COMMIT / _VERSION / _COMPLETED; save_fit_to_file writes
#   them as /meta attributes build_commit, build_version, build_completed
#   (ascii, like the existing input_files).  save_errors_to_file, which
#   appends the error fields into the same file, writes errors_build_* --
#   so a patch that reruns only the error step shows up too.  Unset
#   variables (discover, a local run) mean absent attributes, nothing else.
#   WHY IN THE FILE: job logs are not forever, and the mosaic step reads
#   tiles, not logs.
#   DONE: ATL11_to_ATL15.write_build_provenance(), called from both writers;
#   tests/test_build_provenance.py (4 tests: all fields, the errors_ prefix,
#   unset -> absent, empty -> absent).  With the solver stubbed, the three
#   variables reach it and its argv is unchanged; `conda run` passes them
#   through (checked).  NOTE: h5py 3.16 returns these ascii attributes as
#   str, as it already does input_files.
#
# O12b. [SUGGESTION, NO SOFTWARE]  THE RUN'S LAST STEP: AN ANNOTATED BUILD
#   HISTORY.  Ben, 2026-09-11: "generate a git history for the branch
#   corresponding to the run and annotate it to explain what changes happened
#   during the run.  No software is needed to do this."  Exacting by hand and
#   quick for Claude, so ask Claude to draft it and review the draft:
#     1. The builds the run used: the collector's "Builds that ran these
#        tiles" for each of the run's ledgers, or build_commit /
#        errors_build_commit in the tiles' /meta.
#     2. The branch's history across them, oldest build to newest:
git -C ~/git_repos/ATL1415 log --reverse --date=short \
    --format='%h %ad %s' <oldest build>^..<newest build>
#     3. Annotate it: for each commit, what changed and whether it touches
#        the solve; for each build, which tiles ran on it and why they were
#        rerun.
#     4. Keep it with the run -- e.g. build_history.txt beside the run's args
#        files on the bucket, so writing it never blocks register_algorithm.py.
#   It replaces run_notes.txt and the collector code that read it (built
#   2026-09-10 in 78ed7c3, removed 2026-09-11), which asked for a line per
#   build change DURING the run.
#
# LATER, not now: the mosaic step reading /meta build_*, so a released mosaic
# can list the builds inside it beside the annotated history.
