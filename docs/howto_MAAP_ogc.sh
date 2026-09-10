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
# Steps are numbered O1..O10 so they can be cited ("OGC step 3").  Status
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
#     Outputs are collected by `glob: ./output*`, which run.sh's output/
#     already satisfies.
#
# F9. JOBS ADDRESS A PROCESS BY ID, NOT BY name:version.  submit_job() POSTs
#     to /api/ogc/processes/<process_id>/execution, and each deployed process
#     has a numeric processID beside its string id.  -- maap.py and
#     list_algorithms().  Which of the two submit_job wants is not yet known.
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
# QB. [TEST, O6] Does a SECOND build of the same algorithm_version produce a
#     fresh image?  That is the bug that started today (8aad07d); nothing
#     says the OGC path is immune to it.  The build stamp answers it with one
#     job, and s:commitHash (F10) answers it with none.
#
# QC. [BEN / MAAP] Is the legacy /api/mas path being retired, and when?  It
#     decides whether the ATL14 fallback (F2) is safe to keep.
#
# QD. [TEST, O6] What do get_job_status / get_job_result return, and is the
#     job's _stdout.txt still on the bucket at a prefix we can derive?  The
#     collector and check_build_id both depend on reading that log.
#
# QE. SETTLED 2026-09-10 (Ben): maap-py 5.1.0 in environment.yml, for the
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
# O2. [OK LOCALLY, 2026-09-10; UNTESTED ON DPS]  The CWL calling convention.   [ADE]
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
# O3. [OK AGAINST MOCKS, 2026-09-10; NEVER POSTED]  POST api/build.   [ADE]
# ===========================================================================
# Keep what the script is for -- REFUSE to register work that is not on
# GitHub -- and swap the transport:
#   - read algorithm_config.yml, send it as JSON to POST /api/build using
#     maap-py's _get_api_header() (F6); print the build_id.
#   - print EVERY URL in each response, as f3d5049 already does.
#   - --status <build_id>: GET /api/build/<id>, follow deploymentLink to
#     /api/ogc/deploymentJobs/<id>, and print status plus all three links
#     (pipeline, deployment, deployment pipeline) and any deploymentError.
#   - once deployed, print the processID (F9) that submitters need.
# Works under maap-py 4.2.0 AND 5.x: it uses only the auth header, not
# either version's algorithm methods.
# DONE.  --dry-run now also prints the exact JSON it would POST; identical
# under both maap-py versions.  Tested against mocked responses: accepted
# (202), rejected (400), non-JSON (502), and --status --wait through build
# running -> successful -> deployment running -> deployed -> processID found,
# plus a failed build and an unfinished one.  Exit codes: 0 deployed or
# accepted, 1 invalid config or unpushed work, 2 rejected or failed, 3 not
# finished.  ONE REAL CALL, read-only: --status on a nonexistent id gets HTTP
# 404 {"title": "No build with that build ID found", ...} -- so the path is
# right, and the service speaks RFC 7807 problem documents.


# ===========================================================================
# O4. [BEN] [UNTESTED]  Register, and watch the build.   [ADE]
# ===========================================================================
/srv/conda/envs/notebook/bin/python register_algorithm.py
/srv/conda/envs/notebook/bin/python register_algorithm.py --status $build_id --wait   # id printed above
# Ben runs this one, to open the links in a browser.  RECORD, here: the
# build_id, the three link keys as they really come back, the processID,
# how long the build took, and the s:commitHash in the generated CWL -- the
# first real answer to QB, before any job is spent.


# ===========================================================================
# O5. [NEEDS CODE: check_build_id.py]  Port to the OGC job calls.   [ADE]
# ===========================================================================
# submit_job(process_id, {x0, y0, step: build_id, args_file}, queue), then
# get_job_status / get_job_result (QD).  TWO comparisons instead of one:
# the stamp in the image against origin/<algorithm_version>, AND against the
# s:commitHash in the deployed CWL (F10).  All three agreeing is MATCH.
# The CWL alone would have said "stale" 2026-09-09 without submitting a job.


# ===========================================================================
# O6. [UNTESTED]  One build_id job.   [ADE -> DPS]
# ===========================================================================
/srv/conda/envs/notebook/bin/python scripts/maap/check_build_id.py
# Answers QB (fresh image?) and QD (how logs come back), and QE's remaining
# risk: the BUILD_ID line must say maap_pgt=set, or no tile on this system
# can read NSIDC.  Nothing below runs until this says MATCH and maap_pgt=set.


# ===========================================================================
# O7. [NEEDS CODE: submit_AA_queue.py, collect_AA_queue.py]   [ADE]
# ===========================================================================
# Same two scripts, new calls: submit_job with inputs as a dict and the
# queue as an argument (F1, F9); the collector reads the log wherever O6
# found it (QD).  The half-routing and the ledger format do not change.


# ===========================================================================
# O8. [UNTESTED]  The two pole-hole tiles.   [ADE -> DPS]
# ===========================================================================
# howto_MAAP_AA step 3b-i, unchanged in intent: E220_N20 and E300_N20, and
# N_XO > 0 on both.  E220_N20 read N_XO=0 before the crossover fix.


# ===========================================================================
# O9. [NEEDS CODE: docs]  Point everything that described the old path here.
# ===========================================================================
# howto_MAAP_staging.sh S5/S5b, howto_MAAP_AA 3b-i, Transition_to_maap.md,
# and the four region howtos' step 0.  Replace, do not duplicate, so two
# copies of the procedure cannot drift.


# ===========================================================================
# O10. [BEN]  Retire what is left of the legacy path.
# ===========================================================================
#   - delete the on_s3_v2 branch on GitHub (still at e231ba4, unused).
#   - the legacy registrations ATL1415_tile_solve:on_s3 and :on_s3_v2 still
#     exist on /api/mas/algorithm (F2).  Keep them until QC is answered: they
#     ARE the fallback.
