# howto_MAAP_staging.sh -- stand up a fresh MAAP account for ATL1415
#
# ############################################################################
# ##  TENTATIVE.  Written 2026-09-05 BEFORE any of it has been run end to   ##
# ##  end.  It is the plan, not a record of a successful run: expect steps  ##
# ##  to move, split and change as testing advances.  Revise this file as   ##
# ##  that happens -- it is meant to be edited, not preserved.              ##
# ############################################################################
#
# Every step carries a status tag:
#   [DONE]  already carried out for the ben_smith account; here so a fresh
#           account can repeat it, and so the current state is recorded.
#   [OK]    expected to work as written; the pieces exist.
#   [UNTESTED] the pieces exist but this exact command has not been run.
#   [NEEDS CODE: x] blocked -- x does not exist yet.  See docs/Transition_to_maap.md.
#
# Steps are numbered S1..S7 so they can be referred to from the region howtos
# and from Transition_to_maap.md.  The region howtos assume S1-S6 are done.
#
# What does NOT need staging, and never will:
#   - tide models.  pyTMD publishes zarr stores at s3://pytmd, read by range
#     request; MAAP_dps.txt sets --tide_directory=s3://pytmd.  A 60 km tile
#     reads 1.6 MiB (AA) or 6.3 MiB (GL).  Read ANONYMOUSLY -- that bucket
#     rejects signed requests from other accounts.
#   - ATL11 granules.  Read from NASA Earthdata Cloud via earthaccess.
#   - ATL14/ATL15 previous products.  The published v005 granules ARE the
#     previous product; no copy to our bucket is needed, ever (Q27 F1).
#   - crossover tiling schemas.  Built in memory per Q24; no file, no setup step.


# ===========================================================================
# S1. [DONE 2026-09-06]  Build the ATL1415 conda env in the ADE.   (Q5)
# ===========================================================================
# The ADE notebook env cannot `import pointCollection` or `LSsurf`, so every
# ADE-side stage below and in the region howtos is blocked until this exists.
# build-env.sh is the same script the DPS build runs, so the two environments
# stay in step.  It reads `name:` from environment.yml (currently ATL14).
bash build-env.sh
conda activate ATL14
python -c "import pointCollection, LSsurf, sparseqr; print('env ok')"
#
# RUN 2026-09-06: succeeded, /srv/conda/envs/ATL14, python 3.13.  Every import
# the DPS build checks resolves -- numpy scipy h5py osgeo.gdal pyproj sparseqr
# pointCollection LSsurf ATL1415 earthaccess s3fs -- and sparseqr and LSsurf
# both compiled against the conda toolchain here, as they do on DPS.
# This also unblocks the Q4 measurement, which was waiting on nothing else.
#
# ONE MORE PER-ACCOUNT THING, found the same day: setup_ATL1415_region.py makes
# rel<Release>/<hemi>/<region>/ with os.mkdir, one level at a time, but it does
# NOT make --ATL14_root itself.  On a fresh account that root is absent and the
# first setup call dies with FileNotFoundError on .../ATL14_processing/rel006.
# The guard is deliberate enough to keep -- a typo'd --ATL14_root should fail
# loudly rather than silently build a junk tree -- so make the root by hand:
mkdir -p /home/jovyan/ATL14_processing      # = --ATL14_root in MAAP_dps.txt
#
# AND ONE TRAP, hit on 2026-09-06: build-env.sh runs `pip install .`, which
# COPIES the code into the env.  The console scripts (setup_ATL1415_region.py,
# make_ATL1415_queue.py, ATL11_to_ATL15.py) then run that copy, while
# `python -c "import ATL1415"` from the repo directory runs the repo -- so an
# edit can appear to work and be silently absent from the composed args file.
# For ADE work, reinstall editable after S1, and re-run it after pulling:
conda run -n ATL14 python -m pip install -e . --no-deps --no-build-isolation


# ===========================================================================
# S2. [DONE]  Stage the ice/tide masks from Zenodo to the bucket.
# ===========================================================================
# Canonical source is Zenodo record 22259649 (doi:10.5281/zenodo.22259649),
# "Ice and tide masks for ICESat-2 ATL14/15 data products".  A new version
# appears 3-4x/year with each release, so RE-QUERY THE RECORD rather than
# assuming v4.1 is current.  Do NOT use the repo's git-lfs copies of masks/ --
# the pointers are stale and were never pulled.
#
#   destination: s3://maap-ops-workspace/ben_smith/ATL1415/masks/
#   consumed as: --mask_dir in default_args/MAAP_dps.txt, joined onto the
#                relative --mask_file / --tide_mask_file / --geoid_file /
#                --tide_adjustment_file names in AA_0331.txt / GL_0331.txt.
#
# CAVEAT: not every mask is on Zenodo.  The Antarctic shelf-only variants and
# several Arctic masks are lfs-only; if one is missing from the record, say so
# rather than silently falling back to lfs.
#
# NOTE the bucket is a mountpoint-s3 FUSE mount at ~/my-private-bucket: mkdir,
# cp and rm work, but `mv` fails with "Function not implemented".  Reorganize
# server-side with `aws s3 mv --recursive --dryrun` first.


# ===========================================================================
# S3. [DONE]  Stage the ATL11 per-granule geoIndex, built on discover.  (Q13, Q20)
# ===========================================================================
# The index is built on discover (howto_ATL11.sh) and copied to the bucket.
# THIS IS THE ONLY PER-RELEASE HANDOFF LEFT between discover and MAAP.
#
# Layout MUST be the flat one that pointCollection's
# query_ATL11_cloud.index_path_for_granule() builds:
#   s3://maap-ops-workspace/ben_smith/ATL11_index/ATL11_index_<cycles>_<rel>_<ver>/
# for 0331/007/04 that is ATL11_index_0331_007_04/ -- 8100 files, both
# hemispheres in one flat directory.
#
# The flat tree loses the hemisphere split.  The way back is the manifest at
#   s3://maap-ops-workspace/ben_smith/ATL11_index/hemisphere_manifest/{north,south}.txt
#
# CONSEQUENCE, worth knowing before you plan a local run: there is no
# GeoIndex.h5 for 0331_007_04 any more.  THIS RELEASE IS CLOUD-MODE ONLY.
# Cloud runs pass --ATL11_index=<.../ATL11_index> -- the ROOT, not the subdir.


# ===========================================================================
# S4. [UNTESTED]  Verify from the ADE that each staged input is readable.
# ===========================================================================
# Do this BEFORE registering an algorithm.  A DPS worker reads with its own
# AWS credentials, so an ADE read is necessary but not sufficient -- S7 is what
# actually proves the worker can see these.
conda activate ATL14
python - <<'EOF'
import pointCollection as pc
root = 's3://maap-ops-workspace/ben_smith'
# a mask, through GDAL's /vsis3/
print(pc.grid.data().from_geotif(f'{root}/ATL1415/masks/GreenlandIceMask_100m_v4.1.tif',
                                 bounds=[[-2e5, -1e5], [-2.3e6, -2.2e6]]).shape)
# the ATL11 index root
fs = pc.io_utils.get_s3fs(daac=None)
print(len(fs.ls(f'{root}/ATL11_index/ATL11_index_0331_007_04/')), 'index files')
# a tide store, ANONYMOUSLY -- signed requests are rejected
print(pc.io_utils.get_s3fs(daac=None, anon=True).ls('s3://pytmd')[:3])
EOF
# The arctic .db masks are the one that has never been exercised (Q12): they
# are read by ATL1415/make_mask_from_vector.py through ogr, which as of
# 2026-09-05 routes through pc.io_utils.as_gdal_path() so a raw s3:// URI
# becomes /vsis3/.  Confirm one opens before trusting an arctic run.


# ===========================================================================
# S5. [DONE]  Register the DPS algorithm -- and RE-register after every code change.
# ===========================================================================
# Use the script, which does the four REPL lines plus the push checks below and
# prints the build URL:
#
#   /srv/conda/envs/notebook/bin/python register_algorithm.py
#   /srv/conda/envs/notebook/bin/python register_algorithm.py --dry-run   # checks only
#
# The notebook env is the one that always has maap-py (/srv/conda/envs/notebook,
# maap-py 4.2.0), and it is the default interpreter in a fresh ADE terminal.  An
# ATL14 env built BEFORE 2026-09-08 does not have it, so `conda activate ATL14`
# -- which every region howto starts with -- shadowed it and broke the script;
# environment.yml now installs maap-py==4.2.0, so an ATL14 env rebuilt since
# then works too.  The path above is the one that works either way.  The script
# says so and exits 1 rather than tracebacking.  A notebook is NOT required:
# MAAP_API_HOST and MAAP_PGT are set in the ADE environment, so MAAP()
# authenticates from env with no ~/.maap-py.ini, and a plain python session
# works.  What the script runs is still just:
#
#   from maap.maap import MAAP
#   maap = MAAP(maap_host='api.maap-project.org')
#   response = maap.register_algorithm_from_yaml_file('algorithm_config.yml')
#   response.json()['message']['job_web_url']     # <-- where the build URL is
#
# That last line is the one worth having written down: register_algorithm_from_
# yaml_file() returns a raw requests.Response, and job_web_url is nowhere in
# maap-py -- it is in the server's JSON, one level down under 'message'.
#
# 2026-09-04: ATL1415_tile_solve:on_s3 registered, HTTP 200, build pipeline
# 20059 / job 21091, and describeAlgorithm() answers 200.  Re-register after
# any change to algorithm_config.yml, build-env.sh, run.sh or environment.yml.
#
# Registration makes the signature (x0, y0, step, args_file, queue_name) --
# positionals first, then the file input, then the queue override.
#
# Build logs are BROWSER-ONLY; they are not on this filesystem.  Read them at
# the URL register_algorithm_from_yaml_file() returns.
#
# THE REBUILD RULE, which every region howto's step 0 points here for:
# DPS does not run the ADE working copy.  It clones repository_url at
# algorithm_version (on_s3) FROM GITHUB at build time and bakes the result into
# a container.  Anything uncommitted, unpushed, or committed since the last
# build is NOT on the worker, and no job log will say so -- a stale image
# surfaces as whatever the missing fix was meant to prevent.  So before any
# submission, check both:
#
#   git -C ~/git_repos/ATL1415 status --short                     # nothing uncommitted
#   git -C ~/git_repos/ATL1415 log --oneline origin/on_s3..on_s3  # empty
#
# and if either is non-empty, push and re-register before submitting.
#
# 2026-09-08: THIS HAS ALREADY BITTEN ONCE, before a single job was submitted.
# The 2026-09-04 build predated three commits of fixes that every region needs:
#   0d77630  as_gdal_path, without which the arctic .db masks do not open
#   1023306  pyTMD AWS_NO_SIGN_REQUEST, without which EVERY /vsis3 read on the
#            worker -- mask, geoid, tide mask, in every region -- returns 403
#   b293807  Q27 W1/W3/W4/W5, the previous-product fixes
# on_s3 was pushed to b293807 and re-registration submitted 2026-09-08 by Ben;
# THE BUILD RESULT HAS NOT BEEN SEEN YET, so S7 is still gated on it going
# green.  Had the smoke test gone against the old image it would have failed on
# the Iceland mask read and looked exactly like a credentials problem.


# ===========================================================================
# S6. [UNTESTED]  Ask the MAAP platform team for an organizational DPS queue.
# ===========================================================================
# NOT OPTIONAL FOR PRODUCTION, and it has days of latency -- ask early.
# Public queues are throttled to ~10 jobs/hr, which makes a per-tile fan-out of
# thousands of jobs infeasible.  Queues visible to this account (getQueues(),
# 2026-09-04): maap-dps-sandbox, maap-dps-worker-8gb, -16gb, -32gb, -64gb,
# maap-dps-worker-32vcpu-64gb.
#
# STILL OPEN AS OF 2026-09-08, and the account change did NOT close it.  The
# ben_smith account is now status 'active' and a member of organization
# icesat-2 (id 32) -- see S7 question 1 -- but getQueues() returns THE SAME SIX
# QUEUES, with no icesat-2 queue among them.  Org membership does not carry an
# organizational queue with it; the queue is a separate request to the platform
# team, and there is now an org to attach one to.  Make that request in
# parallel with S7 rather than after it: S7 measures what to ask FOR (cores,
# RAM, walltime per tile), but nothing about S7 produces the queue itself.
#
# Ask for: cores / RAM / disk / walltime per queue, and any max-in-flight limit.
# algorithm_config.yml names -32gb as its default; -32vcpu-64gb is probably the
# better production target, but NOTHING HAS BEEN MEASURED -- S7 is what measures it.


# ===========================================================================
# S7. [UNTESTED]  Smoke-test ONE sandbox job.   <-- THE GATE ON EVERYTHING ELSE
# ===========================================================================
# Steps 5-7 of every region howto are unwritable in final form until this has
# been done once.
#
# GATED ON S5's REBUILD.  Submit only against an image built from b293807 or
# later; against the 2026-09-04 image this test measures a container that has
# neither the pyTMD nor the previous-product fixes, and every answer below
# would be wrong.  Confirm the build is green first.
#
# THE SMOKE-TEST TILE, chosen by Ben 2026-09-05: ICELAND, x0,y0 = 1260, -2620 km.
# 43K points -- a medium-sized dataset, likely on the edge of the ice sheet, so
# it is a real solve rather than a trivial or a saturated one.  IS is also the
# test region for the whole workflow (Q21), so this tile is reused as the mask
# check in howto_MAAP_arctic.sh step 2.  It is on the grid: centers sit at odd
# multiples of tile_spacing/2, and 1260 = 63 x 20 km, -2620 = -131 x 20 km.
#
# Compose and publish the IS args file first -- arctic steps 3 and 4.
#
#   job = maap.submitJob(identifier='ATL1415_smoke',
#                        algo_id='ATL1415_tile_solve', version='on_s3',
#                        queue='maap-dps-sandbox',
#                        x0=1260000, y0=-2620000, step='prelim',
#                        args_file='s3://maap-ops-workspace/ben_smith/ATL1415/run_args/'
#                                  'rel006/north/IS/input_args_IS.txt',
#                        queue_name='maap-dps-sandbox')
#
# THE ARGUMENT IS `queue`, NOT `queue_name`.  Corrected 2026-09-08 by reading
# maap-py 4.2.0: submitJob(identifier, algo_id, version, queue, ...) takes queue
# as a REQUIRED parameter, and everything else in **kwargs is forwarded as a WPS
# algorithm input (DpsHelper._skit, DpsHelper.py:29-46).  The earlier form here
# passed only queue_name, so it would have raised
#   TypeError: submitJob() missing 1 required positional argument: 'queue'
# before sending anything.  Both are passed above, set to the SAME queue, so it
# cannot matter which one the server honours -- registration turns the yaml's
# `queue:` into a queue_name INPUT, while `queue` is maap-py's own job-queue
# field, and which takes precedence is not documented.  A later run that wants a
# per-tile override should set both together, or settle the precedence first.
#
# DO NOT pass username: submitJob overwrites it from profile.account_info()
# (maap.py:345-346), so it is 'ben_smith' whatever is passed.
#
# ALREADY ANSWERED by the successful build -- do not re-ask:
#   - the build container reaches github.com (the git+ deps installed)
#   - sparseqr compiles against the conda suitesparse
#   - input order is positionals then file
#
# WHAT THIS RUN IS FOR -- only a real job settles these:
#   1. does submitJob work at all?  ASKED ORIGINALLY AS "on an account whose
#      status is 'inactive' with no organization" -- THAT PREMISE IS GONE.
#      Checked 2026-09-08: profile.account_info() reports status 'active' and
#      organizations [{'id': 32, 'name': 'icesat-2'}], where it reported
#      inactive and none when this file was written on 2026-09-05.  So the
#      specific fear -- that the account itself would refuse a job -- has no
#      basis now, and it never explained much anyway: listAlgorithms, listJobs
#      and register all returned 200 while it was still inactive.  Submit is
#      nonetheless UNTESTED, so it stays on this list; it is just no longer the
#      question with a named reason to fail.  (account_info() also confirms
#      username 'ben_smith', which is what the submitJob call below passes.)
#      ANSWERED 2026-09-08: YES.  Job 61589c00-7a67-4881-93af-b96d0f0e8c4b was
#      accepted (HTTP 200), queued, and ran on maap-dps-sandbox.  Submission is
#      not a problem on this account.
#   2. is ~/.netrc really bind-mounted into the worker?  earthaccess auth, and
#      therefore every ATL11 read, depends on it.
#      ANSWERED 2026-09-08: NO, AND IT IS THE BLOCKER.  Job
#      1f356303-426f-42aa-9c12-de0f35aeb9ae got all the way into the solve and
#      died on the first ATL11 read:
#        pointCollection/scripts/query_ATL11_cloud.py:52
#            earthaccess.login(strategy='netrc')
#        LoginStrategyUnavailable: ('No .netrc found at /root/.netrc',)
#      The worker runs as uid 0 with HOME=/root, and _docker_params.json lists
#      every mount it gets -- /data/work/{jobs,tasks,workers,cache}, /tmp, the
#      docker socket, and ONE credential-shaped file:
#        /data/work/etc/maap-dps.env -> /home/ops/.maap-dps.env:ro
#      No .netrc anywhere.  NOTE the strategy is HARD-CODED to 'netrc' in
#      pointCollection, so EARTHDATA_USERNAME/PASSWORD in the environment would
#      not help either as the code stands -- earthaccess.login()'s default
#      strategy tries environment first, that call refuses to.
#      RESOLVED 2026-09-08, and NO .netrc IS NEEDED ANYWHERE.  A probe job
#      (b2f411f4-70e1-4c9f-8b2a-848b6810c16a) reported what the worker has:
#        /root/.netrc                 absent
#        /home/ops/.maap-dps.env      present; keys DPS_MACHINE_TOKEN, MAAP_API_HOST
#        set in the container         MAAP_API_HOST, MAAP_PGT
#      MAAP_PGT is the one that matters: it is MAAP's own auth, and
#      maap.aws.earthdata_s3_credentials('https://data.nsidc.earthdatacloud.
#      nasa.gov/s3credentials') exchanges it for NSIDC's short-lived read
#      credentials.  The probe made that call ON THE WORKER and got back
#      accessKeyId/secretAccessKey/sessionToken -- so the DAAC reads need no
#      Earthdata credential of ours at all, and nothing of ours is at rest in
#      the bucket or the DPS cache.  (MAAP docs: science/NISAR/NISAR_access.html
#      cell [2].)  The probe has been removed from run.sh now that it has
#      answered; a broker failure warns from pointCollection instead.
#      THE FIX is in pointCollection, branch maap_s3_credentials, pinned from
#      pyproject.toml until it is merged to main:
#        - get_s3fs(daac=) tries MAAP's broker, falls back to earthaccess
#        - find_ATL11_granules drops login(strategy='netrc'); a CMR search
#          needs no auth, and that call was failing before reaching CMR
#      Verified from the ADE before building: a real read of
#      ATL11_000103_0331_007_04.h5 out of nsidc-cumulus-prod-protected.
#   3. will a `file` input accept an s3://maap-ops-workspace/... URL?
#      ANSWERED 2026-09-08: YES, and note HOW.  DPS does not copy the file into
#      input/ -- it SYMLINKS it into a shared cache:
#        input_args_IS.txt -> /data/work/cache/5/e/4/b/<md5>/input_args_IS.txt
#      That difference is what failed the first smoke test: run.sh selected the
#      args file with `find -type f`, which tests the link (-type l) and matched
#      nothing, so the job exited 2 in its own guard while the ls in that guard
#      printed the symlink it had just refused to find.  Fixed by `find -L`.
#      A LOCALIZED INPUT IS A SYMLINK -- assume it everywhere, including the
#      input/prelim/ tree the --matched step will localize.
#   4. do the s3:// reads in the composed args file work on the worker's own
#      AWS credentials (masks via /vsis3/, index via s3fs, tides anonymously)?
#      ANSWERED 2026-09-08: YES.  Job 08ee68ea got through the Iceland .db
#      mask, the geoid and the ATL11 index on s3://maap-ops-workspace using
#      the worker's own identity, and through the ATL11 granules themselves on
#      MAAP-brokered NSIDC credentials, into three QR iterations of the solve.
#      A 403 on any of them would have stopped it in seconds.  All three
#      credential paths in one job are now exercised.
#      STILL OPEN, but no longer preempted: read_ATL11 is called from
#      ATL11_to_ATL15.py:635, before any mask or geoid read, so until Q2 was
#      fixed the netrc failure hid this question entirely.  The first run on the
#      maap_s3_credentials build is the one that finally tests it.  NOTE the two
#      are different credentials and must not be confused: the DAAC reads use
#      MAAP-brokered NSIDC credentials, while the masks, geoid and ATL11 index
#      on s3://maap-ops-workspace use the WORKER'S OWN AWS identity via
#      get_s3fs(daac=None), and the tide stores are read anonymously.
#   5. how long does one prelim tile take, and how much memory does it need?
#      -> this is what sizes the production queue in S6.
#      PARTIAL, 2026-09-08, job 08ee68ea on maap-dps-sandbox: the Iceland tile
#      read 239613 ATL11 points, decimated to 239613, put 43773 into the fit,
#      and ran QR iterations at ~145 s each on 2 threads (--max_iterations=3).
#      It was killed at 600 s -- see the sandbox limit below -- so the total is
#      still unmeasured.  Use getJobMetrics(), which returns max_mem_usage and
#      job_duration_seconds directly, rather than timing it by hand.
#
# WHAT THE SECOND RUN (1f356303, 2026-09-08) ESTABLISHED BESIDES Q2.  The -L fix
# works: run.sh found the localized args file, printed its banner, and handed
# ATL11_to_ATL15 a valid @argsfile from the cache path.  Every one of the ~35
# composed arguments is echoed in _stdout.txt, which makes that log the easiest
# check that a composed args file is what you think it is.
#
# THE SANDBOX WORKER HAS 2 CORES (run.sh's `threads : 2`, from nproc).  That is
# a data point for Q5 and S6: whatever a tile costs, the sandbox is not where
# its production timing gets measured.
#
# `queue` IS THE ROUTING FIELD, not queue_name.  _context.json's params list
# holds exactly x0, y0, step (positional) and args_file (localize) -- queue_name
# does not appear at all, while `queue` put the job on maap-dps-sandbox.  The
# per-tile queue override that algorithm_config.yml describes therefore needs
# `queue=`; treat queue_name as decorative until something proves otherwise.
#
# 5 IS STILL OPEN.  1, 2, 3 and 4 are answered.
#
# THE SANDBOX QUEUE HAS A 600-SECOND SOFT TIME LIMIT, and nothing in the MAAP
# docs says so -- dps_queues.html says public worker queues have UNLIMITED
# walltime, which is evidently about the worker queues and not about
# maap-dps-sandbox.  Job 08ee68ea died at 600 s with
#   Soft time limit (600s) exceeded for job-ATL1415_tile_solve...
# after ten entirely healthy minutes.  So THE SANDBOX CANNOT COMPLETE A REAL
# TILE: use it to prove a job starts, and move to a worker queue to time one.
#
# AND NOTE THE EXIT CODE: 143.  The docs' one-row error table says 143 means a
# spot interruption, "re-run the job(s) later".  Here it was a walltime kill.
# The two are indistinguishable by exit code alone, which matters for
# check_MAAP_jobs.py: requeueing a job that exceeded its walltime just burns
# the queue again, forever.  Distinguish them on the log -- a walltime kill
# says "Soft time limit (Ns) exceeded" -- or on elapsed time against the
# queue's limit, not on 143 by itself.
#
# WHERE THE LOGS ARE -- and they are NOT browser-only, unlike the BUILD logs.
# A failed job triages itself onto the bucket, readable from the ADE with
# ordinary aws s3:
#
#   aws s3 ls s3://maap-ops-workspace/dataset/triaged_job/v1.4.0/ | tail
#   aws s3 cp s3://maap-ops-workspace/dataset/triaged_job/v1.4.0/<job dir>/_stderr.txt -
#
# The directory is named triaged_job-job-<algo>__<version>-<timestamp>_task-<uuid>
# and holds _stdout.txt, _stderr.txt, _exit_code, _run.sh (the exact command the
# worker ran), _context.json and _docker_params.json (what DPS was told to
# localize).  maap.getJobResult(<job_id>) returns those URLs plus the tail of
# stderr inline, which is usually enough on its own.  Reach for these before the
# browser: _stdout.txt is what showed the symlink in question 3.
