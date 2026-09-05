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
# S1. [OK]  Build the ATL1415 conda env in the ADE.        (Q5)
# ===========================================================================
# The ADE notebook env cannot `import pointCollection` or `LSsurf`, so every
# ADE-side stage below and in the region howtos is blocked until this exists.
# build-env.sh is the same script the DPS build runs, so the two environments
# stay in step.  It reads `name:` from environment.yml (currently ATL14).
bash build-env.sh
conda activate ATL14
python -c "import pointCollection, LSsurf, sparseqr; print('env ok')"


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
# S5. [DONE]  Register the DPS algorithm.
# ===========================================================================
# From a MAAP notebook (the ADE), not the shell:
#
#   from maap.maap import MAAP
#   maap = MAAP(maap_host='api.maap-project.org')
#   maap.register_algorithm_from_yaml_file('algorithm_config.yml')
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


# ===========================================================================
# S6. [UNTESTED]  Ask the MAAP platform team for an organizational DPS queue.
# ===========================================================================
# NOT OPTIONAL FOR PRODUCTION, and it has days of latency -- ask early.
# Public queues are throttled to ~10 jobs/hr, which makes a per-tile fan-out of
# thousands of jobs infeasible.  Queues visible to this account (getQueues(),
# 2026-09-04): maap-dps-sandbox, maap-dps-worker-8gb, -16gb, -32gb, -64gb,
# maap-dps-worker-32vcpu-64gb.
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
# THE SMOKE-TEST TILE, chosen by Ben 2026-09-05: ICELAND, x0,y0 = 1260, -2620 km.
# 43K points -- a medium-sized dataset, likely on the edge of the ice sheet, so
# it is a real solve rather than a trivial or a saturated one.  IS is also the
# test region for the whole workflow (Q21), so this tile is reused as the mask
# check in howto_MAAP_arctic.sh step 2.  It is on the grid: centers sit at odd
# multiples of tile_spacing/2, and 1260 = 63 x 20 km, -2620 = -131 x 20 km.
#
# Compose and publish the IS args file first -- arctic steps 3 and 4.
#
#   job = maap.submitJob(algo_id='ATL1415_tile_solve', version='on_s3',
#                        x0=1260000, y0=-2620000, step='prelim',
#                        args_file='s3://maap-ops-workspace/ben_smith/ATL1415/run_args/'
#                                  'rel006/north/IS/input_args_IS.txt',
#                        queue_name='maap-dps-sandbox',
#                        username='ben_smith', identifier='ATL1415_smoke')
#
# ALREADY ANSWERED by the successful build -- do not re-ask:
#   - the build container reaches github.com (the git+ deps installed)
#   - sparseqr compiles against the conda suitesparse
#   - input order is positionals then file
#
# WHAT THIS RUN IS FOR -- only a real job settles these:
#   1. does submitJob work at all on an account whose status is 'inactive' with
#      no organization?  (listAlgorithms/listJobs/register all return 200, so
#      inactive is NOT the blocker it looked like -- but submit is untested.)
#   2. is ~/.netrc really bind-mounted into the worker?  earthaccess auth, and
#      therefore every ATL11 read, depends on it.
#   3. will a `file` input accept an s3://maap-ops-workspace/... URL?
#   4. do the s3:// reads in the composed args file work on the worker's own
#      AWS credentials (masks via /vsis3/, index via s3fs, tides anonymously)?
#   5. how long does one prelim tile take, and how much memory does it need?
#      -> this is what sizes the production queue in S6.
