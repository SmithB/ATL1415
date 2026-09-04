# Transition to MAAP

## TBDs: 
[ ] Request an organizational DPS queue from the MAAP platform team.
    Public queues are throttled to ~10 jobs/hr, so a per-tile fan-out of thousands of jobs is
    infeasible until this is resolved.  Days of latency -- ask early.
    - queues visible to this account (getQueues(), 2026-09-04): maap-dps-sandbox,
      maap-dps-worker-8gb, -16gb, -32gb, -64gb, maap-dps-worker-32vcpu-64gb.
      algorithm_config.yml currently names -32gb; -32vcpu-64gb is probably the better
      production target, but nothing has been measured.  Ask for cores/RAM/disk/walltime
      per queue and for any max-in-flight-jobs limit.
[ X ] Account status is NOT the blocker it looked like (checked 2026-09-04).
    profile.account_info() still reports status=inactive, organizations=[] (username ben_smith,
    id 1786, created 2026-07-01) -- but that does not stop the API:
    - listAlgorithms() and listJobs() both return HTTP 200
    - register_algorithm_from_yaml_file() returned HTTP 200 and started a real build
    Still UNKNOWN: whether submitJob works on an inactive account with no organization.  That is
    the next thing to find out, and maap-dps-sandbox is where to find it out.
[ ] Smoke-test one sandbox DPS job to settle what the docs do not say.  The three build
    files now exist (build-env.sh, run.sh, algorithm_config.yml -- see Packaging below), so
    this is unblocked apart from the queue question.
    - does the build container reach github.com (for the git+ pip deps)?
    - does sparseqr (PySPQR) actually compile against the conda suitesparse there?
    - is ~/.netrc really bind-mounted into the worker (earthaccess auth depends on it)?
    - will a `file` input accept an s3://maap-ops-workspace/ben_smith/... URL?
    - does DPS pass `file` inputs before or after `positional` inputs?
    - how long does one prelim tile take, and how much memory does it need?

## Registration log
2026-09-04: ATL1415_tile_solve:on_s3 registered from algorithm_config.yml.
  register_algorithm_from_yaml_file() -> HTTP 200, commit 9a9988fd in
  repo.maap-project.org/root/register-job-hysds-v4, build pipeline 20059, job 21091.
  Build log: https://repo.maap-project.org/root/register-job-hysds-v4/-/jobs/21091/raw
  -- that URL needs an OIDC login, so it is readable from a browser but not from a script.
  THAT LOG IS THE ANSWER to two open questions: does the build container reach github.com for
  the git+ pointCollection/LSsurf/PySPQR deps, and does sparseqr compile against the conda
  suitesparse.  build-env.sh fails explicitly at a SuiteSparseQR.hpp check if the headers are
  missing, so look for that message before assuming a generic CFFI failure.
  As of end of day the algorithm had NOT yet appeared in listAlgorithms(); describeAlgorithm()
  returns HTTP 500 ('NoneType' object has no attribute 'get') while a build is still pending,
  which is just "not registered yet" surfacing badly -- not a failure signal in itself.

## Software: 
[ X ] Check build environment for suiteSparse
    - NO MAAP base image ships it.  maap_base is debian:11 + git + Miniforge and nothing else;
      the pangeo ADE image lacks it too.  The DPS build must therefore be conda-based.
[ X ] Commit changes that let pointCollection access ATL11 on earthData
    - merged to pointCollection main (7f1d945); geoindex_bugfix is now fully contained in
      main, and the temporary pin has been dropped from pyproject.toml
[ X ] Update tiling_schema to match true ATL11 configuration
    - ATL11xo bins use _round_ and not _floor_: `ATL1415/scripts/setup_ATL11_xover.py`
      (2ae71a4, 2026-07-03) writes the 200 km schema with mapping_function_name='round',
      which is also pointCollection.tilingSchema's default.  Nothing further to change.
[ ] setup_ATL11_xover must also be able to emit a CLOUD schema.  It currently writes
    directory=<local cycle_src> and no 'source', so resolve_files_for_box() can only find
    local tiles.  Cloud xovers need source={'type':'EarthAccess','short_name':'ATL11XO'}
    and directory=None; the READ side already handles that (read_ATL11_xovers, 299ab26,
    on pointCollection 0e94ef1).
[ X ] Add options for ATL11_to_ATL15 to query earthData for ATL11 granules
    - `--ATL11_earthaccess` (295c181).  Granules for a tile are found by searching
      earthaccess/CMR with the tile's lon/lat bounding box (_lonlat_bounding_box, with
      explicit antimeridian handling -- a near-pole tile straddling +-180 deg otherwise
      produces a bogus near-global box).  Covered by tests/test_read_ATL11.py.
[ X ] Add options to ATL11_to_ATL15 to pass the earthData s3 location to pointCollection
    - there is no separate option, because the S3 location is not something the caller
      supplies: it comes back from the CMR search (granule.data_links(access='direct'))
      and is handed straight to pointCollection's read_ATL11_granule_cloud_items, along
      with one shared s3fs handle reused across granules.  In this mode --ATL11_index is
      reinterpreted as the per-granule geoIndex ROOT.  Crossovers get the equivalent
      treatment through the tiling schema's remote source (299ab26).
    - CAVEAT for DPS: --ATL11_index is os.path.abspath'ed, and the xover schema JSON is
      read with a plain open() inside tilingSchema.from_file, so BOTH still have to name
      a local path.  On a worker with no workspace mount they must either be localized as
      DPS `file` inputs or learn to accept s3:// URIs -- see the MAAP_dps.txt item.
[ ] Create a processing string for ATL11_to_ATL15 to read previous versions of ATL14/15 from the cloud

## Packaging (for the DPS build)
[ X ] pyproject.toml: request the cloud extra -- "pointCollection[cloud] @ git+..."
    - without it earthaccess/s3fs/fsspec are never installed and the cloud read path fails at import
    - the extra is purely additive (cloud = [earthaccess, s3fs, fsspec], no version constraints),
      so it does not perturb the pointCollection that LSsurf resolves to -- one install satisfies both
    - LSsurf/requirements.txt names the same direct URL without the extra; same URL, so pip
      unions the extras rather than reporting a conflict
[ X ] pyproject.toml: SMBcorr is no longer a hard dependency -- moved to a [firn] extra, and
    repointed from the smithb fork to upstream tsutterley/SMBcorr@282e6be (no tags upstream,
    so a sha; a DPS build clones fresh every time).
    - the smithb fork had DIVERGED from upstream: 40 commits behind, 27 ahead, last touched
      2024-02-29 vs upstream 2024-09-19
    - the fork also still used `np.NaN`, removed in numpy 2.0 (this env has numpy 2.4.6), so
      its MERRA2_hybrid path would have raised AttributeError at runtime
    - `ATL1415/assign_firn_variable.py` was a whitespace-only-different copy of SMBcorr's own
      module, carrying the same np.NaN bug; it is now a thin lazy wrapper that delegates, so
      the two cannot drift
    - upstream's `SMBcorr/__init__.py` does NOT re-export assign_firn_variable, so both lazy
      imports name the submodule (`from SMBcorr.assign_firn_variable import ...`); the
      package-level form binds the MODULE and the call raises TypeError
    - committed as f7dd50b
[ ] Confirm the firn correction actually runs against tsutterley/SMBcorr.  Deferred until the
    firn option is next exercised -- the 27 fork-only commits were NOT audited, so if anything
    ATL1415 needs lived only in the fork, it will surface here.  Only `default_args/firn_update.txt`
    uses --firn_model; no production rel_006 string does, and the per-tile DPS solve does not.
    - SMBcorr is not installed in the ADE env, so NOTHING on this path has been executed.  The
      import form above was verified by reading upstream at the pinned sha, not by running it.
    - this branch is deliberately PARKED: do not spend session time auditing it unless the firn
      option is actually being brought back into use.
[ X ] environment.yml: NO scikit-sparse needed -- that earlier note was wrong.  Nothing in
    ATL1415, LSsurf or pointCollection imports `sksparse`; the solver dependency is `sparseqr`
    (PySPQR), which arrives through pip, not conda.  environment.yml is unchanged.
[ ] LSsurf: no change needed.  Its requirements.txt does declare
    `sparseqr @ git+https://github.com/SmithB/PySPQR.git`, so pip will pull PySPQR transitively --
    the open question is only whether it COMPILES in the build container.  It builds CFFI bindings
    against SuiteSparseQR and needs the dev headers, so conda `suitesparse` must be installed
    BEFORE pip runs.  Its bare `gdal` requirement resolves safely only if conda gdal-config is
    already on PATH.  Verify in the smoke-test build; this is the likeliest build failure.

## Algorithm structure
Decided: ONE registered algorithm, not six.  A MAAP algorithm has exactly one run_command, and
the conda+SuiteSparse build is the expensive part -- six registrations means six builds of an
identical environment.  Dispatch the stage inside run.sh instead.

Only the per-tile solve belongs in DPS.  Setup, queue-build, mosaic, to-netcdf and browse stay
in the ADE: they are cheap and serial, and they need the my-private-bucket mount that DPS
workers do not have.  This mirrors the current SLURM split (setup on the login node, queue
lines to sbatch).

There is no chaining/DAG mechanism -- jobs are independent and couple only through S3 paths.
[ ] Have ATL11_to_ATL15 write tiles to a deterministic S3 prefix rather than relying on the
    timestamped dps_output/<algo>/<tag>/YYYY/MM/DD/HH/MM/SS/<us>/ scatter, which the mosaic
    step cannot address.  Keep output/ for logs only.

## MAAP environment
[ ] Identify files that need to be downloaded locally
    [ X ] ATL11 index files 
        - restaged 2026-09-03 to the flat layout index_path_for_granule() builds:
          s3://maap-ops-workspace/ben_smith/ATL11_index/ATL11_index_0331_007_04/ (8100 files)
          old ATL11_007_cycle_03_31_v04/ tree (per-hemisphere GeoIndex.h5 + .tgz) deleted
          hemisphere membership saved under ATL11_index/hemisphere_manifest/
          NOTE: no local-mode GeoIndex.h5 remains for this release -- cloud mode only
    [ X ] pyTMD files
        - CATS2008 staged to s3://maap-ops-workspace/ben_smith/tide_models/CATS2008/
          (grid_CATS2008, hf.CATS2008.out); default_args/MAAP.txt points --tide_directory
          at the ADE mount for that prefix
    [ ] Gr1km-v2 tide model -- NOT staged, and default_args/GL_0331.txt asks for it
        (--tide_model=Gr1km-v2), so a Greenland run has nothing to read.  The only other
        object under tide_models/ is a stray .ipynb_checkpoints/deltat-checkpoint.iers;
        the real deltat.iers is missing too.
    [ X ] Greenland tide mask
        - masks/Arctic/BedMachineGreenland-v6_shelf_edited.tif, 1.6 MiB, staged 2026-09-03;
          that is what GL_0331.txt names as --tide_mask_file
    [ X ] Download Antarctic and Greenland masks from Zenodo
        - SOURCE: Zenodo record 22259649, doi:10.5281/zenodo.22259649, "Ice and tide masks for
          ICESat-2 ATL14/15 data products" (Smith & Sutterley, CC-BY-4.0, published 2026-09-02).
          This is now the canonical distribution point for these masks -- prefer it over the
          repo's lfs copies, and expect a new version 3-4 times a year with each ATL14/15 release.
        - all 5 files downloaded and md5-verified against the Zenodo manifest, then staged to
          s3://maap-ops-workspace/ben_smith/ATL1415/masks/ (2026-09-03) and re-verified by
          md5 after the copy:
            Antarctic/AntarcticIceMask_2018.00_2026.25_240m_v4.1.tif        23.1 MiB
            Antarctic/AntarcticIceMask_..._v4.1_time_stamps.txt
            Antarctic/ATL11_0314_tide_adj_scale_200m.tif                   130.6 MiB
            Arctic/GreenlandIceMask_2018.1_2026.0_100m_v4.1.tif             22.9 MiB
            Arctic/GreenlandIceMask_..._v4.1_time_stamps.txt
        - the ice masks are TIME-VARYING multi-band geotiffs: 40 bands each, one per time slice,
          the decimal year in a per-band 'time' tag (AA 2018.200-2026.250 on EPSG:3031/240 m;
          GL 2018.00-2025.83 on EPSG:3413/100 m; neither has duplicate bands -- that was the
          v4.0 bug v4.1 fixes).  That is exactly what pc.grid.data.from_gdal reads: it builds
          self.t from the 'time' tags and selects bands with t_range, which is how
          ATL11_to_ATL15 already calls from_file for --mask_file.  The sidecar time_stamps.txt
          is informational; nothing in the code reads it.
        - the tide-adjustment grid is now a single-band geotiff (EPSG:3031, 200 m, float32,
          NaN fill), NOT the .h5 the args used to name -- so AA_0331.txt needed BOTH
          --tide_adjustment_file and --tide_adjustment_format changed (see Defaults below).
          Values over Ross are 0-1, mean 0.993, as expected for a flexure scaling.
        - defaults repointed to match (AA_0331.txt, GL_0331.txt) -- see Defaults below.
        - NOT on Zenodo, so still 134-byte git-lfs POINTERS both in the repo and on S3:
          the Antarctic shelf-only variants (Greene_22_shelf_plus_10m_mask{,_1km,_240m,_full},
          scripps_antarctica_IceShelves1km_v1, bedmap2_thickness_gt_50_plus_sio_shelves) and
          Arctic/{U_Texas_ice_mask_2019_100m,_1km, Ice_Ocean_Bed_100m_2019_compress,
          BedMachineGreenland-2021-04-20_shelf_125m}.  No rel_006 args file names any of them,
          so nothing is blocked; they matter only if an older processing string is revived.
        - cause of those pointers: the mask upload copied the repo tree, which had never been
          lfs-pulled, so anything not separately re-uploaded landed as a pointer.  See the lfs
          item under OTHER.
        - LEFT IN PLACE on the bucket, now unreferenced: the superseded Antarctic mask
          Greene_22_shelf_plus_10m_mask_2018_2026.25_v2_240m.h5 (357 MiB), the superseded
          Arctic/GrimpIceMask_2018.1_2026.0_100m.tif, and the dead 134-byte
          Antarctic/ATL11_0314_tide_adj_scale_200m.h5 pointer.  Delete when you are confident
          in the v4.1 run; the .h5 pointer in particular is a trap, since it has the name the
          old args expected but no content.
    [ ] masks/EGM2008_geoid_h.nc -- rel_006_0331.txt passes --geoid_file for it, it is not in
        the repo (only named in .gitattributes) and nothing matching it is on the bucket.
[ X ] Set up a location fromfile for MAAP (ADE)
    - default_args/MAAP.txt written: --mask_dir, --tide_directory, --ATL14_root, all on the
      ~/my-private-bucket mount.
    - --ATL11_root/--ATL11xo_top were deliberately dropped: they are consumed only by
      make_ATL11_index.py / setup_ATL11_xover.py, which are local-filesystem-only.
    - it does NOT yet name the staged index; a cloud run also needs
      --ATL11_earthaccess and --ATL11_index=<my-private-bucket>/ATL11_index/
[ ] MAAP_dps.txt: default_args/MAAP.txt is ADE-only.  DPS workers have no
    /home/jovyan/my-private-bucket, so a variant with s3:// URIs is needed for anything
    submitted to DPS -- and see the abspath/open() caveat under Software.
[ X ] Figure out how the MAAP algorithm definition handles python fromfiles
    - the @argsfile idiom survives unchanged.  Nothing in DPS parses your arguments;
      fromfile_prefix_chars is entirely client-side.  Register the args file as a `file` input;
      DPS localizes it into input/; run.sh does args_file=$(ls -d input/*) and passes
      @${args_file}.  The ~90 argparse options never have to become DPS parameters.
      Only the per-job values need flattening: a tile job is (x0, y0, step, args_file),
      i.e. 3 positionals + 1 file input.
[ X ] define a build command based on maap example repository
    - WRITTEN: `build-env.sh`, `run.sh`, `algorithm_config.yml` at the repo root, modelled on
      MAAP-Project/dps_tutorial (gdal_wrapper) and maap-documentation-examples/gedi-subset.
      No Dockerfile of our own; the base image is chosen by URL.
    - build_command/run_command are paths ONE LEVEL ABOVE the repo root: ATL1415/run.sh
      (DPS clones to /app/<repo>/ and calls /app/dps_wrapper.sh '/app/<repo>/run.sh' args...)
    - output/ is relative to run.sh and is what gets uploaded; input/ is where file inputs land
    - the algorithm repo must be PUBLIC.  algorithm_version is the git ref DPS clones, and is
      currently `on_s3` -- change it to main once this branch merges.
    - docker_container_url = maap_base:v6.0.0, taken from $DOCKERIMAGE_PATH_DEFAULT in the ADE.
      It is the image that already carries conda, which the build needs.

    build-env.sh
    - conda FIRST, pip second, in that order and for a reason: conda supplies suitesparse (the
      dev headers PySPQR compiles its CFFI bindings against) and gdal.  conda-forge's gdal
      installs the python bindings WITH a .dist-info, so pip treats LSsurf's bare `gdal`
      requirement as already satisfied instead of building it against the wrong libgdal.
    - exports CPATH/LIBRARY_PATH/LD_LIBRARY_PATH at the env prefix before pip runs, and hard
      fails if SuiteSparseQR.hpp is not there afterwards -- otherwise the sparseqr failure
      surfaces as an opaque CFFI error.
    - ends with an import check (numpy scipy h5py osgeo.gdal pyproj sparseqr pointCollection
      LSsurf ATL1415 earthaccess s3fs) so a broken env fails at BUILD time, once, rather than
      on every one of thousands of tile jobs.
    - SMBcorr is deliberately not installed; add the [firn] extra here if --firn_model is ever
      brought back into use.

    run.sh
    - `run.sh <x0> <y0> <step>`, step in {prelim, matched}.  One registered algorithm; the stage
      is dispatched here.
    - finds the argsfile as input/*.txt (by extension, so a second file input cannot be picked
      up by mistake) and passes it as @<path>.
    - passes --THREADS=$(nproc) explicitly and BEFORE @argsfile.  This matters: ATL11_to_ATL15
      scrapes THREADS out of sys.argv at IMPORT time to set MKL/OPENBLAS/NUMEXPR/OMP, and it
      does not look inside the @argsfile -- so a --THREADS in the argsfile sets args.THREADS but
      never the BLAS environment.  Before the argsfile so the argsfile can still override.
    - prelim: --base_directory $PWD/output, so tiles land in output/prelim/ and get uploaded.
      Runs the fit and then the --calc_error_for_xy companion, mirroring the single queue line
      make_ATL1415_queue.py writes for SLURM.  It checks for the tile file in between and exits
      0 if it is absent: a tile with too little data is a normal outcome (ATL11_to_ATL15 returns
      0 without writing anything), but --calc_error_for_xy on a missing file exits 1, which at a
      fan-out of thousands of tiles would mark those jobs failed and bury the real failures.
    - matched: --base_directory $PWD/INPUT, not output.  --matched reads the tile's own prelim
      fit and, through prior_edge_include, its NEIGHBOURS' -- so a matched job needs the
      surrounding prelim tiles localized into input/prelim/ (the tile plus its 8 neighbours at
      minimum), and only the result goes to output/.  It also passes --prior_edge_include 1000
      the way the queue does, because ATL11_to_ATL15's own default is None, which silently drops
      the prior-edge constraints.
    - tolerates a leading non-numeric argument, so it reads x0 correctly whether DPS emits the
      file input before or after the positionals -- see the smoke-test TBD.

    Still unconfirmed, and only a real DPS job can settle them (all folded into the smoke test):
    - whether DPS passes `file` inputs before or after `positional` ones (run.sh is defensive
      about it either way);
    - the queue: `maap-dps-worker-32gb` and disk_space 20GB are guesses.  The SLURM runs use 4
      tasks and a 4-hour walltime per tile; nothing here has been measured.
    - whether the build container reaches github.com, which the pointCollection/LSsurf git+ URLs
      require, and whether sparseqr actually compiles there.
[ X ] Figure out how submitted jobs see my_private_bucket and my_public_bucket in a MAAP account
    - THEY DO NOT.  Workspace mounts do not exist on a DPS worker.  Every read of the ATL11
      index, tide models and masks must become an S3 read.
    - ADE mounts are ~/my-private-bucket = s3://maap-ops-workspace/ben_smith/ and
      ~/my-public-bucket = s3://maap-ops-workspace/shared/ben_smith/
      (the docs' s3://maap-ops-workspace/private/<user>/ prefix is empty for this account)
    - the only worker->workspace coupling is one-directional: DPS uploads output/ to
      <workspace-bucket>/<user>/dps_output/..., which then appears under ~/my-private-bucket
    - job payloads show ~/.aws and ~/.netrc ARE bind-mounted into the worker, so earthaccess
      should authenticate there as it does in the ADE -- unconfirmed, see the smoke-test TBD
[ ] Speed tests and costwise benchmarks
    - nothing in the MAAP docs about per-job or per-queue cost; will have to be measured

## Defaults (processing strings)
[ X ] Repoint the rel_006 0331 args at the published Zenodo v4.1 masks (2026-09-03).
    default_args/AA_0331.txt:
      --mask_file             Antarctic/Greene_22_shelf_plus_10m_mask_2018_2026.25_v2_240m.h5
                           -> Antarctic/AntarcticIceMask_2018.00_2026.25_240m_v4.1.tif
      --tide_adjustment_file  Antarctic/ATL11_0314_tide_adj_scale_200m.h5
                           -> Antarctic/ATL11_0314_tide_adj_scale_200m.tif
      --tide_adjustment_format  h5 -> geotif   (argparse choices are geotif/h5/nc)
    default_args/GL_0331.txt:
      --mask_file             Arctic/GrimpIceMask_2018.1_2026.0_100m.tif
                           -> Arctic/GreenlandIceMask_2018.1_2026.0_100m_v4.1.tif
    - the AA tide-adjustment line was BROKEN before this: the .h5 it named is a 134-byte lfs
      pointer both in the repo and on S3, so --tide_adjustment could never have worked on MAAP.
    - the two ice-mask swaps are a DATA VERSION CHANGE, not just a rename: v4.1 is the published
      successor at the same grid and time span, but it has not been run end-to-end here, and
      neither mask has been compared band-for-band against the file it replaces.
    - reading the .tif through the geotif branch also means field_mapping=dict(z='tide_adj_scale')
      in ATL11_to_ATL15.read_ATL11_data is now inert: from_gdal swallows it in **kwargs and puts
      the band in .z anyway, which is the field .interp() defaults to.  Same result, but the
      kwarg no longer documents anything -- worth deleting if the h5 path is retired.
[ ] Run one AA tile and one GL tile against the v4.1 masks before trusting them in production.
    Nothing on this path has been executed: pointCollection and LSsurf are not importable in the
    ADE notebook env, so the checks above are file-format checks (gdal band tags, projections,
    value ranges), not a solve.

## OTHER: 
[ ] There are lfs-tracked files that get overwritten when ATL1415 gets installed.  These should be cleaned up.
    - a DPS build clones fresh every time, so an install that clobbers repo content will bite there
