# Transition to MAAP

## TBDs: 
[ ] Request an organizational DPS queue, and confirm account status, from the MAAP platform team.
    Public queues are throttled to ~10 jobs/hr; this account reports status=inactive and
    organizations=[].  A per-tile fan-out of thousands of jobs is infeasible until this is
    resolved.  Do this first -- it has days of latency and it gates the shape of everything else.
[ ] Smoke-test one sandbox DPS job to settle what the docs do not say:
    - does the build container reach github.com (for the git+ pip deps)?
    - is ~/.netrc really bind-mounted into the worker (earthaccess auth depends on it)?
    - will a `file` input accept an s3://maap-ops-workspace/ben_smith/... URL?

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
    [ ] Download Antarctic and Greenland masks from Zenodo (inc. tide adj mask and Antarctic shelf only)
        - staged and REAL under s3://maap-ops-workspace/ben_smith/ATL1415/masks/ (2026-09-03):
          Arctic/{GrimpIceMask_2018.1_2026.0_100m, GL_d3z_dx2dt_scaling_map, GL_Ed2z0dx2,
          BedMachineGreenland-v6_shelf_edited}.tif, Antarctic/{BedMachineAntarcticaOceanv2,
          AA_Ed2z0dx2}.tif and Antarctic/Greene_22_shelf_plus_10m_mask_2018_2026.25_v2_240m.h5
          (357 MiB).  That covers every mask GL_0331.txt names, and all but one of AA_0331.txt's.
        - STILL MISSING: Antarctic/ATL11_0314_tide_adj_scale_200m.h5, the tide-adjustment grid
          AA_0331.txt requires (--tide_adjustment).  Both the repo copy and the S3 copy are
          134-byte git-lfs POINTERS; the real object is 203 MB.  Same for the shelf-only
          variants (Greene_22_shelf_plus_10m_mask{,_1km,_240m,_full},
          scripps_antarctica_IceShelves1km_v1, bedmap2_thickness_gt_50_plus_sio_shelves) and
          for Arctic/{U_Texas_ice_mask_2019_100m,_1km, Ice_Ocean_Bed_100m_2019_compress,
          BedMachineGreenland-2021-04-20_shelf_125m}.
        - cause: the mask upload copied the repo tree, which had never been lfs-pulled, so
          anything not separately re-uploaded landed as a pointer.  See the lfs item under OTHER.
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
[ ] define a build command based on maap example repository
    - three files in the repo, plus a base image chosen by URL (no Dockerfile of your own):
        build-env.sh  -- conda env update -f environment.yml ; conda run -n <env> pip install .
        run.sh        -- mkdir -p output ; args_file=$(ls -d input/*) ; conda run -n <env> ...
        algorithm_config.yml -- algorithm_name/_version, repository_url, build_command,
                                run_command, docker_container_url, queue, disk_space, inputs
    - build_command/run_command are paths ONE LEVEL ABOVE the repo root: ATL1415/run.sh
      (DPS clones to /app/<repo>/ and calls /app/dps_wrapper.sh '/app/<repo>/run.sh' args...)
    - output/ is relative to run.sh and is what gets uploaded; input/ is where file inputs land
    - the algorithm repo must be PUBLIC
    - model on MAAP-Project/dps_tutorial (gdal_wrapper) and maap-documentation-examples/gedi-subset
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

## OTHER: 
[ ] There are lfs-tracked files that get overwritten when ATL1415 gets installed.  These should be cleaned up.
    - a DPS build clones fresh every time, so an install that clobbers repo content will bite there
