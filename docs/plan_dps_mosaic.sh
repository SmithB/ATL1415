# plan_dps_mosaic.sh -- mosaic and netCDF as DPS jobs (no ADE, no /home)
#
# ############################################################################
# ##  WRITTEN 2026-09-30.  TENTATIVE: nothing below has been run or coded. ##
# ##  Purpose: move N6 (mosaic) and N7 (netCDF) off the ADE, for GL north  ##
# ##  first and then every region.  Blocks plan_GL_north.sh NM1-NM6.       ##
# ############################################################################
#
# Status tags per step: TODO / DONE / BLOCKED.  DECIDED = Ben said so;
# RECOMMENDATION = mine; STATEMENT = checked fact, with where it came from;
# QUESTION = open.
#
# WHY (DECIDED, Ben 2026-09-30):
#   - /home/jovyan has a 150 GB quota (MAAP admin); df does not show it.
#     The GL north fetch filled it (plan_GL_north.sh, REVISED block).
#   - No DPS output is copied to /home.  The authoritative tiles are on
#     $s3_out; a local file of a different size is an error.
#   - Mosaic and netCDF run as DPS jobs: "instances are too unpredictable".
#
# ===========================================================================
# QUESTIONS FOR BEN (answer inline)
# ===========================================================================
#   QD1  How does a job read the 557 + 557 tiles?
#        RECOMMENDATION: in place, by s3:// URI -- no download.
#          pointCollection's grid.data.from_h5 already opens an s3:// path
#          (pc/grid/data.py h5_open -> io_utils.open_remote, default AWS
#          credential chain), and mosaic.from_list hands every str item to
#          grid.data/mosaic.from_file -> from_h5, so it takes URIs as it
#          stands.  It opens each tile at least TWICE (a meta_only pass in
#          setup_bounds_from_list, then add()), which D0's timing includes.  What does NOT: make_mosaic.py finds its inputs with
#          glob.glob(directory + '/' + glob_string) (make_mosaic.py:145), and
#          the netCDF writers' lineage and tile-stats passes glob/listdir
#          --tiles_dir (ATL1415_attrs_meta.set_lineage,
#          make_tile_stats_group.py:52).  Those three listings would learn
#          to list an s3:// prefix.  Whether it is fast enough is D0's job.
#        Other options:
#          B. Download every tile to the worker: ~82 GB prelim + ~80 GB
#             matched per job, x 41 mosaic jobs ~= 6.6 TB moved; worker
#             disk size unknown (QD7).  No code change to the readers.
#          C. Download only what each task reads (its group, from each
#             tile).  Least transfer, most new code.
# AD1: in place
#
#   QD2  Job granularity.
#        RECOMMENDATION: ONE DPS JOB PER MOSAIC TASK, then the netCDF jobs.
#          make_mosaic_jobs.py already splits the work into independent
#          tasks, each writing one .h5 (quarterly GL: 41 = z0 + dz + 9
#          dzdt lags + 3 averaging scales x (dz + 9 lags); lags
#          [1,2,4,8,12,16,20,24,28] from infer_dzdt_lags(0.25, 2018.75-2026.5);
#          monthly IS was 44).  The job would run make_mosaic_jobs.py itself
#          and execute task N, so the task list is never written twice.
#          A failed task is retried alone.  netCDF: one job for ATL14 and
#          one for ATL15, submitted after all mosaic tasks succeed.
#        Other option: one job runs all 41 at -P 4 (the ADE way).  Fewer
#          submissions, but one failure reruns everything, and one worker's
#          walltime and disk hold the whole region.
#  AD2: One DPS job per mosaic task
#
#   QD3  How are the new steps registered?
#        RECOMMENDATION: the SAME algorithm (atl1415_tile_solve), new
#          `step` values `mosaic` and `nc` in run.sh's dispatch, and ONE new
#          input `task` (default "-", never "": the build form drops empty
#          defaults, howto_MAAP_ogc F6).  STATEMENT: one algorithm has one
#          run_command, and the conda+SuiteSparse build is the expensive part
#          (algorithm_config.yml header) -- a second algorithm would be a
#          second identical build.  Reusing x0 as the task number would
#          avoid touching the CWL, but it is a hack a ledger reader would
#          trip on.  One re-registration either way (you register from my
#          checkout; memory: held commits block registration), then
#          check_build_id MATCH.
#
#  AD3: Same algorithm
#
#   QD4  Where do mosaics and netCDFs go?
#        RECOMMENDATION: mosaics at the region prefix, beside the tiles --
#            $s3_out/{z0,dz,dz_10km,...,dzdt_40km_lag28}.h5
#          -- the same layout as a local region dir, so the netCDF writers
#          find them by name; a full-GL run overwrites them.  STATEMENT:
#          this is where IS's netCDFs already are
#          (s3://.../ATL14_processing/rel006/north/IS/ATL14_IS_0332_*.nc).
#        BUT the GL-north netCDFs must NOT take the canonical names at
#          $s3_out: the monthly args would read that as the full-run
#          reference (plan_GL_north.sh N7, QM-A).  So the nc step needs an
#          output location separate from its input tree.  Options:
#          A. a `-` / s3:// `out_prefix` input (a second new input);
#          B. the nc job always writes <tile_prefix>/nc/ and a human copies
#             the product to the canonical name (or the QM-A side key) --
#             RECOMMENDATION, B: no second input, and publishing stays a
#             deliberate act, which N7 already required.
#
#  AD4: Allow additional text in the release directory (e.g. s3://.../ATL14_processing/rel006_0332_v1/north/IS/ATL14_IS_0332_*.nc).  /rel006/ remains the default, but scripting should allow different directory names
#  AD4 (Ben 2026-09-30, follow-up): the nc output of B goes to a
#       rel006_0332_testing release directory:
#         s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006_0332_testing/north/GL/
#       CONSEQUENCE (mine, stated): the nc job reads mosaics + prelim tiles
#       from <tile_prefix> and writes somewhere else, so it needs an output
#       prefix after all -- a second new input `out_prefix` (default "-" =
#       write to <tile_prefix>).  Mosaics stay at the region prefix (QD4
#       recommendation, not objected to).
#
#   QD5  Queue.  RECOMMENDATION: maap-dps-worker-32gb for every mosaic and
#        nc job, as N3/N5, and measure (D0 gives the ADE peak first).
#        ESTIMATE, not measured: the z0 grid for GL north is ~14,600 x 9,400
#        cells at 100 m (~1.1 GB per float64 array); make_mosaic holds one
#        field plus weight and invalid arrays, so a few GB.  Full GL ~3x.
#        IS z0.h5 is 88 MB for 28 tiles -> GL north ~1.8 GB, under
#        outdir_max 20 (algorithm_config.yml).
#   AD5: 32gb for the Antarctic z0 jobs, 16GB for everything else, can revisit if there are problems
#
#   QD6  Where do the checks run?  RECOMMENDATION: on the ADE, READING
#        through ~/my-private-bucket (the mountpoint-s3 view of
#        s3://maap-ops-workspace/ben_smith): check_field_sizes.py,
#        check_mosaic_outputs.py --values and the NM2 rel005 compare only
#        read, and random reads work on that mount; nothing is copied to
#        /home.  (mountpoint-s3 cannot WRITE an HDF5 file: no random writes,
#        no rename -- memory: MAAP bucket is mountpoint-s3.)
#   AD6: on the ADE
#
#   QD7  FOR MAAP (you ask the admin, if we need it): how much local disk
#        does a DPS worker have, and is mountpoint-s3 (or any bucket
#        mount) available on a worker?  Only matters for QD1 B, or if D0
#        says in-place reads are too slow.
#   AD7: Not relevant b/c we're going with the in-place option
# 
# ===========================================================================
# WHAT IS KNOWN (statements, with provenance)
# ===========================================================================
#   - $s3_out for GL north holds 557 prelim + 557 matched tiles (counted
#     2026-09-30); IS holds 28 + 28 plus its five netCDFs.
#   - IS on the ADE (plan_rerun_timing.sh B5): 41 mosaic tasks in 42-59 s
#     at -P 12; the ADE mosaics and netCDFs are still in
#     ~/ATL14_processing/rel006/north/IS (1.4 GB, made 2026-09-25).  They
#     are the reference D6 compares against.
#   - run.sh and algorithm_config.yml both SAY mosaic and to-netcdf stay in
#     the ADE; Transition_to_maap.md Q4/Q18 asked whether the ADE could
#     mosaic AA and never settled it.  All three change with this plan.
#   - A matched job already reads its 3x3 prelim neighbourhood from
#     <tile_prefix>/prelim/ with scripts/s3_tiles.py, and writes its tile
#     back there: the tile tree is addressed by name, not by dps_output.
#   - The netCDF writers read the PRELIM tiles as well as the mosaics:
#     lineage (/meta input_files and granule attrs) and per-tile stats
#     (data/three_sigma_edit, RMS, bias) -- so an nc job reads all 557
#     prelim tiles too, not just ~42 mosaics.
#
# ===========================================================================
# STEPS
# ===========================================================================
# D0. [ADE] DONE 2026-09-30 20:51Z.  Probe in-place reads (QD1, QD5).
#     RESULT, 557 GL-north matched tiles, ADE, one process, sequential
#     (scratchpad d0_probe.py; outputs in /tmp only):
#       run                         wall    net rx   peak RSS  CPU
#       (a) avg_dz_40km  s3, 1 MiB   357 s   2.2 GB   0.29 GiB  0.06 cores
#       (a)              mount       386 s   2.1 GB   0.08 GiB  0.01
#       (b) z0/z0 -w     s3, 1 MiB   802 s   7.6 GB   4.49 GiB  0.14
#       (b)              s3, 256 KiB 841 s   3.3 GB   4.33 GiB  0.13
#       (b)              mount       854 s  13.2 GB   4.06 GiB  0.08
#     Every output IDENTICAL (array_equal, NaN-aware) across s3 1 MiB /
#     256 KiB / mount; z0 grid 9401 x 14601, 91.8 M finite cells.
#     FINDING 1: pointCollection's mosaic.from_list does not pass a
#       block_size to grid.data.from_h5, so s3:// tiles are read with
#       s3fs's 50 MiB default: 4 tiles pulled ~1 GB for a few-kB group,
#       ~3.8 s/tile.  The runs above force it (probe monkeypatch).  D1 must
#       pass it through (1 MiB: fastest here; 256 KiB: half the bytes).
#     FINDING 2: LATENCY-BOUND, not bandwidth- or CPU-bound: 0.6 s/tile for
#       (a), whose data are tiny -- that is per-open cost (each tile is opened
#       twice: setup_bounds_from_list's meta pass, then add()) -- and 1.4
#       s/tile for (b); CPU <= 0.14 cores.  The mount is no faster.
#     WHAT IT MEANS PER TASK (estimate from the above, not measured):
#       make_mosaic_jobs.py writes the z0 task as 7 separate make_mosaic.py
#       calls (6 matched fields + sigma_z0), each re-reading all 557 tiles
#       -> ~7 x 13 min ~= 1.5 h for GL north's z0 task; the other 40 tasks
#       are 2 calls each (matched + prelim sigma) at ~6-13 min -> ~15-25 min.
#       Full GL (~1483 tiles) ~2.7x; AA far more.
#     z0 peak 4.5 GiB for GL north -> full GL ~12 GiB, close to the 16gb
#       queue (AD5 gives 32gb only to AA z0).
#     POSSIBLE SPEEDUPS (none tried): (i) open tiles concurrently (threads;
#       CPU is idle); (ii) read all z0 fields in ONE pass per tile instead of
#       7; (iii) take bounds from the tile names (E<x>_N<y>) instead of the
#       meta_only pass.  Each is a pointCollection/make_mosaic_jobs change.
#     GATE (yours): are these times acceptable for DPS jobs, or do speedups
#       go into D1?
#     GATE ANSWER (Ben 2026-09-30): "Add the functionality to pointCollection
#       and I'll merge.  For Greenland, add the 200-km tile step that is
#       currently in effect for Antarctica to the workflow so that the
#       tiles -> mosaic jobs can be run in parallel as small tasks."  DECIDED.
#     CONCURRENCY MEASURED (40 GL tiles, z0/z0, 256 KiB blocks, ADE): serial
#       32.3 s; 8 THREADS 32.6 s (no gain -- h5py holds its global lock
#       through the file-object reads); 8 processes spawn 14.3 s, FORK 6.6 s.
#       So D1 uses a process pool, not threads.
#     A short script in the scratchpad: pc.grid.mosaic().from_list() over
#     the 557 s3:// matched URIs for (a) a cheap group, avg_dz_40000m, and
#     (b) the heaviest field, z0/z0 with -p 5000 -f 10000; output to /tmp,
#     never /home.  Record wall time, peak RSS (run_with_rusage.py), bytes
#     read.  Compare (a) with the same mosaic read through
#     ~/my-private-bucket (local-path glob) as a second data point.
#     GATE for QD1 A: (b) finishes in a time you accept for a DPS job.
#
# D1. [code, pointCollection] PR OPEN, for Ben to merge:
#     https://github.com/SmithB/pointCollection/pull/60 (branch
#     mosaic_remote_parallel, 72ec08b).  306 tests pass.  Real data, 557 GL
#     tiles on S3: 40 km avg 357 s -> 44.5 s (8 workers), z0 802 s -> 195 s
#     (4 workers), both bit-identical to D0's serial mosaics.  Workers start by
#     FORKSERVER (plain fork after s3fs: "This class is not fork-safe").  Each
#     worker ~0.35 GiB: 8 workers + z0 was OOM-killed on the ADE, whose cgroup
#     memory.max is 7.3 GiB (free shows the 30 GB host).  Also fixes a
#     pre-existing add_to_band bug (in-memory 3-D mosaic + by_band; Ben: add).
#     AFTER MERGE: reinstall pC in the ADE ATL14 env; DPS builds pick it up
#     (LSsurf/pC unpinned on DPS -- memory: reach kernel).
#   D1 as planned:
#     a. io_utils.glob_remote(pattern): the remote glob.glob (sorted URIs).
#     b. make_mosaic.py: --directory may be a URI (listed with a.); -O must
#        then be a LOCAL absolute path (error otherwise -- a relative -O would
#        be joined onto the URI); new --block_size (default for remote:
#        io_utils.DEFAULT_REMOTE_BLOCK_SIZE) and --workers.
#     c. grid.mosaic.from_list(block_size=None, workers=1): block_size is
#        passed to from_file for remote h5/nc items (FINDING 1); workers > 1
#        reads the tiles in a process pool (fork where available; the pool
#        clears io_utils' s3fs session cache in each child), in list order,
#        at most 2 x workers ahead, so the summation order -- and the output
#        -- is the same as the serial read.  Covers the meta pass, the
#        weighted, by_band and replace loops.  workers=1 is today's code path.
#     NOT DONE: (iii) bounds from tile names -- the E<x>_N<y> convention is
#       ATL1415's, not pointCollection's; with the 200 km step each task's
#       meta pass is ~36 tiles and runs in the pool.
#
# D2. [code, ATL1415] DONE 2026-10-01.  paths.list_tiles/open_tile (sorted
#     listing; remote = default AWS chain + pC small block, closed on exit);
#     set_lineage and make_tile_stats_group use them.  Real data: 557 GL
#     prelim tiles listed in 1.2 s; tile-stats reads 1.0 s/tile, identical to
#     the mount read -> ~10 min per pass, two passes per writer (lineage,
#     tile stats), serial.  A pool (as pC's) is the lever if that matters.
#     Also: make_tile_stats_group imports make_nc_projection_variable by full
#     path (`from ATL1415 import` gave the MODULE after __init__'s fallback
#     lookup had imported the submodule -- test-order dependent).
#   D2 as planned: --tiles_dir may be s3://: set_lineage and
#     make_tile_stats_group list the prefix and open tiles by URI.  The
#     writers still read mosaics from -b and write the .nc into -b, so the
#     job keeps -b local.  Tests beside the existing ATL1415 suite (181).
#
# D2b. [code, ATL1415] DONE (code; deploys with the next registration).
#     DECIDED (Ben 2026-09-30, option b): a
#     --no_data_group flag on ATL11_to_ATL15.py, OFF by default, so
#     save_fit_to_file skips /data (80-92% of a matched tile: E80_N-920 192
#     of 208 MiB, E520_N-920 21 of 26 MiB).  Nothing downstream reads a
#     matched tile's /data (matched and error read the PRELIM tile; tile
#     stats and lineage read --tiles_dir = prelim/; mosaics read grids only;
#     only scripts/check_tile_data_vs_DEM.py, a diagnostic, reads any tile's).
#     What a matched fit changes in /data, for the record: z_est and
#     sigma_extra everywhere, three_sigma_edit possibly (up to 6 edit
#     iterations), a few points dropped outside the grids (18 of 131,408 on
#     E520_N-920), `editable` added.
#     CONSTRAINT: prelim and matched share ONE args file on MAAP, and prelim
#     MUST keep /data (matched and error reread it).  So the flag cannot live
#     in input_args_<R>.txt.  RECOMMENDATION: run.sh passes it on the matched
#     command line only (as it does --prior_edge_include), and
#     ATL11_to_ATL15.py exits with an error if it is given without --matched
#     (fail loudly, never a silent prelim tile without /data).
#     The 557 GL-north matched tiles already on S3 keep their /data.
#
# D2c. [code] DONE 2026-10-01.  --release_dir_suffix (e.g. _0332_testing):
#     setup_ATL1415_region.py and make_GL_ATL1415_queue.py name the directory
#     rel<Release><suffix>; --Release (product names) stays 006; the suffix is
#     not written to input_args (-b carries it, like --hemi_suffix).
#     make_ATL1415_queue.py checks the parents of -b instead of a rebuilt
#     rel<Release>/<hemi> (UNTESTED: that script has no tests).
#     HITS LISTED AND LEFT ALONE (discover-only or docs): scripts/run_arctic_*.sh,
#     scripts/run_antarctic_tonc.sh, link_200km_tiles_AA_monthly.py (hard-coded
#     discover path), scripts/old/.  MAAP scripts take their s3 prefixes as
#     given (submit/collect/fetch, s3_tiles.py, run.sh): nothing to change.
#   D2c as planned: AD4: nothing may assume the release
#     directory is exactly rel<NNN>.  Find every place that builds or parses
#     .../ATL14_processing/rel<NNN>/<hemi>/<R> (setup_ATL1415_region.py,
#     default_args, submit/collect/fetch scripts, run.sh comments) and let
#     the release directory carry extra text (rel006_0332_testing), rel006
#     staying the default.  List the hits before changing any.
#
# D3a. [code, ATL1415] TODO.  THE 200 KM STEP FOR GL (Ben, gate answer),
#     as AA does it today:
#       stage 1 -- make_200km_tiles.py: one task per 200 km tile (GL north
#         ~40), each running make_mosaic.py for every group with -r (tiles by
#         name within 10 km of the 200 km square) and -c (crop), matched then
#         prelim sigma, writing <region>/200km_tiles/<group>/<group><bounds>.h5.
#       stage 2 -- make_200km_to_mosaic_jobs.py: one task per group, joining
#         the 200 km tiles into <region>/<group>.h5 (make_mosaic_jobs.py
#         already switches z0 to 200km_tiles/z0 when that directory exists).
#     On DPS: stage 1 = one job per 200 km tile, stage 2 = one job per
#     group; the tiles are read in place (-d s3://..., D1) and each job
#     uploads its outputs.  make_200km_tiles.py finds centers from a local
#     glob of <region>/prelim -- on DPS from glob_remote or a
#     200km_tile_list.txt.  ESTIMATE: stage 1 = 41 groups x 2 calls x ~36
#     tiles per job; stage 2 reads ~40 200 km tiles per field.
#     OPEN (for Ben, later): make_200km_to_mosaic_jobs.py hard-codes
#     'source activate IS2' and one make_mosaic call per field; it has had no
#     tests.  Check it against make_mosaic_jobs.py before relying on it.
#
# D3. [code, ATL1415] TODO.  run.sh: step `mosaic` --
#       make_mosaic_jobs.py -b <local work dir> ... with -d s3://<tile_prefix>
#       for the tile reads, run task $task, upload its one .h5 to
#       <tile_prefix>/ with s3_tiles.py (or a sibling), print the build
#       stamp and worker facts as every job does, keep ./output* non-empty.
#     step `nc` -- fetch the ~42 mosaics (a few GB) into the work dir, run
#       ATL14_write2nc.py or ATL15_write2nc.py (task = ATL14 | ATL15) with
#       --tiles_dir <tile_prefix>/prelim, upload to <out_prefix>/ (AD4
#       follow-up; out_prefix "-" = tile_prefix).
#     algorithm_config.yml: inputs `task` and `out_prefix` (defaults "-"); update the description and the
#     header that says mosaic stays in the ADE.
#     submit_MAAP_jobs.py: a --step mosaic mode that submits task 1..N from
#     the same make_mosaic_jobs.py count, ledgered like a tile run.
#     OPEN DETAIL: make_mosaic_jobs.py puts tile globs and outputs under ONE
#     base (-d {base} ... -O {base}/z0.h5).  Either give it a separate
#     tiles base, or have run.sh rewrite the task's -d.  RECOMMENDATION:
#     a --tiles_base option on make_mosaic_jobs.py (default: base), so the
#     ADE/discover behaviour is unchanged.
#
# D4. [Ben] TODO.  Register from my checkout (clean, pushed); then
#     scripts/maap/check_build_id.py -> MATCH (memory: only MATCH proves the
#     deploy).
#
# D5. [DPS] TODO.  Smoke: ONE mosaic task (avg_dz_40000m, task 2 or so) for
#     IS.  GATE: successful; the .h5 is at the IS prefix; identical to the
#     ADE file in ~/ATL14_processing/rel006/north/IS (same tiles, same code
#     path apart from the listing).
#
# D6. [DPS] TODO.  IS end to end: 41 mosaic jobs, then ATL14 and ATL15 nc
#     jobs, to a TEST prefix (never over IS's canonical products) --
#     RECOMMENDATION: .../ATL14_processing/rel006_0332_testing/north/IS,
#     for BOTH the mosaics and the netCDFs here, since IS's canonical prefix
#     already holds ADE-made products.  GATE:
#     every mosaic identical to the ADE's; netCDFs identical in data and
#     attributes apart from dates/build fields (list the differences, don't
#     assume).  Record per-job time and peak memory.
#
# D7. [DPS] TODO.  GL north: 41 mosaic jobs, ATL14 + ATL15 nc jobs (to
#     out_prefix .../ATL14_processing/rel006_0332_testing/north/GL/, AD4).  ADE checks read through the mount (QD6):
#     check_mosaic_outputs.py --values; quick plot of h and delta_h (N7).
#     Times and memory go to plan_GL_north.sh N8 for MAAP.
#
# D8. [ADE] TODO.  Resume plan_GL_north.sh at NM2 (compare, report only),
#     NM3 (copy the nc S3 -> the QM-A side key; aws s3 cp, S3 to S3), then
#     NM4-NM6 as written, minus every fetch.  The two NM drivers are
#     superseded; a new driver, if any, only submits and waits.
#
# D9. [docs] TODO.  Transition_to_maap.md (Q4/Q18 answered: DPS),
#     run.sh header, howto_MAAP_*.sh: drop every `aws s3 sync ... prelim/`
#     into /home and every fetch_tiles.py step; point checks at the mount.
