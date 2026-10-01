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
# D3 DETAILED DESIGN (written 2026-10-01; TENTATIVE until coded + tested)
#   DECIDED (Ben 2026-10-01): every job reads the solve TILES from
#   <tile_prefix> and reads and writes every DERIVED product -- 200 km tiles,
#   region mosaics, netCDFs -- at <out_prefix> ("-" = tile_prefix).  GL north:
#   tile_prefix .../rel006/north/GL, out_prefix
#   .../rel006_0332_testing/north/GL for all three stages; nothing partial
#   lands at the canonical region prefix.  (Supersedes QD4's "mosaics at the
#   region prefix".)
#   Three new run.sh steps, one DPS job each; new inputs `task`, `out_prefix`:
#     step mosaic200  task = "<x>_<y>" (200 km center, m).  Stage 1.
#       make_200km_tiles.py --center x y --tiles_base <tile_prefix> writes
#       that one tile's task into a local work dir; its ~41 groups run in
#       parallel (-P physical cores), each group's two lines (matched, then
#       prelim sigma) in order; each make_mosaic call reads ~36 tiles by
#       name (-r) straight from S3.  Uploads 200km_tiles/<group>/*.h5 to
#       <out_prefix>.  WHY groups-in-parallel, not make_mosaic -j: every
#       make_mosaic call is its own process, and -j would start a forkserver
#       (~10 s) in each of 82 calls to save ~30 s of reads per call.
#     step mosaic     task = <group> (z0, dz, dzdt_lag4, avg_dz_40000m, ...).
#       Stage 2: make_200km_to_mosaic_jobs.py's commands for that group,
#       reading <out_prefix>/200km_tiles/<group>/ with -j (one call per
#       field, ~40 200 km tiles each), writing <group's file>.h5 locally,
#       uploading to <out_prefix>.
#     step nc         task = ATL14 | ATL15.  Fetches the mosaics from
#       <out_prefix> into -b (a local work dir), --tiles_dir
#       <tile_prefix>/prelim (D2, read in place), writes the .nc files,
#       uploads them to <out_prefix>.
#   Script changes (each with tests):
#     D3a-1 DONE 2026-10-01.  Real data: GL north = 32 200 km tiles (from the
#           S3 listing); task for (100000, -900000) = 82 make_mosaic lines;
#           its z0 line (6 fields, ~36 tiles read in place) 63 s, 0.77 GiB;
#           the 200 km z0 tile is BIT-IDENTICAL to the same 2001 x 2001 square
#           of D0's full-region z0 mosaic.
#     D3a-1 as planned: make_200km_tiles.py: --tiles_base (where the tiles are read, -d,
#           and where the centers are listed from; may be s3://; default
#           region_dir) and --center X Y (write only that tile's task).
#     D3a-2 DONE 2026-10-01.  Restructured around mosaic_commands(); writes
#           the same 41 GL tasks as before (diffed) except z0 now activates the
#           env like the rest.  REAL DATA, stages 1+2 for avg_dz_40000m over
#           all 32 GL-north 200 km tiles from S3: stage 1 129 s (-P 4), stage
#           2 2 s; the region mosaic is BIT-IDENTICAL to D0's direct-path one
#           on the 31 shared bands; the 200 km path drops 2018.75 (31 vs 32
#           bands) -- the t_range STATEMENT above, now measured.
#     D3a-2 as planned: make_200km_to_mosaic_jobs.py: callable per group (--group), reads
#           from --in_base (may be s3://), writes into -b; z0 included by an
#           explicit flag rather than a local isdir() check; the hard-coded
#           'source activate IS2' becomes -e/--environment (as
#           make_mosaic_jobs.py); first tests for it.
#     D3-1  run.sh steps + algorithm_config.yml inputs task, out_prefix.
#     D3-2  submit_MAAP_jobs.py: --step mosaic200 (centers from the tile list
#           on S3), mosaic (groups from make_fields), nc; ledgered as tiles.
#   STATEMENT, to check in D6: stage 1 crops time with --t_range [2019, ...]
#     (dzdt: 2019 + lag*dt/2), so GL's 2018.75 band is dropped from the 200 km
#     mosaics; make_mosaic_jobs.py's direct path keeps it.  The writers crop
#     with --t_crop=2019,... -- whether the netCDFs then agree is D6's to show.
#   CHECKED: make_200km_to_mosaic_jobs.py's groups, fields and output names
#     match make_mosaic_jobs.py's (z0 7 fields; dz 7; dzdt_lagN and avg
#     groups 3 each; dz_40km.h5, dzdt_40km_lag1.h5, ...); stage 1's pad/feather
#     from W=60 km, spacing 40 km = 5000/10000, 0/0 for the 40 and 20 km
#     averages -- as the direct path.
#
# D3-1/D3-2 CODED 2026-10-01 (run.sh steps, algorithm_config inputs task +
#   out_prefix, ATL1415/mosaic_groups.py, submit_MAAP_jobs.py --step
#   mosaic200|mosaic|nc).  LOCAL END-TO-END TEST ON IS (run.sh as DPS would
#   call it, tiles read from S3, products to a scratch prefix): all 4 stage-1,
#   41 stage-2 and 2 nc jobs run (after a manifest fix).  COMPARED WITH THE
#   ADE'S IS PRODUCTS (direct path, 2026-09-25):
#     identical: every averaged mosaic (dz_10/20/40km, dzdt_*km_lag*) on the
#       shared bands; ATL15 20 km and 40 km netCDFs (88/88 variables, and
#       attributes but uuid/date -- so lineage + tile stats read from S3 agree);
#     DIFFERENT: z0, dz, dzdt_lagN (up to 3.5 km in z0, 14 m in dz), ATL14,
#       ATL15 1 km and 10 km.  CAUSE (measured): IS tile centers sit at 20 km
#       offsets from the 200 km grid (all centers mod 40 km = 20).  Stage 1
#       selects tiles by center within 10 km of the 200 km square (-r), but a
#       60 km tile centered 20 km outside still reaches ~5 km inside; 90% of
#       differing cells are within 4.5 km of a 200 km edge, all NaN in the
#       200 km path where the direct path has data.  Also: squares derived from
#       tile CENTERS miss region-edge strips (IS x 990-1000, 1400-1410 km).
#     GL is not affected: GL centers are multiples of 40 km, aligned with the
#       200 km grid (the GL z0 200 km tile was bit-identical).  AA: to check.
#   QD8 (asked 2026-10-01): fix make_200km_tiles.py -- search window W/2 and
#     squares from tile extents?
#   AD8 (Ben 2026-10-01): NO FIX.  "The arctic regions that are not Greenland
#     do not need the 200-km step.  Just use the 200-km step for Antarctica
#     and Greenland."  So make_200km_tiles.py stays as it is, and IS, SV, CN,
#     CS, RA, AK mosaic DIRECTLY from the solve tiles (make_mosaic_jobs.py, as
#     the ADE and discover do).
#
# D3b. [code, ATL1415] THE DIRECT MOSAIC PATH ON DPS, for every region but
#     GL and Antarctica (AD8).  Written 2026-10-01; CODED AND TESTED
#     2026-10-01 (D3b-5).  The step names and inputs do not change: `mosaic` task=<group>
#     is one region mosaic either way; only where it reads from differs.
#     D3b-1 DONE 2026-10-01.  ATL1415/mosaic_groups.py: uses_200km_tiles(region) -- True
#       for GL, AA, A1..A4, False otherwise.  THE ONE PLACE the rule lives;
#       run.sh and the submitter both ask it.
#     D3b-2 DONE 2026-10-01.  make_mosaic_jobs.py: --tiles_base (where the solve tiles
#       are read, default --base_dir; s3:// on DPS -- D3's OPEN DETAIL,
#       RECOMMENDATION taken), --group (that one group's task only, as
#       task_1; names as make_fields: z0, dz, dzdt_lag4, avg_dz_40000m,
#       avg_dzdt_40000m_lag4), -j (make_mosaic.py -j), -e '' (no activate
#       line).  Defaults leave the ADE/discover task files unchanged.
#     D3b-3 DONE 2026-10-01.  run.sh: step mosaic -- uses_200km_tiles(region) ? the
#       200 km stage 2 (as now) : make_mosaic_jobs.py --tiles_base
#       <tile_prefix> --group <task>.  step mosaic200 for any other region:
#       ERROR naming the region and the rule, exit 2, nothing run.
#     D3b-4 DONE 2026-10-01.  submit_MAAP_jobs.py: --step mosaic200 refuses a region
#       without the 200 km step; --step mosaic for such a region lists every
#       group (z0 included: no 200 km z0 tiles to wait for).
#     D3b-5 DONE 2026-10-01.  LOCAL IS END-TO-END run.sh test: 41 mosaic
#       jobs + 2 nc, all exit 0 (943 s in all, longest job 95 s, peak RSS
#       0.82 GiB); scratch prefix deleted.  mosaic200 for IS refused, exit 2.
#       vs the ADE's IS products of 2026-09-25, all bands, all 46 files:
#         netCDFs: every variable of all five identical (ATL14 27/27, ATL15
#           1 km 91/91, 10/20/40 km 88/88); attributes differ only in
#           identifier_file_uuid and creationDate.
#         mosaics: same NaN pattern everywhere; the averaged mosaics, masks
#           and sigmas identical; the WEIGHTED fields of z0, dz, dzdt_lagN
#           differ at rounding level (z0 6.8e-13 m, dz 7.1e-15 m, cell_area
#           2.3e-10 m^2) -- not bit-identical.  CAUSE: not measured; the
#           order tiles enter the weighted sum is the candidate (S3 listing
#           sorted vs local glob).
#           DECIDED (Ben 2026-10-01): these differences are not important;
#           the cause is not pursued.  Gates below read "identical apart
#           from rounding in the weighted fields".
#       as planned: GATE every mosaic and every netCDF identical to the
#       ADE's (netCDFs apart from uuid/creationDate).
#     AA AND GL AGAINST THE TWO IS FAULTS (measured 2026-10-01, from
#       ATL1415/resources/{AA,GL}/40km_tile_list.txt; Ben pointed at AA's):
#       (1) tiles dropped by the center search: NONE.  All 8944 AA and 1483
#         GL centers are multiples of 40 km, so a tile either has its center
#         within the 200 km square (+10 km) or does not reach into it.
#         (IS, centers at 20 km offsets: 14 tile/square overlaps dropped.)
#       (2) squares never made: a tile centered ON a 200 km line reaches
#         25 km (W/2 - pad) into the next square, which exists only if some
#         tile center lies in it.  AA: 9 such squares, reached by 19 tiles
#         (e.g. square (-2100,-300) km by E-2000_N-200..-320).  GL: 5, by 13
#         tiles; GL north as run: 2, by 4 tiles.
#         MEASURED for GL north (the 4 matched tiles on S3): 0 cells with
#         cell_area > 0 strictly inside the unmade squares, z0 and dz -- no
#         data lost.  AA: STATEMENT, not measurable yet (no AA tiles on S3);
#         to repeat on the 19 tiles once AA prelim exists.
#       Also: resources/AA/200km_tile_list.txt has 413 squares, the rule
#         gives 411 from the 40 km list (2 in the file only).
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
# D4. [Ben] DONE 2026-10-01.  Register from my checkout (clean, pushed); then
#     scripts/maap/check_build_id.py -> MATCH (memory: only MATCH proves the
#     deploy).
#     RESULT: 6663090 pushed and registered; check_build_id job a4476420
#     VERDICT MATCH (image, live git and CWL all 6663090; built 15:49:40Z;
#     maap_pgt=set).
#
# D5. [DPS] DONE 2026-10-01 (run 1 failed on credentials, run 2 passed).  Smoke: ONE mosaic task for
#     IS (direct path, D3b: no mosaic200 jobs for IS).  GATE: successful; the
#     .h5 is at the out_prefix; identical to the ADE file in
#     ~/ATL14_processing/rel006/north/IS apart from rounding in the weighted
#     fields (same tiles, same code path apart from the listing).
#     CHANGED from the plan (Ben 2026-10-01: z0 is fine as the smoke): task
#     z0, not avg_dz_40000m -- the submitter can only cut its list with
#     --limit, and z0 is first; out_prefix
#     .../ATL14_processing/rel006_0332_testing/north/IS (D6's test prefix),
#     not IS's canonical prefix.  The ADE reference z0.h5 was made WITH -w
#     (runs/IS_0332_cholmod_w_mosaic, 2026-09-25), so no seam differences
#     are expected.
#     RUN 1: job 3acdb8cc (IS_mosaic_smoke_z0, 16gb, t3a.xlarge
#     i-06407490ca8fee8a7, build 6663090), ledger
#     maap_ledgers/IS_0332_mosaic_smoke_jobs.csv.  FAILED at 97 s; nothing
#     uploaded (test prefix still empty).
#     STATEMENT (triaged_job _stderr.txt): botocore NoCredentialsError,
#       "Unable to locate credentials", raised in pointCollection
#       glob_remote at the start of ONE make_mosaic.py call, 74 s into the
#       task (26 s CPU).  Earlier in the SAME job the same default
#       credential chain worked twice (args file fetch; s3_tiles get_glob
#       for bounds.txt).
#     STATEMENT: the container gets no AWS variables or credential files
#       (_docker_params.json), so the chain ends at the instance metadata
#       service; botocore's default for it is 1 attempt, 1 s timeout
#       (botocore/utils.py), and aiobotocore reads
#       metadata_service_num_attempts / _timeout from the environment.
#     INFERRED (not measured): one metadata lookup timed out or was
#       refused.  The z0 task starts 7 make_mosaic.py processes one after
#       another, each with read workers, and each process looks its
#       credentials up afresh; 74 s is most of the ~95 s the task took
#       locally, so a LATER call failed, not the first.  Which one is not
#       known: the task prints nothing per call.
#     PRECEDENT: plan_GL_north.sh N3, E-40_N-760, NoCredentialsError at
#       8.6 s, 1 job of 555, undiagnosed; its retry succeeded.
#     RUN 2 (Ben: QD9 a), the same job unchanged: a2d48bdb
#       (IS_mosaic_smoke_retry1_z0, 16gb, t3.xlarge i-04444e41015eda95c,
#       build 6663090), ledger IS_0332_mosaic_smoke_retry1_jobs.csv.
#       SUCCESSFUL, 138 s (step mosaic 104 s, steal 17.4 s), peak 0.47 GiB.
#       z0.h5 (91340373 bytes) at the test out_prefix.  GATE PASSES: vs the
#       ADE z0.h5, same datasets and shapes, 0 NaN flips in every field;
#       mask, sigma_z0, x, y identical; the weighted fields differ at
#       rounding only (z0 6.8e-13 m in 6643 cells, cell_area 5.5e-12 m^2,
#       count 7.1e-15, misfit_rms 2.2e-16, misfit_scaled_rms 3.6e-15) --
#       the same z0 figure as the local test (D3b-5).
#     STATEMENT: run 1's failure did not repeat on an unchanged job, so it
#       is intermittent.  1 failure in 2 mosaic jobs says nothing about
#       the rate.
#     QUESTION QD9 (Ben; (a) ANSWERED and done, (b) open): (a) resubmit the same job unchanged -- 1 job, no
#       rebuild; tells intermittent from systematic.  (b) in run.sh export
#       AWS_METADATA_SERVICE_NUM_ATTEMPTS=5 and
#       AWS_METADATA_SERVICE_TIMEOUT=5 for every step, and echo each
#       mosaic command before it runs -- needs a rebuild + register.
#       RECOMMENDATION: (a) now; (b) before D6/D7 whatever (a) shows, since
#       41 mosaic jobs x several processes each multiplies the exposure.
#       QD9 (b) answer (Ben 2026-10-01): NOT NOW.  MAAP admin expects 1-2%
#       of jobs to fail for assorted reasons; the submission API has an
#       auto-resubmit option, for when the on-MAAP algorithms are more
#       mature.  For the time being failed jobs are simply resubmitted;
#       where to harden is decided once there is more data.  So: no run.sh
#       change, no rebuild; every failure's cause is recorded here.
#
# D6. [DPS] DONE 2026-10-01 (result below the step).  IS end to end: 41 mosaic jobs, then ATL14 and ATL15 nc
#     jobs, to a TEST prefix (never over IS's canonical products) --
#     RECOMMENDATION: .../ATL14_processing/rel006_0332_testing/north/IS,
#     for BOTH the mosaics and the netCDFs here, since IS's canonical prefix
#     already holds ADE-made products.  GATE:
#     every mosaic identical to the ADE's; netCDFs identical in data and
#     attributes apart from dates/build fields (list the differences, don't
#     assume).  Record per-job time and peak memory.
#     RESULT 2026-10-01, all on build 6663090, queue 16gb (t3a.xlarge),
#     out_prefix .../ATL14_processing/rel006_0332_testing/north/IS; ledgers
#     and collect output in ~/ATL14_processing/maap_ledgers/
#     IS_0332_D6_{mosaic,mosaic_retry1,nc}_{jobs.csv,collect.txt}.
#       MOSAIC: 41 jobs submitted 16:36Z, all finished by 16:52Z.  40
#         successful, 1 failed (avg_dzdt_20000m_lag20, job 82cfeb69, 15 s):
#         botocore NoCredentialsError at run.sh's args-file fetch, the
#         job's FIRST S3 call, before BUILD_ID.  Resubmitted alone
#         (submit_MAAP_jobs.py --task, new ledger): successful, 49 s.
#         Job time 42-150 s, median 46 s (z0 150 s; ~2080 job-seconds in
#         all, failure and retry included); peak RSS 0.15-0.45 GiB (z0).
#       NC: ATL14 65 s (step 27 s), peak 0.83 GiB; ATL15 125 s (step 87 s,
#         steal 11.7 s), peak 0.66 GiB.  Both successful first time.
#       GATE PASSES, vs the ADE's IS products of 2026-09-25
#       (~/ATL14_processing/rel006/north/IS):
#         mosaics, 41 files, 253 datasets: same datasets and shapes, 0 NaN
#           flips; 225 identical, incl. every averaged mosaic (30 files),
#           masks and sigmas; 28 datasets in the 11 WEIGHTED files (z0, dz,
#           dzdt_lagN) differ at rounding only: z0 6.8e-13 m, dz 7.1e-15 m,
#           dzdt 3.6e-15 m/yr, cell_area 2.3e-10 m^2, count 7.1e-15,
#           misfit 3.6e-15 -- the D3b-5 figures, accepted by Ben.
#         netCDFs, 5 files: same byte sizes as the ADE's; every variable
#           identical (ATL14 28/28, ATL15 1 km 92/92, 10/20/40 km 89/89;
#           counts include coordinate/metadata datasets); attributes differ
#           ONLY in date_created, history (a timestamp),
#           identifier_file_uuid, and METADATA/DatasetIdentification
#           creationDate and uuid.
#       FAILURE DATA so far (Ben 2026-10-01: resubmit, harden later):
#         NoCredentialsError 2 of 44 mosaic-step jobs today (D5 run 1
#         mid-task; this one at the first call), 0 of 2 nc jobs; earlier 1
#         of 555 GL north prelim (plan_GL_north.sh N3).  All three passed
#         on resubmit.
#
# D7. [DPS] IN PROGRESS 2026-10-01 (Ben: go).  GL north: 32 mosaic200 jobs,
#     41 mosaic jobs, ATL14 + ATL15 nc jobs (to
#     out_prefix .../ATL14_processing/rel006_0332_testing/north/GL/, AD4).  ADE checks read through the mount (QD6):
#     check_mosaic_outputs.py --values; quick plot of h and delta_h (N7).
#     Times and memory go to plan_GL_north.sh N8 for MAAP.
#     STAGE 1 (mosaic200), 32 jobs submitted 17:18Z, build 6663090, 16gb,
#     ledger maap_ledgers/GL_0332_north_D7_mosaic200_jobs.csv: 21
#     successful (186-818 s, median 588 s, peak RSS 0.80 GiB), 11 FAILED.
#       STATEMENT (the 11 triaged _stderr.txt): every one is botocore
#         NoCredentialsError in pointCollection glob_remote, at the start
#         of a make_mosaic.py call; 22 such errors in the 11 logs (1-4 per
#         job).  A group that fails ends the job with nothing uploaded.
#       STATEMENT: a mosaic200 job runs 41 groups of ~2 make_mosaic.py
#         calls, 2 groups at a time: ~85 processes, each looking the
#         worker's credentials up afresh.
#       ESTIMATE: >= 22 failed lookups in ~2700 process starts = ~0.8% per
#         lookup; 1 - 0.992^85 = ~50% per job expected, 34% seen.  The
#         same rate fits 1 failure in 41 single-group IS mosaic jobs.
#       RETRY 1: the 11, submitted 17:42Z on 6663090 (--task, ledger
#         ..._mosaic200_retry1_jobs.csv).
#     DECIDED (Ben 2026-10-01, on this data): add the retry settings.
#       run.sh now exports AWS_METADATA_SERVICE_NUM_ATTEMPTS=5 and
#       AWS_METADATA_SERVICE_TIMEOUT=5 (each kept if already set) before
#       any step.  CHECKED on the ADE: botocore and aiobotocore sessions
#       both read 5/5 from the environment, 1/1 without.  INFERRED, not yet
#       measured: that the failed lookups are the metadata service's
#       (the container has no other credential source).  The measurement
#       is the failure count of the remaining D7 jobs on the new build.
#       NEEDS: push, Ben registers when nothing is in flight,
#       check_build_id MATCH.
#
# D8. [ADE] TODO.  Resume plan_GL_north.sh at NM2 (compare, report only),
#     NM3 (copy the nc S3 -> the QM-A side key; aws s3 cp, S3 to S3), then
#     NM4-NM6 as written, minus every fetch.  The two NM drivers are
#     superseded; a new driver, if any, only submits and waits.
#
# D9. [docs] TODO.  Transition_to_maap.md (Q4/Q18 answered: DPS),
#     run.sh header, howto_MAAP_*.sh: drop every `aws s3 sync ... prelim/`
#     into /home and every fetch_tiles.py step; point checks at the mount.
