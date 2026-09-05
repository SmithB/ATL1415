# Transition to MAAP

## QUESTIONNAIRE -- input needed before the MAAP howtos can be written
Raised 2026-09-04 while planning docs/howto_MAAP_<region>.sh (see "Porting the howto
workflows to MAAP" at the end of this file for the plan these questions come out of).
Each item has a recommendation; answer on the `A:` line and the plan follows from it.

STILL OPEN as of the 2026-09-05 pass, and only this one:
  Q4  -- A: is still "TBD", but it is NOT waiting on you: Q18 turned it into a measurement,
        and that measurement is blocked on building the ATL1415 env in the ADE (Q5).  It
        resolves itself as soon as that env exists and one z0 field group can be timed.
Everything else (Q1-Q3, Q5-Q27) is answered.  What remains in this file is work items, not
questions.

Longer entries are sectioned so the kinds of claim are distinguishable: FINDINGS and
MEASUREMENTS are verified statements with their provenance, WORK ITEMS are consequences with a
status, RECOMMEND marks a proposal of mine rather than a fact, and THE QUESTION is whatever
needed your decision.  Q27 is the worked example; say the word and the other long entries get
the same treatment.

### Scope and naming
Q1. Which howtos do you want?  Proposed: howto_MAAP_GL, howto_MAAP_AA,
    howto_MAAP_arctic (RA IS CN CS SV), and howto_MAAP_ATL11 (indexing/staging).
    Recommend all four, with ATL11 written last -- see Q13.
    A: The howto_MAAP_ATL11 is not needed b/c that will be done on discover.  All 
    others are needed

Q2. MAAP howtos alongside the SLURM ones, or replacing them?  Recommend ALONGSIDE:
    the discover runs are not going away this release, and the two schedulers diverge
    only in the submit step.  The existing four files would get a one-line header
    saying they are the discover/SLURM variants.
    A: Alongside

Q3. You asked for `.txt`; the existing four are `.sh`.  They are copy-paste command
    logs, not executable scripts, so `.txt` is arguably the more honest extension.
    Say if you would rather keep `.sh` for consistency.
    A: Keep .sh.  Otherwise when I open these in emacs they get line-wrapped in unhelpful ways.

### What runs where
Q4. Confirm the split the algorithm structure section already assumes: ONLY the
    per-tile prelim and matched solves go to DPS; setup, queue-build, 200 km tiles,
    mosaic, to-netcdf and browse all stay in the ADE.  Specifically -- is the ADE
    instance big enough (cores, RAM, walltime) to mosaic Antarctica?  On discover
    that is a 4-hour, 4-task SLURM job per field group.
    A: TBD

Q5. Is there an ADE conda env that can `import pointCollection` and `LSsurf`?  The
    Cloud-paths section records that neither is importable in the ADE notebook env.
    setup_ATL1415_region.py, make_ATL1415_queue.py, make_200km_tiles.py and the
    mosaic scripts ALL need pointCollection, so every ADE-side stage is blocked until
    one exists.  build-env.sh already knows how to build exactly this env; recommend
    running it once in the ADE to make `ATL1415` there.
    A:  Yes.  Add that build to the workflow.

### The tile grid
Q6. THE v4.1 MASKS HAVE NO 1 km VERSION, and make_ATL1415_queue.py needs one to lay
    out tile centers.  It derives the grid from a 1 km sibling of --mask_file:
      AA  AntarcticIceMask_..._240m_v4.1.tif -> neither '100m' nor '125m' in the name,
          so it raises ValueError outright
      GL  GreenlandIceMask_..._100m_v4.1.tif -> GreenlandIceMask_..._1km_v4.1.tif,
          which is not in Zenodo record 22259649 and is not on the bucket
    So no cloud tile list can be built for either ice sheet today.  Which fix?
      (a) publish 1 km decimations alongside the masks (Zenodo + bucket), keep the
          sibling-name convention;
      (b) add an explicit --grid_mask_file to make_ATL1415_queue.py and stop guessing
          the name;
      (c) freeze the tile list once per release into region_files/ and drive every
          later run from --xy_list_file.
    Recommend (b) plus (c): (b) removes a fragile string substitution, (c) makes the
    tile set reproducible and auditable, which matters more once jobs are billed.
    A: (b) + (c), and a 1-km mask should be calculated on the fly (if needed) with gdal_translate from band 1 of each file when the setup is run, then cached in the masks directory.

Q7. AA is currently submitted as two halves (--min_xy 360000 north, --max_xy 440000
    south) purely so the two could get different SLURM walltimes.  On DPS walltime is
    a property of the queue, not the submission, so the split can collapse to one
    fan-out -- unless you want the big near-pole tiles sent to a larger queue.  One
    submission, or keep two?
    A: keep two

### DPS mechanics
Q8. HOW DOES A `matched` JOB GET ITS NEIGHBOURS?  run.sh expects input/prelim/ to hold
    the tile plus its 8 neighbours (prior_edge_include reads them), but
    algorithm_config.yml declares a single `file` input and DPS localizes one URL per
    input.  Options:
      (a) give run.sh a `prelim_prefix` input and let it `aws s3 cp` the 9 keys itself;
      (b) build a per-tile tarball of the 3x3 neighbourhood in the ADE and register it
          as a second file input (thousands of tarballs; also collides with run.sh's
          `find input -name '*.txt' | head -1` if it were a .txt);
      (c) run matched in the ADE and keep only prelim in DPS.
    Recommend (a): one small change to run.sh, one config input, nothing staged.  It
    depends on Q9.
    A: (a)

Q9. DPS scatters output to dps_output/<algo>/<tag>/YYYY/MM/DD/HH/MM/SS/<us>/, which
    the mosaic step cannot address (this is already an open TBD).  Proposed
    deterministic prefix, written by the job itself:
      s3://maap-ops-workspace/ben_smith/ATL14_processing/rel<NNN>/<hemi>/<REGION>/{prelim,matched}/E%d_N%d.h5
    i.e. the same tree the ADE already builds, mirrored on the bucket, with output/
    kept for logs only.  Good, or do you want a different layout?
    A: Yes, but the directory structure will need a quarterly / monthly specifier as well

Q10. The monitor.  Proposed `check_MAAP_jobs.py <ledger.csv>` mirroring
     slurm_run_status.py's interface -- counts by state, failures grouped by exception
     parsed from the job's stderr, and `--requeue failed|timedout|both [--dry_run]`.
     The ledger would be a CSV written by the submitter: x0,y0,step,job_id,submitted_at.
     Is that the shape you want, or should job state be recovered from listJobs() and
     the tag rather than from a local file?
     A: local file sounds good.

Q11. Submission rate.  Blocked on the organizational-queue TBD at the top of this file,
     but the submitter needs a policy: max jobs in flight, sleep between submissions,
     and what to do on a submitJob error (retry with backoff vs record and continue).
     Recommend: --max_in_flight (poll and top up), --rate, and record-and-continue with
     the failure in the ledger.  Any limits you already know of?
     A:  No limits known.  For testing purposes, the current limits are not yet a problem

### Data staging
Q12. The arctic five (RA IS CN CS SV) use RGI `.db` masks plus a `_40km.tif` sibling,
     so they DODGE the Q6 problem -- but two things are unverified:
       - is masks/RGI_reduced/ actually on the bucket as real files?  The mask upload
         copied a repo tree that had never been lfs-pulled, and the .db files were added
         later (35d86f4) with the lfs filter commented out (136c224).
       - can ATL11_to_ATL15 read a `.db` mask over s3 at all?  The cloud read path
         covers geotiff (/vsis3) and h5/nc (s3fs); a sqlite .db is not obviously covered.
     Should I verify both before writing howto_MAAP_arctic, or is the arctic set out of
     scope for the first MAAP release?
     A: The arctic set is a good place for testing.  Please verify the db mask- if this fails, consider a simpler file format (i.e. geojson)
     VERIFIED 2026-09-04, both halves.  NO FORMAT CHANGE IS NEEDED -- but one line of code is.
       - masks/RGI_reduced/ IS on the bucket as REAL FILES, not lfs pointers: all 5 .db, all
         5 _40km.tif and all 5 _80km.tif are there at byte sizes identical to the repo copies
         (e.g. Iceland .db 995328 both places).  The lfs-pointer worry does not apply here.
       - GDAL READS THE .db OVER S3 FINE.  The driver is SQLite (not GeoPackage), and
         ogr.Open('/vsis3/maap-ops-workspace/ben_smith/ATL1415/masks/RGI_reduced/
         06_rgi60_Iceland_reduced.db') opens it, finds the layer and counts its 14 features,
         on gdal 3.12.3.  geojson is not needed.
       - WHAT FAILS is the raw URI: ogr.Open('s3://...') raises "No such file or directory".
         ATL1415/make_mask_from_vector.py line 22 calls ogr.Open(mask_file, 0) with the path
         exactly as given, bypassing pc.io_utils.as_gdal_path(), which is what maps
         s3:// -> /vsis3/ for every raster read.  Routing that one call through as_gdal_path()
         is the whole fix, and it is the only thing standing between the arctic five and a
         cloud run.  (from_geotif already routes correctly, so the _40km.tif the queue builder
         reads is fine as is.)

Q13. THE ATL11 INDEX FOR THE NEXT RELEASE.  make_ATL11_index.py globs a local ATL11
     archive and symlinks it; there is no such archive on MAAP, and the current
     ATL11_index_0331_007_04 tree was staged by hand.  For 0332 and after, does the
     per-granule index get built from earthaccess granules by a new script, or do you
     keep staging it by hand from discover?  howto_MAAP_ATL11 says something quite
     different depending on the answer.
     A: I keep building the ATL11 indexes on discover.

### Output differences to confirm before they are baked into a howto
Q14. CROSSOVERS ARE OFF in cloud mode: setup_ATL11_xover.py still cannot write a schema
     with a remote `source`, so MAAP_dps.txt deliberately omits --ATL11_xover_dir.  A
     MAAP run therefore does not match a discover run of the same processing string.
     Acceptable for rel006 production, or does the cloud schema have to land first?
     A: Haven't we added this to the read_ATL11 mechanics?  If not, that's a TBD
     CHECKED 2026-09-04 -- YOU ARE MOSTLY RIGHT, and the remaining gap is smaller than this
     question implied.  Three separate things were being conflated:
       1. ALONG-TRACK ATL11 over CMR+S3: done (read_ATL11_at, --ATL11_earthaccess, 295c181).
       2. THE CROSSOVER READ over S3: also done.  read_ATL11_xovers() passes a shared fs=
          through to every group it opens, and tilingSchema.resolve_files_for_box() has a full
          EarthAccess branch -- one batched search_data(granule_name=[...]) over the candidate
          tile names, returning direct S3 URLs (299ab26).  So the read mechanics ARE there.
       3. WRITING THE SCHEMA THAT POINTS AT THEM: still missing, and that is the only gap.
          setup_ATL11_xover.py writes directory=<local cycle_src> with no 'source', and it
          refuses to run unless that local directory exists -- so it cannot be run in the ADE.
     NEWLY VERIFIED, and it makes the fix worth doing rather than speculative:
       - earthaccess.search_datasets(short_name='ATL11XO') returns one collection, version 007.
         ATL11XO is genuinely published, so the EarthAccess source resolves.
       - real granule names -- ATL11XO_AA_E600_N-1400_c01_007_03.h5,
         ATL11XO_AR_E-2200_N-1400_c01_007_03.h5 -- match setup_ATL11_xover.py's format_str
         character for character, including the 007/03 that --ATL11xo_version=007_cycle_03_30_v03
         parses to.  Nothing needs renaming, converting or staging.
     So crossovers do not have to stay off.  See Q19 below.
     A: If we're doing S3 reads, then if a schema is needed it can be generated internally by a call to pointCollection and kept in memory, never written to disk.  Follow-up question: is a schema needed at all?  Second follow-up point to point out: the current ATL11xo products have filenames that point to the centers of the tiles, not the corners as described in the ATL11 ATBD.  
     ANSWERED 2026-09-04, by reading the real product.  See Q24/Q25/Q26 below; in short:
     YOU ARE RIGHT ON BOTH COUNTS, setup_ATL11_xover.py is NOT needed on MAAP, and reading
     the product turned up a third thing that matters more than either.
       - CENTERS CONFIRMED, not corners.  ATL11XO_AR_E1200_N-2400_c05_007_03.h5 holds
         x in [1100018, 1297561] and y in [-2499800, -2300013] -- i.e. the tile spans
         1100-1300 km by -2500..-2300 km and the label (1200, -2400) is its CENTER.  So the
         ATBD is wrong about the shipped product, and mapping_function_name='round' in
         setup_ATL11_xover.py is right.  Recorded here because the next person to read the
         ATBD will get this backwards, and getting it backwards silently returns the wrong
         tiles rather than failing.
       - NO SCHEMA FILE IS NEEDED IN THE CLOUD.  Everything tilingSchema contributes is a
         box -> list of names mapping; nothing about it has to be persisted.  Two ways to
         drop the file, both viable -- Q24.
       - AND THE THING NEITHER OF US ASKED: the product has 31 cycles per tile and the code
         reads 2.  Q26.  This one is not a MAAP question at all.

Q15. AA now actually gets --tide_adjustment (the bare-flag bug fix), and AA_0331.txt now
     names CATS2008-v2023 rather than CATS2008.  Both are real output changes that a
     howto would bake in.  Confirm both before I write howto_MAAP_AA.
     A: Confirmed.

### Follow-ups, raised 2026-09-04 after the answers above

Q16. THE ON-THE-FLY 1 km MASK (your Q6 answer) needs three details settled before it is code:
     (a) WHERE DOES THE CACHE GO?  "cached in the masks directory" works on discover, but on
         MAAP --mask_dir is s3://maap-ops-workspace/ben_smith/ATL1415/masks/ and the queue
         builder runs in the ADE.  Write the cached 1 km tif back to the bucket (so it is
         built once ever and shared), or to a local cache dir, or next to the composed args
         file under <ATL14_root>/rel<NNN>/<hemi>/<REGION>/?  Recommend writing it back to the
         bucket beside the source mask and keying it on the source's ETag, so a re-published
         mask invalidates it.
     (b) HOW IS IT RESAMPLED?  gdal_translate's default is nearest neighbour, which on a
         240 m ice mask decimated 4.2x will drop thin margins and outlet glaciers -- exactly
         the tiles you least want to lose.  For a 0/1 mask the meaningful rule is "any ice in
         the cell": -r max (gdal >= 3.3), or -r average with a >0 threshold.  Recommend -r max.
     (c) WHICH BAND?  You said band 1, which is the EARLIEST time slice (AA 2018.200,
         GL 2018.00), so the tile layout would follow start-of-record extent.  Since the tile
         grid is dilated by ~1.5 tile widths anyway that is probably harmless, but a max over
         all 40 bands costs nothing at 1 km and removes the question.  Band 1, or all-band max?
     A: (a): write the cache back to the bucket.  (b): use max.  (c): Use band 1.  The dilation renders the question moot.


Q17. QUARTERLY/MONTHLY IN THE S3 LAYOUT (your Q9 answer).  On discover this is already
     expressed as a hemisphere suffix -- --hemi_suffix=_monthly in monthly.txt gives
     rel006/north_monthly/GL/ -- and setup_ATL1415_region.py builds exactly that tree.
     Mirroring it needs nothing new:
       s3://maap-ops-workspace/ben_smith/ATL14_processing/rel006/north_monthly/GL/prelim/
     Is reusing hemi_suffix what you meant, or do you want an explicit separate segment
     (.../rel006/north/GL/monthly/prelim/)?  Recommend reusing hemi_suffix -- the ADE and the
     bucket then have literally the same tree, which is what makes the gather step trivial.
     A: reuse the hemi suffix

Q18. Q4 IS STILL "TBD" AND IT GATES A WHOLE STAGE.  If the ADE cannot mosaic Antarctica, the
     mosaic step is not an ADE step and the plan changes shape.  Want me to settle it by
     measurement -- report the ADE instance's cores/RAM/disk, then time one field-group mosaic
     for a small region (IS or GL z0) and extrapolate to AA?  That turns a TBD into a number
     without committing to anything.
     A: Yes.  test.
     HARDWARE, measured in the ADE 2026-09-04: 16 cores (Xeon Platinum 8259CL @ 2.50 GHz),
     124.4 GiB RAM (118.8 available), 782 GiB free on /home/jovyan and 177 GiB on /.
     That is FOUR TIMES the tasks a discover mosaic job asks for (4) and no walltime limit at
     all, so the shape of the answer is already encouraging -- an AA mosaic looks like a
     long-running ADE job rather than an impossible one.
     The TIMING half still has to wait on Q5's env: pointCollection is not importable in the
     ADE notebook env today (confirmed -- `import pointCollection` raises ModuleNotFoundError),
     so no mosaic can actually be run yet.  Sequence: build-env.sh in the ADE, then time one
     z0 field group for IS, then extrapolate by tile count to AA.

Q19. GIVEN THE Q14 FINDINGS, do crossovers go back on for MAAP before the first production
     run?  The change is: setup_ATL11_xover.py gains a cloud mode (source=EarthAccess,
     directory=None, no local-directory check), MAAP_dps.txt sets --ATL11_xover_dir, and
     setup_ATL1415_region.py stops suppressing it in cloud mode.  Recommend yes -- otherwise
     every MAAP tile differs from its discover twin for a reason that is now cheap to remove.
     If yes: where should the two cycle_NN/200km_tiling_{AA,AR}.json files live?  They are a
     few hundred bytes each and tilingSchema.from_file() reads them over s3, so
     s3://maap-ops-workspace/ben_smith/ATL11_xover_schema/<release>/ would work.
     A:  Crossovers should be included.  There might be a code change as described in the answer to Q14

Q20. WITH NO howto_MAAP_ATL11 (your Q1/Q13 answers), THE discover -> S3 HANDOFF IS NOWHERE.
     You build the index on discover; the tile jobs read it from
     s3://maap-ops-workspace/ben_smith/ATL11_index/ATL11_index_0331_007_04/ (8100 files).
     Nothing documents or scripts the copy between those two facts, and it has to happen once
     per release.  Where should it live -- a short "staging" section at the top of each MAAP
     howto, or one howto_MAAP_staging.sh that all three reference?  And was the 0331 upload a
     hand-rolled aws s3 cp, or is there a script to point at?
     A: There should be a howto_MAAP_staging.sh that descibes how to setup a fresh maap account with everything that ATL1415 needs.
     NOTED, and that is wider than what I asked -- good.  Scope taken as: build the ADE env
     (Q5), stage the masks from Zenodo, stage the ATL11 index built on discover, register the
     DPS algorithm, and verify each of those is readable.  Tide models need no staging
     (s3://pytmd) and, per Q24, neither do crossover schemas.

Q21. THE TEST PLAN, given the ~10 jobs/hr public-queue throttle.  Your Q11 answer says current
     limits are not yet a problem, which is true only for a small fan-out: a full GL prelim run
     is well into the thousands of tiles, i.e. weeks at that rate.  Suggest the first MAAP runs
     be explicitly bounded -- Iceland (the smallest arctic region) end to end, then a
     --bounds-limited GL patch -- rather than a whole region, and that the howtos say so.
     Agree, and is Iceland the right first target?
     A: Iceland

Q22. AA STAYS SPLIT IN TWO (your Q7 answer) -- but on DPS the split only buys something if the
     two halves go to different queues, since walltime is a queue property.  Should the near-pole
     half target maap-dps-worker-32vcpu-64gb (or whatever the org queue turns out to be) and the
     outer half a smaller one?  That is worth naming in the queue request, which is still open at
     the top of this file.
     A:  The south partition should target a larger instance than the north partition.  It may need to be established in testing how large an instance is needed for each partition.

Q23. NEW, FOUND WHILE AUDITING: default_args/MAAP.txt (the ADE variant) still says
     --tide_directory=/home/jovyan/my-private-bucket/tide_models, and CATS2008-v2023 is NOT
     staged there -- only the old CATS2008, which is itself missing hf.CATS2008.out.  So an
     ADE-side local run of AA_0331 composed with MAAP.txt fails on a missing tide model, while
     the DPS run composed with MAAP_dps.txt succeeds.  Point MAAP.txt at s3://pytmd as well
     (recommended -- ATL1415.tides handles it identically and nothing needs staging), or stage
     the TMD3 netCDF to the mount?
     A: Point MAAP.txt at s3://pytmd 

Q24. IS A SCHEMA NEEDED AT ALL?  No -- not as a FILE, and on MAAP not as a setup step either.
     There are two ways to get there and they differ more than they look.
       OPTION A -- BUILD THE SCHEMA IN MEMORY (what you proposed).  In cloud mode
       read_ATL11_xovers constructs, per cycle,
         pc.tilingSchema(tile_spacing=200e3, mapping_function_name='round', scale=1000,
                         extension='.h5', format_str=f'ATL11XO_{hemi}_E%d_N%d_c{cc}_{rel}_{ver}',
                         source={'type':'EarthAccess','short_name':'ATL11XO'})
       and calls the SAME resolve_files_for_box() that already works.  ~8 lines, no file, no
       setup_ATL11_xover.py on MAAP, no --ATL11_xover_dir.  It depends on the granule-name
       convention -- which is now verified exact, including the centers question above.
       OPTION B -- DROP THE SCHEMA CONCEPT AND SEARCH CMR SPATIALLY, mirroring read_ATL11_at:
       earthaccess.search_data(short_name='ATL11XO', bounding_box=_lonlat_bounding_box(...)).
       MEASURED: over an Iceland box that returns 13 distinct tiles, all AR, all on the 200 km
       grid -- precise, not over-broad, so this genuinely works.  It is immune to any future
       filename change, and it reuses the antimeridian handling _lonlat_bounding_box() already
       has.  Its cost: it returns EVERY cycle (368 granules for those 13 tiles), so cycle and
       version filtering has to happen client-side before anything is opened.
     RECOMMEND OPTION A.  It is smaller, it reuses a code path that is already tested, and the
     one risk that argued for B -- trusting the filename convention -- is the thing I just
     verified against the shipped product.  Which do you want?
     A: schema in memory.  (Option A, confirmed 2026-09-05.) Build the schema in memory.

Q25. OPTION A NEEDS RELEASE/VERSION TO REACH THE SOLVE, and there is a free way to do it.
     read_ATL11_xovers currently learns '007'/'03' only from the format_str that
     setup_ATL11_xover.py baked into the schema file.  Build the schema in memory and that
     information has to arrive some other way.
     THE GOOD NEWS: --ATL11xo_version=007_cycle_03_30_v03 is ALREADY in every composed args
     file (it comes from rel_006_0331.txt), and ATL11_to_ATL15 parses with parse_known_args(),
     so today it is silently swallowed as an unknown argument.  Promoting it to a real
     add_argument() costs nothing and changes no defaults file.  setup_ATL11_xover.py already
     has the regex that turns it into (release, version) = ('007','03').
     Also needed: a switch for "read crossovers from the cloud".  Today the switch IS
     --ATL11_xover_dir (read_ATL11 returns early when it is None), which in cloud mode has
     nothing to point at.  Cleanest is to let --ATL11_earthaccess cover crossovers too, with
     --ATL11_xover_dir remaining the local-mode switch.  Agreed?
     A:  Agreed.

Q27. WHERE DO THE PREVIOUS ATL14/15 PRODUCTS COME FROM?  (your question, 2026-09-04)
     SHORT ANSWER: yes, CMR + S3 works, nothing has to be staged, and the numbers are good --
     but there are five code changes, one of which is a silent-failure bug.

     HOW TO READ THIS ENTRY.  Most of it is findings, not questions.  FINDINGS and
     MEASUREMENTS are verified statements, with how each was verified.  WORK ITEMS are
     consequences of those findings, each with a status; none of them needs your input.
     Anything that is my proposal rather than a fact is marked RECOMMEND.  The single thing
     that needed a decision is under THE QUESTION, which you have answered.

     ---- THE QUESTION (answered) ----
     How should the release/version to search for reach the solve?
       (a) reinterpret --previous_product under a flag, the way --ATL11_index is reinterpreted
           under --ATL11_earthaccess -- it becomes "the release to search for" rather than a
           directory; or
       (b) add a separate --previous_product_earthaccess plus --previous_product_release.
     RECOMMEND (a), for symmetry with the ATL11 path.
     A: reinterpret --previous_product under a flag the way --ATL11_index was
        reinterpreted under --ATL11_earthaccess

     ---- FINDINGS (statements; checked against CMR or the code, 2026-09-04) ----
     F1. THE PUBLISHED PRODUCT IS THE PREVIOUS PRODUCT.  ATL14 short_name='ATL14' version 005
         has 11 granules; ATL15 has 40 (4 averaging scales x 10 regions).  Both are exactly the
         0329 products that --previous_product_top=/discover/.../rel005_0329/ already names.
         CONSEQUENCE: no copy to our bucket is needed, ever.
     F2. THE COLLECTION MIXES CYCLE RANGES.  ATL14_A1_0328_100m_005_01.nc sits alongside
         ATL14_A1_0329_100m_005_02.nc under the same short_name.  CONSEQUENCE: a search must
         filter on the cycles string AND version, not just short_name, or 0328 and 0329 get
         mosaicked together -- wrong data, no error.
     F3. THE FILES ARE BIG: ATL14_GL 1.4 GiB, ATL14_A1 3.4 GiB, ATL15_A1_01km 2.9 GiB.  All of
         AA (A1-A4, ATL14 + ATL15 1 km) is about 20 GiB.  This is why the read path matters.
     F4. CMR HAS v005 ONLY TODAY.  When rel006 publishes, rel007's "previous product" becomes
         v006.  CONSEQUENCE: release/version must be an argument, not a constant.

     ---- MEASUREMENTS (a 60 km GL tile, read straight from the NSIDC bucket) ----
       ATL14  h        chunks (2491,1401) float32 gzip -> 721x721 window,  5.7-35 MiB, 0.3-1.3 s
       ATL15  delta_h  chunks (8,686,386)  gzip       -> 29x73x73 window, 16.2 MiB,   2.2 s
     i.e. roughly 22 MiB and a couple of seconds per tile for GL, against 2.9 GiB of files;
     about 90 MiB for AA across four sectors.  Entirely practical.
     THE 5.7-35 MiB RANGE IS THE s3fs BLOCK SIZE, and it is worth setting explicitly: 5 MiB
     blocks (the default) fetch 35 MiB, 256 KiB blocks fetch 6.6 MiB, and an untuned first
     attempt pulled 320 MiB.  Same lesson as the tide stores -- the access pattern costs more
     than the data.

     ---- WORK ITEMS (consequences, with status; not questions) ----
     W1 [OPEN -- THE ONE THAT MATTERS MOST].  THE SILENT-SKIP BUG.
        set_three_sigma_edit_from_previous_product finds its files with
        glob.glob(os.path.join(directory,'ATL14_*.nc')).  glob over an s3:// URI returns [],
        so the function finds no coverage, prints "no previous-product coverage for tile,
        skipping", and RETURNS 0.  Point --previous_product at a URI today and the three-sigma
        pre-edit is silently disabled on every tile, with no error and a message that looks
        like a normal edge-of-domain outcome.  --previous_product is already typed
        path_or_uri, so nothing stops someone doing exactly that.
        RECOMMEND: it should fail loudly.
     W2 [DONE UPSTREAM 2026-09-05].  THE READ PATH.  pointCollection PR #53 (fa66500, merged
        as 96ea0c3).  from_nc() on a URI used to read the WHOLE OBJECT into memory, so a
        single AA tile job would have pulled ~20 GiB.  It now opens the granule through h5py
        and range requests -- new module grid/nc_h5.py presents an h5py.File through the
        netCDF4 API from_nc() uses, so from_nc()'s body is unchanged -- and a windowed read
        costs only the chunks the window touches.  nc_open()/from_nc() gained
        engine=('auto','h5py','netcdf4'), block_size and rdcc_nbytes; NETCDF3/classic falls
        back to netCDF4, and local uncompressed reads take exactly the path they did before.
        Nothing in ATL1415 has to change to get this -- pyproject.toml tracks pointCollection
        main unpinned, so the next DPS build picks it up.  One thing to know: netCDF4 silently
        applies scale_factor/add_offset and h5py does not, so the adapter replicates that
        arithmetic.
     W3 [OPEN].  DISCOVERY.  The current code globs A1-A4 and mosaics them by filling NaN
        left-to-right.  RECOMMEND replacing the glob with an earthaccess search by bounding
        box, filtered per F2 to the right cycles/version.  That is a simplification as well as
        a cloud fix: a spatial search returns only the sectors that actually intersect the tile.
     W4 [OPEN].  THE SETUP SIDE.  setup_ATL1415_region.py resolves --previous_product_top by
        local glob (glob.glob(<top>/south/A?) plus os.path.isdir), which finds nothing in the
        ADE.  RECOMMEND: in cloud mode emit a cloud marker plus the release/version to search
        for, rather than a list of directories that do not exist.
     W5 [OPEN, NOT ON THE QUARTERLY CRITICAL PATH].  --ATL14_reference_file has the SAME shape
        of problem by a different route: line 682 calls glob.glob(ATL14_reference_file).  It is
        used only for monthly runs, but it needs the same treatment as W1/W3.

Q26. NOT A MAAP QUESTION, AND PROBABLY THE MOST IMPORTANT THING IN THIS FILE: THE CROSSOVER
     READ COVERS 2 CYCLES OUT OF 31, ON DISCOVER AS WELL AS ON MAAP.
     - read_ATL11_xovers(..., xover_cycles=[1,2]) is a default that ATL11_to_ATL15 never
       overrides -- nothing in the repo passes xover_cycles.  setup_ATL11_xover(cycles=['01','02'])
       matches it, and has since the file was first written (2ae71a4, 2026-07-03); it has never
       been revisited.
     - the shipped product has cycles 1..31 per tile.  Over the Iceland box: 13 tiles x 31
       cycles = 368 granules.
     - the cycle label IS the ICESat-2 mission cycle, confirmed from the granules' own
       ancillary_data: c01 covers 2018-10-14 to 2018-12-27, c05 covers 2019-09-26 to
       2019-12-25.
     So every fit spanning -t=2018.75,2026.5 is currently constrained by crossovers from
     October 2018 to about March 2019 only -- the first ~6 months of an 8-year record -- while
     select_best_xover_index() is written to pick the best measurement per cycle across
     however many cycles it is given.
     IS THAT DELIBERATE?  If it is, it needs a comment saying why, because it reads as an
     early-mission default that outlived its data.  If it is not, then discover runs have the
     same gap and this wants fixing before rel006 production rather than after -- and it is a
     bigger deal than anything else on this list, since it changes the fit everywhere, not
     just on MAAP.
     A: The only data present in cycles 1 and 2 is crossover data.  ATL11 along-track begins in cycle 3.  
     That is why only cycles 1 and 2 are read from the crossovers.  There are no plans to include crossover 
     data from later cycles.

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
[ ] Smoke-test one sandbox DPS job.  THE BUILD IS NO LONGER THE OPEN PART -- the algorithm
    is registered and describeAlgorithm answers 200 (see Registration log).  What remains is
    purely a run-time question, and listJobs() is still empty.
    ANSWERED by the successful build / registration, do not re-ask:
    - the build container DOES reach github.com (the git+ pointCollection/LSsurf/PySPQR deps
      installed)
    - sparseqr DOES compile against the conda suitesparse there (once environment.yml
      supplied a toolchain -- f9543fb)
    - input order is positionals then file: (x0, y0, step, args_file, queue_name)
    THE TILE, chosen by Ben 2026-09-05: ICELAND, x0,y0 = 1260, -2620 km (43K points, likely
    on the ice-sheet edge, so a real solve rather than a trivial or saturated one).  IS is
    also the test region for the whole workflow (Q21).  On the grid: centers are odd
    multiples of tile_spacing/2, and 1260 = 63 x 20 km, -2620 = -131 x 20 km.  It appears in
    howto_MAAP_staging.sh S7 and, as the mask check, howto_MAAP_arctic.sh step 2.
    STILL OPEN, and only a real job settles them:
    - is ~/.netrc really bind-mounted into the worker (earthaccess auth depends on it)?
    - will a `file` input accept an s3://maap-ops-workspace/ben_smith/... URL?
    - do the s3:// reads in the composed args file work on the worker's own AWS credentials?
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

2026-09-04, ~13 h later: THE BUILD FAILED.  ATL1415_tile_solve is still absent from
  listAlgorithms() (200, 649 algorithms) and describeAlgorithm still returns the same 500.
  At that age it is a failure, not a slow build.  The log is the only diagnostic and it is
  browser-only -- repo.maap-project.org/api/v4 answers 401/404 unauthenticated and the raw
  job URL 302s to users/sign_in, so no script can fetch it.
  Cheap hypotheses ELIMINATED locally, so do not re-check them when reading the log:
    - the on_s3 branch IS pushed, at exactly the registered commit (ab37947)
    - all four git+ dep repos (pointCollection, LSsurf, PySPQR, ATL1415) answer 200 to an
      unauthenticated git-upload-pack probe, so they are public
    - `conda env update --name X --file Y --prune` DOES create a missing env (tested on
      conda 26.3.2), so build-env.sh line 27 is not the failure
    - build-env.sh and run.sh are mode 100755 in git, environment.yml does list suitesparse,
      and ATL11_to_ATL15.py is a real [project.scripts] entry point
  So look at the log for the SuiteSparseQR.hpp check first, then pip/network.

2026-09-04, from the build log: CAUSE FOUND AND FIXED -- IT WAS NOT SuiteSparse.  conda and
  the SuiteSparseQR.hpp check both PASSED; the build died later, at `pip install .`, with
  "No such file or directory: 'gcc'" / "'g++'" from LSsurf (Cython), sparseqr (CFFI) and
  cartopy (C++).  maap_base ships no toolchain and environment.yml asked for cython without
  compilers.  Fixed in f9543fb: c-compiler/cxx-compiler added (with a build-env.sh guard next
  to the SuiteSparse one), python pinned to 3.13 to match the ADE (unpinned the solver took
  3.14.7, which has no cartopy wheel -- hence the C++ compile), numpy pinned to 2.4, cartopy
  taken from conda, matplotlib -> matplotlib-base (dropping ~140 MB of Qt).
  Then 89ff11a: the build verification caught LSsurf and ATL1415 failing on
  ModuleNotFoundError: threadpoolctl (LSsurf/smooth_fit.py imports it and LSsurf's
  requirements.txt does not declare it -- that would have failed on EVERY tile job, not just
  at import; the real fix belongs upstream in SmithB/LSsurf).  pyTMD and timescale were
  likewise undeclared.
  Re-registered after those fixes.

2026-09-04, confirmed against the API: THE BUILD SUCCEEDED AND THE ALGORITHM IS LIVE.
  - listAlgorithms() -> HTTP 200, 650 algorithms, including
    {'type': 'ATL1415_tile_solve', 'version': 'on_s3'}.  (It was 649 and absent when the
    first build failed, so the count moving is the build landing.)
  - describeAlgorithm('ATL1415_tile_solve:on_s3') -> HTTP 200 with a real WPS
    ProcessOffering, no longer the 500 that "still building / never registered" produced.
  So the conda+SuiteSparse+compiler build works in the DPS container, sparseqr and LSsurf
  compile there, and the build container does reach github.com for the git+ deps.  Three of
  the smoke-test questions are answered by that alone.
  TWO THINGS THE REGISTERED SIGNATURE SETTLES:
  - the input order is (x0, y0, step, args_file, queue_name) -- declaration order survived
    registration, positionals then file.  run.sh's defensive skip of a leading non-numeric
    argument is still worth keeping until a real job shows what the wrapper builds for argv,
    but it is no longer the coin-flip it was.
  - `queue` in algorithm_config.yml became a queue_name INPUT carrying that value as its
    default, so a submitJob can override the queue per tile without re-registering.  That is
    what makes Q22 (AA's two halves to two different queues) free.
  STILL UNRUN: listJobs() returns [], so no tile has ever been solved on DPS.  Everything
  about runtime, memory, walltime and cost remains unmeasured, and the smoke test is still
  the next thing to do -- but it is now a job-submission question, not a build question.

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
[ X ] Chunk-wise remote reads for netCDF4 rasters (the ATL14/ATL15 previous-product path)
    - pointCollection PR #53, fa66500 -> main at 96ea0c3, 2026-09-05.  See Q27 W2.
    - full suite green on the merge: 265 passed, 3 skipped (the S3 tests are network-gated).
    - unpinned git+ dependency, so no ATL1415 change is needed to consume it.

[ ] setup_ATL11_xover must also be able to emit a CLOUD schema.  THIS IS THE ONLY THING
    KEEPING CROSSOVERS OFF IN CLOUD MODE, and it is a small change -- verified 2026-09-04:
    - the READ side is complete.  tilingSchema.resolve_files_for_box() has an EarthAccess
      branch that batches one earthaccess.search_data(granule_name=[...]) over the candidate
      tile names and returns direct S3 URLs, and read_ATL11_xovers() already plumbs the shared
      fs= through pc.io_utils.open_h5 for every group it reads (299ab26).
    - ATL11XO IS A REAL CMR COLLECTION: earthaccess.search_datasets(short_name='ATL11XO')
      returns exactly one, version 007.  This was NOT previously confirmed and it is what makes
      the EarthAccess source viable rather than hypothetical.
    - THE GRANULE NAMES MATCH THE SCHEMA EXACTLY.  Sampled granules look like
      ATL11XO_AA_E600_N-1400_c01_007_03.h5 / ATL11XO_AR_E-2200_N-1400_c01_007_03.h5, which is
      setup_ATL11_xover.py's format_str f"ATL11XO_{AA|AR}_E%d_N%d_c{cycle}_{release}_{version}"
      plus extension='.h5', with release/version 007/03 -- i.e. exactly what
      --ATL11xo_version=007_cycle_03_30_v03 in rel_006_0331.txt parses to.  No renaming, no
      staging, and no local xover tree is needed on the cloud side.
    - what is missing: setup_ATL11_xover.py writes directory=<local cycle_src> and no 'source',
      and it requires the local cycle_src directory to EXIST before it will write anything
      (it raises otherwise), so it cannot be run in the ADE at all today.  A cloud mode wants
      source={'type':'EarthAccess','short_name':'ATL11XO'}, directory=None, and no local check.
    - then MAAP_dps.txt can set --ATL11_xover_dir at wherever the two cycle_NN/200km_tiling_*.json
      files live (tilingSchema.from_file reads a remote schema, so an s3:// prefix is fine), and
      setup_ATL1415_region.py's cloud-mode branch must stop suppressing the option.
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
    - [ X ] the abspath/open() caveat is FIXED (2026-09-04): the code now accepts s3:// URIs
      for these.  See the "Cloud paths" section below for what changed and what was verified.
[ ] Create a processing string for ATL11_to_ATL15 to read previous versions of ATL14/15 from
    the cloud.  SCOPED 2026-09-04 -- see Q27 in the questionnaire.  ATL14/ATL15 v005 are
    published on CMR and ARE the 0329 previous product, so nothing needs staging; a 60 km GL
    tile costs ~22 MiB and a couple of seconds to read windowed.  Four code changes, of which
    one is a silent-failure bug: glob.glob() over an s3:// URI returns [], so an s3
    --previous_product disables the three-sigma pre-edit without any error.

## Cloud paths (s3://)
Decided 2026-09-04: rather than localizing the static data as DPS `file` inputs, the CODE
learned to read s3:// URIs.  A DPS worker has no ~/my-private-bucket mount, so masks, ancillary
grids and the ATL11 geoIndex are now named by URI in the args file and read straight from the
bucket.  Credentials come from the worker's ordinary AWS chain (~/.aws is bind-mounted), which
is deliberately NOT the earthaccess/NSIDC session: those buckets are ours, not a DAAC's.

TIDES ARE NO LONGER AN EXCEPTION (superseded 2026-09-04, commit 2535dfb).  This paragraph used
to say that pyTMD's plain file I/O forced --tide_directory to be a local path and that run.sh
had to stage ~1-1.5 GiB of OTIS model per tile job.  That is all gone: pyTMD publishes chunked
zarr stores at s3://pytmd, ATL1415/tides.py reads them by range request (1.6 MiB for an AA
tile, 6.3 MiB for GL), MAAP_dps.txt sets --tide_directory=s3://pytmd, and run.sh's staging
block was deleted.  Nothing is staged and nothing needs baking into the image.  See "Tide
models in the cloud" below for the measurements and the caveats that DO remain.

pointCollection (branch s3_static_paths, off origin/main):
    - io_utils.get_s3fs(daac=None) returns a plain s3fs.S3FileSystem on the default AWS
      credential chain, alongside the existing earthaccess DAAC sessions.  Everything reading
      OUR data defaults to it.
    - io_utils.as_gdal_path() maps s3:// -> /vsis3/, gs:// -> /vsigs/, http(s):// -> /vsicurl/;
      local and already-/vsi paths pass through.  grid.data.from_geotif() routes through it.
    - grid.data.h5_open() opens a remote file through s3fs and hands the object to
      h5py/bz2/gzip.  grid.data.nc_open() DID have to read the whole object into memory and use
      netCDF4's memory=; as of PR #53 (96ea0c3, 2026-09-05) it goes through h5py and range
      requests instead, so remote netCDF4 rasters are now windowed like .h5 rather than
      downloaded whole.  External whole-file compression (a gzip/bzip2 .nc.gz) still
      decompresses whole -- that is inherent, not a gap.  Both take fs=, plumbed through
      from_h5/from_nc, and from_nc also takes engine/block_size/rdcc_nbytes.
    - tilingSchema.from_file() reads a remote .json/.h5 schema; its directory default then
      correctly resolves to the remote directory holding the tiles.
    - geoIndex.from_file() takes fs= and opens a remote index.
    - query_ATL11_cloud: TWO BUGS FIXED.  (1) the index existence check was os.path.isfile(),
      which is False for ANY URI -- so with an s3 index root every granule would have been
      warned-and-skipped and every tile would have come back EMPTY rather than failing.
      (2) the index needs different credentials from the granule, so there is now a separate
      index_fs= parameter; the granule keeps the earthaccess session in fs=.

ATL1415:
    - new ATL1415/paths.py: path_or_uri() (the argparse type -- abspath() would rewrite
      's3://bucket/key' as '<cwd>/s3:/bucket/key'), join_path_or_uri() (os.path.join treats a
      URI as relative, since it does not start with '/'), and exists().
    - ATL11_to_ATL15: twelve input path arguments now use path_or_uri --- ATL11_index,
      ATL14_reference_file, E_d2z0dx2_file, DEM_file, mask_file, rock_mask_file, geoid_file,
      tide_mask_file, tide_adjustment_file, data_file, calc_error_file, previous_product.
      base_directory, out_name and map_dir stay local (they are outputs), and so does
      tide_directory (pyTMD).
    - setup_ATL1415_region.py, three fixes:
      * BARE FLAGS WERE BEING SILENTLY DROPPED.  The key=value regex does not match a line like
        '--tide_adjustment', so any store_true option asked for in a defaults file never reached
        the composed args file and the run proceeded with it OFF.  This bit --tide_adjustment in
        AA_0331.txt, and would have made --ATL11_earthaccess impossible to compose.  NOTE this
        means AA runs now actually get --tide_adjustment, which they did not before: that is a
        REAL OUTPUT CHANGE, and wants a look before it goes to production.
      * the key is now captured non-greedily, so only the first '=' splits key from value, and
        both sides are stripped.
      * --ATL11_earthaccess makes --ATL11_index the ROOT of the per-granule geoIndex tree, i.e.
        a directory.  The old os.path.isfile() check rejected that outright, so no cloud args
        file could be composed at all; the index is also no longer synthesized from
        --ATL11_release, and --ATL11_xover_dir is no longer derived, in that mode.
    - run.sh: --base_directory now goes AFTER @argsfile.  setup writes '-b=<region_dir>' as the
      composed file's last line and -b is the same dest, so argparse's last-wins rule meant the
      worker would have written its tiles to an ADE path that does not exist on it.

Verified against the real bucket from the ADE (not simulated):
    - windowed, band-selected read of the 23 MiB Antarctic v4.1 ice mask through /vsis3
      (250x250x1 over Ross, t=2020.2, all ice) in ~0.7 s
    - the tide-adjustment grid (130 MiB) and the tide mask, likewise windowed
    - nc_open() over s3fs on EGM2008_geoid_h.nc
    - geoIndex.from_file() on a per-granule index in ATL11_index_0331_007_04/ (34 bins)
    - every s3:// path in a composed AA args file resolves, and none is an lfs pointer
Not verified: the ATL11 granule read itself (needs the earthaccess session, unchanged here),
and anything requiring pyTMD or sparseqr, neither of which installs in the ADE env.
Tests: pointCollection tests/test_cloud_paths.py (11) and ATL1415 tests/test_paths.py (7) +
tests/test_setup_region.py (8); suites are 211 and 18 green.  Four of the setup_region tests
fail against the pre-change code, which is the point of them.

## Tide models in the cloud -- measured, and the plan
ANSWER: use the zarr stores pyTMD already publishes on its public bucket, read with s3fs, and
predict through the .tmd accessors -- the "Predict tides from model hosted on s3" recipe in
https://pytmd.readthedocs.io/en/latest/user_guide/Cloud-Access.html .  s3://pytmd/ already
holds CATS2008-v2023.zarr, Gr1km-v2.zarr and FES2022.zarr, publicly readable (anon=True), so
BOTH hemispheres are covered and nothing has to be converted or staged by us.

Measured on one 50k-point 60 km Ross tile, against pyTMD v3.0.9:

    route                                        time      bytes pulled
    netCDF over https, compute, chunks='auto'    27.5 s    1809.8 MiB   (362 range requests)
    netCDF over https, compute, chunks=None      14.3 s     582.4 MiB   (7 requests)
    netCDF over https, compute, crop=True        20.9 s     605.2 MiB
    ZARR ON S3, .tmd accessors                    2.6 s       1.56 MiB  (19 object reads)

    Gr1km-v2 zarr, SE Greenland tile              1.8 s       6.32 MiB  (17 object reads)

That is a ~370x reduction in bytes and ~6x in time against the best netCDF configuration.
The two routes AGREE: over the same 50000 points, max|A-B| = 7.2e-05 m and rms = 3.5e-05 m on
tides spanning +/-1.07 m.  The residual is consistent with the netCDF storing hRe/hIm as scaled
int16 while the zarr store holds complex64 -- i.e. the zarr is the more precise of the two, and
the difference is four orders of magnitude below --sigma_geo (4.5 m).  Treat that comparison as
the acceptance test; it has been run.

Why the netCDF route is so expensive, for the record: every variable in CATS2008_v2023.nc is
CONTIGUOUS AND UNCOMPRESSED (h5py reports chunks=None for all of them), so no spatial subset can
avoid reading through a variable's full extent -- 582 MiB = hRe (256.9) + hIm (256.9) + coords
and mask, exactly.  chunks='auto' makes it worse, not better: dask issues overlapping ranged
reads totalling MORE than the whole 1695.7 MiB file, and it fails outright when combined with
either extrapolate=True or crop=True ("inconsistent chunks along dimension y").  The zarr
stores are chunked (CATS2008-v2023 254x416, Gr1km-v2 442x274), which is the whole difference.
Also worth knowing: reaching the netCDF at all needs two pyTMD one-line fixes -- open_tmd3_dataset
hands a pyTMD.utilities.URL straight to xr.open_dataset (no inferable backend, needs
fsspec-wrapping and an explicit engine), and the next line's pathlib.Path(input_file).name
raises on the same URL.  The zarr route needs NEITHER; it does not go through model.open_dataset.

MODEL NAMES: Greenland is Gr1km-v2 (NOT greenlandTMD_v2, which is the directory and is not a
model name pyTMD accepts).  default_args/GL_0321.txt, GL_0329.txt and GL_0331.txt already say
Gr1km-v2 and need no change.  Antarctica currently says CATS2008; the zarr store is for
CATS2008-v2023, so pointing AA at the cloud model is a MODEL VERSION CHANGE and wants a
deliberate decision, not a silent rename.

DONE (2026-09-04) -- the cloud tide path is implemented and wired up:
  [ X ] ATL1415/tides.py.  tide_elevations() dispatches on --tide_directory: a local directory
      goes to pyTMD.compute.tide_elevations with exactly the arguments it always got, so local
      and SLURM runs are untouched; an s3:// prefix opens s3://<bucket>/<model>.zarr and runs
      the same post-open sequence compute does (coords_as -> interp -> predict + infer).  No
      tide computation is reimplemented -- only the dataset is substituted.
      - THE STORE IS READ ANONYMOUSLY, and that is not cosmetic: pyTMD's bucket is public but
        REJECTS SIGNED REQUESTS from other accounts.  s3fs.S3FileSystem() with the ADE's own
        credentials gets AccessDenied; anon=True works.  A DPS worker has its own credentials,
        so it would have hit exactly this.  Override with ATL1415_TIDE_ANON=0 for a private store.
  [ X ] ATL11_to_ATL15 calls ATL1415.tides.tide_elevations instead of pyTMD.compute directly,
      and --tide_directory is back to path_or_uri (it was briefly local-only).
  [ X ] AA moved to CATS2008-v2023 (default_args/AA_0331.txt), which is the model the zarr store
      is published for.  NOTE this is a MODEL VERSION CHANGE for Antarctica, not a rename, and
      it applies to LOCAL runs of the 0331 string too -- a local run now needs the TMD3 netCDF
      (CATS2008_v2023/CATS2008_v2023.nc) in its tide_directory, or a cloud --tide_directory.
      AA.txt and AA_0329.txt are older processing strings and were deliberately left on CATS2008.
      Greenland needed no change: GL_0321/0329/0331.txt already say Gr1km-v2.
  [ X ] MAAP_dps.txt --tide_directory=s3://pytmd, and run.sh's 39-line per-job staging block is
      gone along with the ATL1415_TIDE_S3 escape hatch.  There is nothing left to stage.
  [ X ] pyproject.toml declares pyTMD (pulls timescale) and zarr.  pyTMD was imported at module
      level by ATL11_to_ATL15 but declared NOWHERE, so `import ATL1415` could not have worked in
      a clean install -- build-env.sh's import check would have failed the DPS build there.
  [ X ] tests/test_tides.py: 5 offline tests (dispatch, store naming, that a local directory
      still reaches pyTMD with the same arguments, and the anonymous-access default) plus an
      opt-in live read of both stores, enabled with ATL1415_TIDE_NETWORK_TESTS=1.  Composing a
      real AA DPS args file end to end was checked too.

Still open:
  [ ] Ask upstream for a zarr branch in pyTMD's model.open_dataset, or for compute.tide_elevations
      to accept an already-open dataset.  Either would let ATL1415/tides.py's cloud branch be
      deleted and --tide_directory passed straight through.  Worth doing: the accessor sequence
      we now carry is short, but it is a copy of pyTMD internals and will drift.
  [ ] The two URL bugs in the netCDF path are still worth a PR even though we no longer depend
      on it: open_tmd3_dataset hands a pyTMD.utilities.URL to xr.open_dataset (no inferable
      backend), and the next line's pathlib.Path(input_file).name raises on it.
  [ ] --tide_adjustment is unaffected but unverified end to end: ATL1415 applies its own
      tide_adj_scale grid after the prediction and never passes pyTMD's apply_flexure, so the
      CATS2008-v2023 store's own `flexure` variable is not double-counted.  Confirm on the first
      real AA tile rather than assuming.
  [ ] CATS2008 (the OTIS model) still has no hf.CATS2008.out on our bucket.  Now that AA_0331
      uses CATS2008-v2023 this blocks nothing current, but AA.txt/AA_0329.txt still name it, so
      copy hf.CATS2008.v1.1 to hf.CATS2008.out if an older string is ever revived.

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
    [ X ] pyTMD files -- MOOT FOR DPS as of 2026-09-04 (commit 2535dfb).  A DPS run reads
        s3://pytmd/<model>.zarr and stages nothing, so none of the staging below is on the
        critical path any more.  What is on our bucket still matters only for a LOCAL or
        ADE-side run that points --tide_directory at a directory.
        - CATS2008 was staged to s3://maap-ops-workspace/ben_smith/tide_models/CATS2008/
          (grid_CATS2008, hf.CATS2008.old, hf.CATS2008.v1.1)
        - Gr1km-v2 IS complete there -- CONFIRMED against pyTMD's own database (v3.0.9),
          which resolves it to greenlandTMD_v2/{h,grid}_Greenland8.v2, both present.  (The
          model NAME and the DIRECTORY name differ; that used to matter because run.sh's
          staging step keyed off it, and that step no longer exists.)
        - z-only file sizes, if anyone ever stages these again:
          CATS2008     grid 25.8 MiB + hf 257.1 MiB = 283 MiB   (uv.CATS2008.out, 514 MiB,
                       and the duplicate hf.*.old, 257 MiB, are not needed for elevations)
          Gr1km-v2     grid 59.0 MiB + h  471.8 MiB = 531 MiB   (u_Greenland8_rot.v2, 944 MiB,
                       is not needed for elevations)
        - the real deltat.iers is still missing.
    [ X ] "CATS2008 IS BROKEN" -- RESOLVED BY OBSOLESCENCE, not by the fix it proposed.
        The old finding was real: pyTMD resolves CATS2008 to CATS2008/hf.CATS2008.out, the
        bucket has hf.CATS2008.old and hf.CATS2008.v1.1 and no .out, and .v1.1 is
        byte-identical to the canonical file on the public pyTMD bucket (same 269539196 bytes,
        same md5 over the first and last 4 MiB) while .old differs in the tail.  But AA_0331.txt
        now names CATS2008-v2023 and reads it from s3://pytmd, so no current processing string
        touches hf.CATS2008.out at all.
        STILL TRUE, and the reason not to delete this item: AA.txt and AA_0329.txt still name
        plain CATS2008, so reviving either older string -- locally or on MAAP -- needs
        hf.CATS2008.v1.1 copied to hf.CATS2008.out first.
        ALSO STILL TRUE, and newly relevant: default_args/MAAP.txt (the ADE variant) still says
        --tide_directory=/home/jovyan/my-private-bucket/tide_models, and CATS2008-v2023 is NOT
        staged there.  An ADE-side local run of AA_0331 composed with MAAP.txt will fail on a
        missing model.  Either stage the TMD3 netCDF there or point MAAP.txt at s3://pytmd too.
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
    [ X ] masks/EGM2008_geoid_h.nc -- IS on the bucket, at
        s3://maap-ops-workspace/ben_smith/ATL1415/masks/EGM2008_geoid_h.nc, 81.7 MiB, a real
        netCDF (dims lon=8641, lat=4321; vars lat, lon, geoid_h, geoid_free2mean) and not an
        lfs pointer.  The earlier "nothing matching it is on the bucket" note was wrong.
        Because it sits directly under masks/, --mask_dir joins the bare --geoid_file name
        onto it correctly, and it was read over s3 successfully (see Cloud paths below).
        NOTE its coordinates are lon/lat, not x/y, so whatever reads it must say so.
[ X ] Set up a location fromfile for MAAP (ADE)
    - default_args/MAAP.txt written: --mask_dir, --tide_directory, --ATL14_root, all on the
      ~/my-private-bucket mount.
    - --ATL11_root/--ATL11xo_top were deliberately dropped: they are consumed only by
      make_ATL11_index.py / setup_ATL11_xover.py, which are local-filesystem-only.
    - it does NOT yet name the staged index; a cloud run also needs
      --ATL11_earthaccess and --ATL11_index=<my-private-bucket>/ATL11_index/
[ X ] MAAP_dps.txt: WRITTEN (2026-09-04).  default_args/MAAP.txt stays the ADE variant;
    default_args/MAAP_dps.txt is the one to compose with for anything submitted to DPS.
    Use it INSTEAD OF MAAP.txt, not alongside it.  It sets:
      --mask_dir=s3://maap-ops-workspace/ben_smith/ATL1415/masks/   (read over /vsis3 + s3fs)
      --ATL11_earthaccess and --ATL11_index=s3://maap-ops-workspace/ben_smith/ATL11_index/
      --tide_directory=/tmp/ATL1415_static/tide_models  -- LOCAL, staged by run.sh
      --ATL14_root=/home/jovyan/ATL14_processing        -- ADE-side only, read by setup
    It deliberately does NOT set --ATL11_xover_dir: setup no longer derives it in cloud mode
    (that derivation assumes the local '<...>/index/GeoIndex.h5' layout) and
    setup_ATL11_xover.py still cannot write a schema with a remote 'source', so there is no
    cloud crossover tree to point at.  Crossovers are skipped in this mode.
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
    - the queue: `maap-dps-worker-32gb` and disk_space 20GB are guesses.  The SLURM runs use 4
      tasks and a 4-hour walltime per tile; nothing here has been measured.  Note the queue is
      overridable per submitJob via the queue_name input, so this is a default, not a commitment.
    - whether ~/.netrc reaches the worker, and whether the args file's s3:// paths resolve on
      the worker's own AWS credentials.
    (Resolved by the successful build: github.com reachability, sparseqr compiling, and the
    positional-vs-file input order.  See the Registration log.)
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

## Porting the howto workflows to MAAP
Plan written 2026-09-04; revised the same day against the answers in the QUESTIONNAIRE.
Target: docs/howto_MAAP_GL.sh, howto_MAAP_AA.sh and howto_MAAP_arctic.sh, alongside (not
replacing) the existing howto_GL.sh / howto_AA.sh / howto_arctic.sh, which get a one-line
header naming them the discover/SLURM variants.  `.sh` per Q3; no howto_MAAP_ATL11 per Q1 and
Q13 -- indexing stays on discover -- but the discover-to-S3 handoff still needs a home (Q20).

### What changes, and what does not
The layered @argsfile idiom, the region directory tree, the queue-file-of-shell-commands
and every script from make_200km_tiles.py onward all survive unchanged.  Three things
change:

  1. THE LOCATION LAYER.  default_args/discover.txt -> default_args/MAAP_dps.txt.
  2. THE HEMISPHERE LAYER GOES AWAY.  north.txt and south.txt exist only to supply
     --ATL11_index, --ATL11_xover_dir and --Hemisphere.  In cloud mode MAAP_dps.txt
     carries --ATL11_earthaccess plus the s3 index root, and there is no xover dir, so
     all that is left in them is --Hemisphere -- which is already a command-line option
     on setup_ATL1415_region.py.  The MAAP howtos therefore compose
       MAAP_dps.txt + latest_release.txt + <REGION>_latest.txt + quarterly|monthly.txt
     and pass --Hemisphere=1 / --Hemisphere=-1 directly.  (howto_GL.sh and howto_AA.sh
     already do this; it is the north.txt/south.txt pair that becomes dead weight.)
  3. THE SCHEDULER.  `setup_slurm_run.py ... ; sbatch slurm_run.sh` becomes either a DPS
     fan-out (prelim, matched) or a local parallel run of the same queue directory
     (everything else).

### The howtos themselves
WRITTEN 2026-09-05, and marked TENTATIVE in their own headers -- they were written BEFORE
any of the workflow has been run on DPS, deliberately, so that there is a plan to revise
rather than a blank page.  Revise them as testing advances; they are meant to be edited.

  docs/howto_MAAP_staging.sh   steps S1-S7, once per MAAP account.  S7 (the smoke test)
                               is the gate on everything else.
  docs/howto_MAAP_GL.sh        steps 1-12.  THE REFERENCE WORKFLOW -- the other two mark
                               only their differences and refer back to "GL step N".
  docs/howto_MAAP_AA.sh        steps 1-14.  Two halves, two queues, plus 200 km tiles,
                               sectors and per-sector mosaics.
  docs/howto_MAAP_arctic.sh    steps 1-11.  Five regions; IS alone first (Q21).

Every step carries a status tag -- [OK], [UNTESTED], [NEEDS CODE: x] -- and an [ADE]/[DPS]
marker.  The [NEEDS CODE] tags are the same list as the section below, so the two stay in
step: if you add an item there, tag the step that needs it.

The existing four discover howtos now carry a one-line header naming them the SLURM
variants (Q2).  howto_ATL11.sh says explicitly that it has no MAAP counterpart (Q1/Q13).

### Code that has to exist before the howtos are runnable
Listed in the order they block the sequence above.

[ ] make_ATL1415_queue.py HAS NOT HAD THE CLOUD-PATH TREATMENT setup_ATL1415_region.py
    got, and it is step 3 of every howto.  Four separate problems:
    - line 134, `os.path.isfile(defaults['--ATL11_index'])`: always False for a URI, and
      in cloud mode the index is a directory anyway.  It then re-joins the URI onto
      --ATL14_root and exits 1.  Blocks every cloud queue build outright.
    - lines 107-111, the --mask_dir join: plain os.path.join plus os.path.isfile, i.e.
      exactly the bugs join_path_or_uri()/paths.exists() were written to fix.
    - line 86, `defaults_re=re.compile(r'(.*)\s*=\s*(.*)')`: the greedy regex that drops
      bare store_true flags.  Fixed in setup_ATL1415_region.py, still live here -- so
      make_ATL1415_queue.py silently loses --ATL11_earthaccess and --tide_adjustment.
    - it needs a --xy_out (or an equivalent emitter) so a DPS fan-out gets tile centers
      rather than shell command lines.  The parse already exists as
      setup_ATL1415_run.py's xy0_re if a separate converter is preferred.

[ ] THE 1 km MASK.  See Q6 -- as things stand neither ice sheet can produce a tile list
    from the v4.1 masks.  This gates step 3 just as hard as the item above.  Q6 chooses
    (b) an explicit --grid_mask_file, (c) a frozen tile list in region_files/, and building the
    1 km decimation on the fly with gdal_translate and caching it.  Cache location, resampling
    rule and band selection are still open -- see Q16, all three of which change the tile set.

[ ] scripts/submit_MAAP_jobs.py and scripts/check_MAAP_jobs.py.  Nothing in the repo
    talks to maap.submitJob today.  The submitter loops the xy list, applies a rate /
    max-in-flight policy (Q11), and writes the ledger; the monitor is the
    slurm_run_status.py analogue (Q10).  This is the biggest new piece and the one that
    makes the howtos executable rather than aspirational.

[ ] scripts/run_queue_local.sh.  make_mosaic_jobs.py, make_200km_tiles.py,
    make_200km_to_mosaic_jobs.py and the to-nc wrappers all still emit
    <run_dir>/queue/task_N plus a slurm_run.sh, and there is no sbatch in the ADE.  A
    runner that reproduces packable_job.txt's bookkeeping -- mv queue->running->done,
    log to active_logs/ then logs/, copy to error_logs/ on nonzero exit or a
    ##TASK_LINE_FAILED## marker -- with `xargs -P` for concurrency lets
    slurm_run_status.py keep working unchanged on those directories.  One small script
    that unblocks steps 8 and 9 of all three howtos.

[ ] Deterministic tile output on S3 (Q9), then the matched neighbourhood (Q8).  These
    two are one problem: matched cannot be fanned out until a job can name its
    neighbours' prelim tiles by key.

[ X ] The arctic mask question is SETTLED (Q12), and it reduced to one line of code:
    masks/RGI_reduced/ is on the bucket as real files, and gdal's SQLite driver reads a .db
    over /vsis3 -- but ATL1415/make_mask_from_vector.py called ogr.Open(mask_file, 0) directly
    instead of through pc.io_utils.as_gdal_path(), so a raw s3:// URI raised "No such file or
    directory".  DONE 2026-09-05: that call now goes through as_gdal_path().  No format change,
    no re-staging.  Not yet exercised against the bucket -- the five regions Q21 wants to test
    on are the check.

[ ] Cloud crossovers (Q19: yes, include them).  REVISED after Q24: setup_ATL11_xover.py is
    NOT needed on MAAP at all, and neither is a schema file or --ATL11_xover_dir.  The work is
    entirely inside read_ATL11_xovers: in cloud mode build the per-cycle tilingSchema in
    memory with source={'type':'EarthAccess','short_name':'ATL11XO'} and hand it to the same
    resolve_files_for_box() that already works.  Supporting changes: promote --ATL11xo_version
    to a real ATL11_to_ATL15 argument (it is already in every composed args file and is being
    silently discarded by parse_known_args), and let --ATL11_earthaccess be the switch that
    turns crossovers on, since --ATL11_xover_dir is the local-mode switch and has nothing to
    point at in the cloud.  setup_ATL11_xover.py stays exactly as it is for discover.
    Verified prerequisites: ATL11XO is a real CMR collection (v007), granule names match the
    format_str character for character, and the labels are tile CENTERS as 'round' assumes.

[ ] SEPARATELY, AND NOT A MAAP ITEM: settle the xover_cycles=[1,2] question (Q26) before
    rel006 production.  The product has 31 cycles; the code reads 2, on discover too.

[ X ] docs/howto_MAAP_staging.sh (Q20) WRITTEN 2026-09-05, tentative: standing up a fresh MAAP account with everything
    ATL1415 needs -- build the ADE env, stage the Zenodo masks, stage the ATL11 index built on
    discover, register the DPS algorithm, and verify each is readable.  Tide models need no
    staging (s3://pytmd), and per Q24 neither do crossover schemas, so the only per-release
    handoff left is the ATL11 index itself.

### Suggested order
DECIDED 2026-09-05: WRITE THE HOWTOS FIRST, marked tentative -- "you have to have a plan
before you can change it".  They are the spec and the acceptance criteria for the code
items below, rather than trailing them.  Item 0 is done.

  0. [DONE 2026-09-05] Write the four howtos, tentative and numbered.
  1. Smoke-test one sandbox DPS job (staging S7).  It settles five things that no amount
     of reading can, including whether submitJob works on this account at all, and it is
     what sizes the production queue.        -> unblocks GL step 5, AA 6, arctic 6
  2. make_ATL1415_queue.py cloud fixes + --xy_out, and settle the 1 km mask (Q6/Q16).
                                             -> GL step 4, AA 4-5, arctic 5
  3. submit_MAAP_jobs.py + check_MAAP_jobs.py.  The biggest new piece.
                                             -> GL steps 5-6, AA 6-7, arctic 6-7
  4. Deterministic output prefix, then the matched neighbourhood (Q9 then Q8).  These are
     one problem: matched cannot fan out until a job can name its neighbours by key.
                                             -> GL steps 7 and 9
  5. run_queue_local.sh.                     -> GL steps 10-11, AA 10/12/13, arctic 10
  6. The arctic .db mask check (Q12) and the previous-product items (Q27 W1/W3/W4).
     W1 is a silent-failure bug and should not wait for a production run to surface it.

The howtos are now the thing to revise as each of these lands, rather than the thing to
write at the end.
