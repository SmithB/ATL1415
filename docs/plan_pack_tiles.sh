# plan_pack_tiles.sh -- several tiles per DPS job, and retries at the MAAP
# API calls that fail at job start
#
# ############################################################################
# ##  WRITTEN 2026-10-06.  TENTATIVE: nothing implemented yet.              ##
# ##  Ben 2026-10-06: "It seems like it will turn out to make sense to pack ##
# ##  multiple tiles into single jobs, so that we can have more than ~100   ##
# ##  tiles running at once.  Small tiles can run in parallel on a single   ##
# ##  32gb node, but larger tiles will need to run serially (or in their    ##
# ##  own jobs.)  We should also add a catch-and-try-again to the call to   ##
# ##  MAAP() and to any other points of failure that are catching us        ##
# ##  (although we don't want to pay for nodes that are waiting around)."   ##
# ############################################################################
#
# Status tags per step: TODO / DONE / BLOCKED.  DECIDED = Ben said so;
# RECOMMENDATION = mine; QUESTION = open.
#
# QUESTIONS FOR BEN -- ALL DECIDED 2026-10-06 (Ben: "Go with your
# recommendations on QK1-QK4, start K1"):
#   QK1 DECIDED (recommendation taken). Job size.  RECOMMENDATION: about 1 h of wall time per job, 2 tiles
#        at a time on a 32gb node (lanes, below) -> ~10 GL prelim tiles per
#        job, so the 1112 GL prelim tiles are ~110 jobs, all in flight at
#        once under today's 100-in-flight cap.  Longer jobs lose more to a
#        mid-job failure; shorter ones bring back the start-up load.
#   QK2 DECIDED (recommendation taken). The 401 case (our credentials step after the runner's
#        get_maap_pgt_token.py timed out): retrying with the same MAAP_PGT
#        cannot succeed.  RECOMMENDATION: stop at the first 401 instead of
#        retrying 5 times (saves ~1 min of node per such job), unless K1
#        finds that a job can fetch a fresh token itself.
#   QK3 DECIDED (recommendation taken). pointCollection's broker retry (K3) is in your repo.
#        RECOMMENDATION: I write it on a branch, you merge, as with
#        blockcache_reads.
#   QK4 DECIDED (recommendation taken). Scope: prelim and matched only.  The mosaic steps already batch
#        (one job per mosaic task).  RECOMMENDATION: leave them alone.
#
# WHAT IS KNOWN (statements, with provenance):
#   - FAILURES ARE ALL MAAP API CALLS AT JOB START (plan_GL_maskv5.sh V4b,
#     all 540 failed jobs of rounds 0-1 classified).  Four call sites:
#       a. MAAP's stage_in.py: `maap = MAAP()` at import ->
#          /api/environment/config, no retry.  MAAP's code; we cannot
#          change it.  326 of 534 in round 0, 4 of 6 in round 1.
#       b. MAAP's runner get_maap_pgt_token.py timing out, then HTTP 401 on
#          every call that uses MAAP_PGT.  MAAP's code.  81 + 1.
#       c. Our scripts/workspace_credentials.py: 5 attempts, 10 s apart,
#          30 s timeout each (line 40); a new MAAP() per attempt (line 54).
#          76 refused/timeout (plus the 81 of b).
#       d. pointCollection io_utils._s3fs_from_maap (NSIDC credentials): 5
#          attempts, 10 s apart; a new MAAP() per attempt, so each attempt
#          is two API calls (config + credentials).  49 in round 0, inside
#          the fit (~50 s, 0.7 GiB in).
#     a and b happen once per JOB, so packing N tiles per job divides them
#     by N.  Retries (K2-K3) can only cover c and d.
#   - THE 100-IN-FLIGHT CAP IS MINE, NOT MEASURED.  Round 0 (1112 starting
#     within minutes) failed 48%; round 1 (100 in flight) failed 1.1%.
#     Nothing between those has been tried.
#   - GL PRELIM TILE SIZES (578 successful round-0 jobs, collect_jobs,
#     e6d7051): peak memory p10 4.7, p50 9.1, p90 10.4, max 13.2 GiB --
#     none above 14 GiB, so any two fit in 32 GiB.  Wall time p10 338, p50
#     568, p90 731, max 907 s.  Average cores in use: fit ~1.1, error step
#     ~1.7 on r5.xlarge (4 vCPU = 2 physical cores).  The 32gb queue also
#     lands on m5.2xlarge (127 of 578) and c5.4xlarge (61).
#   - NOT YET MEASURED: GL matched, IS, AA tile sizes (AA tiles are larger;
#     the 28 GiB gate in plan_GL_maskv5.sh V4a came from AA).
#
# DESIGN (RECOMMENDATION):
#   - New optional job input `tiles`: "x0,y0;x0,y0;..." (default '-').  With
#     '-', run.sh behaves exactly as now (x0/y0), so nothing already
#     registered or scripted changes meaning.
#   - Each job runs its tiles in L LANES: L tiles at a time, each lane
#     working through its share in series.  The SUBMITTER picks L and the
#     split from each tile's predicted peak memory (from the previous run of
#     that tile, else from N_ATL11), with the sum of L tiles' predicted
#     peaks <= ~26 GiB.  Big tiles get L=1 or a job of their own: Ben's
#     "larger tiles will need to run serially (or in their own jobs)".
#   - Credentials (c) and the NSIDC session (d) are fetched once per job and
#     refreshed by the existing expiry logic, not once per tile.
#   - Each tile is a separate process, so one tile's crash does not kill
#     the others; run.sh prints one `TILE_STATUS <tile> ok|nodata|failed`
#     line per tile and exits 1 if any tile failed.
#   - The DRIVER decides what to resubmit PER TILE, not per job: a tile is
#     done when its .h5 is under tile_prefix with LastModified after the job
#     started, or it is listed nodata.  A job that dies part-way loses only
#     its unfinished tiles.
#   - Waiting costs node time: retries at c and d use backoff with jitter
#     (e.g. 10, 20, 40, 60, 60 s) and a cap of ~4 min per call site.  An r5
#     node waiting 4 min costs about $0.02, against losing the whole job.
#
# ===========================================================================
# K1. [ADE] DONE 2026-10-06.  NO: a job cannot fetch a fresh MAAP_PGT.
#   - The runner (verdi side, before cwltool) runs /app/get_maap_pgt_token.py:
#     GET https://api.maap-project.org/api/members/<username> (traceback in
#     the round-0 logs; requests.get with no timeout, no retry).
#   - maap-api-nasa (7f5bd6e) api/endpoints/members.py Member.get: the
#     session_key (the PGT) is returned only if valid_dps_request() -- header
#     `dps-token` equal to settings.DPS_MACHINE_TOKEN (api/auth/security.py
#     269) -- or to the user's own session, which needs a PGT already.
#   - cwltool runs our container with --preserve-environment MAAP_PGT and
#     MAAP_API_HOST only, so the machine token never reaches our code (and
#     using it would be going around MAAP's design if it did).
#   - In all 80 round-0 jobs with the token timeout + 401, MAAP_PGT was still
#     non-empty (run.sh did not print "MAAP_PGT is not set"), so it held a
#     stale or invalid value; every one of 5 attempts got 401.
#   => QK2 as decided: stop at the first 401, naming the runner's token
#     fetch as the likely cause.  Only MAAP can fix the fetch itself (a
#     timeout + retry in get_maap_pgt_token.py): for the admin note.
# K2. [repo] DONE 2026-10-06 (uncommitted; not deployed -- needs
#   registration).  scripts/workspace_credentials.py: pauses 10, 20, 40, 60,
#   60 s (up to 6 tries), each x uniform(0.5, 1.5); no pause that would end
#   past 240 s from the first try (worst case ~270 s with the 30 s per-try
#   timeout; was 5 tries 10 s apart, ~3 min); one MAAP client kept across
#   tries (built again only if building it failed); HTTP 401 stops at the
#   first try with the runner's token fetch named.  tests/
#   test_workspace_credentials.py 26 passed (8 new or rewritten); whole
#   suite 257 passed, 2 skipped.
# K3. [pointCollection] DONE 2026-10-07 on branch maap_broker_backoff
#   (f7c3ca4, PR SmithB/pointCollection#62; Ben merges).  _s3fs_from_maap: K2's pauses,
#   jitter and 240 s budget; one MAAP() across tries (built again only if
#   building failed); 401 stops at once.  No per-try timeout (no SIGALRM off
#   the main thread).  test_maap_broker.py 16 passed (5 new); suite 323
#   passed, 3 skipped.
# K4. [repo] TODO.  run.sh `tiles` input + lanes + TILE_STATUS;
#   algorithm_config.yml input `tiles` (string, default '-').  Needs
#   registration (Ben).
# K5. [repo] TODO.  submit_MAAP_jobs.py --pack: bin tiles into jobs from a
#   prediction table (tile -> peak GiB, secs); the driver's per-tile
#   resubmit.
# K6. [DPS] TODO.  Smoke: 3 jobs of 4-6 already-solved GL prelim tiles to a
#   TEST prefix; gates: every tile identical to its single-tile run
#   (compare_products), peak memory under 30 GiB, and per-tile wall time vs
#   the single-tile runs (the 2 lanes share 2 physical cores on r5.xlarge,
#   so tiles may run slower; the gain is fewer job starts).
# K7. [DPS] TODO.  Scale test on a real step: the in-flight cap raised in
#   steps (100 jobs, then 200) with failure rates recorded, to measure
#   what the API takes instead of guessing.
