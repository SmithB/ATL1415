# plan_merge_to_main.sh -- MERGING on_s3 INTO main, AND MOVING DISCOVER ONTO IT
#
# ############################################################################
# ##  TENTATIVE, AND DEFERRED.  Written 2026-09-12.  NOTHING HAS BEEN DONE. ##
# ##  Ben, 2026-09-12: "I'm not ready to move on_s3 to main until it's      ##
# ##  closer to final."  Saved so the strategy is not re-derived.  Its      ##
# ##  FINDINGS are a snapshot and WILL go stale as on_s3 moves -- run M0    ##
# ##  again before acting on any of them.                                   ##
# ############################################################################
#
# BRANCH ROLES (Ben, 2026-09-12): pre_rel006 is the ACTIVE branch on discover;
# on_s3 is the development branch on MAAP; the intent is to merge on_s3 into
# main, assuming on_s3 contains everything from pre_rel006.
#
# Tags: [READY] runnable as written.  [DEFERRED] waits for on_s3 to be closer
#       to final.  [DECISION] Ben's call.  [ADE] / [DISCOVER] / [GITHUB] where
#       it runs.
# EVERY CLAIM IS LABELLED.  STATEMENT = verified, and how.  RECOMMEND = my
# proposal.  QUESTION = I do not know and did not guess.
# Split out of plan_remove_lfs.sh L7, which now points here.


# ===========================================================================
# THE STRATEGY, in one paragraph
# ===========================================================================
# RECOMMEND: merge everything INTO on_s3 and resolve it there, where the code
# runs.  After that, main moves forward to on_s3 with nothing to resolve, and
# discover's pre_rel006 is already part of main's history, so discover
# switches branches without a merge.  NEVER squash or rebase any of it: build
# stamps and tiles record commit shas, and those commits have to stay in
# main's history (plan_remove_lfs QL2).


# ===========================================================================
# FINDINGS -- SNAPSHOT of 2026-09-12, after `git fetch`
#   on_s3 = 41d5bd3   origin/main = 048b66f   origin/pre_rel006 = d5ef8f5
#   Trial merges were run with `git merge-tree`, which does not touch the
#   working tree.
# ===========================================================================
# F1. STATEMENT (`git cherry`, trial merge): origin/on_s3 has every patch on
#     origin/pre_rel006 (0 unique).  The one commit not in its history,
#     d5ef8f5 (the .db line in .gitattributes), matches a change on_s3 already
#     has.  Trial merge: clean.
#
# F2. STATEMENT (`git cherry`, then reading on_s3's code -- READ, NOT RUN):
#     origin/main has 10 commits not on on_s3; 4 of them are changes on_s3 does
#     not have, and each is already done on on_s3 another way:
#       7a204ef  pyTMD v3 (PR #25)        on_s3: same tide_elevations arguments
#       57fb4de  ATL15_write2nc main()    on_s3: main() at ATL15_write2nc.py:285
#       b77926d  attrs csv name           on_s3: reads ATL15_output_attrs.csv
#       20c8adf  datetime / numpy bug     on_s3: rewrote the block with int()
#
# F3. STATEMENT (trial merge of origin/main into on_s3): 5 conflicts --
#     ATL11_to_ATL15.py (3 hunks), ATL1415_attrs_meta.py (2; main's copy still
#     has the old ATL14_attrs_meta.py name), scripts/ATL15_write2nc.py (3),
#     pyproject.toml (1), and resources/ATL15_monthly_output_attrs.csv (on_s3
#     deleted it in dfe4b7c; main modified it).  In every hunk, on_s3's side is
#     the one to keep.
#
# F4. STATEMENT -- THE TRAP (measured with `git merge-tree -X ours`, and git
#     grep): taking on_s3's side only in the CONFLICTING HUNKS is not enough.
#     Git also applies main's non-conflicting changes.  One of them makes
#     get_SRS_info() return only the proj4 string, while on_s3's only caller,
#     ATL11_to_ATL15.py:562, does `SRS_proj4, EPSG = get_SRS_info(hemisphere)`
#     -- the first tile would fail with a ValueError.  The same merge also adds a
#     duplicate "crs:" docstring line and a commented-out block, and brings the
#     deleted monthly CSV back.
#     So do NOT use `-X ours`.  `git checkout --ours <file>` takes on_s3's WHOLE
#     file and avoids all of it; the merged tree should then be exactly
#     on_s3's (M3 checks this).
#
# F5. STATEMENT (gh api): main has no branch protection; merge commits,
#     squash and rebase are all enabled; main took #23 and #25 as merge commits.


# ===========================================================================
# STEPS
# ===========================================================================
# M0. [ADE] [READY -- RUN THIS FIRST, EVERY TIME]  Refresh the snapshot.
#     Anything on on_s3 or main after 2026-09-12 can change F1-F4.
cd ~/git_repos/ATL1415 && git fetch origin
git cherry origin/on_s3 origin/pre_rel006 | grep -c '^+'    # F1: expect 0
git cherry origin/on_s3 origin/main | grep '^+'             # F2: re-check each
git merge-tree --write-tree --name-only --messages origin/on_s3 origin/main   # F3
T=$(git merge-tree --write-tree -X ours origin/on_s3 origin/main | head -1)
git diff --stat origin/on_s3 $T                             # F4: what main adds on top of on_s3's side
git grep -n 'get_SRS_info' origin/on_s3 origin/main         # F4's caller
#     If the conflict list or F4's diff has changed, revise this file before M2.

# M1. [ADE] [READY]  Land plan_remove_lfs L4 (the LFS commit) on on_s3 before
#     M3, so the merge carries it to main.  Independent of this plan's timing.

# M2. [ADE] [DEFERRED]  Merge pre_rel006.  Nothing changes in the files; git
#     only records that pre_rel006 is merged, which is what makes M5 a plain
#     switch.  REDO it if discover pushes to pre_rel006 before M4.
git switch on_s3
git merge origin/pre_rel006 -m "Merge pre_rel006 (discover) into on_s3"

# M3. [ADE] [DEFERRED]  Merge main; take on_s3's version of every conflicted file.
git merge origin/main                      # stops with F3's conflicts
git checkout --ours -- ATL1415/ATL11_to_ATL15.py ATL1415/ATL1415_attrs_meta.py \
                       ATL1415/scripts/ATL15_write2nc.py pyproject.toml
git rm -q ATL1415/resources/ATL15_monthly_output_attrs.csv
git diff --cached --stat HEAD              # EXPECT EMPTY: the result equals on_s3's files.
#                                            If not empty, stop and read the diff.
git commit -m "Merge main into on_s3; main's patches already exist on on_s3 (plan_merge_to_main M3)"
git push origin on_s3
#     Because no file changes, whatever has already been tested on on_s3 still
#     holds; the merge changes history only.  Push before the next DPS
#     registration (held commits block it).
#     RECOMMEND: M2 and M3 can happen well before M4 if you want main's
#     history on on_s3 early.  Doing it early keeps the later conflicts from
#     growing, but M0 still applies.

# M4. [GITHUB] [DEFERRED -- "closer to final"]  Move main.  After M3, main
#     is in on_s3's history, so there is nothing to resolve.  RECOMMEND a PR
#     on_s3 -> main merged with "Create a merge commit" (F5).  Without a PR,
#     the equivalent is:
#       git switch main && git merge --ff-only origin/on_s3 && git push origin main
#     NOT "Squash and merge" or "Rebase and merge".

# M5. [DISCOVER] [DEFERRED; DECISION: timing -- production runs from pre_rel006]
#     Move discover to main:
#       1. Copy real masks/ files out of the tree (plan_remove_lfs L3).  If
#          main no longer tracks the 16 LFS pointer paths, git deletes any that
#          are still pointers, and refuses to switch if a real file sits on one.
#       2. git fetch origin && git switch main   (pre_rel006 is in main's
#          history, so no merge is needed)
#       3. copy the masks back; reinstall ATL1415 if it is not an editable install.
#     QUESTION: I have no record of on_s3's code running in discover's local
#     (SLURM, non-cloud) mode since the cloud-read work.  RECOMMEND one tile
#     there before production uses main.

# M6. [ADE] [DEFERRED; DECISION: Ben]  After M4: algorithm_config.yml says to
#     change algorithm_version to main "once the cloud-read work is merged".
#     Whether development moves to main or stays on on_s3 is your call;
#     changing it means re-registering.
