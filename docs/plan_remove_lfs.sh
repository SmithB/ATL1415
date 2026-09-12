# plan_remove_lfs.sh -- TAKE GIT-LFS OUT OF THE ATL1415 REPO
#
# ############################################################################
# ##  TENTATIVE.  Written 2026-09-12.  NOTHING BELOW HAS BEEN DONE: no      ##
# ##  file untracked, no attribute removed, nothing rescued.  This is the   ##
# ##  diagnosis and the sequence.  QL1-QL4 ANSWERED by Ben, 2026-09-12.   ##
# ############################################################################
#
# THE GOAL (Ben, 2026-09-12): no LFS anywhere in the repo.  Masks are copied
# into masks/ as needed, and copied-in files must stop being overwritten with
# LFS pointers.
#
# Tags: [READY] runnable as written.  [DECISION: QLn] waits on a question.
#       [ADE] / [DISCOVER] / [GITHUB] where it runs.
# EVERY CLAIM IS LABELLED.  STATEMENT = verified, and how.  RECOMMEND = my
# proposal, yours to take or drop.  QUESTION = I do not know and did not guess.
# The open item at docs/Transition_to_maap.md "OTHER" points here.


# ===========================================================================
# QUESTIONS -- answer on the A: lines; the steps that depend on each say so
# ===========================================================================
# QL1. Twelve LFS-tracked files have NO real copy on the bucket or in this
#      checkout (F5): the Greene_22 shelf masks (.h5, _full.h5, _1km, _240m),
#      scripps_antarctica_IceShelves1km_v1, bedmap2_..._sio_shelves{,_edited},
#      the four Arctic U_Texas / Ice_Ocean_Bed / BedMachineGreenland shelf
#      files, and ATL11_0314_tide_adj_scale_200m.h5 -- plus notebooks/RGT_001.h5.
#      As far as I can tell, their only real copies are on GitHub's LFS store
#      and possibly on discover.  Removing the pointers loses nothing that is
#      there now, but once the repo stops using LFS nobody will think to go
#      back for them.  What to do with them?
#        (a) copy them from discover, if the originals are still there;
#        (b) fetch just those objects from GitHub LFS, ONCE, into scratch,
#            check sha256 against the pointer oid, and upload them to the bucket
#            (a one-off archival exception to "masks come from Zenodo, never lfs");
#        (c) abandon them -- no rel_006 args file names any of them.
#      RECOMMEND (a), with (b) as the fallback for anything (a) cannot find.
#      About 700 MB, or ~905 MB counting the tide .h5, which the Zenodo .tif
#      has replaced.  (b) counts against the owner's LFS bandwidth quota.
#   A: (c) ANSWERED 2026-09-12 (Ben): let them go.  L1 is dropped.
#
# QL2. Rewrite history to take the pointers out of old commits and tags, or
#      leave history alone?
#        (a) leave it: only the live branch tips change;
#        (b) rewrite (git lfs migrate export / git filter-repo) on every
#            branch and on tags v1.0.0-v2.1.0.
#      RECOMMEND (a).  (b) changes every commit sha, which breaks the commit
#      recorded in each DPS build stamp and tile, the tags, and every other
#      clone.  The cost of (a) is residual risk L9: checking out an OLD commit
#      still writes pointers.
#   A: (a) ANSWERED 2026-09-12 (Ben): leave the history alone.  L9 is accepted.
#
# QL3. Which branches get the fix?  on_s3 for sure.  main by merging on_s3,
#      or with its own commit now?  And is pre_rel006 still alive, or retired?
#      RECOMMEND on_s3 now, main through the on_s3 merge you already plan,
#      and pre_rel006 left alone if it is retired.
#   A: ANSWERED 2026-09-12 (Ben): pre_rel006 is the ACTIVE branch on
#      discover; on_s3 is the development branch on MAAP; the intent is to
#      merge on_s3 into main.  The fix lands on on_s3 and reaches main and
#      discover through that merge -- L7, deferred; the strategy is in
#      docs/plan_merge_to_main.sh.
#
# QL4. The Transition note says files "get overwritten when ATL1415 gets
#      installed".  I could NOT find an install step that writes to the
#      working tree (F4, last paragraph).  Do you remember what you ran just
#      before it happened -- a pull, a branch switch, a stash, or the pip
#      install itself?  If it really was the install, there is a mechanism
#      I have not found, and L8 needs a test for it.
#   A: ANSWERED 2026-09-12 (Ben): misspoke -- it is certain git actions, not
#      the install.  That is F4's M1/M2; no install test is needed.


# ===========================================================================
# FINDINGS -- how the repo is configured today (on_s3 @ 41d5bd3)
# ===========================================================================
# F1. STATEMENT (read the files): the LFS configuration is ATTRIBUTES ONLY.
#       .gitattributes (root):
#         *.tif filter=lfs ...          <- REPO-WIDE, not just masks/
#         *.h5  filter=lfs ...          <- REPO-WIDE
#         ATL1415/resources/BRW_template.h5 !filter      (exempts it; real, 27756 B)
#         #*.db filter=lfs ...          <- commented out on on_s3, ACTIVE on main
#         masks/EGM2008_geoid_h.nc -filter=lfs ...
#       masks/Antarctic/.gitattributes: five Greene_22 names (duplicates the
#         root *.tif / *.h5 rules; one names a _2022_update.h5 that is not in
#         git).
#       build/lib/masks/Antarctic/.gitattributes: a stale copy in the ignored
#         build/ directory, which git does not track.
#
# F2. STATEMENT (git config --show-origin, .git/hooks, which git-lfs): NO LFS
#     FILTER IS INSTALLED, so nothing converts pointers to or from real files.
#       - no filter.lfs.* in the global or local config; there is no
#         /etc/gitconfig.
#       - .git/hooks holds only the stock *.sample files, so no LFS hooks.
#       - local config holds just lfs.repositoryformatversion=0, and there is
#         a .git/lfs/ directory created 2026-09-03 23:08 holding one fetched
#         object (GL_Ed2z0dx2.tif, 7514 B).  Both are left over from a one-off
#         lfs command.
#       - git-lfs 3.7.1 is on PATH from /srv/conda/envs/notebook (the platform
#         image, not ours).  ATL14_notebook.yml also lists git-lfs; no other
#         file references that yml.
#       - DPS: the maap_base image is debian + git + Miniforge
#         (algorithm_config.yml), and environment.yml has no git-lfs, so a
#         DPS clone gets the pointer TEXT.  run.sh reads masks from s3
#         (MAAP_dps.txt --mask_dir), so DPS never reads repo masks/.
#     CONSEQUENCE: LFS is not "running" anywhere.  The overwrites are ordinary
#     git putting committed 130-byte pointer blobs back on disk (F4).
#
# F3. STATEMENT (git lfs ls-files per ref, plus a scan of HEAD blobs for the
#     pointer header): the committed pointer blobs.
#       on_s3, origin/on_s3, pre_rel006: 16 files, and every working-tree
#       copy in this checkout is a pointer (129-134 B):
#         masks/Antarctic/  AA_Ed2z0dx2.tif (42 KB real)
#                           ATL11_0314_tide_adj_scale_200m.h5 (203 MB)
#                           BedMachineAntarcticaOceanv2.tif (3.3 MB)
#                           Greene_22_shelf_plus_10m_mask.h5 (221 MB)
#                           Greene_22_shelf_plus_10m_mask_1km.tif (685 KB)
#                           Greene_22_shelf_plus_10m_mask_240m.tif (26 MB)
#                           Greene_22_shelf_plus_10m_mask_full.h5 (405 MB)
#                           bedmap2_thickness_gt_50_plus_sio_shelves.tif (3.4 MB)
#                           scripps_antarctica_IceShelves1km_v1.tif (1.6 MB)
#                           updates/bedmap2_..._sio_shelves_edited.tif (918 KB)
#         masks/Arctic/     BedMachineGreenland-2021-04-20_shelf_125m.tif (6.3 MB)
#                           GL_Ed2z0dx2.tif (7.5 KB)
#                           Ice_Ocean_Bed_100m_2019_compress.tif (10 MB)
#                           U_Texas_ice_mask_2019_100m.tif (30 MB)
#                           U_Texas_ice_mask_2019_1km.tif (161 KB)
#         notebooks/RGT_001.h5 (133 KB)
#       main / origin/main: those 16, PLUS the five RGI_reduced/*.db files and
#         masks/EGM2008_geoid_h.nc as pointers.  on_s3 fixed those in
#         35d86f4 / 136c224 / d5ef8f5 / c8add6e; main has not had that merge.
#       tags v1.0.0 ... v2.1.0: the .db files, EGM2008 and RGT_001 are pointers.
#       history only: RGI_reduced/*_40km.tif / *_80km.tif (cc4c3d5, 3a6800b),
#         GrimpIceMask_2018.1_2023.9_100m.tif (untracked in a9c7ac3).
#
# F4. STATEMENT (reproduced in a throwaway repo, 2026-09-12): the ways a real
#     file gets replaced by a pointer.
#     M1. TRACKED PATH.  Copying a real mask onto one of F3's paths only
#         makes git see that file as modified.  `git stash`, `git checkout -- .`,
#         `git reset --hard`, or any pull, merge or switch that touches the
#         path puts the pointer back.  `git commit -a` would do the opposite
#         and commit the full binary.
#     M2. IGNORED PATH THAT ANOTHER COMMIT TRACKS.  If a real file sits at a
#         path that .gitignore covers (masks/* does) and you check out a commit
#         that tracks that path, git overwrites the file SILENTLY, with no
#         "would be overwritten" error, because git treats ignored files as
#         disposable.  .gitignore does not protect these files; it is the
#         reason there is no warning.
#         Concrete, today: from on_s3, `git switch main` replaces the real .db
#         files with pointers, and would do the same to a real
#         EGM2008_geoid_h.nc; checking out any tag does the same.
#     M3. (a stopgap, not the fix) `git update-index --skip-worktree <path>`
#         makes git REFUSE instead of overwrite -- the demo's switch failed
#         with "local changes would be overwritten", and the real file
#         survived.
#     QL4 confirmed it: the overwrites come from git actions, not the install
#     (build-env.sh's `pip install .` does not write the source tree).
#
# F5. STATEMENT (stat on ~/my-private-bucket/ATL1415/masks, compared against
#     the pointers' size fields -- SIZE only, not checksum):
#       real on the bucket, size matches: AA_Ed2z0dx2.tif, GL_Ed2z0dx2.tif,
#         BedMachineAntarcticaOceanv2.tif -- the three that AA_0331.txt /
#         GL_0331.txt use -- and EGM2008_geoid_h.nc and all five .db files.
#       130-byte POINTERS on the bucket: the other twelve masks in F3.  That
#         is every pointer on the bucket (find -size -300c).
#       not on the bucket: GrimpIceMask_2018.1_2023.9, 01_Alaska_80km,
#         05_GreenlandPeriphery_80km (history only; nothing current uses them).
#       RGI_reduced 40/80 km tifs: real on the bucket; they are not tracked at
#         any branch tip.
#
# F6. STATEMENT (git grep): what still NAMES F3's files.
#       current release: AA_0331.txt (AA_Ed2z0dx2, BedMachineAntarcticaOceanv2),
#         GL_0331.txt (GL_Ed2z0dx2); scripts/maap/{find_tide_tiles,make_AA_queue}.py
#         (BedMachineAntarcticaOceanv2).  All three are real on the bucket (F5).
#       older processing strings only: AA.txt, AA_0329.txt, GL_0321.txt,
#         GL_0329.txt, old/GL_0319.txt, masks/Antarctic/make_time_varying_mask.ipynb.
#       The code never reads repo masks/ through a relative path: --mask_dir
#       supplies the directory (MAAP.txt: the bucket mount; MAAP_dps.txt: s3).


# ===========================================================================
# STEPS
# ===========================================================================
# L0. [ADE] [READY]  Start clean.
#     on_s3 clean and level with origin/on_s3, because Ben registers from this
#     checkout and held commits block registration.
git -C ~/git_repos/ATL1415 status --short && git -C ~/git_repos/ATL1415 rev-list --left-right --count origin/on_s3...on_s3

# L1. [DROPPED -- QL1 (c)]  No rescue.  The twelve lfs-only files are
#     abandoned.  Their 130-byte pointer files on the bucket are handled in L10.

# L2. [ADE] [READY]  Checksum the three in-use masks on the bucket.  F5 matched
#     SIZE only; run this before the pointer oids are no longer on the branch tip.
cd ~/git_repos/ATL1415
for p in masks/Antarctic/AA_Ed2z0dx2.tif masks/Antarctic/BedMachineAntarcticaOceanv2.tif masks/Arctic/GL_Ed2z0dx2.tif; do
  want=$(git show HEAD:$p | sed -n 's/^oid sha256://p')
  have=$(sha256sum ~/my-private-bucket/ATL1415/$p | cut -d' ' -f1)
  [ "$want" = "$have" ] && echo "OK   $p" || echo "DIFF $p"
done

# L3. [EVERY CLONE] [READY]  Back up real files before pulling L4 anywhere else.
#     STATEMENT (git behaviour): when a clone pulls a commit that untracks a
#     path, git deletes that path from the working tree if it is unmodified
#     (a pointer -- harmless).  If a real file was copied over it, the pull
#     refuses instead, and `git stash` would put the pointer back (M1).  So
#     in each other clone (discover included):
#       copy real masks/ files to a directory outside the repo;
#       git checkout -- masks/ notebooks/ ; git pull ; copy the files back.

# L4. [ADE] [READY]  THE FIX -- one commit on on_s3.
cd ~/git_repos/ATL1415
git rm --cached -- $(git lfs ls-files -n)          # the 16 in F3; leaves the files on disk
git rm -- .gitattributes masks/Antarctic/.gitattributes
#     Removing the whole root .gitattributes is RECOMMENDED: every line in it
#     is LFS-related, and the BRW_template / EGM2008 exemptions mean nothing
#     without the LFS rules.
#     STATEMENT (.gitignore): masks/* and notebooks/ already ignore the 16
#     files once they are untracked, so no .gitignore change is needed.
#     Then delete the leftover pointer files on disk -- only files that are
#     tiny AND start with the LFS header, so a real copied-in mask is never
#     removed:
find masks notebooks -type f -size -300c -exec grep -lq '^version https://git-lfs' {} \; -print   # review first
#     ...then run the same command with -delete instead of -print.
git commit -m "Remove git-lfs: untrack the pointer files, drop the lfs attributes"

# L5. [ADE] [READY]  Remove local LFS leftovers from this clone (F2).
#     Do NOT run `git lfs uninstall`: it edits the GLOBAL config, and F2
#     found nothing there to remove.
git -C ~/git_repos/ATL1415 config --local --unset lfs.repositoryformatversion
rm -rf ~/git_repos/ATL1415/.git/lfs
rm -rf ~/git_repos/ATL1415/build      # ignored stale build copy, incl. its .gitattributes (optional)

# L6. [ADE] [READY]  RECOMMEND: delete `- git-lfs` from ATL14_notebook.yml
#     (nothing references that file).  Leave git-lfs in the /srv notebook env
#     alone -- the platform image owns it.

# L7. [DEFERRED -- Ben 2026-09-12: "not ready to move on_s3 to main until it's
#     closer to final"]  Getting the fix to main and discover.  The strategy
#     lives in docs/plan_merge_to_main.sh (steps M0-M6); it is not repeated here.
#     Until then: L4 lives on on_s3 only, discover's pre_rel006 still tracks
#     the pointer files, and L3 applies whenever discover finally moves.

# L8. [ADE] [READY after L4]  Verify.
cd ~/git_repos/ATL1415
git lfs ls-files                                       # expect: nothing
git grep -lI '^version https://git-lfs' -- ':!docs/'   # expect: nothing
git check-attr filter -- masks/Arctic/x.tif            # expect: unspecified
#     And the test that matters -- M1 and M2 no longer happen at the tip:
#     in a fresh clone in scratch, copy a real mask into masks/Arctic/, then
#     `git stash`, `git checkout -- .`, and pull on on_s3.  The file must survive
#     every step.  Switching to main is part of the test only after L7.

# L9. [RESIDUAL, ACCEPTED -- QL2 (a)]
#     Old commits and tags still hold pointer blobs, so M2 still happens when
#     one is checked out: `git checkout v2.1.0` silently replaces real .db
#     files, EGM2008 and RGT_001.
#     RECOMMEND: before checking out an old ref, keep real masks outside the
#     work tree (MAAP.txt already reads from the bucket mount), or copy them
#     aside first.
#     QUESTION, not verified: whether GitHub keeps (and bills for) the LFS
#     objects after the pointers are gone, and whether they can be deleted
#     without deleting the repository.  Check GitHub's documentation before
#     acting on it.

# L10. [ADE] [READY after L4]  Update the docs and notes that describe LFS.
#     - Transition_to_maap.md "OTHER": close the item, pointing to L4's commit.
#     - howto_MAAP_staging.sh S2: "Do NOT use the repo's git-lfs copies" --
#       once L4 lands there are no such copies; reword it.
#     - the Claude memory note "masks-come-from-zenodo-not-lfs" says masks/ is
#       lfs-tracked; revise it.
#     - Transition_to_maap.md, the "NOT on Zenodo, so still 134-byte git-lfs
#       POINTERS" note: record that the twelve were abandoned (QL1).
#     - RECOMMEND, confirm before doing it: delete the twelve pointer files from
#       s3://maap-ops-workspace/ben_smith/ATL1415/masks/ (F5's list).  Each one
#       is a trap -- the right name, no content -- as Transition already says
#       of the tide .h5.  Use `aws s3 rm` per file, never --recursive.
