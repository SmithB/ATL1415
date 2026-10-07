# howto_MAAP_ogc.sh -- register ATL1415 as a MAAP OGC process and verify the build.
# Reasoning and history: ~/ATL14_processing/dev/howto_long/howto_MAAP_ogc.sh.
# Register only when solver code changed (ATL1415/, run.sh, scripts/ used by
# run.sh, environment, or a dependency such as pointCollection on GitHub main).

cd ~/git_repos/ATL1415

# 1. Nothing in flight on DPS (a build mid-run splits the run across images),
#    and the checkout clean and pushed -- register_algorithm.py refuses otherwise.
git status -sb
/srv/conda/envs/notebook/bin/python register_algorithm.py --dry-run   # "push check: OK"

# 2. Register.  The build and deploy are followed on the web pages it prints.
/srv/conda/envs/notebook/bin/python register_algorithm.py

# 3. Verify: one DPS job reports the commit it actually runs.  Only
#    "VERDICT: MATCH" proves the deploy -- the CWL link and lastModifiedTime
#    change before the image does.  Want also: workspace_credentials=ok.
conda activate ATL14
scripts/maap/check_build_id.py \
    s3://maap-ops-workspace/ben_smith/ATL1415/run_args/rel006/north/GL/input_args_GL.txt \
    maap-dps-worker-16gb --expect $(git rev-parse HEAD)
# After a timeout, re-read the same job:  check_build_id.py ... --job <job_id>
