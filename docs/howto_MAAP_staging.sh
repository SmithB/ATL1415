# howto_MAAP_staging.sh -- one-time setup of the MAAP ADE and the bucket.
# MAAP/DPS counterpart of the discover setup.  Reasoning and history:
# ~/ATL14_processing/dev/howto_long/howto_MAAP_staging.sh (on the MAAP ADE).
# Bucket root s3://maap-ops-workspace/ben_smith is mounted read-write at
# ~/my-private-bucket (mountpoint-s3: no rename; reorganize with `aws s3 mv`).

# 1. The ATL14 conda env, in the persistent home.  [ADE]
CONDA_ENVS_PATH=/home/jovyan/.conda/envs bash build-env.sh
conda activate ATL14
python -c "import pointCollection, LSsurf, sparseqr; print('env ok')"
mkdir -p /home/jovyan/ATL14_processing      # = --ATL14_root in default_args/MAAP_dps.txt
conda run -n ATL14 python -m pip install -e . --no-deps --no-build-isolation
# Home is an NFS volume with a 150 GB quota (df does not show it): keep DPS
# output on the bucket, never in /home.  The ADE itself is capped at ~7.3 GiB.

# 2. Masks.  [ADE]
# Copy from Zenodo record 22259649 (doi:10.5281/zenodo.22259649) to
# s3://maap-ops-workspace/ben_smith/ATL1415/masks/{Arctic,Antarctic}/.
# Never fetch masks with git-lfs.

# 3. ATL11 per-granule geoIndex.  [DISCOVER -> ADE]
# Build on discover (howto_ATL11.sh); copy to
# s3://maap-ops-workspace/ben_smith/ATL11_index/ATL11_index_<cycles>_<release>_<version>/
# FLAT: both hemispheres in one directory.  The generation must be one CMR
# still serves (0332_007_05 as of 2026-09-17).

# 4. Check that each staged input is readable.  [ADE]
python - <<'EOF'
import pointCollection as pc
root = 's3://maap-ops-workspace/ben_smith'
fs = pc.io_utils.get_s3fs(daac=None)
print(len(fs.ls(f'{root}/ATL11_index/ATL11_index_0332_007_05/')), 'index files')
print(fs.ls(f'{root}/ATL1415/masks/Arctic/')[:3])
print(pc.io_utils.get_s3fs(daac=None, anon=True).ls('s3://pytmd')[:3])
EOF

# 5. Register the algorithm and verify the build: howto_MAAP_ogc.sh.
