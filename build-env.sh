#!/usr/bin/env bash
# DPS build command for the ATL1415 per-tile solve.
#
# Registered as build_command: ATL1415/build-env.sh  (DPS clones this repo to
# /app/<repo>/ and paths in algorithm_config.yml are relative to /app).
#
# The base image (maap_base) is debian:11 + git + Miniforge and nothing else, so
# the environment has to be built here.  It MUST be conda-based: LSsurf pulls
# sparseqr (PySPQR), which compiles CFFI bindings against SuiteSparseQR and needs
# the SuiteSparse dev headers at pip time.  No MAAP base image ships SuiteSparse.
set -euo pipefail

repo_dir=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
cd "$repo_dir"

# environment.yml declares `name: ATL14`; keep the two in sync with run.sh.
env_name=$(sed -n 's/^name:[[:space:]]*//p' environment.yml | head -1)
: "${env_name:?could not read 'name:' from environment.yml}"

echo "=== ATL1415 DPS build: env '${env_name}' from ${repo_dir} ==="
conda --version

# 1. conda first.  This is what supplies suitesparse (for PySPQR) and gdal
#    (conda-forge's gdal installs the python bindings with a .dist-info, so pip
#    treats LSsurf's bare `gdal` requirement as already satisfied instead of
#    trying to build it from source against the wrong libgdal).
conda env update --name "$env_name" --file environment.yml --prune

env_prefix=$(conda run -n "$env_name" python -c 'import sys; print(sys.prefix)')
echo "=== conda env prefix: ${env_prefix} ==="

# 2. Point the compiler at the conda env before pip builds anything.  PySPQR's
#    CFFI build looks for SuiteSparseQR.hpp / libspqr in the default search
#    paths, which do not include $CONDA_PREFIX.
export CPATH="${env_prefix}/include:${env_prefix}/include/suitesparse${CPATH:+:$CPATH}"
export LIBRARY_PATH="${env_prefix}/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
export LD_LIBRARY_PATH="${env_prefix}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# Fail early and legibly if the conda step did not actually deliver the headers.
if [ ! -f "${env_prefix}/include/suitesparse/SuiteSparseQR.hpp" ] \
   && [ ! -f "${env_prefix}/include/SuiteSparseQR.hpp" ]; then
    echo "ERROR: SuiteSparseQR headers not found under ${env_prefix}/include." >&2
    echo "       sparseqr (PySPQR) cannot build; check the suitesparse entry in environment.yml." >&2
    exit 1
fi

# 2a. Same for the compilers.  The base image has no toolchain at all, so the
#     c-compiler/cxx-compiler entries in environment.yml are the only source of
#     one; without them pip gets all the way to the compile step and then dies
#     with "No such file or directory: 'gcc'" on LSsurf, sparseqr and cartopy,
#     which is exactly how the first DPS build failed.
#
#     Nothing needs to set CC/CXX.  These packages ship no activate.d scripts,
#     but they do put plain gcc/g++ in the env's bin, and `conda run` puts that
#     bin on PATH -- so distutils resolves the bare 'gcc' from its own sysconfig
#     against the conda toolchain, keeping the -B python_compiler_compat flag
#     that the conda interpreter expects.  Overriding CC would drop it.
for tool in gcc g++; do
    if [ ! -x "${env_prefix}/bin/${tool}" ]; then
        echo "ERROR: ${tool} not found at ${env_prefix}/bin/${tool}." >&2
        echo "       LSsurf (Cython) and sparseqr (CFFI) cannot compile; check the" >&2
        echo "       c-compiler / cxx-compiler entries in environment.yml." >&2
        exit 1
    fi
done
echo "=== toolchain: $(${env_prefix}/bin/gcc --version | head -1) ==="

# 3. pip, inside the conda env.  This pulls pointCollection[cloud] (earthaccess,
#    s3fs, fsspec -- required by the ATL11 cloud read path) and LSsurf, both as
#    git+ URLs, so the build container must be able to reach github.com.
#    SMBcorr is deliberately NOT installed: it is only needed for --firn_model,
#    which no production rel_006 string and no per-tile solve uses.  Add the
#    [firn] extra here if that ever changes.
conda run --no-capture-output -n "$env_name" python -m pip install --no-cache-dir .

# 4. Prove the pieces that matter are importable, so a broken build fails at
#    build time rather than on every one of thousands of tile jobs.
conda run --no-capture-output -n "$env_name" python - <<'PYEOF'
import importlib
missing = []
for mod in ("numpy", "scipy", "h5py", "osgeo.gdal", "pyproj",
            "sparseqr", "pointCollection", "LSsurf", "ATL1415",
            "earthaccess", "s3fs"):
    try:
        importlib.import_module(mod)
        print(f"  OK   {mod}")
    except Exception as exc:
        print(f"  FAIL {mod}: {type(exc).__name__}: {exc}")
        missing.append(mod)
if missing:
    raise SystemExit("build verification failed: " + ", ".join(missing))
PYEOF

conda run --no-capture-output -n "$env_name" which ATL11_to_ATL15.py
echo "=== ATL1415 DPS build complete ==="
