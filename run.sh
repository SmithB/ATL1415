#!/usr/bin/env bash
# DPS run command for the ATL1415 per-tile solve.
#
# Registered as run_command: ATL1415/run.sh.  Under MAAP's OGC system the
# generated CWL's baseCommand is /app/<repo>/run.sh, called with every input
# as --<name> <value>, from a working directory whose ./output* is collected
# when the job finishes.  The legacy system called it through
# /app/dps_wrapper.sh with positionals and localized file inputs into input/;
# both conventions are accepted (see ARGUMENTS below).
#
# Usage: run.sh --x0 <m> --y0 <m> --step <step> --args_file <uri|path>   (OGC)
#        run.sh <x0> <y0> <step>        (legacy DPS, local runs; args file in input/)
#   x0, y0     tile center, in meters (polar stereographic; may be negative)
#   step       prelim | matched | build_id
#   args_file  s3:// URI or local path of the composed input_args_<REGION>.txt
#
# step=build_id prints the build stamp and exits 0 without solving anything, so
# ONE cheap job says which commit the image was built from.  --build-id does the
# same from a shell.  See "BUILD ID" below.
#
# There is ONE registered algorithm rather than one per stage (see
# docs/Transition_to_maap.md): a MAAP algorithm has a single run_command, and
# the conda+SuiteSparse build is the expensive part.  Only the per-tile solve
# belongs in DPS; setup, queue-build, mosaic, to-netcdf and browse stay in the
# ADE, where the my-private-bucket mount exists.
set -euo pipefail

repo_dir=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
env_name=$(sed -n 's/^name:[[:space:]]*//p' "${repo_dir}/environment.yml" | head -1)
: "${env_name:?could not read 'name:' from environment.yml}"

# ===========================================================================
# BUILD ID -- answer "what is actually in this image?" in one job.
# ===========================================================================
# A build clones code_repository at algorithm_version and bakes the result into a
# container, so the working copy in the ADE has nothing to do with what runs on
# a worker, and no job log has ever said which commit it carries.  That made a
# stale image indistinguishable from a fix that did not work: on 2026-09-09 an
# algorithm_version that had been built once already ran the OLD code, which
# MAAP support has since confirmed is not expected behaviour.
#
# build-env.sh writes ${repo_dir}/.atl1415_build_id at BUILD time; this prints
# it and exits 0 without touching input/, the args file or the solver, so a
# single build_id job is a complete answer (scripts/maap/check_build_id.py).
#
# Checked BEFORE the non-numeric skip loop below, which would otherwise shift
# '--build-id' away as a leading non-numeric argument.
build_id_file="${repo_dir}/.atl1415_build_id"
BUILD_ID_PY='import importlib.metadata as md
try:
    import ATL1415
    print("import_ok=True")
    print("file=%s" % getattr(ATL1415, "__file__", "unknown"))
except Exception as exc:
    print("import_ok=False")
    print("error=%s: %s" % (type(exc).__name__, exc))
try:
    print("version=%s" % md.version("ATL1415"))
except Exception as exc:
    print("version=unknown (%s)" % type(exc).__name__)
try:
    print("maap_py=%s" % md.version("maap-py"))
except Exception as exc:
    print("maap_py=absent (%s)" % type(exc).__name__)
'
git_q () { git -c safe.directory='*' -C "$repo_dir" "$@" 2>/dev/null; }
config_version () {
    sed -n 's/^algorithm_version:[[:space:]]*//p' "${repo_dir}/algorithm_config.yml" 2>/dev/null | head -1
}
# Reads one key out of the stamp, and succeeds (printing nothing) when there is
# no stamp to read.  It has to succeed: under `set -e` a bare
# `x=$(sed ... missing_file)` exits the script with sed's status 2, which is
# how the first version of this flag exited 2 after printing everything but its
# summary line.
stamp_field () {
    [ -f "$build_id_file" ] || return 0
    sed -n "s/^$1=//p" "$build_id_file" | head -1
}

# THE ONE-LINE SUMMARY, shared by --build-id and the header of EVERY tile job.
# Every field but the last two comes from the stamp, so the line describes one
# build rather than pairing a build-time commit with a run-time config (the
# first version did that, printing e231ba4 beside a version it never had);
# only with no stamp at all does it fall back to live.  maap_py (the argument)
# and maap_pgt are run-time facts about this machine.
build_id_summary () {
    local commit built version pgt
    commit=$(stamp_field commit)
    built=$(stamp_field build_completed)
    version=$(stamp_field algorithm_version)
    if [ -n "${MAAP_PGT:-}" ]; then pgt=set; else pgt=unset; fi
    echo "BUILD_ID: commit=${commit:-unknown} built=${built:-INCOMPLETE_OR_ABSENT} algorithm_version=${version:-$(config_version)} maap_py=${1:-unknown} maap_pgt=${pgt}"
}

print_build_id () {
    echo "=========================================================="
    echo "  ATL1415 build id"
    echo "=========================================================="
    echo "repo dir    : ${repo_dir}"

    if [ -f "$build_id_file" ]; then
        echo "--- build stamp (written by build-env.sh at build time) ---"
        sed 's/^/  /' "$build_id_file"
    else
        # An image built before this flag existed, or a build that never got
        # past its first step.  Say which, rather than printing nothing.
        echo "--- NO BUILD STAMP at ${build_id_file} ---"
        echo "  This image predates the build-id flag, or build-env.sh did not"
        echo "  reach its first step.  Falling back to live git below."
    fi

    # Cross-check: what the clone in the image says NOW.  It should agree with
    # the stamp; a disagreement means the image was modified after its build.
    echo "--- live git in the image ---"
    live_commit=
    if git_q rev-parse --git-dir >/dev/null; then
        live_commit=$(git_q rev-parse HEAD || echo unknown)
        echo "  live_commit=${live_commit}"
        echo "  live_ref=$(git_q describe --all --always HEAD || echo unknown)"
    else
        echo "  (no git metadata in the image -- the stamp above is the only record)"
    fi
    # Say it, rather than leaving the reader to diff two 40-character hashes.
    stamp_commit=$(stamp_field commit)
    if [ -n "$stamp_commit" ] && [ -n "$live_commit" ] && [ "$stamp_commit" != "$live_commit" ]; then
        echo "  WARNING: STAMP AND LIVE GIT DISAGREE -- the tree moved after it was built."
        echo "           The stamp is what build-env.sh built; on DPS (pip install .) that is what runs."
    fi

    # What this image THINKS it was registered as.  If it disagrees with the
    # algorithm_version you submitted to, the image is not the one you meant.
    echo "  image_algorithm_version=$(config_version)"

    # The installed copy, which is what actually runs: build-env.sh does
    # `pip install .`, so the console scripts run a COPY of the repo, not the
    # repo itself.  A mismatch here is the stale-code failure in miniature.
    # Guarded -- a broken env must not stop the stamp above from being printed.
    echo "--- installed ATL1415 and maap-py (the copies the solve imports) ---"
    py_out=$(conda run --no-capture-output -n "$env_name" python -c "$BUILD_ID_PY" 2>&1) \
        || py_out="${py_out}
(could not run python in conda env '${env_name}')"
    printf '%s\n' "$py_out" | sed 's/^/  /'
    maap_py=$(printf '%s\n' "$py_out" | sed -n 's/^maap_py=//p' | head -1)

    # MAAP_PGT decides whether this machine can read NSIDC at all.  A worker
    # has no Earthdata login of its own; pointCollection gets the DAAC's
    # temporary S3 credentials from maap.aws.earthdata_s3_credentials(), and
    # ONLY TRIES when MAAP_PGT is set -- it treats its absence as "not on
    # MAAP", says nothing, and falls back to earthaccess, which has no
    # credentials on a worker either.  So a missing MAAP_PGT surfaces much
    # later as an ATL11 read failure that looks like a permissions or data
    # error.  Legacy workers set it; whether OGC workers do is unverified
    # (howto_MAAP_ogc QE), and this answers it without reading any ATL11.
    # The VALUE is a credential and is never printed -- only whether it is set.
    echo "--- MAAP credentials on this machine ---"
    if [ -n "${MAAP_PGT:-}" ]; then maap_pgt=set; else maap_pgt=unset; fi
    echo "  MAAP_PGT=${maap_pgt}   (value never printed)"
    echo "  MAAP_API_HOST=${MAAP_API_HOST:-<unset: maap-py defaults to api.maap-project.org>}"
    if [ "$maap_pgt" = unset ]; then
        echo "  WARNING: MAAP_PGT IS NOT SET -- NSIDC credentials will not be brokered"
        echo "           here, and every ATL11 read from NSIDC will fail."
    fi

    # ONE greppable line, so a collector need not parse the block above.
    echo "=========================================================="
    build_id_summary "${maap_py:-unknown}"
    echo "=========================================================="
}

# OUTPUT/ IS NOT OPTIONAL, EVEN HERE.  The generated CWL collects the job's
# products with `glob: ./output*`, and a job with nothing matching is a
# permanentFail -- which is exactly how the first OGC build_id job ended
# (db93c7f3..., 2026-09-10): it printed a complete, correct report and was
# then failed with "Did not find output file with glob pattern:
# ['./output*']".  So the report also goes to output/build_id.txt: that
# satisfies the glob, and makes the answer a product uploaded with the job,
# readable from its output prefix as well as from the log.
build_id_and_exit () {
    mkdir -p output
    print_build_id | tee output/build_id.txt
    exit 0
}

for arg in "$@"; do
    case "$arg" in
        --build-id|--build_id|build_id|build-id) build_id_and_exit ;;
    esac
done

# ===========================================================================
# ARGUMENTS -- two conventions, because there are two callers.
# ===========================================================================
# OGC (since 2026-09-10, docs/howto_MAAP_ogc.sh O2): the generated CWL binds
# every input as a PREFIXED option, so a job arrives as
#     run.sh --x0 220000 --y0 20000 --step prelim --args_file s3://.../input_args_AA.txt
# LEGACY, and every local run: positionals, with the args file found in input/
#     run.sh 220000 20000 prelim
# Any --x0/--y0/--step/--args_file anywhere in argv selects the first.
x0= ; y0= ; step= ; args_src=
prefixed=false
for arg in "$@"; do
    case "$arg" in
        --x0|--x0=*|--y0|--y0=*|--step|--step=*|--args_file|--args_file=*) prefixed=true ;;
    esac
done

if $prefixed; then
    while [ "$#" -gt 0 ]; do
        case "$1" in
            --x0|--y0|--step|--args_file)
                if [ "$#" -lt 2 ]; then
                    echo "ERROR: $1 needs a value" >&2; exit 2
                fi
                # $2 is taken verbatim, so a negative coordinate is a value,
                # not mistaken for another option.
                name=${1#--}; value=$2; shift 2 ;;
            --x0=*|--y0=*|--step=*|--args_file=*)
                name=${1%%=*}; name=${name#--}; value=${1#*=}; shift ;;
            *)
                echo "run.sh: ignoring unexpected argument '$1'"; shift; continue ;;
        esac
        case "$name" in
            x0) x0=$value ;;
            y0) y0=$value ;;
            step) step=$value ;;
            args_file) args_src=$value ;;
        esac
    done
else
    # Legacy DPS passed the declared inputs as positionals, and which group it
    # emitted first was never settled, so tolerate a leading localized-path
    # argument rather than mis-reading it as x0.
    while [ "$#" -gt 0 ] && ! [[ $1 =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; do
        echo "run.sh: skipping non-numeric leading argument '$1'"
        shift
    done
    if [ "$#" -ge 3 ]; then
        x0=$1; y0=$2; step=$3
    fi
fi

if [ -z "$x0" ] || [ -z "$y0" ] || [ -z "$step" ]; then
    echo "usage: run.sh --x0 <m> --y0 <m> --step <prelim|matched> --args_file <uri|path>" >&2
    echo "       run.sh <x0> <y0> <prelim|matched>      (args file found in input/)" >&2
    echo "       run.sh --build-id" >&2
    exit 2
fi
for v in "$x0" "$y0"; do
    if ! [[ $v =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
        echo "ERROR: tile center must be numeric meters, got '${v}'" >&2; exit 2
    fi
done

case "$step" in
    prelim|matched) ;;
    # The pre-scan above catches build_id as a bare token; this catches every
    # other spelling that parses to it, e.g. --step=build_id, which the first
    # version rejected with "must be ... 'build_id', got 'build_id'".
    build_id|build-id) build_id_and_exit ;;
    *) echo "ERROR: step must be 'prelim', 'matched' or 'build_id', got '${step}'" >&2; exit 2 ;;
esac

mkdir -p output

# ---------------------------------------------------------------------------
# THE ARGS FILE.  The ~90 argparse options never become DPS parameters: the
# @argsfile idiom is entirely client-side (fromfile_prefix_chars), so the
# composed input_args_<REGION>.txt is one input, read here and passed as @path.
# ---------------------------------------------------------------------------
if [ -n "$args_src" ]; then
    case "$args_src" in
        s3://*)
            # A STRING input that run.sh fetches, not a CWL File the runner
            # stages -- whether the runner can stage an s3:// File is untested
            # (howto_MAAP_ogc QA).  s3fs rather than the aws CLI: build-env.sh
            # proves s3fs importable, and nothing guarantees an aws binary on
            # maap_base.  It uses the worker's own credential chain, as every
            # other bucket read in the solve does.
            mkdir -p input
            args_file="${PWD}/input/$(basename "$args_src")"
            echo "run.sh: fetching ${args_src} -> ${args_file}"
            conda run --no-capture-output -n "$env_name" python -c \
                'import sys, s3fs; s3fs.S3FileSystem().get(sys.argv[1], sys.argv[2])' \
                "$args_src" "$args_file"
            ;;
        *)
            args_file=$args_src ;;
    esac
    if [ ! -f "$args_file" ]; then
        echo "ERROR: args file '${args_src}' is not a readable file here." >&2
        exit 2
    fi
else
    # LEGACY: DPS localized the `file` input into input/.  Select it by
    # extension so a second file input (the prelim tile set, for --matched)
    # cannot be picked up by mistake.
    #
    # -L IS LOAD-BEARING.  Legacy DPS did not copy a localized input into
    # input/ -- it SYMLINKED it into a shared cache, e.g.
    #   input_args_IS.txt -> /data/work/cache/5/e/4/b/<md5>/input_args_IS.txt
    # and `find -type f` tests the LINK, which is -type l, so without -L this
    # matches nothing and the job dies in the guard below on a file that
    # localized perfectly.  That is exactly how the first smoke test failed
    # (job 61589c00-7a67-4881-93af-b96d0f0e8c4b, 2026-09-08).  -L follows the
    # link, so a symlink to a regular file tests as -type f, and a dangling
    # one is correctly still skipped.
    # `|| true`: with no input/ at all, find fails, pipefail fails the pipeline,
    # and set -e would end the job right here with status 1 and no message --
    # instead of in the guard below, which says what is missing.
    args_file=$(find -L input -maxdepth 1 -type f -name '*.txt' 2>/dev/null | sort | head -1) || true
    if [ -z "${args_file:-}" ]; then
        echo "ERROR: no --args_file given, and no *.txt args file found in input/." >&2
        ls -la input 2>&1 >&2 || true
        exit 2
    fi
fi
args_file=$(readlink -f "$args_file")

# One thread per available core.  ATL11_to_ATL15 sets MKL/OPENBLAS/NUMEXPR/OMP
# from a --THREADS= it scrapes out of sys.argv at import time -- it does not look
# inside the @argsfile -- so this has to be an explicit command-line argument to
# take effect.  It goes BEFORE @${args_file} so the args file can still override.
threads=${ATL1415_THREADS:-$(nproc)}

echo "=========================================================="
echo "  ATL1415 DPS tile job"
echo "  step        : ${step}"
echo "  xy0         : ${x0} ${y0}"
echo "  args file   : ${args_file}"
echo "  threads     : ${threads}"
echo "  conda env   : ${env_name}"
echo "  working dir : ${PWD}"
echo "=========================================================="
# EVERY TILE RECORDS ITS OWN BUILD, so a run can be audited after the fact
# (howto_MAAP_ogc O11): a build_id job only vouches for the worker it lands
# on.  collect_jobs.py reports it per tile.  maap_py=unchecked: reading it
# costs a conda start, and the build_id job reports it.
build_id_summary unchecked
echo "=========================================================="
# ...and hands the same facts to the solver, which writes them into the tile's
# /meta as build_commit / build_version / build_completed (howto_MAAP_ogc
# O12a) -- the record of which build made a tile that outlives this log.
export ATL1415_BUILD_COMMIT="$(stamp_field commit)"
export ATL1415_BUILD_VERSION="$(stamp_field algorithm_version)"
export ATL1415_BUILD_COMPLETED="$(stamp_field build_completed)"
grep -v '^[[:space:]]*$' "$args_file" | sed 's/^/  arg: /'
echo "=========================================================="

# Every solve is wrapped so the job reports its own peak memory.  DPS will not
# tell us: get_job_metrics() returns null memory fields, and on the legacy path it came back
# an empty dict for the first tile that succeeded (fdc4d767, 2026-09-08), and
# retrieve_attributes() populated only `status`.  Sizing the production queue
# needs peak RSS per tile, so the job measures itself.  The wrapper is
# transparent -- it exits with the solver's own status -- and costs nothing.
run_solve () {
    "${repo_dir}/scripts/run_with_rusage.py" "$1" \
        conda run --no-capture-output -n "$env_name" ATL11_to_ATL15.py "${@:2}"
}

if [ "$step" = "prelim" ]; then
    # Fit, then the error-calculation companion, mirroring the single queue line
    # that make_ATL1415_queue.py writes for SLURM.  Tiles land in output/prelim/,
    # which is what DPS uploads.
    base_directory="${PWD}/output"
    tile_name=$(awk -v x="$x0" -v y="$y0" 'BEGIN{printf "E%d_N%d.h5", int(x/1000), int(y/1000)}')

    # --base_directory goes AFTER @${args_file}, unlike --THREADS: argparse takes
    # the last occurrence, and setup_ATL1415_region.py writes '-b=<region_dir>'
    # (the same dest as --base_directory) as the final line of the composed args
    # file.  Passed before the args file it would be overridden by that ADE path,
    # which does not exist on a worker.
    run_solve fit --THREADS="${threads}" --xy0 "$x0" "$y0" --prelim \
              "@${args_file}" --base_directory "$base_directory"

    # A tile with too little data is a normal outcome: ATL11_to_ATL15 returns 0
    # without writing a file.  Running the error calculation on it would then
    # exit 1 and mark the whole DPS job failed, which at a fan-out of thousands
    # of tiles would bury the real failures.  Stop cleanly instead.
    if [ ! -f "${base_directory}/prelim/${tile_name}" ]; then
        echo "no fit written for ${tile_name} (insufficient data); skipping error calculation"
        exit 0
    fi

    run_solve error --THREADS="${threads}" --xy0 "$x0" "$y0" --prelim \
              "@${args_file}" --base_directory "$base_directory" --calc_error_for_xy
else
    # --matched reads the tile's own prelim fit AND its neighbours', through
    # prior_edge_include, so a matched job needs the surrounding prelim tiles
    # localized into input/prelim/ (the tile itself plus its 8 neighbours at
    # minimum).  base_directory therefore points at input/, not output/: that is
    # where ATL11_to_ATL15 looks for <base>/prelim/E*_N*.h5.  Only the result is
    # written to output/.
    if [ ! -d input/prelim ]; then
        echo "ERROR: --matched needs the prelim tiles for this tile and its" >&2
        echo "       neighbours localized into input/prelim/ ." >&2
        exit 2
    fi
    base_directory="${PWD}/input"
    # Same name ATL11_to_ATL15 builds: 'E%d_N%d.h5' % (x0/1e3, y0/1e3), i.e.
    # kilometers truncated toward zero.  awk int() truncates the same way.
    tile_name=$(awk -v x="$x0" -v y="$y0" 'BEGIN{printf "E%d_N%d.h5", int(x/1000), int(y/1000)}')
    prelim_file="${base_directory}/prelim/${tile_name}"
    if [ ! -f "$prelim_file" ]; then
        echo "ERROR: prelim tile ${prelim_file} not found; localized files are:" >&2
        ls -la input/prelim >&2 || true
        exit 2
    fi

    # make_ATL1415_queue.py passes --prior_edge_include on every matched line
    # (default 1000); ATL11_to_ATL15's own default is None, which silently drops
    # the prior-edge constraints, so pass it here too.  Before @${args_file}, so
    # an args file that sets it still wins.
    prior_edge_include=${ATL1415_PRIOR_EDGE_INCLUDE:-1000}

    run_solve matched --THREADS="${threads}" --matched \
              --prior_edge_include "$prior_edge_include" \
              --data_file "$prelim_file" \
              "@${args_file}" \
              --out_name "${PWD}/output/${tile_name}" \
              --base_directory "$base_directory"
fi

echo "=== tile job complete; output/ contains: ==="
find output -type f | sed 's/^/  /'
