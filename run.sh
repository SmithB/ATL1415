#!/usr/bin/env bash
# DPS run command for the ATL1415 per-tile solve.
#
# Registered as run_command: ATL1415/run.sh.  DPS invokes
#   /app/dps_wrapper.sh '/app/<repo>/run.sh' <x0> <y0> <step>
# from the job's working directory, where `input/` holds the localized `file`
# inputs and `output/` is what gets uploaded when the job finishes.
#
# Usage: run.sh <x0> <y0> <step>
#   x0, y0  tile center, in meters (polar stereographic; may be negative)
#   step    prelim | matched
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

# DPS passes the declared inputs as positional arguments.  algorithm_config.yml
# lists the positionals before the file input, but which group DPS emits first is
# one of the things the sandbox smoke test still has to settle, so tolerate a
# leading localized-path argument rather than mis-reading it as x0.
while [ "$#" -gt 0 ] && ! [[ $1 =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; do
    echo "run.sh: skipping non-numeric leading argument '$1'"
    shift
done

if [ "$#" -lt 3 ]; then
    echo "usage: run.sh <x0> <y0> <prelim|matched>" >&2
    exit 2
fi
x0=$1
y0=$2
step=$3

case "$step" in
    prelim|matched) ;;
    *) echo "ERROR: step must be 'prelim' or 'matched', got '${step}'" >&2; exit 2 ;;
esac

mkdir -p output

# The ~90 argparse options never become DPS parameters: the @argsfile idiom is
# entirely client-side (fromfile_prefix_chars), so the composed args file that
# setup_ATL1415_region.py writes (input_args_<REGION>.txt) is registered as a
# `file` input and DPS localizes it into input/.  Select it by extension so a
# second file input (the prelim tile set, for --matched) cannot be picked up
# by mistake.
args_file=$(find input -maxdepth 1 -type f -name '*.txt' | sort | head -1)
if [ -z "${args_file:-}" ]; then
    echo "ERROR: no *.txt args file found in input/ -- register the composed" >&2
    echo "       input_args_<REGION>.txt as a DPS 'file' input." >&2
    ls -la input 2>&1 >&2 || true
    exit 2
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
grep -v '^[[:space:]]*$' "$args_file" | sed 's/^/  arg: /'
echo "=========================================================="

run_solve () {
    conda run --no-capture-output -n "$env_name" ATL11_to_ATL15.py "$@"
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
    run_solve --THREADS="${threads}" --xy0 "$x0" "$y0" --prelim \
              "@${args_file}" --base_directory "$base_directory"

    # A tile with too little data is a normal outcome: ATL11_to_ATL15 returns 0
    # without writing a file.  Running the error calculation on it would then
    # exit 1 and mark the whole DPS job failed, which at a fan-out of thousands
    # of tiles would bury the real failures.  Stop cleanly instead.
    if [ ! -f "${base_directory}/prelim/${tile_name}" ]; then
        echo "no fit written for ${tile_name} (insufficient data); skipping error calculation"
        exit 0
    fi

    run_solve --THREADS="${threads}" --xy0 "$x0" "$y0" --prelim \
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

    run_solve --THREADS="${threads}" --matched \
              --prior_edge_include "$prior_edge_include" \
              --data_file "$prelim_file" \
              "@${args_file}" \
              --out_name "${PWD}/output/${tile_name}" \
              --base_directory "$base_directory"
fi

echo "=== tile job complete; output/ contains: ==="
find output -type f | sed 's/^/  /'
