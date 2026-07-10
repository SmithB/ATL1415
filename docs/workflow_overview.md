# ATL1415 Processing Workflow

This document summarizes the end-to-end pipeline used to generate the ATL14
(reference DEM) and ATL15 (surface-height-change) products, as encoded in the
driver scripts `docs/howto_AA.sh` (Antarctica), `docs/howto_arctic.sh` /
`docs/howto_GL.sh` (Arctic regions). It is meant as a map of the moving parts,
not a replacement for the ATBD.

## Big picture

Processing happens on a per-**tile** basis (100 km x 100 km, or other
`tile_spacing`/width combos), in two passes:

1. **`prelim`** — initial least-squares surface fit (z0 + dz/dt at several
   scales) to all ATL11 data in a tile.
2. **`matched`** — re-fit using crossover-matched constraints and the prior
   `prelim` tiles, so tile edges agree with their neighbors.

Per-tile HDF5 outputs are then **mosaicked** into seamless regional grids,
optionally via an intermediate 200 km tiling step for large domains (e.g.
Antarctica), and finally converted to NASA-format **NetCDF** products
(`ATL14_*.nc`, `ATL15_*.nc`).

Every heavy step is farmed out as a SLURM job array built from a flat text
"queue file" of shell commands, one per tile/task.

## Configuration: layered `default_args/*.txt` files

Scripts take CLI args, but also accept `@path/to/file.txt` argparse syntax to
load arguments from a file. The repo composes config in layers, each
contributing `--key=value` lines, combined positionally:

- **Location** (`discover.txt`): cluster paths — `--mask_dir`,
  `--tide_directory`, `--ATL14_root`.
- **Release** (`latest_release.txt`, a symlink to the current release file, e.g. `rel_005_0329.txt`): science/release parameters — tile
  spacing/width, time range (`-t`), error covariance terms
  (`--E_d3zdx2dt`, `--E_d2z0dx2`, `--E_d2zdt2`), `--cycles`, `--Release`,
  `--version`, `--ATL11_release`, DEM/geoid tolerances, bias parameters.
- **Hemisphere**: passed directly as `--Hemisphere=1` (north) or `--Hemisphere=-1`
  (south) on the `setup_ATL1415_region.py` command line. `--ATL11_index` is
  derived from `--ATL11_release` (in the release file) and hemisphere;
  `--ATL11_xover_dir` is derived as the `xover/` sibling of the index directory.
- **Period** (`quarterly.txt` / monthly variant): `--dzdt_lags`, output grid
  spacing (`-g`).
- **Region** (`GL_0329.txt`, `AA_0329.txt`, etc.): `--region`, mask/tide
  files specific to that ice sheet/region (e.g. `--mask_file`,
  `--tide_mask_file`, `--tide_model`, error-scaling grids).

`setup_ATL1415_region.py` merges all of these into one composite
`input_args_<REGION>.txt`, which every later step (`make_ATL1415_queue.py`,
`ATL14_write2nc.py`, `ATL15_write2nc.py`, ...) consumes via `@input_args_*.txt`.

## Step-by-step

### 1. Region setup — `setup_ATL1415_region.py`
```
setup_ATL1415_region.py discover.txt latest_release.txt GL_0329.txt --Hemisphere=1
```
Merges the layered defaults files, resolves relative mask/index paths
against `--ATL14_root`/`--mask_dir`, creates the
`<ATL14_root>/rel<NNN>/<north|south>/<REGION>/` directory tree, and writes
the composite `input_args_<REGION>.txt`.

### 2. Tile queue generation — `make_ATL1415_queue.py`
```
make_ATL1415_queue.py prelim <region_dir>/input_args_GL.txt
```
Reads the region's mask raster, finds the valid tile-center grid (aligned
to half the tile width), and writes one `ATL11_to_ATL15.py ...` command per
tile (plus, by default, a second `--calc_error_for_xy` command) into a queue
file, e.g. `1415_queue_GL_prelim.txt`. Step argument selects the processing
stage:
- `prelim` — full tile grid.
- `centers` / `edges` / `corners` — diagnostic single-tile / offset subsets.
- `matched` — re-reads tiles produced by `prelim`, adds
  `--matched`/`--prior_edge_include` so edges are constrained by neighbors.

Supports bounding by `--min_xy`/`--max_xy` (used to split Antarctica into
north/south halves with different SLURM time limits) or explicit
`--xy_list_file`/`--tile_list_file`.

### 3. SLURM job setup — `setup_slurm_run.py`
```
setup_slurm_run.py --run_name GL_prelim -q 1415_queue_GL_prelim.txt --time 04:00:00 -j 7 -e ATL14
```
Converts the queue file into a SLURM-ready run directory: per-task shell
scripts under `queue/`, `running/`/`done/`/`logs/` bookkeeping
subdirectories, and a `slurm_run.sh` (or chunked `slurm_run_part*.sh`/
`slurm_run_*cpu.sh` if the queue is large). If `-j`/`--jobs_per_task` isn't
given explicitly, CPU count per task is chosen automatically from tile
radius (more CPUs for the costlier near-pole tiles). The user submits with
`sbatch slurm_run.sh`.

### 4. Per-tile fit — `ATL1415/ATL11_to_ATL15.py` (run by the queue, not called directly)
Reads ATL11 data for one tile (`--xy0 X Y`), applies tide/DEM/bias
corrections, fits the surface (z0) and elevation-change fields (dz, dzdt at
multiple lags and averaging scales) by least squares, and writes one HDF5
file per tile (e.g. `E480_N-1360.h5`) under `prelim/` or `matched/`.

### 5. Sanity check — `make_field_size_report.py`
```
make_field_size_report.py <region_dir>/prelim GL prelim
```
Summarizes field array shapes for each tile HDF5 into JSON reports, used
with `notebooks/check_tiles.ipynb` to spot malformed/short tiles before
mosaicking.

### 6. Repeat for `matched` step
Same `make_ATL1415_queue.py matched ...` → `setup_slurm_run.py` →
`sbatch` sequence, now using the `matched` step and prior `prelim` outputs.

### 7. Mosaicking 100 km tiles into a regional grid
Two paths, chosen by region size:

- **Direct** (`scripts/make_mosaic_jobs`, or `make_mosaic_jobs.py`): builds a
  SLURM job that runs `make_mosaic.py` (from the `pointCollection`
  dependency) directly over all 100 km tiles per field group (z0, dz,
  `dzdt_lag<N>`, `avg_dz_<scale>`, ...), producing regional HDF5s like
  `z0.h5`, `dz.h5`, `dzdt_lag1.h5`.
- **Two-stage, for large domains like Antarctica** (`make_200km_tiles.py`
  then `make_200km_to_mosaic_jobs(.py|.sh)`): first blends 100 km tiles into
  intermediate 200 km tiles (with computed pad/feather based on tile
  overlap), written under `200km_tiles/<field>/`, then a second mosaic pass
  combines those 200 km tiles into the final regional HDF5s. Sigma
  (uncertainty) fields are always sourced from `prelim/` tiles for
  consistency, even when the regular fields come from `matched/`.

Antarctica's `howto_AA.sh` additionally splits the continent into
`min_xy`/`max_xy`-bounded north/south halves and four `A1..A4` sectors for
mosaicking, due to its size; see `notebooks/make_and_check_links.ipynb` for
verifying that the sector tiles link up correctly before the final mosaic.

### 8. NetCDF conversion
```
ATL14_write2nc.py @<region_dir>/input_args_GL.txt
ATL15_write2nc.py @<region_dir>/input_args_GL.txt
```
- `ATL14_write2nc.py` turns the regional `z0.h5` mosaic into
  `ATL14_<region>_<cycles>_100m_<release>_<version>.nc` (CF-compliant,
  polar-stereographic projection, ATL14 metadata template).
- `ATL15_write2nc.py` turns `dz.h5`/`dzdt_lag<N>.h5`/averaged-scale HDF5s
  into `ATL15_<region>_<cycles>_<xxkm>_<release>_<version>.nc`, one file per
  averaging scale (1, 10, 20, 40 km by default), with `delta_h` and
  `dhdt_<N>mo` groups. The same script handles monthly vs. quarterly output
  (it picks the cadence from `--delta_t`/`--grid_spacing`); the standalone
  `ATL15_write2nc_monthly.py` script was merged into `ATL15_write2nc.py`.

Both scripts are run as SLURM jobs themselves for the Arctic/Antarctic
regions, queued by `scripts/run_arctic_to_nc.sh` / `scripts/run_antarctic_tonc.sh`.

### 9. Error recovery
The slurm job template (`ATL1415/resources/slurm_templates/packable_job.txt`,
used by `setup_slurm_run.py`, `make_mosaic_jobs.py`, and
`make_200km_to_mosaic_jobs.py`) writes each task's log to `active_logs/`
while it runs, then moves it to `logs/` once the task script returns. A
task killed by the scheduler (time/memory limit, node failure, `scancel`)
never reaches that step, so its log is left behind in `active_logs/` and
its task file is left behind in `running/`. A task that completes but
exits nonzero additionally has its log copied to `error_logs/`. These two
directories distinguish "killed" from "completed but failed" without
having to grep log text for words like "killed" or "error".

`slurm_run_status.py <run_dir>` reports queue/running/done counts plus
error counts (grouped by exception type, parsed from tracebacks in
`error_logs/`) and killed-task counts (from `active_logs/`), and can
requeue either or both groups by moving their task files straight back
into `queue/`:
```
slurm_run_status.py GL_prelim
slurm_run_status.py GL_prelim --requeue both --dry_run
slurm_run_status.py GL_prelim --requeue both
```
`howto_arctic.sh` shows the typical check after each stage:
```
for j in RA IS CN CS SV; do echo $j; slurm_run_status.py $j"_prelim"; done
```

The older `requeue_slurm_errors.py`/`requeue_slurm_killed.py`/
`requeue_slurm_running.py` scripts (which classified failures by grepping
log text) have been deprecated in favor of `slurm_run_status.py` and moved
to `ATL1415/scripts/old/`; they are no longer installed as console scripts.

## Region-specific wrapper scripts (Arctic)

`scripts/run_arctic_prelim.sh`, `run_arctic_matched.sh`, `run_arctic_mosaic.sh`,
and `run_arctic_to_nc.sh` each take `(release_file, location_file,
period_file, [region list])` and loop over the Arctic regions
(`RA IS CN CS SV` by default — Russian Arctic, Iceland, Canada North/South,
Svalbard) running the corresponding setup/queue/slurm steps for each region,
submitting jobs with `sbatch` at the end of each stage. They detect monthly
vs. quarterly processing by grepping the period file for "monthly" and
adjust output paths/scripts accordingly (`*_monthly` suffix,
`make_mosaic_jobs_monthly`; `ATL15_write2nc.py` itself picks monthly vs.
quarterly from `--delta_t`/`--grid_spacing`, not from a separate script).

Antarctica has no single all-in-one wrapper; `docs/howto_AA.sh` runs the
setup/queue/slurm/mosaic steps directly, split into north/south by
`--min_xy`/`--max_xy` and into four sectors (A1-A4) for final mosaicking.

## Typical command sequence (Greenland, from `docs/howto_GL.sh`)

```bash
setup_ATL1415_region.py default_args/discover.txt default_args/latest_release.txt default_args/GL_0329.txt --Hemisphere=1

make_ATL1415_queue.py prelim <ATL14_root>/rel005/north/GL/input_args_GL.txt
setup_slurm_run.py --run_name GL_prelim -q 1415_queue_GL_prelim.txt --time 04:00:00 -j 7 -e ATL14
# sbatch GL_prelim/slurm_run.sh, then inspect with check_tiles.ipynb

make_ATL1415_queue.py matched <ATL14_root>/rel005/north/GL/input_args_GL.txt
setup_slurm_run.py --run_name GL_matched -q 1415_queue_GL_matched.txt --time 04:00:00 -j 8 -e ATL14

scripts/make_mosaic_jobs <ATL14_root>/rel005/north/GL

ATL14_write2nc.py @<ATL14_root>/rel005/north/GL/input_args_GL.txt
ATL15_write2nc.py @<ATL14_root>/rel005/north/GL/input_args_GL.txt
```
