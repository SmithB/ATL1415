#!/usr/bin/env python3
"""
What a region's mosaic is made of, without pointCollection: the mosaic groups
(make_fields), the 200 km tiles (read_200km_centers: the region's canonical
list; footprint_centers_200km: the cells its tiles reach; centers_200km: the
cell each tile's center is in), and which regions take the 200 km step at all
(uses_200km_tiles).

Light on purpose -- numpy only -- so the ADE-side DPS submitter, which runs in
an interpreter without pointCollection, lists a mosaic run's jobs from the
same definitions the workers use (docs/plan_dps_mosaic.sh D3-2).
make_200km_tiles.py imports make_fields from here.
"""
import os
import re

import numpy as np


# Ben 2026-10-01 (docs/plan_dps_mosaic.sh AD8): "Just use the 200-km step for
# Antarctica and Greenland."  Every other region mosaics directly from its
# solve tiles (make_mosaic_jobs.py).
REGIONS_200KM = ('GL', 'AA', 'A1', 'A2', 'A3', 'A4')


def uses_200km_tiles(region):
    """
    True if the region's mosaics are made by way of 200 km tiles
    (make_200km_tiles.py, then make_200km_to_mosaic_jobs.py); False if they
    are made directly from the solve tiles (make_mosaic_jobs.py).
    """
    return region in REGIONS_200KM


def make_fields(dzdt_lags, t_res=0.25, skip_z0=False ):
    """
    Build the field lists and time ranges to mosaic for each output group.

    Parameters
    ----------
    dzdt_lags : iterable
        dzdt lags (in grid-spacing units) to generate dzdt/avg_dzdt groups for.
    t_res : float, optional
        dz/dt grid time resolution, used to convert lags to years. The default is 0.25.
    skip_z0 : bool, optional
        if true, omit the z0 group. The default is False.

    Returns
    -------
    fields : dict
        mapping of group name to list of fields to mosaic.
    time_ranges : dict
        mapping of group name to [start, end] year range for the mosaic.

    """

    fields={}
    if not skip_z0:
        fields['z0']="z0 sigma_z0 misfit_rms misfit_scaled_rms mask cell_area count".split(' ')

    fields['dz']="dz sigma_dz count misfit_rms misfit_scaled_rms mask cell_area".split(' ')

    time_ranges={}
    time_ranges['dz']=[2019, 2050]

    lags = [ f'_lag{lag}' for lag in dzdt_lags ]

    for lag in lags:
        field_str='dzdt'+lag
        fields[field_str] = ["dzdt"+lag, "sigma_dzdt"+lag, "cell_area"]
        time_ranges[field_str] = [2019 + t_res * int(lag.replace('_lag',''))/2, 2050]
    for res in ["_40000m", "_20000m", "_10000m"]:
        fields['avg_dz'+res] = ["avg_dz"+res, "sigma_avg_dz"+res,'cell_area']
        time_ranges['avg_dz'+res] = [2019, 2050]
        for lag in lags:
            field_str='avg_dzdt'+res+lag
            fields[field_str]=[field_str, 'sigma_'+field_str, 'cell_area']
            time_ranges[field_str] = [2019 + t_res * int(lag.replace('_lag',''))/2, 2050]
    #for key, item in fields.items():
    #print(key+" : "+str(item))
    #print(fields)
    return fields, time_ranges


def centers_200km(tile_names, tile_W=200e3):
    """
    The 200 km cells holding the tiles' CENTERS, SORTED (by x, then y).

    Each tile is in exactly one of these, so this is the rule that gives a
    tile to one 200 km tile (ATL1415.tile_meta).  It is NOT the list of 200 km
    tiles to make: a tile reaches W/2 past its center, into cells that may
    hold no tile center -- footprint_centers_200km (plan_200km_footprint.sh).

    Parameters
    ----------
    tile_names : iterable of str
        tile file names or paths, E<x>_N<y>.h5 (km).
    tile_W : float
        200 km tile width, m.

    Returns
    -------
    list of [x, y] floats, m
    """
    tile_re = re.compile(r'E(-?[0-9]+)_N(-?[0-9]+)\.h5$')
    xy0 = [[1000.*int(v) for v in tile_re.search(os.path.basename(name)).groups()]
           for name in tile_names if tile_re.search(os.path.basename(name))]
    if not xy0:
        return []
    xyc = np.floor(np.array(xy0)/tile_W)*tile_W + tile_W/2
    return [list(map(float, row)) for row in np.unique(xyc, axis=0)]


def footprint_centers_200km(tile_names, half_width=30e3, tile_W=200e3):
    """
    Every 200 km cell that some tile's square (center +- half_width, m) reaches
    with positive area, as [x, y] centers SORTED (by x, then y).

    The 200 km tiles a region needs.  centers_200km gives only the cells holding
    a tile center, and the parts of tiles reaching past a cell edge into a cell
    with no tile center were never mosaicked: southern Greenland lost the strip
    y -3200..-3230 km (Ben, 2026-10-08; plan_200km_footprint.sh).  A square
    whose edge only touches a cell edge does not count: the node on the edge is
    in the neighbouring 200 km tile, which keeps its bounds.
    """
    tile_re = re.compile(r'E(-?[0-9]+)_N(-?[0-9]+)\.h5$')
    cells = set()
    for name in tile_names:
        m = tile_re.search(os.path.basename(name))
        if not m:
            continue
        x, y = (1000. * int(v) for v in m.groups())
        ix = range(int(np.floor((x - half_width + 1) / tile_W)), int(np.floor((x + half_width - 1) / tile_W)) + 1)
        iy = range(int(np.floor((y - half_width + 1) / tile_W)), int(np.floor((y + half_width - 1) / tile_W)) + 1)
        cells.update((i * tile_W + tile_W / 2, j * tile_W + tile_W / 2) for i in ix for j in iy)
    return [list(map(float, c)) for c in sorted(cells)]


def region_200km_list(region):
    """The path of the region's canonical 200 km list (a package resource), or None."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources', str(region),
                        '200km_tile_list.txt')
    return path if os.path.isfile(path) else None


def read_200km_centers(region):
    """The region's canonical 200 km tile centers ([x, y], m, file order), or None."""
    path = region_200km_list(region)
    if path is None:
        return None
    centers = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.strip():
                parts = line.split()
                if len(parts) != 2:
                    raise ValueError(f'{path}:{n}: not an "<x> <y>" pair: {line.strip()!r}')
                centers.append([float(parts[0]), float(parts[1])])
    return centers

