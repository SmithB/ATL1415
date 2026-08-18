#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Split Antarctic AA / AA_44km region outputs into four polar-stereographic
sector quadrants (A1..A4): symlinks 200km-tile and per-tile outputs into
each sector, and generates each sector's bounds.txt / input_args_A{n}.txt.

Converted from notebooks/make_and_check_links.ipynb (cells before the
"check the output" markdown cell); the notebook's post-"check the output"
diagnostic/plotting cells are not part of this script.

@author: ben
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

N_SECTORS = 4

xy_re = re.compile('E(.*)_N(.*).h5')


def get_200km_loc(filename, field):
    """
    Parse a 200km-tile filename into its bounds and center.

    Parameters
    ----------
    filename : str
        200km-tile file path, of the form '.../<field><xmin>_<xmax>_<ymin>_<ymax>.h5'.
    field : str
        field name to strip from the filename before parsing the bounds.

    Returns
    -------
    bounds : list
        [[xmin, xmax], [ymin, ymax]], in km.
    ctr : numpy.ndarray
        [x_center, y_center], in km.

    """
    file = os.path.basename(filename)
    bounds = np.array([*map(int, file.replace(field, '').replace('.h5', '').split('_'))])
    return [bounds[0:2], bounds[2:]], np.array([np.mean(bounds[0:2]), np.mean(bounds[2:])])


def get_quadrant(xy):
    """
    Classify point(s) into a polar-stereographic sector quadrant.

    Parameters
    ----------
    xy : numpy.ndarray
        [x, y] pair, or Nx2 array of [x, y] pairs.

    Returns
    -------
    quadrant : numpy.ndarray
        integer sector number (1-4): 1=(x>=0,y>=0), 2=(x<0,y>=0),
        3=(x<0,y<0), 4=(x>=0,y<0).

    """
    if xy.ndim == 1:
        quadrant = np.zeros([1], dtype=int)
    else:
        quadrant = np.zeros(xy.shape[0], dtype=int)
    ii = (xy[0] >= 0) & (xy[1] >= 0)
    quadrant[ii] = 1
    ii = (xy[0] < 0) & (xy[1] >= 0)
    quadrant[ii] = 2
    ii = (xy[0] < 0) & (xy[1] < 0)
    quadrant[ii] = 3
    ii = (xy[0] >= 0) & (xy[1] < 0)
    quadrant[ii] = 4
    return quadrant


def check_source_dirs(top_N, top_S):
    """
    Verify that the source region directories exist, exiting with an error if not.

    Parameters
    ----------
    top_N : str
        outer-domain source region directory.
    top_S : str
        near-pole/fine-grid source region directory.

    Returns
    -------
    None.

    """
    for thedir in [top_N, top_S]:
        if thedir is None:
            continue
        if not os.path.isdir(thedir):
            print(f"setup_AA_sectors.py:\n\tError: source directory not found: {thedir}")
            sys.exit(1)


def make_sector_dirs(top_dir):
    """
    Create the A1..A{N_SECTORS} sector directory trees under top_dir.

    Parameters
    ----------
    top_dir : str
        directory in which to create the sector directories.

    Returns
    -------
    None.

    """
    def mkdir(thedir):
        if not os.path.isdir(thedir):
            os.mkdir(thedir)

    for quadrant in range(1, N_SECTORS+1):
        sector_dir = os.path.join(top_dir, f'A{quadrant}')
        mkdir(sector_dir)
        mkdir(os.path.join(sector_dir, '200km_tiles'))
        mkdir(os.path.join(sector_dir, 'prelim'))
        mkdir(os.path.join(sector_dir, 'matched'))


def mirror_200km_field_dirs(top_dir, top_N):
    """
    Create one 200km_tiles/<field> subdirectory per sector, matching top_N's field list.

    Parameters
    ----------
    top_dir : str
        directory containing the sector directories.
    top_N : str
        outer-domain source region directory, whose 200km_tiles field-name
        subdirectories are mirrored into each sector.

    Returns
    -------
    None.

    """
    subs_200 = glob.glob(os.path.join(top_N, '200km_tiles', '*'))
    for quadrant in range(1, N_SECTORS+1):
        quad_200km = os.path.join(top_dir, f'A{quadrant}', '200km_tiles')
        for sub in subs_200:
            new_dir = os.path.join(quad_200km, os.path.basename(sub))
            if not os.path.isdir(new_dir):
                os.mkdir(new_dir)


def link_200km_tiles(top_dir, top_S, top_N, near_pole_radius):
    """
    Symlink 200km-tile outputs from top_S/top_N into each sector.

    Tiles whose center lies within near_pole_radius meters of the origin in
    both x and y are drawn from top_S; all other tiles are drawn from top_N.

    Parameters
    ----------
    top_dir : str
        directory containing the sector directories.
    top_S : str
        near-pole/fine-grid source region directory.
    top_N : str
        outer-domain source region directory.
    near_pole_radius : float
        |x|,|y| threshold, in meters, below which tiles come from top_S.

    Returns
    -------
    None.

    """
    def is_near_pole(ctr_m):
        return np.all(np.abs(ctr_m) < near_pole_radius)

    for src_top, use_south in zip([top_S, top_N], [True, False]):
        if src_top is None:
            continue
        src_subs = sorted(glob.glob(os.path.join(src_top, '200km_tiles', '*')))
        print(f"linking 200km tiles from {src_top}: {len(src_subs)} field(s)")
        for src_sub in src_subs:
            field = os.path.basename(src_sub)
            tile_files = glob.glob(os.path.join(src_sub, '*'))
            for tile_file in tile_files:
                bds, ctr = get_200km_loc(tile_file, field)
                if is_near_pole(ctr*1000) != use_south:
                    continue
                quad = get_quadrant(ctr)
                dst = os.path.join(top_dir, f'A{int(quad[0])}', '200km_tiles',
                                    field, os.path.basename(tile_file))
                if not os.path.islink(dst):
                    os.symlink(tile_file, dst)


def link_tiles(top_dir, top_S, top_N, near_pole_radius, steps):
    """
    Symlink per-tile outputs from top_S/top_N into each sector.

    Tiles with |x|,|y| both below near_pole_radius (in km) are drawn from
    top_S; all other tiles are drawn from top_N.

    Parameters
    ----------
    top_dir : str
        directory containing the sector directories.
    top_S : str
        near-pole/fine-grid source region directory.
    top_N : str
        outer-domain source region directory.
    near_pole_radius : float
        |x|,|y| threshold, in meters, below which tiles come from top_S.
    steps : iterable of str
        per-tile subdirectories to link (e.g. ['prelim', 'matched']).

    Returns
    -------
    None.

    """
    near_pole_radius_km = near_pole_radius / 1000

    def south_test(xy):
        return np.all(np.abs(xy) < near_pole_radius_km)

    def north_test(xy):
        return np.any(np.abs(xy) >= near_pole_radius_km)

    for src, test in zip([top_S, top_N], [south_test, north_test]):
        for step in steps:
            if src is None:
                continue
            tiles = glob.glob(os.path.join(src, step, 'E*.h5'))
            print(f"linking {len(tiles)} {step} tiles from {src}")
            for tile in tiles:
                xy = np.array([*map(int, xy_re.search(tile).groups())])
                if not test(xy):
                    continue
                quad = get_quadrant(xy)
                dst = os.path.join(top_dir, f'A{int(quad[0])}', step, os.path.basename(tile))
                if not os.path.islink(dst):
                    os.symlink(tile, dst)


def compute_full_bounds(top_N, pad):
    """
    Compute the padded [[xmin,xmax],[ymin,ymax]] extent of top_N's prelim tiles.

    Parameters
    ----------
    top_N : str
        outer-domain source region directory.
    pad : float
        padding, in meters, added to the tile-center extent on each side.

    Returns
    -------
    full_bounds : list
        [[xmin, xmax], [ymin, ymax]], in meters.

    """
    files = glob.glob(os.path.join(top_N, 'prelim', 'E*'))
    xy = np.c_[[[*map(int, xy_re.search(file).groups())] for file in files]]
    pad_km = pad/1000
    full_bounds = [[np.min(xy[:, 0])-pad_km, np.max(xy[:, 0])+pad_km],
                   [np.min(xy[:, 1])-pad_km, np.max(xy[:, 1])+pad_km]]
    return [1000*np.array(jj) for jj in full_bounds]


def write_sector_bounds(top_dir, full_bounds):
    """
    Write each sector's bounds.txt from the padded full-domain extent.

    Parameters
    ----------
    top_dir : str
        directory containing the sector directories.
    full_bounds : list
        [[xmin, xmax], [ymin, ymax]] full-domain extent, in meters, as
        returned by compute_full_bounds.

    Returns
    -------
    None.

    """
    quadrant_bounds = {
        1: (0, full_bounds[0][1], 0, full_bounds[1][1]),
        2: (full_bounds[0][0], 0, 0, full_bounds[1][1]),
        3: (full_bounds[0][0], 0, full_bounds[1][0], 0),
        4: (0, full_bounds[0][1], full_bounds[1][0], 0),
    }
    for quadrant in range(1, N_SECTORS+1):
        with open(os.path.join(top_dir, f'A{quadrant}', 'bounds.txt'), 'w') as fh:
            fh.write('%d %d %d %d\n' % quadrant_bounds[quadrant])


def write_sector_input_args(top_dir, north_name):
    """
    Generate each sector's input_args_A{n}.txt from north_name's input_args file.

    Every line is copied verbatim except lines starting with '-b=' or
    '--region=', where north_name is replaced with the sector name.

    Parameters
    ----------
    top_dir : str
        directory containing the sector directories.
    north_name : str
        name of the outer-domain source region, whose input_args_<north_name>.txt
        is used as the template.

    Returns
    -------
    None.

    """
    src = os.path.join(top_dir, north_name, f'input_args_{north_name}.txt')
    for quadrant in range(1, N_SECTORS+1):
        sector = f'A{quadrant}'
        dst = os.path.join(top_dir, sector, f'input_args_{sector}.txt')
        with open(src, 'r') as fh_in, open(dst, 'w') as fh_out:
            for line in fh_in:
                if line.startswith('-b=') or line.startswith('--region='):
                    fh_out.write(line.replace(north_name, sector))
                else:
                    fh_out.write(line)


def main():
    parser = argparse.ArgumentParser(fromfile_prefix_chars="@")
    parser.add_argument('top_dir', type=str,
                         help='directory containing north_name/south_name region dirs; '
                              'A1..A4 sector dirs are created here')
    parser.add_argument('--north_name', type=str, default='AA',
                         help='name of the outer-domain source region (subdir of top_dir)')
    parser.add_argument('--south_name', type=str, default='AA_44km',
                         help='name of the near-pole/fine-grid source region (subdir of top_dir)')
    parser.add_argument('--near_pole_radius', type=float, default=4.e5,
                         help='|x|,|y| threshold (m) below which tiles are drawn from '
                              '--south_name instead of --north_name')
    parser.add_argument('--bounds_pad', type=float, default=4.e4,
                         help='padding (m) added to the full prelim-tile extent when '
                              'computing each sector bounds.txt')
    parser.add_argument('--steps', type=str, nargs='+', default=['prelim', 'matched'],
                         help='per-tile subdirectories to symlink (in addition to 200km_tiles)')
    args, _ = parser.parse_known_args()

    top_N = os.path.join(args.top_dir, args.north_name)
    
    if args.near_pole_radius > 0:
        top_S = os.path.join(args.top_dir, args.south_name)
    else:
        top_S = None

    check_source_dirs(top_N, top_S)

    print(f"making sector directories under {args.top_dir}")
    make_sector_dirs(args.top_dir)

    print("mirroring 200km_tiles field directories into sectors")
    mirror_200km_field_dirs(args.top_dir, top_N)

    print("linking 200km tiles into sectors")
    link_200km_tiles(args.top_dir, top_S, top_N, args.near_pole_radius)

    print(f"linking {args.steps} tiles into sectors")
    link_tiles(args.top_dir, top_S, top_N, args.near_pole_radius, args.steps)

    print("computing sector bounds")
    full_bounds = compute_full_bounds(top_N, args.bounds_pad)
    write_sector_bounds(args.top_dir, full_bounds)

    print("writing sector input_args files")
    write_sector_input_args(args.top_dir, args.north_name)

    print("done.")


if __name__ == '__main__':
    main()
