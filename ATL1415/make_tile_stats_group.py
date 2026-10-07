#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 09:39:20 2025

@author: ben
"""
import importlib
import numpy as np
import csv
import os
import re
import h5py
# by full path: `from ATL1415 import` returns the MODULE once anything has
# imported the submodule (the package attribute is rebound on import)
from ATL1415.make_nc_projection_variable import make_nc_projection_variable


def make_tile_stats_group(nc, args, tile_spacing = 40):
    """
    make the tile_stats group for an ATL14 or ATL15 file

    Parameters
    ----------
    nc: netCDF4 file handle
        file handle for ouput file
    args: dict
        dictionary giving input arguments
    tile_spacing: float
        tile-center spacing in km

    Returns:
    -------
    tilegrp: netCDF4 group handle
        group handle for tile_stats group

    """
    tilegrp = nc.createGroup('tile_stats')
    with importlib.resources.open_text('ATL1415.resources','tile_stats_output_attrs.csv', encoding='utf-8-sig') as fh:
        tile_reader = list(csv.DictReader(fh))

    tile_attr_names=[x for x in tile_reader[0].keys() if x != 'field' and x != 'group']

    tile_field_attrs_by_name = {row['field']: {tile_attr_names[ii]:row[tile_attr_names[ii]] for ii in range(len(tile_attr_names))} for row in tile_reader}

    tile_field_names = [row['field'] for row in tile_reader]

    tile_stats={}        # dict for appending data from the tile files
    for field in tile_field_names:
        if field not in tile_stats:
            tile_stats[field] = { 'data': [], 'mapped':np.array(())}


    # Each tile's values come from ATL1415.tile_meta: saved by the 200 km jobs
    # (args.tile_meta_dir), or read from the tiles -- which may be an s3://
    # prefix on MAAP -- once per run (read_tile_stats has the field list).
    from ATL1415.tile_meta import tile_records
    for file_path, rec in tile_records(args):
        file = os.path.basename(file_path)
        if rec['stats'] is None:
            print(f"ATL15_write2nc: problem in write_tile_stats with [ {file} ], skipping")
            continue
        tile_stats['x']['data'].append(int(re.match(r'^.*E(.*)\_.*$',file).group(1)))
        tile_stats['y']['data'].append(int(re.match(r'^.*N(.*)\..*$',file).group(1)))
        for field, value in rec['stats'].items():
            tile_stats[field]['data'].append(value)

    # establish output grids from min/max of x and y
    for key in tile_stats.keys():
        if key == 'N_data' or key == 'N_bias':  # key == 'x' or key == 'y' or
            tile_stats[key]['mapped'] = np.zeros( [len(np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+tile_spacing,tile_spacing)),
                                                    len(np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+tile_spacing,tile_spacing))],
                                                    dtype=int)
        else:
            tile_stats[key]['mapped'] = np.zeros( [len(np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+tile_spacing,tile_spacing)),
                                                    len(np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+tile_spacing,tile_spacing))],
                                                    dtype=float)
    # put data into grids
    for key in tile_stats.keys():
        # fact helps convert x,y in km to m
        if key == 'x' or key == 'y':
            continue
        for (yt, xt, dt) in zip(tile_stats['y']['data'], tile_stats['x']['data'], tile_stats[key]['data']):
            if not np.isfinite(dt):
                print(f"ATL14_write2nc: found bad tile_stats value in field {key} : {dt} at x={xt}, y={yt}")
                continue
            row=int((yt-np.min(tile_stats['y']['data']))/tile_spacing)
            col=int((xt-np.min(tile_stats['x']['data']))/tile_spacing)
            tile_stats[key]['mapped'][row,col] = dt
        tile_stats[key]['mapped'] = np.ma.masked_where(tile_stats[key]['mapped'] == 0, tile_stats[key]['mapped'])

    # make dimensions, fill them as variables
    tilegrp.createDimension('y',len(np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+tile_spacing,tile_spacing)))
    tilegrp.createDimension('x',len(np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+tile_spacing,tile_spacing)))

    # create tile_stats/ variables in .nc file
    for field in tile_field_names:
        tile_field_attrs = {field: tile_field_attrs_by_name[field]}
        if field == 'x':
            dsetvar = tilegrp.createVariable('x', tile_field_attrs[field]['datatype'], ('x',), fill_value=np.finfo(tile_field_attrs[field]['datatype']).max, zlib=True)
            dsetvar[:] = np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+tile_spacing,tile_spacing) * 1000 # convert from km to meter
            dsetvar.setncattr('standard_name','projection_x_coordinate')
        elif field == 'y':
            dsetvar = tilegrp.createVariable('y', tile_field_attrs[field]['datatype'], ('y',), fill_value=np.finfo(tile_field_attrs[field]['datatype']).max, zlib=True)
            dsetvar[:] = np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+tile_spacing,tile_spacing) * 1000 # convert from km to meter
            dsetvar.setncattr('standard_name','projection_y_coordinate')
        elif field == 'N_data' or field == 'N_bias':
            dsetvar = tilegrp.createVariable(field, tile_field_attrs[field]['datatype'],('y','x'),fill_value=np.iinfo(tile_field_attrs[field]['datatype']).max, zlib=True)
        else:
            dsetvar = tilegrp.createVariable(field, tile_field_attrs[field]['datatype'],('y','x'),fill_value=np.finfo(tile_field_attrs[field]['datatype']).max, zlib=True)

        if field != 'x' and field != 'y':
            dsetvar[:] = tile_stats[field]['mapped'][:]
            dsetvar.setncattr('coordinates', "tile_stats/y tile_stats/x")
        else:
            dsetvar.setncattr('coordinates', f"tile_stats/{field}")

        for attr in ['units','dimensions','datatype','description','long_name','source']:
            dsetvar.setncattr(attr,tile_field_attrs[field][attr])
        dsetvar.setncattr('grid_mapping','Polar_Stereographic')


    crs_var = make_nc_projection_variable(args.region, tilegrp)
    crs_var.GeoTransform = str(tilegrp['x'][0])+" "+str(tilegrp['x'][1]-tilegrp['x'][0])+" 0.0 "+str(tilegrp['y'][0])+" 0.0 "+str(tilegrp['y'][1]-tilegrp['y'][0])

    return tilegrp
