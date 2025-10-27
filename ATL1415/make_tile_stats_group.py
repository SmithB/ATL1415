#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 09:39:20 2025

@author: ben
"""
import pkg_resources
import numpy as np
import csv
import os
import re
import h5py
from ATL1415 import  make_nc_projection_variable


def make_tile_stats_group(nc, args):
    """
    make the tile_stats group for an ATL14 or ATL15 file

    Parameters
    ----------
    nc: netCDF4 file handle
        file handle for ouput file
    args: dict
        dictionary giving input arguments

    Returns:
    -------
    tilegrp: netCDF4 group handle
        group handle for tile_stats group

    """
    tilegrp = nc.createGroup('tile_stats')
    tileFile = pkg_resources.resource_filename('ATL1415','resources/tile_stats_output_attrs.csv')
    with open(tileFile,'r', encoding='utf-8-sig') as tilefile:
        tile_reader=list(csv.DictReader(tilefile))

    tile_attr_names=[x for x in tile_reader[0].keys() if x != 'field' and x != 'group']

    tile_field_names = [row['field'] for row in tile_reader]

    tile_stats={}        # dict for appending data from the tile files
    for field in tile_field_names:
        if field not in tile_stats:
            tile_stats[field] = { 'data': [], 'mapped':np.array(())}


    files = os.listdir(args.tiles_dir)
    files = [f for f in files if f.endswith('.h5')]
    for file in files:
        try:
            tile_stats['x']['data'].append(int(re.match(r'^.*E(.*)\_.*$',file).group(1)))
        except Exception:
            print(f"ATL15_write2nc problem in write_tile_stats with [ {file} ], skipping")
            continue
        tile_stats['y']['data'].append(int(re.match(r'^.*N(.*)\..*$',file).group(1)))

        # with h5py.File(os.path.join(args.base_dir,sub,file),'r') as h5: #used for sub in ['centers','edges','corners']
        with h5py.File(os.path.join(args.tiles_dir,file),'r') as h5:
            tile_stats['N_data']['data'].append( np.sum(h5['data']['three_sigma_edit'][:]) )
            tile_stats['RMS_data']['data'].append( h5['RMS']['data'][()] )  # use () for getting a scalar.
            tile_stats['RMS_bias']['data'].append( np.sqrt(np.mean((h5['bias']['val'][:]/h5['bias']['expected'][:])**2)) )
            tile_stats['N_bias']['data'].append( len(h5['bias']['val'][:]) )  #### or all BUT the zeros.
            tile_stats['RMS_d2z0dx2']['data'].append( h5['RMS']['grad2_z0'][()] )
            tile_stats['RMS_d2zdt2']['data'].append( h5['RMS']['d2z_dt2'][()] )
            tile_stats['RMS_d2zdx2dt']['data'].append( h5['RMS']['grad2_dzdt'][()] )
            tile_stats['sigma_xx0']['data'].append( h5['E_RMS']['d2z0_dx2'][()] )
            tile_stats['sigma_tt']['data'].append( h5['E_RMS']['d2z_dt2'][()] )
            tile_stats['sigma_xxt']['data'].append( h5['E_RMS']['d3z_dx2dt'][()] )

    # establish output grids from min/max of x and y
    for key in tile_stats.keys():
        if key == 'N_data' or key == 'N_bias':  # key == 'x' or key == 'y' or
            tile_stats[key]['mapped'] = np.zeros( [len(np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+40,40)),
                                                    len(np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+40,40))],
                                                    dtype=int)
        else:
            tile_stats[key]['mapped'] = np.zeros( [len(np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+40,40)),
                                                    len(np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+40,40))],
                                                    dtype=float)
    # put data into grids
    for key in tile_stats.keys():
        # fact helps convert x,y in km to m
        if key == 'x' or key == 'y':
            continue
        for (yt, xt, dt) in zip(tile_stats['y']['data'], tile_stats['x']['data'], tile_stats[key]['data']):
            if not np.isfinite(dt):
                print(f"ATL14_write2nc: found bad tile_stats value in field {key} at x={xt}, y={yt}")
                continue
            row=int((yt-np.min(tile_stats['y']['data']))/40)
            col=int((xt-np.min(tile_stats['x']['data']))/40)
            tile_stats[key]['mapped'][row,col] = dt
        tile_stats[key]['mapped'] = np.ma.masked_where(tile_stats[key]['mapped'] == 0, tile_stats[key]['mapped'])

    # make dimensions, fill them as variables
    tilegrp.createDimension('y',len(np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+40,40)))
    tilegrp.createDimension('x',len(np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+40,40)))

    # create tile_stats/ variables in .nc file
    for field in tile_field_names:
        tile_field_attrs = {row['field']: {tile_attr_names[ii]:row[tile_attr_names[ii]] for ii in range(len(tile_attr_names))} for row in tile_reader if field in row['field']}
        if field == 'x':
            dsetvar = tilegrp.createVariable('x', tile_field_attrs[field]['datatype'], ('x',), fill_value=np.finfo(tile_field_attrs[field]['datatype']).max, zlib=True)
            dsetvar[:] = np.arange(np.min(tile_stats['x']['data']),np.max(tile_stats['x']['data'])+40,40.) * 1000 # convert from km to meter
            dsetvar.setncattr('standard_name','projection_x_coordinate')
        elif field == 'y':
            dsetvar = tilegrp.createVariable('y', tile_field_attrs[field]['datatype'], ('y',), fill_value=np.finfo(tile_field_attrs[field]['datatype']).max, zlib=True)
            dsetvar[:] = np.arange(np.min(tile_stats['y']['data']),np.max(tile_stats['y']['data'])+40,40.) * 1000 # convert from km to meter
            dsetvar.setncattr('standard_name','projection_y_coordinate')
        elif field == 'N_data' or field == 'N_bias':
            dsetvar = tilegrp.createVariable(field, tile_field_attrs[field]['datatype'],('y','x'),fill_value=np.iinfo(tile_field_attrs[field]['datatype']).max, zlib=True)
        else:
            dsetvar = tilegrp.createVariable(field, tile_field_attrs[field]['datatype'],('y','x'),fill_value=np.finfo(tile_field_attrs[field]['datatype']).max, zlib=True)

        if field != 'x' and field != 'y':
            dsetvar[:] = tile_stats[field]['mapped'][:]

        for attr in ['units','dimensions','datatype','coordinates','description','coordinates','long_name','source']:
            dsetvar.setncattr(attr,tile_field_attrs[field][attr])
        dsetvar.setncattr('grid_mapping','Polar_Stereographic')


    crs_var = make_nc_projection_variable(args.region, tilegrp)
    crs_var.GeoTransform = str(tilegrp['x'][0])+" "+str(tilegrp['x'][1]-tilegrp['x'][0])+" 0.0 "+str(tilegrp['y'][0])+" 0.0 "+str(tilegrp['y'][1]-tilegrp['y'][0])

    return tilegrp
