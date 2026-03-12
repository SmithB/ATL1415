#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:45:47 2020

@author: ben05
"""

import numpy as np
from scipy import stats
import sys, os, h5py, glob, csv
import io, re
import ast
#import pointCollection as pc
import importlib
from netCDF4 import Dataset
#import matplotlib.pyplot as plt
#from ATL1415 import ATL14_attrs_meta
from ATL1415 import ATL14_attrs_meta, make_nc_projection_variable, make_tile_stats_group


#from ATL11.h5util import create_attribute

def ATL14_write2nc(args):
    dz_dict ={'x':'x',   # ATL14 .nc varname : z0.h5 varname
              'y':'y',
              'h':'z0',
              'h_sigma':'sigma_z0',
              'ice_area':'cell_area',
              'data_count':'count',
              'misfit_rms':'misfit_rms',
              'misfit_scaled_rms':'misfit_scaled_rms',
              }
    nctype = {'float64':'f8',
              'float32':'f4',
              'int8':'i1'}

    # establish output file
    fileout = args.base_dir.rstrip('/') + '/ATL14_' + args.region + '_' + args.cycles + '_100m_' + args.Release + '_' + args.version +'.nc'
    print('output file:',fileout)

    with Dataset(fileout,'w',clobber=True) as nc:
        nc.setncattr('GDAL_AREA_OR_POINT','Area')
        nc.setncattr('Conventions','CF-1.6')
        crs_var_root = make_nc_projection_variable(args.region, nc)
        tilegrp = make_tile_stats_group(nc, args)

        # get handle for input file with ROOT and height_change variables.
        FH = h5py.File(args.base_dir.rstrip('/')+'/z0.h5','r')
        if 'z0' not in FH:
            print('no z0.h5 file')
            FH.close()
            exit(-1)
        else:
            print('Reading file:',args.base_dir.rstrip('/')+'/z0.h5')

        with importlib.resources.open_text('ATL1415.resources', 'ATL14_output_attrs.csv', encoding='utf-8-sig') as attrfile:
            reader=list(csv.DictReader(attrfile))

        attr_names=[x for x in reader[0].keys() if x != 'field' and x != 'group']

        field_names = [row['field'] for row in reader if 'ROOT' in row['group']]

        # create dimensions
        for field in ['x', 'y', 'ice_area']:
            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
            dimensions = field_attrs[field]['dimensions'].split(',')
            dimensions = tuple(x.strip() for x in dimensions)
            # if field != 'ice_mask':
            fill_value = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
            # else:
                # fill_value = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
            data = np.array(FH['z0'][dz_dict[field]])

            if field == 'ice_area':
                data[data==0.0] = np.nan
                ice_area_mask = data  # where ice_area is invalid, so are h and h_sigma
            if field == 'x':
                nc.createDimension(field_attrs[field]['dimensions'],data.shape[0])
                x = data
                xll = np.min(x)
                dx = x[1]-x[0]
            if field == 'y':
                nc.createDimension(field_attrs[field]['dimensions'],data.shape[0])
                y = data
                yll = np.max(y)
                dy = y[0]-y[1]

            data = np.nan_to_num(data,nan=fill_value)
            dsetvar = nc.createVariable(field,
                                        nctype[field_attrs[field]['datatype']],
                                        dimensions, zlib=True,
                                        least_significant_digit=ast.literal_eval(field_attrs[field]['least_significant_digit']),
                                        fill_value=fill_value)
            dsetvar[:] = data

            for attr in attr_names:
                dsetvar.setncattr(attr,field_attrs[field][attr])
            # add attributes for projection
            if field == 'x':
                dsetvar.standard_name = 'projection_x_coordinate'
            if field == 'y':
                dsetvar.standard_name = 'projection_y_coordinate'
            if field == 'ice_area':  # 'ice_mask'
                dsetvar.setncattr('grid_mapping','Polar_Stereographic')

        crs_var_root.GeoTransform = (xll,dx,0,yll,0,dy)

        for field in [item for item in field_names if item != 'x' and item != 'y' and item != 'ice_area']:
            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
            dimensions = field_attrs[field]['dimensions'].split(',')
            dimensions = tuple(x.strip() for x in dimensions)
            data = np.array(FH['z0'][dz_dict[field]])
            if field.startswith('h'):
                data[np.isnan(ice_area_mask)] = np.nan
            if field_attrs[field]['datatype'].startswith('int'):
                fill_value = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
            elif field_attrs[field]['datatype'].startswith('float'):
                fill_value = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
            data = np.nan_to_num(data,nan=fill_value)
            dsetvar = nc.createVariable(field,
                                        nctype[field_attrs[field]['datatype']],
                                        dimensions, zlib=True, least_significant_digit=None,
                                        fill_value=fill_value)
            dsetvar[:] = data
            for attr in attr_names:
                dsetvar.setncattr(attr,field_attrs[field][attr])
            dsetvar.setncattr('grid_mapping','Polar_Stereographic')
        with importlib.resources.path('ATL1415.resources', 'atl14_metadata_template.nc') as ncTemplate:
            ATL14_attrs_meta.write_atl14meta(nc, fileout, ncTemplate, args)

        FH.close()

    return fileout

def main():

    import argparse
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, fromfile_prefix_chars='@')
    parser.add_argument('-b','--base_dir', type=str, default=os.getcwd(), help='directory in which to look for mosaicked .h5 files')
    parser.add_argument('-rr','--region', type=str, help='2-letter region indicator \n'
                                                         '\t A[1-4]: Antarctica, by quadrant \n'
                                                         '\t AK: Alaska \n'
                                                         '\t CN: Arctic Canada North \n'
                                                         '\t CS: Arctic Canada South \n'
                                                         '\t GL: Greeland and peripheral ice caps \n'
                                                         '\t IS: Iceland \n'
                                                         '\t SV: Svalbard \n'
                                                         '\t RA: Russian Arctic')
    parser.add_argument('-c','--cycles', type=str, help="4-digit number specifying first/last cycles for output filename")
    parser.add_argument('-R','--Release', type=str, help="3-digit release number for output filename")
    parser.add_argument('-v','--version', type=str, help="2-digit version number for output filename")
    parser.add_argument('-list11','--ATL11_lineage_dir', type=str, help='directory in which to look for ATL11 .h5 filenames')
    parser.add_argument('--tiles_dir', type=str, help='directory in which to look for tile .h5 files, defaults to [base_dir]/prelim')
    parser.add_argument('--ATL11_xover_dir', type=str, help="directory in which to look for ATL11 crossover files")
    parser.add_argument('--ATL11_index', type=str, help='GeoIndex file pointing to ATL11 data')
    args, _=parser.parse_known_args()

    if args.ATL11_lineage_dir is None:
        # if ATL11 lineage_dir is not specified, assume that the grandparent of the ATL11_index works
        args.ATL11_lineage_dir=os.path.dirname(os.path.dirname(args.ATL11_index))

    if args.tiles_dir is None:
        args.tiles_dir=os.path.join(args.base_dir, 'prelim')

    print('args:',args)



    fileout = ATL14_write2nc(args)

if __name__=='__main__':
    main()
