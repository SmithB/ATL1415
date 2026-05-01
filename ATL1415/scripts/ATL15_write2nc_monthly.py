#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:45:47 2020

@author: ben05
"""
import numpy as np
import  os, h5py, csv, re, glob
import ast
import importlib
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#import matplotlib as mpl
#import cartopy.crs as ccrs
#import cartopy.feature
from scipy import stats
from ATL1415 import ATL14_attrs_meta, make_nc_projection_variable, make_tile_stats_group

def update_attr_dict(attr_template, args, res_m, res_km):
    '''
    update the field attributes for the current resolution and lags

    This function updates a template attribute dictionary by selecting
        appropriate to downsampled or non-downsampled output files, and
        replacing occurences of {res} and {res_m} with the resolution in km
        or meters, respectively. It also duplicates lines containing 'lag'
        to match the required number of lagged measurements specified in
        input arguments, and translates the duration of each lag into an
        adjective

    Parameters
    ----------
    attr_template : iterable
        iterable of dicts containing field attribute templates
    args : namespace
        input arguments.
    res_m : str
        resolution in meters.
    res_km : str
        resolution in km.

    Returns
    -------
    new_attrs : iterable
        iterable of dicts containig updated field attribtues.

    '''
    # define the words used to describe the lagged dh/dt resolutions
    lag_adjectives = {1/12: 'monthly',
                1/4: 'quarterly',
                1/2: 'semiannual',
                1: 'annual',
                2: 'biennial',
                3: 'triennial',
                4: 'quadrennial',
                5: 'pentennial',
                6: 'hexennial',
                7: 'heptennial',
                8:  'octennial'}
    lag_names = {}
    lag_months = {}
    lag_monthly_names = {}
    lag_adjective_key_list = list(lag_adjectives.keys())
    for lag in args.dzdt_lags:
        months = int(round(lag * args.delta_t*12))
        this = np.argmin(np.abs(args.delta_t*lag - np.array(lag_adjective_key_list)))
        lag_monthly_names[lag] = f'{months:03d}mo'
        lag_names[lag] = lag_adjectives[lag_adjective_key_list[this]]
        lag_months[lag] = months
    if args.verbose:
        print(lag_names)
    lags = args.dzdt_lags
    new_attrs=[]
    average_key='None' if args.avg_scale is None else '{res}'
    for row in attr_template:
        if row['average'] != average_key:
            continue
        if row['lag'] is None:
            new_attrs.append({key:val.replace('{lag_name}', lag_names[1]) for key, val in row})
        else:
            for lag in lags:
                new_row = {attr_name: attr.replace('{lag}', str(lag))\
                                           .replace('{lag_name}', lag_names[lag])\
                                           .replace('{lag_mo}', lag_monthly_names[lag])\
                                           .replace('{lag_mo_num}', str(lag_months[lag]))
                           for attr_name, attr in row.items()}
                new_attrs.append(new_row)
    for rcount, row in enumerate(new_attrs):
        for attr_name, attr in row.items():
            if '{res}' in attr:
                attr=attr.replace('{res}', res_km)
            if '{res_m}' in attr:
                attr=attr.replace('{res_m}', res_m)
            row[attr_name] = attr
    return new_attrs

def make_dataset(field, data, attrs, file_obj, group_obj, nctype, dimScale=False):
    '''
    make a dataset in an netcdf output file

    Parameters
    ----------
    field : str
        field to be created.
    data : numpy.array
        data.
    attrs : dict
        attributes for the field.
    file_obj : netcdf4 file object
        file object to be written to .
    group_obj : netcdf4 group object
        group object to be written to.
    nctype : dict
        definitions for netcdf data types.
    dimScale : book, optional
        If true, define the field as a dimension. The default is False.

    Returns
    -------
    file_obj : netcdf4 file object
        file object that was written.

    '''
    dimensions = attrs['dimensions'].split(',')
    dimensions = tuple(x.strip() for x in dimensions)
    if attrs['datatype'].startswith('int'):
        fill_value = np.iinfo(np.dtype(attrs['datatype'])).max
    elif attrs['datatype'].startswith('float'):
        fill_value = np.finfo(np.dtype(attrs['datatype'])).max
    data = np.nan_to_num(data,nan=fill_value)
    if dimScale:
        group_obj.createDimension(attrs['dimensions'],data.shape[0])

    dsetvar = group_obj.createVariable(field,
                                       nctype[attrs['datatype']],
                                       dimensions,
                                       fill_value=fill_value,
                                       zlib=True,
                                       least_significant_digit=ast.literal_eval(attrs['least_significant_digit']))

    dsetvar[:] = data
    for attr in ['units','dimensions','datatype','coordinates','description','long_name','source']:
        dsetvar.setncattr(attr,attrs[attr])
    # add attributes for projection
    if not field.startswith('time'):
        dsetvar.setncattr('grid_mapping','Polar_Stereographic')
    if field == 'x':
        dsetvar.standard_name = 'projection_x_coordinate'
    if field == 'y':
        dsetvar.standard_name = 'projection_y_coordinate'

    return file_obj

def get_group_attrs(group,  all_attrs):
    group_attrs={}
    for row in all_attrs:
        if row['out_group']==group:
            if row['field'] is None or row['field'] == '':
                group_description = row['description']
                continue
            group_attrs[row['field']] = {attr_name: attr for attr_name, attr in row.items()}
    return group_attrs, group_description

def ATL15_write2nc_monthly(args):

    with importlib.resources.open_text('ATL1415.resources', 'ATL15_monthly_output_attrs.csv', encoding='utf-8-sig') as attrfile:
        all_attrs = list(csv.DictReader(attrfile))

    res = args.grid_spacing[1]
    if args.avg_scale is not None:
        res = args.avg_scale
    res_m = f'{int(res)}'
    if np.mod(res, 1000) < 1:
        res_km = f'{int(round(res/1000))}'
    else:
        res_km = f'{res/1000:0.1f}'
    if args.verbose:
        print(f"ATL15_write2nc_monthly: making nc for resolution {res_km} km, for lags {args.dzdt_lags}")

    # establish output file
    avg_name = f'{res_km}km'
    if np.abs(args.delta_t-1/12)<0.001:
        delta_t_str = '1mo'
    elif np.abs(args.delta_t-1/4) < 0.001:
        delta_t_str = '3mo'
    else:
        raise ValueError(f'time resolution of {args.delta_t} not recognized')

    fileout = os.path.join(args.base_dir , '_'.join(['ATL15', args.region , args.cycles, delta_t_str,  avg_name, args.Release, args.version]) + '.nc')

    if args.verbose:
        print('output file:',fileout)

    all_attrs = update_attr_dict(all_attrs, args, res_m, res_km)

    nctype = {'float64':'f8',
              'float32':'f4',
              'int8':'i1'}

    fh_in={}

    with Dataset(fileout,'w',clobber=True) as nc:
        nc.setncattr('GDAL_AREA_OR_POINT','Area')
        nc.setncattr('Conventions','CF-1.6')
        _ = make_tile_stats_group(nc, args)

        for lag in [None] + args.dzdt_lags:
            if lag is None:
                group = 'delta_h'
            else:
                months = int(round(lag * args.delta_t*12))
                group = f'dhdt_{months:03d}mo'
            if args.verbose:
                print('-'*20 + '\n'+ group+'\n'+'-'*20)
            try:
                group_attrs, group_description = get_group_attrs(group,  all_attrs)
            except Exception as e:
                print("group:" + group)
                print(all_attrs)
                raise(e)
            nc_group = nc.createGroup(group)
            nc_group.setncattr('description', group_description)
            crs_var = make_nc_projection_variable(args.region, nc.groups[group])

            for field, attrs in group_attrs.items():
                dimscale = False
                # out field
                in_file = attrs['source_file']
                # cache the file handle
                if in_file not in fh_in:
                    fh_in[in_file] = h5py.File(os.path.join(args.base_dir, in_file),'r')
                data = np.array(fh_in[in_file][attrs['source_group']][attrs['source_var']])

                if field == 'x':
                    xll = np.min(data)
                    dx = data[1]-data[0]
                    dimscale=True
                elif field == 'y':
                    yll = np.min(data)
                    dy = data[1]-data[0]
                    dimscale = True
                elif field == 'time':
                    data = (data-2018.)*365.25
                    dimscale = True
                if data.ndim==3:
                    # hdf5 files store data as y, x, t
                    data = np.moveaxis(data,2,0)  # t, y, x
                if field == 'ice_area':
                    data[data==0.0] = np.nan
                    mask_lag = np.isfinite(data)
                if 'delta_h' in field or 'dhdt' in field:
                    data[~mask_lag]=np.nan
                # make the dataset
                make_dataset(field, data, attrs, fh_in[in_file], nc_group, nctype, dimScale=dimscale)

            # we should have an x and a y defined by now
            crs_var.GeoTransform = (xll,dx,0,yll,0,dy)

        #ncTemplate=pkg_resources.resource_filename('ATL1415','resources/atl15_metadata_template.nc')
        with importlib.resources.path('ATL1415.resources', 'atl15_metadata_template.nc') as ncTemplate: 
            ATL14_attrs_meta.write_atl14meta(nc, fileout, ncTemplate, args)
    for fh in fh_in.values():
        try:
            fh.close()
        except:
            pass
    return fileout


def main():
    import argparse
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,  fromfile_prefix_chars='@')
    parser.add_argument('-b','--base_dir', type=str, default=os.getcwd(), help='directory in which to look for mosaicked .h5 files')
    parser.add_argument('-rr','--region', type=str, help='2-letter region indicator \n'
                                                         '\t A(1-4): Antarctica, by quadrant \n'
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
    parser.add_argument('--ATL11_index', type=str, help='GeoIndex file pointing to ATL11 data')
    parser.add_argument('--ATL11_xover_dir', type=str, help="directory containing ATL11 crossover cycle directories")
    parser.add_argument('-list11','--ATL11_lineage_dir', type=str, help='directory in which to look for ATL11 .h5 filenames')
    parser.add_argument('--tiles_dir', type=str, help='directory in which to look for tile .h5 files')
    parser.add_argument('--avg_scale', type=str, default=None, help='average scale')
    parser.add_argument('--avg_scales', type=str, default=None, help='list of average scales, comma separated')
    parser.add_argument('--dzdt_lags', type=str, default='1,4', help='lags for which to calculate dz/dt, comma-separated list, no spaces')
    parser.add_argument('--grid_spacing','-g', type=str, help='DEM (meters),dh maps xy (meters),dh_maps time (years): comma-separated, no spaces', default='100.,1000.,0.25')
    parser.add_argument('--delta_t', type=str, help='time-step spacing, yr')
    parser.add_argument('--verbose', action='store_true')

    args, unknown = parser.parse_known_args()
    args.grid_spacing = args.grid_spacing.split(',')
    for ind, spacing in enumerate(args.grid_spacing):
        if '/' in spacing:
            temp=[*map(float, spacing.split('/'))]
            args.grid_spacing[ind] = temp[0]/temp[1]
        else:
            args.grid_spacing[ind] = float(args.grid_spacing[ind])

    if args.delta_t is None:
        args.delta_t = args.grid_spacing[2]

    if args.delta_t is not None:
        if isinstance(args.delta_t,str) and '/' in args.delta_t:
            args.delta_t = [*map(float, args.delta_t.split('/'))]
            args.delta_t = args.delta_t[0] / args.delta_t[1]

    args.dzdt_lags =  [*map(int, args.dzdt_lags.split(','))]

    if args.ATL11_lineage_dir is None:
        # if ATL11 lineage_dir is not specified, assume that the grandparent of the ATL11_index works
        args.ATL11_lineage_dir=os.path.dirname(os.path.dirname(args.ATL11_index))

    if args.tiles_dir is None:
        args.tiles_dir=os.path.join(args.base_dir, 'prelim')

    print('args',args)
    if args.avg_scales is not None:
        for avg_scale in [None]+[*map(int, args.avg_scales.split(','))]:
            args.avg_scale=avg_scale
            fileout = ATL15_write2nc_monthly(args)
    else:
        fileout = ATL15_write2nc_monthly(args)

if __name__=='__main__':
    main()
