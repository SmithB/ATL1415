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



def ATL15_write2nc(args):

    def make_dataset(field,fieldout,data,field_attrs,file_obj,group_obj,nctype,dimScale=False):
        # where field is the name from ATL15_output_attrs.csv file
        # where fieldout is the name of the output variable in the .nc file
        dimensions = field_attrs[field]['dimensions'].split(',')
        dimensions = tuple(x.strip() for x in dimensions)
        if field_attrs[field]['datatype'].startswith('int'):
            fill_value = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
        elif field_attrs[field]['datatype'].startswith('float'):
            fill_value = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
        data = np.nan_to_num(data,nan=fill_value)

        if dimScale:
            group_obj.createDimension(field_attrs[field]['dimensions'],data.shape[0])

        dsetvar = group_obj.createVariable(fieldout,
                                           nctype[field_attrs[field]['datatype']],
                                           dimensions,
                                           fill_value=fill_value, zlib=True,
                                           least_significant_digit=ast.literal_eval(field_attrs[field]['least_significant_digit']))
                                           # significant_digits=ast.literal_eval(field_attrs[field]['least_significant_digit'])) DOESN'T WORK

        dsetvar[:] = data
        for attr in attr_names:
            if attr != 'group description':
                dsetvar.setncattr(attr,field_attrs[field][attr])
        # add attributes for projection
        if not field.startswith('time'):
            dsetvar.setncattr('grid_mapping','Polar_Stereographic')
        if field == 'x':
            dsetvar.standard_name = 'projection_x_coordinate'
        if field == 'y':
            dsetvar.standard_name = 'projection_y_coordinate'

        return file_obj

    # find which lags are in attributes file
    with importlib.resources.open_text('ATL1415.resources', 'ATL15_monthly_output_attrs.csv', encoding='utf-8-sig') as attrfile:
        with open(attrFile,'r',encoding='utf-8-sig') as attrfile:
            reader=list(csv.DictReader(attrfile))
    qtrs = []
    for row in reader:
        if row['group'] == 'height_change' and row['field'].startswith('time'):
            qtrs.append(re.search('(time)(.*)',row['field']).group(2))
    # print(qtrs)

    # fileslag = glob.glob(os.path.join(args.base_dir,'*lag*.h5'))
    # qtrs=[]
    # for filesl in fileslag:
    #     qtrs.append(re.search('lag(\d+).',filesl).group(1))
    # ll = sorted([int(x) for x in list(set(qtrs))])
    # print('lags available for this region',ll)
    # print()
    #
    # lagkeys = [f'_lag{x}' for x in ll]
    # print(lagkeys.insert(0,''))
    # exit(-1)

    dz_dict = {}
    for qtr in qtrs:
        dz_dict[f'time{qtr}'] = 't'   # {ATL15 outgoing var name: hdf5 incoming var name}
    dz_dict2 = {'x':'x',                # {ATL15 outgoing var name: hdf5 incoming var name}
                'y':'y',
                'ice_area':'cell_area',
                'data_count':'count',
                'misfit_rms':'misfit_rms',
                'misfit_scaled_rms':'misfit_scaled_rms',
                'delta_h':'dz',
                'delta_h_sigma':'sigma_dz',
                'delta_h_10km':'avg_dz_10000m',
                'delta_h_sigma_10km':'sigma_avg_dz_10000m',
                'delta_h_20km':'avg_dz_20000m',
                'delta_h_sigma_20km':'sigma_avg_dz_20000m',
                'delta_h_40km':'avg_dz_40000m',
                'delta_h_sigma_40km':'sigma_avg_dz_40000m',
                }
    dz_dict.update(dz_dict2)
    # print(dz_dict)

    nctype = {'float64':'f8',
              'float32':'f4',
              'int8':'i1'}
    lags = {'file': [f'FH{qtr}' for qtr in qtrs]}
    lags['vari'] = [qtr for qtr in qtrs]
    lags['varigrp'] = ['delta_h' if qtr=='' else 'dhdt'+qtr for qtr in qtrs]
    # print(lags)

    # lags = {
    #         'file' : ['FH','FH_lag1','FH_lag4','FH_lag8','FH_lag12'],
    #         'vari' : ['','_lag1','_lag4','_lag8','_lag12'],
    #         'varigrp' : ['delta_h','dhdt_lag1','dhdt_lag4','dhdt_lag8','dhdt_lag12']
    #        }
    avgs = ['','_10km','_20km','_40km']
    # # open data attributes file
    # attrFile = pkg_resources.resource_filename('ATL1415','resources/ATL15_output_attrs.csv')
    # with open(attrFile,'r',encoding='utf-8-sig') as attrfile:
    #     reader=list(csv.DictReader(attrfile))

    attr_names=[x for x in reader[0].keys() if x != 'field' and x != 'group']

    for ave in avgs:
            # establish output file, one per average
            if ave=='':
                fileout = args.base_dir.rstrip('/') + '/ATL15_' + args.region + '_' + args.cycles + '_01km_' + args.Release + '_' + args.version + '.nc'
            else:
                fileout = args.base_dir.rstrip('/') + '/ATL15_' + args.region + '_' + args.cycles + ave + '_' + args.Release + '_' + args.version + '.nc'
            print('output file:',fileout)

            with Dataset(fileout,'w',clobber=True) as nc:
                nc.setncattr('GDAL_AREA_OR_POINT','Area')
                nc.setncattr('Conventions','CF-1.6')                
                tilegrp = make_tile_stats_group(nc, args)

                ice_area_mask=None
                # loop over dz*.h5 files for one ave
                for jj in range(len(lags['file'])):
                    if jj==0:
                        filein = args.base_dir.rstrip('/')+'/dz'+ave+lags['vari'][jj]+'.h5'
                    else:
                        filein = args.base_dir.rstrip('/')+'/dzdt'+ave+lags['vari'][jj]+'.h5'

                    if not os.path.isfile(filein):
                        print('No file:',args.base_dir.rstrip('/')+'/'+os.path.basename(filein))
                        continue
                    else:
                        print('Reading file:',args.base_dir.rstrip('/')+'/'+os.path.basename(filein))
                    lags['file'][jj] = h5py.File(filein,'r')  # file object
                    dzg=list(lags['file'][jj].keys())[0]      # dzg is group in input file

                    nc.createGroup(lags['varigrp'][jj])

                    # make projection variable for each group
                    crs_var = make_nc_projection_variable(args.region,nc.groups[lags['varigrp'][jj]])

                    # dimension scales for each group
                    for field in ['x','y']:
                        data = np.array(lags['file'][jj][dzg][dz_dict[field]])
                        if field == 'x':
                            x = data
                            xll = np.min(x)
                            dx = x[1]-x[0]
                        if field == 'y':
                            y = data
                            yll = np.max(y)
                            dy = y[0]-y[1]
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                        make_dataset(field,field,data,field_attrs,nc,nc.groups[lags['varigrp'][jj]],nctype,dimScale=True)
                    crs_var.GeoTransform = (xll,dx,0,yll,0,dy)

                    if jj==0:  # no lag
                        field = 'time'
                        data = np.array(lags['file'][jj][dzg]['t'])
                        # convert to decimal days from 1/1/2018
                        data = (data-2018.)*365.25
                        ntime = data.shape[0]
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                        make_dataset(field,field,data,field_attrs,nc,nc.groups[lags['varigrp'][jj]],nctype,dimScale=True)

                        for fld in ['ice_area','delta_h','delta_h_sigma','data_count','misfit_rms','misfit_scaled_rms']:  # fields that can be ave'd but not lagged

                            if (len(ave) > 0) and (fld.startswith('misfit') or fld=='data_count'): # not in ave'd groups  or fld=='ice_mask'
                                break
                            field = fld+ave
                            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field == row['field'] if row['group']=='height_change'+ave}

                            # get data from .h5
                            if fld.startswith('delta_h'):  # fields with complicated name changes
                                #print("from:" + lags['file'][jj])
                                print("\t reading:" + str([dzg, dz_dict[field]]))
                                data = np.array(lags['file'][jj][dzg][dz_dict[field]])
                                data[np.isnan(ice_area_mask)] = np.nan
                                if fld=='delta_h':  # add group description
                                    nc.groups[lags['varigrp'][jj]].setncattr('description',field_attrs[field]['group description'])
                            else:
                                data = np.array(lags['file'][jj][dzg][dz_dict[fld]])
                                if fld == 'ice_area':
                                    data[data==0.0] = np.nan
                                    ice_area_mask = data # where ice_area is invalid, so are delta_h and dhdt variables.
                            if len(data.shape)==3:
                                data = np.moveaxis(data,2,0)  # results in t, y, x

                            make_dataset(field,fld,data,field_attrs,nc,nc.groups[lags['varigrp'][jj]],nctype,dimScale=False)

                    else:  # one of the lags
                        field = 'time'+lags['vari'][jj]
                        data = np.array(lags['file'][jj][dzg]['t'])
                        # convert to decimal days from 1/1/2018
                        data = (data-2018.)*365.25
                        ntime = data.shape[0]
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                        make_dataset(field,'time',data,field_attrs,nc,nc.groups[lags['varigrp'][jj]],nctype,dimScale=True)

                        # get ice_area first, because that's the mask for dhdt and _sigma.
                        field = 'ice_area'+lags['vari'][jj]+ave
                        data = np.array(lags['file'][jj][dzg]['cell_area'])
                        data[data==0.0] = np.nan
                        data = np.moveaxis(data,2,0)  # t, y, x
                        mask_lag = data;
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                        make_dataset(field,'ice_area',data,field_attrs,nc,nc.groups[lags['varigrp'][jj]],nctype,dimScale=False)

                        field = 'dhdt'+lags['vari'][jj]+ave
                        data = np.array(lags['file'][jj][dzg][dzg])
                        data = np.moveaxis(data,2,0)  # t, y, x
                        data[np.isnan(mask_lag)] = np.nan
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                        make_dataset(field,'dhdt',data,field_attrs,nc,nc.groups[lags['varigrp'][jj]],nctype,dimScale=False)
                        # add group description
                        nc.groups[lags['varigrp'][jj]].setncattr('description',field_attrs[field]['group description'])

                        field = 'dhdt'+lags['vari'][jj]+'_sigma'+ave
                        data = np.array(lags['file'][jj][dzg]['sigma_'+dzg])
                        data = np.moveaxis(data,2,0)  # t, y, x
                        data[np.isnan(mask_lag)] = np.nan
                        field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field'] if row['group']=='height_change'+ave}
                        make_dataset(field,'dhdt_sigma',data,field_attrs,nc,nc.groups[lags['varigrp'][jj]],nctype,dimScale=False)

                for jj in range(len(lags['file'])):
                    try:
                        lags['file'][jj].close()
                    except:
                        pass
                #BS: converted from pkg_resources
                with importlib.resources.path('ATL1415.resources', 'atl15_metadata_template.nc') as ncTemplate:
                    ATL14_attrs_meta.write_atl14meta(nc, fileout, ncTemplate, args)

    return fileout


if __name__=='__main__':
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
    args, unknown = parser.parse_known_args()

    if args.ATL11_lineage_dir is None:
        # if ATL11 lineage_dir is not specified, assume that the grandparent of the ATL11_index works
        args.ATL11_lineage_dir=os.path.dirname(os.path.dirname(args.ATL11_index))

    if args.tiles_dir is None:
        args.tiles_dir=os.path.join(args.base_dir, 'prelim')

    print('args',args)

    fileout = ATL15_write2nc(args)
