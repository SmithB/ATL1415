#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:45:47 2020

@author: ben05
"""
#import ATL11
import numpy as np
from scipy import stats
import sys, os, h5py, glob, csv
import io
import pointCollection as pc


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import cartopy.crs as ccrs
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature
import osgeo.gdal
#import imageio
import datetime as dt
from ATL11.h5util import create_attribute


def ATL14_write():
    dz_dict ={'x':'x',   # ATL14 varname : z0.h5 varname
              'y':'y',
              'h':'z0',
              'h_sigma':'sigma_z0',
#              'cell_area':'area',
#              'data_count':'count',
              'misfit_rms':'misfit_rms',
              'misfit_scaled_rms':'misfit_scaled_rms',
              }
    scale = {'Nx':'x',
             'Ny':'y',
             }

    # establish output file
    fileout = 'ATL14_tile01.h5'
    print('output file',fileout)
    if os.path.isfile(fileout):
        os.remove(fileout)
    with h5py.File(fileout.encode('ASCII'),'w') as fo:
        # get handle for input file with ROOT and height_change variables.
        FH = h5py.File('z0.h5','r')
        if 'z0' not in FH:
            print('no zo')
            FH.close()
            exit(-1)
        
        with open('ATL14_output_attrs.csv','r', encoding='utf-8-sig') as attrfile:
            reader=list(csv.DictReader(attrfile))
#        group_names = set([row['group'] for row in reader])
    
        attr_names=[x for x in reader[0].keys() if x != 'field' and x != 'group']
        
        # work ROOT group first
        field_names = [row['field'] for row in reader if 'ROOT' in row['group']]
        print(field_names)
        #establish variables that are dimension scales first
        for field in ['x', 'y']: 
            print('field',field)
            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
            dimensions = field_attrs[field]['dimensions'].split(',')
            data = np.array(FH['z0'][dz_dict[field]])
            data = np.nan_to_num(data,nan=np.finfo(np.dtype(field_attrs[field]['datatype'])).max)
            fillvalue = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
            dset = fo.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
            dset.make_scale(field)
            for ii,dim in enumerate(dimensions):
                print(field,dim)
                dset.dims[ii].label = scale[dim.strip()]
            for attr in attr_names:
                 if 'dimensions' not in attr and 'datatype' not in attr:
                     create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
            if field_attrs[field]['datatype'].startswith('int'):
                dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
            elif field_attrs[field]['datatype'].startswith('float'):
                dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
        
        for field in [item for item in field_names if item != 'x' and item != 'y']:
            # read attrs from .csv
            field_attrs = {row['field']: {attr_names[ii]:row[attr_names[ii]] for ii in range(len(attr_names))} for row in reader if field in row['field']}
            dimensions = field_attrs[field]['dimensions'].split(',')
            if dz_dict.get(field)!=None:   # if key in dz_dict
                data = np.squeeze(np.array(FH['z0'][dz_dict[field]]))
            else:
                data = np.ndarray(shape=tuple([ii+1 for ii in range(len(dimensions))]),dtype=float)
            print('line 93',field, dimensions)
            data = np.nan_to_num(data,nan=np.finfo(np.dtype(field_attrs[field]['datatype'])).max)
            print('shape',data.shape)
            fillvalue = np.finfo(np.dtype(field_attrs[field]['datatype'])).max
            dset = fo.create_dataset(field.encode('ASCII'),data=data,fillvalue=fillvalue,chunks=True,compression=6,dtype=field_attrs[field]['datatype'])
            for ii,dim in enumerate(dimensions):
                print('line 98',ii,dim)
                dset.dims[ii].label = scale[dim.strip()]
                dset.dims[ii].attach_scale(fo[scale[dim.strip()]])
            for attr in attr_names:
                 if 'dimensions' not in attr and 'datatype' not in attr:
                     create_attribute(dset.id, attr, [], str(field_attrs[field][attr]))
            if field_attrs[field]['datatype'].startswith('int'):
                dset.attrs['_FillValue'.encode('ASCII')] = np.iinfo(np.dtype(field_attrs[field]['datatype'])).max
            elif field_attrs[field]['datatype'].startswith('float'):
                dset.attrs['_FillValue'.encode('ASCII')] = np.finfo(np.dtype(field_attrs[field]['datatype'])).max

        FH.close()

#   
    return fileout
    
if __name__=='__main__':
    import argparse
#    parser=argparse.ArgumentParser()
#    parser.add_argument('ATL11_file', type=str)
#    parser.add_argument('--Hemisphere','-H', type=int, default=1, help='1 for Norhtern, -1 for Southern')
#    parser.add_argument('--mosaic', '-m', type=str)
#    parser.add_argument('--out_path', '-o', type=str, help='default is ATL11_file path')
#    parser.add_argument('--pdf', action='store_true', default=False, help='write images to .pdf file')
#    parser.add_argument('--nolog', action='store_true', default=False, help='no writing errors to .log file')
#    args=parser.parse_args()
    fileout = ATL14_write()



