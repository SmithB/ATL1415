#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 15:28:08 2022

@author: ben
"""

import numpy as np
import pointCollection as pc
from osgeo import ogr, gdal

def make_mask_from_vector(mask_file, W, ctr, spacing, srs_proj4=None):
    
    x,y = [np.arange(ctr[dim]-W[dim]/2, ctr[dim]+W[dim]/2+spacing, spacing) \
           for dim in ['x','y']]

    z=np.zeros((len(y), len(x)))
    mask_ds=pc.grid.data().from_dict({'x':x,'y':y,'z':z}).to_gdal(srs_proj4=srs_proj4)
    
    #with ogr.Open(mask_file) as vector_mask:
    vector_mask=ogr.Open(mask_file, 0) # readonly
    mask_layer=vector_mask.GetLayer()
    gdal.RasterizeLayer(mask_ds, [1], mask_layer, burn_values=[1])
    #vector_mask=None  # this statement crashes the python kernel, at least in debug mode.  Perhaps we don't need it.  ??
    
    mask=pc.grid.data().from_dict({
        'x' : x,\
        'y' : y,\
        'z' : np.flipud(mask_ds.GetRasterBand(1).ReadAsArray())})
    
    return mask