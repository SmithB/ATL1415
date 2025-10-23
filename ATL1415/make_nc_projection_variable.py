#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 09:12:09 2025

@author: ben
"""
import numpy as np

def make_nc_projection_variable(region, group):
    """
    create a variable in a netcdf file holding crs information

    Parameters
    ----------
    region : str
        string identifying the region for which the projection will be
        created.  Options are AK, CN, CS, GL, IS, SV, RA, A1-A4
    group : netcdf group
        group in which the variable is to be created.

    Returns
    -------
    crs_var : netcdf variable
        handle for variable that was created

    """
    if region in ['AK','CN','CS','GL','IS','SV','RA']:
        crs_var = group.createVariable('Polar_Stereographic',np.byte,())
        crs_var.standard_name = 'Polar_Stereographic'
        crs_var.grid_mapping_name = 'polar_stereographic'
        crs_var.straight_vertical_longitude_from_pole = -45.0
        crs_var.latitude_of_projection_origin = 90.0
        crs_var.standard_parallel = 70.0
        crs_var.scale_factor_at_projection_origin = 1.
        crs_var.false_easting = 0.0
        crs_var.false_northing = 0.0
        crs_var.semi_major_axis = 6378.137
        crs_var.semi_minor_axis = 6356.752
        crs_var.inverse_flattening = 298.257223563
        crs_var.spatial_epsg = '3413'
        crs_var.spatial_ref = 'PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3413"]]'
        crs_var.crs_wkt = ('PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3413"]]')
    elif region in ['A1','A2','A3','A4']:
        crs_var = group.createVariable('Polar_Stereographic',np.byte,())
        crs_var.standard_name = 'Polar_Stereographic'
        crs_var.grid_mapping_name = 'polar_stereographic'
        crs_var.straight_vertical_longitude_from_pole = 0.0
        crs_var.latitude_of_projection_origin = -90.0
        crs_var.standard_parallel = -71.0
        crs_var.scale_factor_at_projection_origin = 1.
        crs_var.false_easting = 0.0
        crs_var.false_northing = 0.0
        crs_var.semi_major_axis = 6378.137
        crs_var.semi_minor_axis = 6356.752
        crs_var.inverse_flattening = 298.257223563
        crs_var.spatial_epsg = '3031'
        crs_var.spatial_ref = 'PROJCS["WGS 84 / Antarctic Polar Stereographic",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",-71],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","3031"]]'
        crs_var.crs_wkt = ('PROJCS["WGS 84 / Antarctic Polar Stereographic",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",-71],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","3031"]]')
    return crs_var
