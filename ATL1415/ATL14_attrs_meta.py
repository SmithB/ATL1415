#!/usr/bin/env python3

import numpy as np
import numpy.ma as ma
import sys, os
from netCDF4 import Dataset
import netCDF4
import h5py
from osgeo import osr, ogr
import csv
import json
import re
import glob
import uuid
from importlib import resources
import warnings
from datetime import datetime
from ATL1415.version import softwareVersion,softwareDate,softwareTitle,identifier,series_version

def write_atl14meta(dst,fileout,ncTemplate,args):

    # setup basic dictionary of attributes to touch
    root_info={'date_created':'', 'fileName':'', 'geospatial_lat_max':0., \
        'geospatial_lat_min':0., 'geospatial_lon_max':0., 'geospatial_lon_min':0., \
        'netcdfversion':'', 'history':'SET_BY_PGE', \
        'identifier_product_format_version':'SET_BY_PGE', 'time_coverage_duration':0., \
        'time_coverage_end':'', 'time_coverage_start':'', 'identifier_file_uuid':''}

    # copy attributes, dimensions, variables, and groups from template
    if 'ATL15' in os.path.basename(fileout):
        ncTemplate = ncTemplate.replace('atl14','atl15')
    with Dataset(ncTemplate,'r') as src:
    # copy attributes
        for name in src.ncattrs():
            dst.setncattr(name, src.getncattr(name))
    # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited else None))
    # copy variables
        for name, variable in src.variables.items():
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst.variables[name][:] = src.variables[name][:]
            for attribute in src.variables[name].ncattrs():
                dst.variables[name].setncattr(attribute, src.variables[name].getncattr(attribute))
    # copy groups, recursively
        for grp in walktree(src):
            for child in grp:
                dg = dst.createGroup(child.path)
                for name in child.ncattrs():
                    dg.setncattr(name,child.getncattr(name))
                for name, dimension in child.dimensions.items():
                    dg.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
                for name, variable in child.variables.items():
                    x = dg.createVariable(name, variable.datatype, variable.dimensions)
                    dg.variables[name][:] = child.variables[name][:]
                    for attribute in child.variables[name].ncattrs():
                        dg.variables[name].setncattr(attribute, child.variables[name].getncattr(attribute))
    # build ATL11 lineage
    set_lineage(dst,root_info,args)
    # lat/lon bounds
    set_geobounds(dst,fileout,root_info)

    # set file and date attributes
    root_info.update({'netcdfversion': netCDF4.__netcdf4libversion__})
    root_info.update({'identifier_file_uuid': str(uuid.uuid4())})
    dst['METADATA/DatasetIdentification'].setncattr('uuid', str(uuid.uuid4()).encode('ASCII'))
    dateval = str(datetime.now().date())
    dateval = dateval+'T'+str(datetime.now().time())+'Z'
    root_info.update({'date_created': dateval})
    root_info.update({'history': dateval})
    dst['METADATA/DatasetIdentification'].setncattr('creationDate', str(datetime.now().date()))
    root_info.update({'fileName': os.path.basename(fileout)})
    dst['METADATA/DatasetIdentification'].setncattr('fileName', os.path.basename(fileout))
    dst['METADATA/DatasetIdentification'].setncattr('VersionID', os.path.basename(fileout).split('_')[4])
    root_info.update({'identifier_product_format_version': series_version()})
    dst['METADATA/SeriesIdentification'].setncattr('VersionID', series_version())
    dst['METADATA/ProcessStep/PGE'].setncattr('softwareDate', softwareDate())
    dst['METADATA/ProcessStep/PGE'].setncattr('softwareTitle', softwareTitle())
    dst['METADATA/ProcessStep/PGE'].setncattr('softwareVersion', softwareVersion())
    # apply dict of root level attributes
    for key, keyval in root_info.items():
        dst.setncattr(key, keyval)

# To recursively step through groups
def walktree(top):
    yield top.groups.values()
    for value in top.groups.values():
        yield from walktree(value)

def set_lineage(dst,root_info,args):
    tilepath = args.tiles_dir
    atl11path = args.ATL11_lineage_dir
# list of lineage attributes
    lineage = []
# regular expression for extracting ATL11 parameters
    rx = re.compile(r'(ATL\d{2})_(\d{4})(\d{2})_(\d{2})(\d{2})_(\d{3})_(\d{2})(.*?).h5$')
# For each tile:
    min_start_delta_time = np.finfo(np.float64()).max
    max_end_delta_time = np.finfo(np.float64()).tiny
    ATL11_files=set()
    for tile in glob.iglob(os.path.join(tilepath,'*.h5')):
        with h5py.File(tile,'r') as h5f:
            inputs=str(h5f['/meta/'].attrs['input_files'])
            if inputs[0]=='b':
                inputs=inputs[1:]
            ATL11_files.update(inputs.replace("'",'').split(','))
    for FILE in ATL11_files:
        # extract parameters from filename
        PRD,TRK,GRAN,SCYC,ECYC,RL,VERS,AUX = rx.findall(FILE).pop()
        with h5py.File(os.path.join(atl11path,FILE),'r') as fileID:
            # extract ATL11 attributes from files
            UUID = fileID['METADATA']['DatasetIdentification'].attrs['uuid'].decode('utf-8')
            SGEOSEG = fileID['ancillary_data/start_geoseg'][0]
            EGEOSEG = fileID['ancillary_data/end_geoseg'][0]
            SORBIT = fileID['ancillary_data/start_orbit'][0]
            EORBIT = fileID['ancillary_data/end_orbit'][0]
            sdeltatime = fileID['ancillary_data/start_delta_time'][0]
            edeltatime = fileID['ancillary_data/end_delta_time'][0]
            # track earliest and latest delta time and UTC
            if sdeltatime < min_start_delta_time:
                sUTCtime = fileID['ancillary_data/data_start_utc'][0].decode('utf-8')
                min_start_delta_time = sdeltatime
            if edeltatime > max_end_delta_time:
                eUTCtime = fileID['ancillary_data/data_end_utc'][0].decode('utf-8')
                max_end_delta_time = edeltatime

        # merge attributes as a tuple
        attrs = (FILE,PRD,int(TRK),int(GRAN),int(SCYC),int(ECYC),int(VERS),UUID,int(SGEOSEG),
                int(EGEOSEG),int(SORBIT),int(EORBIT))
        # add attributes to list, if not already present
        if attrs not in lineage:
            lineage.append(attrs)
    # reduce to unique lineage attributes (no repeat files)
    #    sorted(set(lineage))

# sort and set lineage attributes
    slineage = sorted(lineage,key=lambda x: (x[0]))
    dst['METADATA/Lineage/ATL11'].setncattr('fileName',list(zip(*slineage))[0])
    dst['METADATA/Lineage/ATL11'].setncattr('shortName',list(zip(*slineage))[1])
    dst['METADATA/Lineage/ATL11'].setncattr('start_rgt',list(zip(*slineage))[2])
    dst['METADATA/Lineage/ATL11'].setncattr('end_rgt',list(zip(*slineage))[2])
    dst['METADATA/Lineage/ATL11'].setncattr('start_region',list(zip(*slineage))[3])
    dst['METADATA/Lineage/ATL11'].setncattr('end_region',list(zip(*slineage))[3])
    dst['METADATA/Lineage/ATL11'].setncattr('start_cycle',list(zip(*slineage))[4])
    dst['METADATA/Lineage/ATL11'].setncattr('end_cycle',list(zip(*slineage))[5])
    dst['METADATA/Lineage/ATL11'].setncattr('version',list(zip(*slineage))[6])
    dst['METADATA/Lineage/ATL11'].setncattr('uuid',list(zip(*slineage))[7])
    dst['METADATA/Lineage/ATL11'].setncattr('start_geoseg',list(zip(*slineage))[8])
    dst['METADATA/Lineage/ATL11'].setncattr('end_geoseg',list(zip(*slineage))[9])
    dst['METADATA/Lineage/ATL11'].setncattr('start_orbit',list(zip(*slineage))[10])
    dst['METADATA/Lineage/ATL11'].setncattr('end_orbit',list(zip(*slineage))[11])

# set time attributes
    root_info.update({'time_coverage_start': sUTCtime})
    root_info.update({'time_coverage_end': eUTCtime})
    root_info.update({'time_coverage_duration': max_end_delta_time-min_start_delta_time})
    dst['/METADATA/Extent'].setncattr('rangeBeginningDateTime',sUTCtime)
    dst['/METADATA/Extent'].setncattr('rangeEndingDateTime',eUTCtime)

# buuild lat/lon geo boundaries
def set_geobounds(dst,fileout,root_info):
    if 'ATL14' in os.path.basename(fileout):
        georoot = ''
    else:
        georoot = '/delta_h'
    polar_srs=osr.SpatialReference()
    polar_srs.ImportFromEPSG(int(dst[georoot+'/Polar_Stereographic'].getncattr('spatial_epsg')))
    ll_srs=osr.SpatialReference()
    ll_srs.ImportFromEPSG(4326)
    if hasattr(osr,'OAMS_TRADITIONAL_GIS_ORDER'):
        ll_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        polar_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ct=osr.CoordinateTransformation(polar_srs, ll_srs)

    xmin,xmax = (np.min(dst[georoot+'/x']),np.max(dst[georoot+'/x']))
    ymin,ymax = (np.min(dst[georoot+'/y']),np.max(dst[georoot+'/y']))
    N = 2
    dx = (xmax-xmin)/N
    dy = (ymax-ymin)/N

    multipoint = ogr.Geometry(ogr.wkbMultiPoint)
    for x in range(N+1):
        for y in range(N+1):
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(ymax - y*dy,xmin + x*dx)
            multipoint.AddGeometry(point)

    multipoint.Transform(ct)
    lonmin,lonmax,latmin,latmax = multipoint.GetEnvelope()
    if (lonmin == -180.0) | (lonmax == 180.0):
        lonmin,lonmax = (-180.0,180.0)
# set variables and attributes, from JSON polygons, if present
    try:
      region = os.path.basename(fileout).split("_")[1]
      polyfile = os.path.join(resources.files('ATL1415'),'resources','region_extent_polygons.json')
      with open (polyfile) as poly_f:
        poly_data = poly_f.read()
      reg_poly = region+'_poly'
      poly = json.loads(poly_data)
      x = [row[0] for row in poly[reg_poly]]
      y = [row[1] for row in poly[reg_poly]]
      dst['/orbit_info'].variables['bounding_polygon_dim1'][:] = np.arange(1,np.size(x)+1)
      dst['/orbit_info'].variables['bounding_polygon_lon1'][:] = np.array(x)[:]
      dst['/orbit_info'].variables['bounding_polygon_lat1'][:] = np.array(y)[:]
      latmin = min(np.array(y))
      latmax = max(np.array(y))
      lonmin=ma.min(ma.masked_where(abs(np.array(y)) > 88.0, np.array(x)))
      lonmax=ma.max(ma.masked_where(abs(np.array(y)) > 88.0, np.array(x)))
    except:
      warnings.filterwarnings("always")
      warnings.warn("Deprecated. Use polygon from json file", DeprecationWarning)
      dst['/orbit_info'].variables['bounding_polygon_dim1'][:] = np.arange(1,4+1)
      dst['/orbit_info'].variables['bounding_polygon_lon1'][:] = np.array([lonmin,lonmax,lonmax,lonmin])[:]
      dst['/orbit_info'].variables['bounding_polygon_lat1'][:] = np.array([latmax,latmax,latmin,latmin])[:]
    dst['/METADATA/Extent'].setncattr('westBoundLongitude',lonmin)
    dst['/METADATA/Extent'].setncattr('eastBoundLongitude',lonmax)
    dst['/METADATA/Extent'].setncattr('northBoundLatitude',latmax)
    dst['/METADATA/Extent'].setncattr('southBoundLatitude',latmin)
    root_info.update({'geospatial_lon_min': lonmin})
    root_info.update({'geospatial_lon_max': lonmax})
    root_info.update({'geospatial_lat_min': latmin})
    root_info.update({'geospatial_lat_max': latmax})

