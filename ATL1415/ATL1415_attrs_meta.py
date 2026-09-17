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
import timescale
from importlib import resources
import warnings
from datetime import datetime, timedelta
from ATL1415.version import softwareVersion,softwareDate,softwareTitle,identifier,series_version

def write_atl1415meta(dst,fileout,ncTemplate,args):

    # setup basic dictionary of attributes to touch
    root_info={'date_created':'', 'fileName':'', \
        'geospatial_bounds':'SET_BY_PGE', 'geospatial_bounds_crs':'SET_BY_PGE', 'geospatial_lat_max':0., \
        'geospatial_lat_min':0., 'geospatial_lon_max':0., 'geospatial_lon_min':0., \
        'netcdfversion':'', 'history':'SET_BY_PGE', \
        'identifier_product_format_version':'SET_BY_PGE', 'time_coverage_duration':0., \
        'time_coverage_end':'', 'time_coverage_start':'', 'identifier_file_uuid':''}

    # copy attributes, dimensions, variables, and groups from template
    #if 'ATL15' in os.path.basename(fileout):
    #    ncTemplate = ncTemplate.replace('atl14','atl15')
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
    # set the starting time and time span:
    set_time_range(dst, root_info, args)
    # build ATL11 lineage
    set_lineage(dst, root_info, args)
    # lat/lon bounds
    set_geobounds(dst,fileout,root_info)
    print(fileout)
    if os.path.basename(fileout).lower().startswith('atl14'):
        # output file format:
        # ATL14_IS_0331_100m_006_01.nc
        out_file_match = re.compile('(ATL\d\d)_(\D.)_(\d\d)(\d\d)_(.*)_(\d\d\d)_(\d\d).nc')\
                           .search(os.path.basename(fileout)).groups()
        out_file_attrs ={'shortname':out_file_match[0],
                     'region':out_file_match[1],
                     'c0':out_file_match[2],
                     'c1':out_file_match[3],
                     'resolution':out_file_match[4],
                     'release':out_file_match[5],
                     'version':out_file_match[6]}
    elif os.path.basename(fileout).lower().startswith('atl15'):
        # output file format:
        #ATL15_IS_0331_3mo_10km_006_01.nc
        out_file_match = re.compile('(ATL\d\d)_(\D.)_(\d\d)(\d\d)_(.*mo)_(.*m)_(\d\d\d)_(\d\d).nc')\
                           .search(os.path.basename(fileout)).groups()
        out_file_attrs ={'shortname':out_file_match[0],
                     'region':out_file_match[1],
                     'c0':out_file_match[2],
                     'c1':out_file_match[3],
                     'time_res':out_file_match[4],
                     'resolution':out_file_match[5],
                     'release':out_file_match[6],
                     'version':out_file_match[7]}
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
    dst['METADATA/DatasetIdentification'].setncattr('shortName', out_file_attrs['shortname'])
    dst['METADATA/DatasetIdentification'].setncattr('VersionID', out_file_attrs['release'])
    # Add RevisionID
    dst['METADATA/DatasetIdentification'].setncattr('RevisionID', out_file_attrs['version'])
    root_info.update({'identifier_product_format_version': series_version()})
    dst['METADATA/SeriesIdentification'].setncattr('VersionID', series_version())
    dst['METADATA/ProcessStep/PGE'].setncattr('softwareDate', softwareDate())
    dst['METADATA/ProcessStep/PGE'].setncattr('softwareTitle', softwareTitle())
    dst['METADATA/ProcessStep/PGE'].setncattr('softwareVersion', softwareVersion())
    # apply dict of root level attributes
    for key, keyval in root_info.items():
        dst.setncattr(key, keyval)

# Lineage attributes that only the granule itself can supply.  TEMPORARY
# (docs/plan_IS_run.sh I9g2): the netCDF step no longer opens ATL11 granules.
# The prelim step is to record these in the tile metadata at solve time; until
# it does, they stay 'NOT_SET', marking them invalid in the product.
FILE_ONLY_LINEAGE_ATTRS = {
    'along-track': ['uuid', 'start_geoseg', 'end_geoseg', 'start_orbit', 'end_orbit'],
    'xo': ['uuid', 'start_geoseg', 'end_geoseg', 'start_rgt', 'end_rgt']}

def attributes_for_ATL11_file(file):
    """
    Lineage attributes for one ATL11 or ATL11XO file, from its NAME alone.

    Attributes that need the granule opened (FILE_ONLY_LINEAGE_ATTRS) are left
    'NOT_SET'.

    inputs:
        file: basename of the ATL11 or ATL11XO file, as in a tile's
            meta/input_files
    outputs:
        fa: dict of lineage attributes
        this_format: 'along-track' or 'xo'
    """
    # regular expression for extracting ATL11 parameters
    rx = re.compile(r'(ATL\d{2})_(\d{4})(\d{2})_(\d{2})(\d{2})_(\d{3})_(\d{2}).*?.h5$')
    rx_xo = re.compile(r'(ATL11XO)_.._E.*_N.*_c(\d\d)_(\d\d\d)_(\d\d).h5$', flags=re.I)
    lineage_attrs=['end_cycle', 'end_geoseg', 'end_orbit', 'end_region', 'end_rgt',
                    'fileName', 'shortName', 'start_cycle', 'start_geoseg',
                    'start_orbit', 'start_region', 'start_rgt',
                    'uuid', 'version', 'release']
    # initialize the file_attribute (fa) dict:
    fa= {attr : 'NOT_SET' for attr in lineage_attrs}
    fa['fileName'] = os.path.basename(file)
    # extract attributes from filename
    m = rx.search(file)
    if m is not None:
        fa['shortName'], \
        fa['start_rgt'], \
        fa['start_region'],\
        fa['start_cycle'],\
        fa['end_cycle'],\
        fa['release'],\
        fa['version'] = m.groups()
        #start_region, end_region, start_orbit, end_orbit are not defined for an ATL11xo file
        fa['end_region'] = fa['start_region']
        # an along-track granule covers one rgt, the one in its name.  NOT for
        # ATL11XO, whose start_rgt and end_rgt differ (e.g. 238 and 1381)
        fa['end_rgt'] = fa['start_rgt']
        this_format='along-track'
    else:
        m = rx_xo.search(file)
        if m is None:
            raise ValueError(f'attributes_for_ATL11_file: {file} is neither an ATL11 '
                             'nor an ATL11XO file name')
        fa['shortName'],\
        fa['start_cycle'],\
        fa['release'],\
        fa['version'] = m.groups()
        fa['end_cycle'] = fa['start_cycle']
        this_format='xo'

    return fa, this_format

# To recursively step through groups
def walktree(top):
    yield top.groups.values()
    for value in top.groups.values():
        yield from walktree(value)

def set_lineage(dst,root_info,args):
    tilepath = args.tiles_dir
# list of lineage attributes
    lineage = []
    ATL11_files={}
    for tile in glob.iglob(os.path.join(tilepath,'*.h5')):
        try:
            with h5py.File(tile,'r') as h5f:
                inputs=str(h5f['/meta/'].attrs['input_files'])
                if inputs[:1]=='b':
                    inputs=inputs[1:]
                inputs=inputs.replace("'",'')
        except Exception:
            print("ATL14_attrs_meta.py: failed to open tile file : "+tile)
            continue
        # a tile that read no ATL11 (a matched tile) has input_files == ''
        for file in filter(None, inputs.split(',')):
            ATL11_files.setdefault(file, tile)
    invalid={}
    for file, tile in ATL11_files.items():
        try:
            fa, this_format = attributes_for_ATL11_file(file)
        except ValueError as e:
            raise ValueError(f'{e} (listed in {tile})') from e
        invalid.setdefault(this_format, 0)
        invalid[this_format] += 1
        # add attributes to list, if not already present
        if fa not in lineage:
            lineage.append(fa)
    for this_format, count in invalid.items():
        print(f'set_lineage: WARNING: lineage is INVALID for {count} {this_format} files: '
              f'{", ".join(FILE_ONLY_LINEAGE_ATTRS[this_format])} are NOT_SET '
              '(not yet recorded in the tiles; plan_IS_run.sh I9g2)')

    # reduce to unique lineage attributes (no repeat files)
    #    sorted(set(lineage))
    slineage={ key:[] for key in lineage[0] }
    for l_i in sorted(lineage, key=lambda x: (x['fileName'])):
        for key, val in l_i.items():
            slineage[key].append(val)
    for field, val in slineage.items():
        dst['METADATA/Lineage/ATL11'].setncattr(field, val)

# set time range
def set_time_range(dst, root_info, args):
    # set the nominal time range of the product.
    # N.B.  This is coded explicitly to satisfy requirements from SIPS and NSIDC
    t_span = args.t_crop
    datetime_start = datetime(2019, 1, 1, 0, 0, 0) + timedelta(days = (t_span[0]-2019.0)*365.25)
    datetime_end = datetime(2019, 1, 1, 0, 0, 0) + timedelta(days = (t_span[1]-2019.0)*365.25)

    # Give each region a unique time offset:
    unicode_vals=[]
    for char in args.region:
        unicode_vals.append(ord(char))
    datetime_start += timedelta( seconds = int(unicode_vals[0]*3+unicode_vals[1]*2) )
    # convert starting and ending delta times to UTC
    sUTCtime = (str(datetime_start.date())+'T'+
                    datetime_start.strftime("%H:%M:%S.%f")+'Z')
    eUTCtime = (str(datetime_end.date())+'T'+
                    datetime_end.strftime("%H:%M:%S.%f")+'Z')
    if args.verbose:
        print(f"ATL1415_attrs_meta: UTC time range: {sUTCtime} - {eUTCtime}")

    # set time attributes
    root_info.update({'time_coverage_start': sUTCtime})
    root_info.update({'time_coverage_end': eUTCtime})
    root_info.update({'time_coverage_duration': int((datetime_start-datetime_end).seconds)})
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
      # Build csv style geospatial bounds for root attribute (lat lon, lat lon)
      geo_bounds_str = ','.join(f"{a},{b}" for a, b in zip(y, x))
      root_info.update({'geospatial_bounds':'POLYGON(('+geo_bounds_str+'))'})

    except (FileNotFoundError, KeyError):
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
    root_info.update({'geospatial_bounds_crs': 'EPSG:4326'})

