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

def attributes_for_ATL11_file(file, args):

    # regular expression for extracting ATL11 parameters
    rx = re.compile(r'(ATL\d{2})_(\d{4})(\d{2})_(\d{2})(\d{2})_(\d{3})_(\d{2}).*?.h5$')
    lineage_attrs=['end_cycle', 'end_geoseg', 'end_orbit', 'end_region', 'end_rgt',
                    'fileName', 'shortName', 'start_cycle', 'start_geoseg',
                    'start_orbit', 'start_region', 'start_rgt',
                    'uuid', 'version', 'release']
    # initialize the file_attribute (fa) dict:
    fa= {attr : 'NOT_SET' for attr in lineage_attrs}
    fa['fileName'] = os.path.basename(file)
    # extract attributes from filename
    try:
        fa['shortName'], \
        fa['start_rgt'], \
        fa['start_region'],\
        fa['start_cycle'],\
        fa['end_cycle'],\
        fa['release'],\
        fa['version'] = rx.search(file).groups()
        atl11path = args.ATL11_lineage_dir
        this_format='along-track'
    except Exception:
        rx=re.compile('(ATL11XO)_.._E.*_N.*_c(\d\d)_(\d\d\d)_(\d\d).h5$', flags=re.I)
        fa['shortName'],\
        fa['start_cycle'],\
        fa['release'],\
        fa['version'] = rx.search(file).groups()
        fa['end_cycle'] = fa['start_cycle']
        atl11path = os.path.join(args.ATL11_xover_dir, f'cycle_{fa["start_cycle"]}')
        if not os.path.isfile( os.path.join(atl11path,file) ):
            # we may be using a tiling schema to point to the file location
            schema_files = glob.glob(os.path.join(atl11path, '*til*.json'))
            if not schema_files:
                raise FileNotFoundError(f"could not find {file}")
            with open(schema_files[0],'r') as fh:
                atl11path = json.load(fh)['directory']
        this_format='xo'

    with h5py.File(os.path.join(atl11path,file),'r') as fileID:
        # extract ATL11 attributes from files
        fa['uuid'] = fileID['METADATA']['DatasetIdentification'].attrs['uuid'].decode('utf-8')
        fa['start_geoseg'] = fileID['ancillary_data/start_geoseg'][0]
        fa['end_geoseg'] = fileID['ancillary_data/end_geoseg'][0]
        if this_format=='xo':
            # start_rgt and end_rgt are not in the filename, read them from the file
            fa['start_rgt'] = fileID['ancillary_data/start_rgt'][0]
            fa['end_rgt'] = fileID['ancillary_data/end_rgt'][0]
        else:
            #start_region, end_region, start_orbit, end_orbit are not defined for an ATL11xo file
            fa['start_orbit'] = fileID['ancillary_data/start_orbit'][0]
            fa['end_orbit'] = fileID['ancillary_data/end_orbit'][0]
            fa['end_region'] = fa['start_region']
        sdeltatime = fileID['ancillary_data/start_delta_time'][0]
        edeltatime = fileID['ancillary_data/end_delta_time'][0]

    fa['end_rgt'] = fa['start_rgt']

    return fa

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
    ATL11_files=set()
    for tile in glob.iglob(os.path.join(tilepath,'*.h5')):
        try:
            with h5py.File(tile,'r') as h5f:
                inputs=str(h5f['/meta/'].attrs['input_files'])
                if inputs[0]=='b':
                    inputs=inputs[1:]
                ATL11_files.update(inputs.replace("'",'').split(','))
        except Exception:
            print("ATL14_attrs_meta.py: failed to open tile file : "+tile)
    for file in ATL11_files:
        fa = attributes_for_ATL11_file(file, args)
        # add attributes to list, if not already present
        if fa not in lineage:
            lineage.append(fa)

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

