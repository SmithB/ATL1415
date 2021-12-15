#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:45:47 2020

@author: ben05
"""
import numpy as np
from scipy import stats
from netCDF4 import Dataset
import shutil
import uuid
import io, re, os, glob
import h5py

import matplotlib.pyplot as plt
#import warnings
#warnings.filterwarnings('ignore',category=DeprecationWarning)
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import osgeo
import imageio
import datetime as dt

def ATL15_browse_plots(args):
        
    # get projection
    if args.Hemisphere==1:
        # set cartopy projection as EPSG 3413
        projection = ccrs.Stereographic(central_longitude=-45.0, 
                                        central_latitude=+90.0,
                                        true_scale_latitude=+70.0)   
#        DEM_file = '/Volumes/ice3/suzanne/arcticdem_mosaic_100m_v3.0.tif'
    else:
        # set cartopy projection as EPSG 3031
        projection = ccrs.Stereographic(central_longitude=+0.0,
                                        central_latitude=-90.0,
                                        true_scale_latitude=-71.0)  
#        DEM_file = '/Volumes/ice3/suzanne/REMA_100m_dem.tif'
        
    lagtext = {'lag1':'quarterly','lag4':'annual', 'lag8':'biennial', 'lag12':'triennial', 'lag16':'quadriennial'}

    # make a log file for errors 
    if not args.nolog:        
        log_file = '{}/ATL15_BrowsePlots_{}.log'.format(args.base_dir.rstrip('/'), dt.datetime.now().date())
        fhlog = open(log_file,'a')
     
    # list of spatial averaging ATL15 files
    avgs = ['_01km','_10km','_20km','_40km']
    for ii, ave in enumerate(avgs):
        filein = args.base_dir.rstrip('/') + '/ATL15_' + args.region + '_' + args.cycles + ave + '_' + args.Release + '_' + args.version + '.nc'
        print('Making browse figures from ',filein)
        pngfile = args.base_dir.rstrip('/') + '/ATL15_' + args.region + '_' + args.cycles + ave + '_' + args.Release + '_' + args.version + '_BRW'
    
        ds = Dataset(filein)
    #    # find group of largest lag 
    #    lag=1
    #    for grp in ds.groups:
    #        if grp.startswith('dhdt_lag'):
    #            if int(grp[8:])>lag:
    #                lag=int(grp[8:])
    #    grp = 'dhdt_lag'+str(lag)
        
        x = ds.groups['dhdt_lag1']['x']
        y = ds.groups['dhdt_lag1']['y']
        extent=[np.min(x),np.max(x),np.min(y),np.max(y)]
        
        dhdt = ds.groups['dhdt_lag1']['dhdt'] 
        dhdt[:][dhdt[:]==ds.groups['dhdt_lag1']['dhdt']._FillValue] = np.nan  
        dhdtmn = np.nanmean(dhdt,axis=0)
        dhdtstd = np.nanstd(dhdt,axis=0)
    
        if np.any(dhdtmn.ravel()):
            # get limits for colorbar
            h05mn = stats.mstats.scoreatpercentile(dhdtmn[~np.isnan(dhdtmn)].ravel(),5)    # dhdt is a masked array what?
            h95mn = stats.mstats.scoreatpercentile(dhdtmn[~np.isnan(dhdtmn)].ravel(),95)
            h05std = stats.mstats.scoreatpercentile(dhdtstd[~np.isnan(dhdtstd)].ravel(),5)    # dhdt is a masked array what?
            h95std = stats.mstats.scoreatpercentile(dhdtstd[~np.isnan(dhdtstd)].ravel(),95)
        else:
            fhlog.write('{}: No valid dhdt data, no browse plot written.\n'.format(filein))
            exit(-1)
    
        fig,ax = plt.subplots(1,1)
        ax = plt.subplot(1,1,1,projection=projection)
        ax.add_feature(cfeature.LAND,facecolor='0.8')
        ax.coastlines(resolution='50m',linewidth=0.5)
        ax.gridlines(crs=ccrs.PlateCarree())
        h = ax.imshow(dhdtmn, extent=extent, cmap='Spectral', vmin=h05mn, vmax=h95mn, origin='lower', interpolation='nearest')
        fig.colorbar(h,ax=ax,label='dh/dt, m',shrink=1/2, extend='both')
        ax.set_title(f'Mean quarterly dh/dt: {os.path.basename(filein)}',wrap=True)
        if args.Hemisphere==1:
            plt.figtext(0.1,0.01,f'Figure 1. Average quarterly rate of height change (dhdt_lag1/dhdt) at {ave[1:]}-resolution, in meters, from cycle {args.cycles[0:2]} to cycle {args.cycles[2:4]}. Map is plotted in a polar-stereographic projection with a central longitude of 45W and a standard latitude of 70N.',wrap=True)
        elif args.Hemisphere==-1:
            plt.figtext(0.1,0.01,f'Figure 1. Average quarterly rate of height change (dhdt_lag1/dhdt) at {ave[1:]}-resolution, in meters, from cycle {args.cycles[0:2]} to cycle {args.cycles[2:4]}. Map is plotted in a polar-stereographic projection with a central longitude of 0W and a standard latitude of 71S.',wrap=True)
        fig.savefig(f'{pngfile}_default1.png')
        
        fig,ax = plt.subplots(1,1)
        ax = plt.subplot(1,1,1,projection=projection)
        ax.add_feature(cfeature.LAND,facecolor='0.8')
        ax.coastlines(resolution='50m',linewidth=0.5)
        ax.gridlines(crs=ccrs.PlateCarree())
        h = ax.imshow(dhdtstd, extent=extent, cmap='viridis', vmin=h05std, vmax=h95std, origin='lower',interpolation='nearest')
        fig.colorbar(h,ax=ax,label='dh/dt standard deviation, m',shrink=1/2, extend='both')
        ax.set_title(f'Standard deviation of quarterly dh/dt: {os.path.basename(filein)}',wrap=True)
        if args.Hemisphere==1:
            plt.figtext(0.1,0.01,f'Figure 2. Standard deviation of quarterly rate of height change (dhdt_lag1/dhdt) at {ave[1:]}-resolution, in meters, from cycle {args.cycles[0:2]} to cycle {args.cycles[2:4]}. Map is plotted in a polar-stereographic projection with a central longitude of 45W and a standard latitude of 70N.',wrap=True)
        elif args.Hemisphere==-1:
            plt.figtext(0.1,0.01,f'Figure 2. Standard deviation of quarterly rate of height change (dhdt_lag1/dhdt) at {ave[1:]}-resolution, in meters, from cycle {args.cycles[0:2]} to cycle {args.cycles[2:4]}. Map is plotted in a polar-stereographic projection with a central longitude of 0W and a standard latitude of 71S.',wrap=True)
        fig.savefig(f'{pngfile}_default2.png')
        
    #    print(glob.glob(f'{args.base_dir.rstrip("/")}/ATL15_{args.region}_{args.cycles}{ave}_{args.Release}_{args.version}_BRW_default*.png'))
    
        # write images to browse .h5 file
        brwfile = args.base_dir.rstrip('/') + '/ATL15_' + args.region + '_' + args.cycles + ave + '_' + args.Release + '_' + args.version + '_BRW.h5'
        print(f'Making file {brwfile}') 
        if os.path.isfile(brwfile):
            os.remove(brwfile)
        shutil.copyfile('surfaceChange/resources/BRW_template.h5',brwfile)
        with h5py.File(brwfile,'r+') as hf:   
            hf.require_group('/default')
            for ii, name in enumerate(sorted(glob.glob(f'{args.base_dir.rstrip("/")}/ATL15_{args.region}_{args.cycles}{ave}_{args.Release}_{args.version}_BRW_default*.png'))):
                img = imageio.imread(name, pilmode='RGB')
                #print(ii,name)
                
    #            ave = os.path.basename(name).split('_')[3]
                dset = hf.create_dataset(f'default/default{ii+1}', \
                                         img.shape, data=img.data, \
                                         chunks=img.shape, \
                                         compression='gzip',compression_opts=6)
                dset.attrs['CLASS'] = np.string_('IMAGE')
                dset.attrs['IMAGE_VERSION'] = np.string_('1.2')
                dset.attrs['IMAGE_SUBCLASS'] = np.string_('IMAGE_TRUECOLOR')
                dset.attrs['INTERLACE_MODE'] = np.string_('INTERLACE_PIXEL')
            
    #    plt.show(block=False)
    #    plt.pause(0.001)
    #    input('Press enter to end.')
    #    plt.close('all')
    #    exit(-1)
    
        fhlog.close()
    
#    plt.show()
#
    
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,  fromfile_prefix_chars='@')
    parser.add_argument('-b','--base_dir', type=str, default=os.getcwd(), help='directory in which to look for dz .h5 files')
    parser.add_argument('-rr','--region', type=str, help='2-letter region indicator \n'
                                                         '\t AA: Antarctica \n'
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
    parser.add_argument('--Hemisphere','-H', type=int, default=1, help='1 for Northern, -1 for Southern')
    parser.add_argument('--mosaic', '-m', type=str)
    parser.add_argument('--out_path', '-o', type=str, help='default is ATL15_file path')
    parser.add_argument('--pdf', action='store_true', default=False, help='write images to .pdf file')
    parser.add_argument('--nolog', action='store_true', default=False, help='no writing errors to .log file')
    args, unknown = parser.parse_known_args()
    print(args)

    #-- digital elevation model
    elevation_dir = {}
    elevation_tile_index = {}
    #-- ArcticDEM
    resstr = '32m'
    res=float(resstr.strip('m'))

    elevation_dir['ArcticDEM'] = 'https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/' + resstr + '/'
    elevation_tile_index['ArcticDEM'] = 'https://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/indexes/ArcticDEM_Tile_Index_Rel7.zip'
    #-- REMA DEM
    elevation_dir['REMA'] = 'https://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v1.1/100m/'
    elevation_tile_index['REMA'] = 'https://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/REMA_Tile_Index_Rel1.1.zip'

    
    # from github.com/tsutterley/read-ICESat-2/blob/main/scripts/MPI_DEM_ICESat2_ATL11.py
    #-- PURPOSE: set the DEM model based on region 
#    def set_DEM_model(region):
#        if region == 'AA': 
#            DEM_MODEL = 'REMA'
#        else:
#            DEM_MODEL = 'ArcticDEM'
#        return DEM_MODEL
    
    #-- PURPOSE: read zip file containing index shapefiles for finding DEM tiles
    def read_DEM_index(shape, epsg, DEM_MODEL):
        #-- read the compressed shapefile and extract entities
#        shape = fiona.open('zip://{0}'.format(os.path.expanduser(index_file)))
#        epsg = shape.crs['init']
        #-- extract attribute indice for DEM tile (REMA,GIMP) or name (ArcticDEM)
        if (DEM_MODEL == 'REMA'):
            #-- REMA index file attributes:
            #-- name: DEM mosaic name for tile (file name without suffix)
            #-- tile: DEM tile identifier (IMy_IMx)
            #-- nd_value: fill value for elements with no data
            #-- resolution: DEM horizontal spatial resolution (meters)
            #-- creationda: creation date
            #-- raster: (empty)
            #-- fileurl: link to file on PGC server
            #-- spec_type: specific type (DEM)
            #-- qual: density of scenes within tile (0 to 1)
            #-- reg_src: DEM registration source (ICESat or neighbor align)
            #-- num_gcps: number of ground control points
            #-- meanresz: mean vertical residual (meters)
            #-- active: (1)
            #-- qc: (2)
            #-- rel_ver: release version
            #-- num_comp: number of components
            #-- st_area_sh: tile area (meters^2)
            #-- st_length_: perimeter length of tile (meters)
            field = 'tile'
        elif (DEM_MODEL == 'GIMP'):
            #-- GIMP index file attributes (from make_GIMP_tile_shapefile.py):
            #-- name: DEM mosaic name for tile (file name without suffix)
            #-- tile: DEM tile identifier (IMy_IMx)
            #-- nd_value: fill value for elements with no data
            #-- resolution: DEM horizontal spatial resolution (meters)
            #-- fileurl: link to file on NSIDC server
            #-- spec_type: specific type (DEM)
            #-- reg_src: DEM registration source (ICESat or neighbor align)
            #-- rel_ver: release version
            #-- num_comp: number of components
            #-- st_area_sh: tile area (meters^2)
            #-- st_length_: perimeter length of tile (meters)
            field = 'tile'
        elif (DEM_MODEL == 'ArcticDEM'):
            #-- ArcticDEM index file attributes:
            #-- objectid: DEM tile object identifier for sub-tile
            #-- name: DEM mosaic name for sub-tile (file name without suffix)
            #-- tile: DEM tile identifier (IMy_IMx) (non-unique for sub-tiles)
            #-- nd_value: fill value for elements with no data
            #-- resolution: DEM horizontal spatial resolution (meters)
            #-- creationda: creation date
            #-- raster: (empty)
            #-- fileurl: link to file on PGC server
            #-- spec_type: specific type (DEM)
            #-- qual: density of scenes within tile (0 to 1)
            #-- reg_src: DEM registration source (ICESat or neighbor align)
            #-- num_gcps: number of ground control points
            #-- meanresz: mean vertical residual (meters)
            #-- active: (1)
            #-- qc: (2)
            #-- rel_ver: release version
            #-- num_comp: number of components
            #-- st_area_sh: tile area (meters^2)
            #-- st_length_: perimeter length of tile (meters)
            field = 'name'
        #-- create python dictionary for each polygon object
        poly_dict = {}
        attrs_dict = {}
        #-- extract the entities and assign by tile name
        for i,ent in enumerate(shape.values()):
            #-- tile or name attributes
            if DEM_MODEL in ('REMA','GIMP'):
                tile = str(ent['properties'][field])
            else:
                tile, = re.findall(r'^(\d+_\d+_\d+_\d+)',ent['properties'][field])
                print('line 539',tile)
                
            #-- extract attributes and assign by tile
            attrs_dict[tile] = {}
            for key,val in ent['properties'].items():
                attrs_dict[tile][key] = val
                if key=='fileurl':
                    print('line 546', val)
            #-- upper-left, upper-right, lower-right, lower-left, upper-left
            ul,ur,lr,ll,ul2 = ent['geometry']['coordinates'].pop()
            #-- tile boundaries
            attrs_dict[tile]['xmin'] = ul[0]
            attrs_dict[tile]['xmax'] = lr[0]
            attrs_dict[tile]['ymin'] = lr[1]
            attrs_dict[tile]['ymax'] = ul[1]
            #-- extract Polar Stereographic coordinates for entity
            x = [ul[0],ur[0],lr[0],ll[0],ul2[0]]
            y = [ul[1],ur[1],lr[1],ll[1],ul2[1]]
            poly_obj = Polygon(list(zip(x,y)))
            #-- Valid Polygon may not possess overlapping exterior or interior rings
            if (not poly_obj.is_valid):
                poly_obj = poly_obj.buffer(0)
            poly_dict[tile] = poly_obj
        #-- close the file
        shape.close()
        #-- return the dictionaries of polygon objects and attributes
        return (poly_dict,attrs_dict,epsg)
    
    #-- PURPOSE: read DEM tile file from gzipped tar files
    def read_DEM_file(elevation_file, nd_value):
        #-- open file with tarfile (read)
        tar = tarfile.open(fileobj=io.BytesIO(elevation_file.read()), mode='r:gz')  # a mod from TS (io.BytesIO)
        #-- find dem geotiff file within tar file
        member, = [m for m in tar.getmembers() if re.search(r'dem\.tif',m.name)]
        print('line 547',member.name)
        #-- use GDAL memory-mapped file to read dem
        mmap_name = "/vsimem/{0}".format(uuid.uuid4().hex)
        osgeo.gdal.FileFromMemBuffer(mmap_name, tar.extractfile(member).read())
        ds = osgeo.gdal.Open(mmap_name)
        #-- read data matrix
        im = ds.GetRasterBand(1).ReadAsArray()
        fill_value = ds.GetRasterBand(1).GetNoDataValue()
        fill_value = 0.0 if (fill_value is None) else fill_value
        #-- get dimensions
        xsize = ds.RasterXSize
        ysize = ds.RasterYSize
        #-- create mask for finding invalid values
        mask = np.zeros((ysize,xsize),dtype=bool)
        indy,indx = np.nonzero((im == fill_value) | (~np.isfinite(im)) |
            (np.ceil(im) == np.ceil(fill_value)))
        mask[indy,indx] = True
        #-- verify that values are finite by replacing with nd_value
        im[indy,indx] = nd_value
        #-- get geotiff info
        info_geotiff = ds.GetGeoTransform()
        #-- calculate image extents
        xmin = info_geotiff[0]
        ymax = info_geotiff[3]
        xmax = xmin + (xsize-1)*info_geotiff[1]
        ymin = ymax + (ysize-1)*info_geotiff[5]
        #-- close files
        ds = None
        osgeo.gdal.Unlink(mmap_name)
        tar.close()
        #-- create image x and y arrays
        xi = np.arange(xmin,xmax+info_geotiff[1],info_geotiff[1])
        yi = np.arange(ymax,ymin+info_geotiff[5],info_geotiff[5])
        #-- return values (flip y values to be monotonically increasing)
        return (im[::-1,:],mask[::-1,:],xi,yi[::-1],member.name)
    
    
    
    
    ATL15_browse_plots(args) #.ATL15_file, hemisphere=args.Hemisphere, mosaic=args.mosaic, out_path=args.out_path, pdf=args.pdf)



