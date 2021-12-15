#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:45:47 2020

@author: ben05
"""
import numpy as np
from scipy import stats
import os, glob, sys
from netCDF4 import Dataset
import shutil
import h5py
#import pointCollection as pc
#from PointDatabase.mapData import mapData

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from _usgs import _usgs
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import imageio
import datetime as dt

def ATL14_browse_plots(args):
    # set the color map
    cmap = _usgs()  

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

    # make a log file for errors 
    if not args.nolog:        
        log_file = '{}/ATL14_BrowsePlots_{}.log'.format(args.base_dir.rstrip('/'), dt.datetime.now().date())
        fhlog = open(log_file,'a')
     
    filein = args.base_dir.rstrip('/') + '/ATL14_' + args.region + '_' + args.cycles + '_100m_' + args.Release + '_' + args.version + '.nc'
    print('Making browse figures from ',filein)
    pngfile = args.base_dir.rstrip('/') + '/ATL14_' + args.region + '_' + args.cycles + '_100m_' + args.Release + '_' + args.version + '_BRW'

    ds = Dataset(filein)
    x = ds['x']
    y = ds['y']
    extent=[np.min(x),np.max(x),np.min(y),np.max(y)]

    h = np.array(ds['h'])
    h[h==ds['h']._FillValue] = np.nan
    h_sigma = np.array(ds['h_sigma'])
    h_sigma[h_sigma==ds['h_sigma']._FillValue] = np.nan
    
    if ~np.any(h.ravel()):   # convert to numpy array so can use ravel and other .methods()
        fhlog.write('{}: No valid h data, no browse plot written.\n'.format(filein))
        exit(-1)
        
    fig,ax = plt.subplots(1,1)    
    ax = plt.subplot(1,1,1,projection=projection)
    ax.add_feature(cfeature.LAND,facecolor='0.8')
    ax.coastlines(resolution='50m',linewidth=0.5)
    ax.gridlines(crs=ccrs.PlateCarree())
    h05 = stats.mstats.scoreatpercentile(h[~np.isnan(h)].ravel(),5)    
    h95 = stats.mstats.scoreatpercentile(h[~np.isnan(h)].ravel(),95)
    handle = ax.imshow(h, extent=extent, cmap=cmap, vmin=h05, vmax=h95, origin='lower', interpolation='nearest')
    fig.colorbar(handle,ax=ax,label='DEM, m',shrink=1/2, extend='both')
    ax.set_title(f'DEM: {os.path.basename(filein)}')
    if args.Hemisphere==1:
        plt.figtext(0.1,0.01,f'Figure 1. DEM surface height (h), referenced to WGS84, in meters, from cycle {args.cycles[0:2]} to cycle {args.cycles[2:4]}. Map is plotted in a polar-stereographic projection with a central longitude of 45W and a standard latitude of 70N.',wrap=True)
    elif args.Hemisphere==-1:
        plt.figtext(0.1,0.01,f'Figure 1. DEM surface height (h), referenced to WGS84, in meters, from cycle {args.cycles[0:2]} to cycle {args.cycles[2:4]}. Map is plotted in a polar-stereographic projection with a central longitude of 0W and a standard latitude of 71S.',wrap=True)
    fig.savefig(f'{pngfile}_default1.png')

    fig,ax = plt.subplots(1,1)    
    ax = plt.subplot(1,1,1,projection=projection)
    ax.add_feature(cfeature.LAND,facecolor='0.8')
    ax.coastlines(resolution='50m',linewidth=0.5)
    ax.gridlines(crs=ccrs.PlateCarree())
    h01 = stats.mstats.scoreatpercentile(h_sigma[~np.isnan(h_sigma)].ravel(),1)  
    h99 = stats.mstats.scoreatpercentile(h_sigma[~np.isnan(h_sigma)].ravel(),99)
    handle = ax.imshow(h_sigma, extent=extent, cmap='viridis', norm=LogNorm(vmin=h01, vmax=h99), origin='lower', interpolation='nearest') 
    fig.colorbar(handle,ax=ax,label='DEM uncertainty, m',shrink=1/2, extend='both')
    ax.set_title(f'DEM Uncertainty: {os.path.basename(filein)}')
    if args.Hemisphere==1:
        plt.figtext(0.1,0.01,f'Figure 2. Uncertainty in the DEM surface height (h_sigma), in meters, from cycle {args.cycles[0:2]} to cycle {args.cycles[2:4]}. Map is plotted in a polar-stereographic projection with a central longitude of 45W and a standard latitude of 70N.',wrap=True)
    elif args.Hemisphere==-1:
        plt.figtext(0.1,0.01,f'Figure 2. Uncertainty in the DEM surface height (h_sigma), in meters, from cycle {args.cycles[0:2]} to cycle {args.cycles[2:4]}. Map is plotted in a polar-stereographic projection with a central longitude of 0W and a standard latitude of 71S.',wrap=True)
    fig.savefig(f'{pngfile}_default2.png')

    # write images to browse .h5 file
    brwfile = args.base_dir.rstrip('/') + '/ATL14_' + args.region + '_' + args.cycles + '_100m_' + args.Release + '_' + args.version + '_BRW.h5'
    print('making file', brwfile)
    if os.path.isfile(brwfile):
        os.remove(brwfile)
    shutil.copyfile('surfaceChange/resources/BRW_template.h5',brwfile)
    with h5py.File(brwfile,'r+') as hf:  
        hf.require_group('/default')
        for ii, name in enumerate(sorted(glob.glob(f'{args.base_dir.rstrip("/")}/ATL14_{args.region}_{args.cycles}_100m_{args.Release}_{args.version}_BRW_default*.png'))):
            img = imageio.imread(name, pilmode='RGB')
            dset = hf.create_dataset('default/default'+str(ii+1), \
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
    
    ATL14_browse_plots(args) 



