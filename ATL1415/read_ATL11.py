#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 09:14:51 2025

@author: ben
"""
import pointCollection as pc
import numpy as np
import os


def select_best_xovers(D):
    _, i_pts = pc.unique_by_rows(np.c_[D.rgt, D.cycle_number, D.beam_pair], return_dict=True)
    ii = np.zeros(len(i_pts), dtype=int)
    for count, (pt, i_pt) in enumerate(i_pts.items()):
        if len(i_pt)==0:
            ii[count]=i_pt
        else:
            ii[count]=i_pt[np.argsort(D.h_corr_sigma[i_pt])[0]]
    D.index(ii)


def read_ATL11(xy0, Wxy, index_file, SRS_proj4, xover_tile_root=None, sigma_geo=6.5, sigma_radial=0.03, xover_cycles=[1,2]):


    bounds = [xy0[0]+np.array([-Wxy/2, Wxy/2]), xy0[1]+np.array([-Wxy/2, Wxy/2])]

    D_at, ATL11_file_list = read_ATL11_at(bounds, index_file, SRS_proj4,
                  sigma_geo=sigma_geo,
                  sigma_radial=sigma_radial)

    if xover_tile_root is None:
        return D_at, ATL11_file_list
    
    # Otherwise, read the crossover tiles
    D_xo, xover_file_list = read_ATL11_xovers(bounds, D_at, SRS_proj4,
                                              xover_tile_dir = xover_tile_root, 
                                              xover_cycles = xover_cycles)
    return pc.data().from_list(D_at, D_xo), ATL11_file_list + xover_file_list


def read_ATL11_at(bounds, index_file, SRS_proj4,
               sigma_geo=6.5, sigma_radial=0.03):
    '''
    read ATL11 data from an index file

    inputs:
        xy0 : 2-element iterable specifying the domain center
        Wxy : Width of the domain
        index_file : file made by pointCollection.geoindex pointing at ATL11 data
        SRS_proj4: projection information for the data

    output:
        D: data structure
        file_list: list of ATL11 files read
    '''

    field_dict_11={None:['latitude','longitude','delta_time',\
                        'h_corr','h_corr_sigma','h_corr_sigma_systematic', 'ref_pt'],\
                        '__calc_internal__' : ['rgt'],
                        'cycle_stats' : {'tide_ocean','dac'},
                        'ref_surf':['e_slope','n_slope', 'x_atc', 'fit_quality', 'dem_h', 'geoid_h']}

    try:
        # catch empty data
        D11_list=pc.geoIndex().from_file(index_file).query_xy_box(
            *bounds, fields=field_dict_11)
    except ValueError:
        return None, []
    if D11_list is None:
        return None, []
    D_list=[]

    D11_files=[]
    for D11 in D11_list:
        D11.get_xy(proj4_string=SRS_proj4)
        # select the subset of the data within the domain
        keep = (D11.x[:,0] >= bounds[0][0]) & (D11.x[:,0] <= bounds[0][1]) &\
             (D11.y[:,0] >= bounds[1][0]) & (D11.y[:,0] <= bounds[1][1])
        D11.index(keep)
        if D11.size==0:
            continue
        D11_files += [D11.filename]
        sigma_corr=np.sqrt((sigma_geo*np.abs(np.median(D11.n_slope)))**2+\
                           (sigma_geo*np.abs(np.median(D11.e_slope)))**2+sigma_radial**2)

        n_cycles=np.sum(np.isfinite(D11.h_corr), axis=1)
        n_cycles=np.reshape(n_cycles, (D11.shape[0],1))
        n_cycles=np.tile(n_cycles, [1, D11.shape[1]])

        D_list += [pc.data().from_dict({'z':D11.h_corr,
           'sigma_corr':sigma_corr+np.zeros_like(D11.h_corr),
           'sigma':D11.h_corr_sigma,
           'x':D11.x,
           'y':D11.y,
           'x_atc': D11.x_atc,
           'latitude':D11.latitude,
           'longitude':D11.longitude,
           'rgt':D11.rgt,
           'pair':np.zeros_like(D11.x)+D11.pair,
           'ref_pt':D11.ref_pt,
           'cycle':D11.cycle_number,
           'n_cycles': n_cycles,
           'fit_quality': D11.fit_quality,
           'dem_h': D11.dem_h,
           'tide_ocean': D11.tide_ocean,
           'dac': D11.dac,
           'geoid_h':D11.geoid_h,
           'delta_time': D11.delta_time,
           'time':D11.delta_time/24/3600/365.25+2018,
           'n_slope':D11.n_slope,
           'e_slope':D11.e_slope,
           'along_track':np.ones_like(D11.x, dtype=bool)})]

    return pc.data().from_list(D_list), D11_files

def read_ATL11_xovers(bounds, D_at, SRS_proj4, xover_tile_dir=None, xover_cycles=[1,2]):

    D_x=[]
    D_r=[]
    for x_cycle in xover_cycles:
        schema_file = os.path.join(xover_tile_dir,
                                   f'cycle_{x_cycle:02d}',
                                   '200km_tiling.json')
        xover_files = pc.tilingSchema().from_file(schema_file).filenames_for_box(bounds)
        for xover_file in xover_files:
            D_xi = pc.data().from_h5(xover_file, group='crossing_track').get_xy(proj4_string=SRS_proj4)
            keep = (D_xi.x >= bounds[0][0]) & (D_xi.x <= bounds[0][1]) &\
                 (D_xi.y >= bounds[1][0]) & (D_xi.y <= bounds[1][1])
            if not np.any(keep):
                continue
            D_xi.index(keep)
            D_ri = pc.data().from_h5(xover_file, group='datum_track', fields=['rgt','ref_pt','pair_track','cycle_number'])
            D_ri.index(keep)
            D_x += [D_xi]
            D_r += [D_ri]

    D_x = pc.data().from_list(D_x)
    D_r = pc.data().from_list(D_r)
    # Make sortable indexes for the along-track data
    rtp_at = np.core.records.fromarrays([D_at.rgt, D_at.pair, D_at.ref_pt])
    rtp_xo = np.core.records.fromarrays([D_r.rgt, D_r.pair_track, D_r.ref_pt])
    i_at = rtp_at.argsort()
    i_xo = i_at[np.searchsorted(rtp_xo, rtp_at[i_at])]
    # map parameters from the along-track data to the across-track data
    for field in ['geoid_h', 'dem_h','n_slope', 'e_slope','sigma_corr']:
        D_x.assign({field : getattr(D_at, field)[i_xo]})

    # choose the smallest_sigma xover for each rgt and pair
    select_best_xovers(D_x)

    blank = np.zeros_like(D_x.h_corr) + np.NaN
    D_xo = [pc.data().from_dict({
        'z':D_x.h_corr,
        'sigma':D_x.h_corr_sigma,
        'sigma_corr': D_x.h_corr,
        'x':D_x.x,
        'y':D_x.y,
        'latitude':D_x.latitude,
        'longitude':D_x.longitude,
        'dem_h':D_x.dem_h,
        'geoid_h':D_x.geoid_h,
        'rgt':D_x.rgt,
        'pair':D_x.beam_pair,
        'ref_pt':blank,
        'cycle':D_x.cycle_number,
        'n_cycles':blank,
        'fit_quality':D_r.fit_quality,
        'tide_ocean':D_x.tide_ocean,
        'dac':D_x.dac,
        'delta_time':D_x.delta_time,
        'n_slope':D_x.n_slope,
        'e_slope':D_x.e_slope,
        'time': D_x.delta_time/24/3600/365.25+2018,
        'along_track':np.zeros_like(D_x.x, dtype=bool)})]
    return D_xo, xover_files
