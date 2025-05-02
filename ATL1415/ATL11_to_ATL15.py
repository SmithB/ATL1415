#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 17:15:08 2019

@author: ben
"""
import os
import sys
import re

threads_re=re.compile("THREADS=(\S+)")
n_threads="1"
for arg in sys.argv:
    try:
        n_threads=str(threads_re.search(arg).group(1))
    except Exception:
        pass

# if THREADS was not specified as an input argument, check if it's set by slurm
if n_threads=="1" and "SLURM_NTASKS" in os.environ:
    n_threads=os.environ['SLURM_NTASKS']

os.environ["MKL_NUM_THREADS"]=n_threads
os.environ["OPENBLAS_NUM_THREADS"]=n_threads

import numpy as np
from LSsurf.smooth_fit import smooth_fit
from SMBcorr import assign_firn_variable
import pointCollection as pc

import re
import sys
import h5py
import traceback
from ATL1415.reread_data_from_fits import reread_data_from_fits
from ATL1415.make_mask_from_vector import make_mask_from_vector
from ATL1415.SMB_corr_from_grid import SMB_corr_from_grid
import pyTMD
import scipy.optimize

def get_SRS_info(hemisphere):
    if hemisphere==1:
        return '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs', 3413
    else:
        return '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs', 3031

def manual_edits(D):
    '''
    Remove known problematic data from a data structure

    inputs:
        D: data structure
    outputs:
        None (modifies input structure in place)
    '''
    bad=(D.rgt == 741)  & (D.cycle == 7)

    #N.B.  This fixes a glitch on the east coast of Greenland:
    if np.max(np.abs(np.array([np.mean(D.x), np.mean(D.y)])-np.array([480000, -1360000])))<1.e5:
        print("EDITING TRACK 1092 CYCLE 1")
        bad |= ((D.rgt==1092) & (D.cycle==1))

    D.index(~bad)
    return

def select_best_xovers(D):
    _, i_pts = pc.unique_by_rows(np.c_[D.rgt, D.cycle_number, 1+np.floor((D.spot_crossing-1)/2)], return_dict=True)
    ii = np.zeros(len(i_pts), dtype=int)
    for count, (pt, i_pt) in enumerate(i_pts.items()):
        if len(i_pt)==0:
            ii[count]=i_pt
        else:
            ii[count]=i_pt[np.argsort(D.h_corr_sigma[i_pt])[0]]
    D.index(ii)

def read_ATL11(xy0, Wxy, index_file, SRS_proj4, \
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
    xover_fields = pc.ATL11.crossover_data().__default_XO_field_dict__()
    xover_fields = xover_fields[list(xover_fields.keys())[0]] + ['spot_crossing']

    try:
        # catch empty data
        D11_list=pc.geoIndex().from_file(index_file).query_xy_box(
            xy0[0]+np.array([-Wxy/2, Wxy/2]), \
            xy0[1]+np.array([-Wxy/2, Wxy/2]), fields=field_dict_11)
    except ValueError:
        return None, []
    if D11_list is None:
        return None, []
    D_list=[]
    XO_list=[]
    xover_count=[0, 0]
    file_list= [ Di.filename for Di in D11_list ]
    for D11 in D11_list:
        D11.get_xy(proj4_string=SRS_proj4)
        # select the subset of the data within the domain
        D11.index((np.abs(D11.x[:,0]-xy0[0]) <= Wxy/2) &\
                 (np.abs(D11.y[:,0]-xy0[1]) <= Wxy/2))
        if D11.size==0:
            continue
        sigma_corr=np.sqrt((sigma_geo*np.abs(np.median(D11.n_slope)))**2+\
                           (sigma_geo*np.abs(np.median(D11.e_slope)))**2+sigma_radial**2)
        # fix for early ATL11 versions that had some zero error estimates.
        bad=np.any(D11.h_corr_sigma==0, axis=1)
        D11.h_corr_sigma[bad,:]=np.NaN

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

        if len(D11.ref_pt) == 0:
            continue
        # N.B.  D11 is getting indexed in this step, and it's leading to the warning in
        # line 76.  Can fix by making crossover_data.from_h5 copy D11 on input
        D_x = pc.ATL11.crossover_data().from_h5(D11.filename, pair=D11.pair, D_at=D11,\
                                                crossover_fields=xover_fields)
        if D_x is None:
            continue
        # constant fields in D11 apply to all cycles, but are mapped
        # only to some specific cycle of the reference cycle
        constant_fields = ['fit_quality', 'dem_h', 'geoid_h']
        for field in constant_fields:
            temp=getattr(D_x, field)
            val=np.nanmax(temp[:,:,0], axis=1)
            for cycle in range(D_x.shape[1]):
                temp[:, cycle, 1] = val
        #temp={'fit_quality':np.nanmax(D_x.fit_quality[:,:,0], axis=1),
        #      'dem_h':np.nanmax(D_x.dem_h[:,:,0], axis=1)}
        #for cycle in range(D_x.shape[1]):
        #    D_x.fit_quality[:,cycle,1]=temp['fit_quality']
        #    D_x.dem_h[:,cycle,1]=temp['dem_h']

        D_x.get_xy(proj4_string=SRS_proj4)

        good=np.isfinite(D_x.h_corr)[:,0:2,1].ravel()
        for field in D_x.fields:
            # Pull out only cycles 1 and 2
            temp=getattr(D_x, field)[:,0:2,1]
            setattr(D_x, field, temp.ravel()[good])
        # select the subset of the data within the domain
        D_x.index((np.abs(D_x.x-xy0[0]) <= Wxy/2) &\
                 (np.abs(D_x.y-xy0[1]) <= Wxy/2))
        if D_x.size==0:
            continue

        # choose the smallest_sigma xover for each rgt and pair
        xover_count[0]+=D_x.size
        select_best_xovers(D_x)
        xover_count[1]+=D_x.size

        #N.B.  Check whether n_slope and e_slope are set correctly.
        zero = np.zeros_like(D_x.h_corr)
        blank = zero+np.NaN
        XO_list += [pc.data().from_dict({'z':D_x.h_corr,
            'sigma':D_x.h_corr_sigma,
            'sigma_corr': sigma_corr+np.zeros_like(D_x.h_corr),
            'x':D_x.x,
            'y':D_x.y,
            'latitude':D_x.latitude,
            'longitude':D_x.longitude,
            'dem_h':D_x.dem_h,
            'geoid_h':D_x.geoid_h,
            'rgt':D_x.rgt,
            'pair':np.zeros_like(D_x.x)+D_x.pair,
            'ref_pt':blank,
            'cycle':D_x.cycle_number,
            'n_cycles':blank,
            'fit_quality':D_x.fit_quality,
            'tide_ocean':D_x.tide_ocean,
            'dac':D_x.dac,
            'delta_time':D_x.delta_time,
            'n_slope':D_x.n_slope,
            'e_slope':D_x.e_slope,
            'time':D_x.delta_time/24/3600/365.25+2018,
            'along_track':np.zeros_like(D_x.x, dtype=bool)})]
    try:
        D=pc.data().from_list(D_list+XO_list).ravel_fields()
    except ValueError:
        # catch empty data
        return None, file_list
    if hasattr(D,'z'):
        D.index(np.isfinite(D.z))
    else:
        return None, file_list
    # accept crossover points, along-track points that include at least 5 cycles, and along-track-points that have good fit_quality stats
    D.index( (D.along_track & (D.n_cycles>5)) |
            (( D.fit_quality ==0 ) | ( D.fit_quality == 2 )))
    print(f'xover_count={xover_count}')
    return D, file_list

def apply_tides(D, xy0, W,
                tide_mask_file=None,
                tide_mask_data=None,
                tide_directory=None,
                tide_model=None,
                tide_adjustment=False,
                tide_adjustment_file=None,
                tide_adjustment_format='h5',
                extrapolate=True,
                cutoff=200,
                EPSG=None,
                verbose=False):

    '''
    read in the tide mask, calculate ocean tide elevations, and
    apply dynamic atmospheric correction (dac) and tide to ice-shelf elements

    inputs:
        D: data structure
        xy0: 2-element iterable specifying the domain center
        W: Width of the domain


    keyword arguments:
        tide_mask_file: geotiff file for masking to ice shelf elements
        tide_mask_data: pc.grid.data() object containing the tide mask (alternative to tide_mask_file)
        tide_directory: path to tide models
        tide_model: the name of the tide model to use
        tide_adjustment: adjust amplitudes of tide model amplitudes to account for ice flexure
        tide_adjustment_file: File for adjusting tide and dac values for ice shelf flexure
        tide_adjustment_format: file format of the scaling factor grid
        extrapolate: extrapolate outside tide model bounds with nearest-neighbors
        cutoff: extrapolation cutoff in kilometers

    output:
        D: data structure corrected for ocean tides and dac
    '''

    # the tide mask should be 1 for non-grounded points (ice shelves), zero otherwise
    if tide_mask_file is not None and tide_mask_data is None:
        try:
            tide_mask = pc.grid.data().from_geotif(tide_mask_file,
                        bounds=[np.array([-0.6, 0.6])*W+xy0[0], np.array([-0.6, 0.6])*W+xy0[1]])
        except IndexError:
            return None
        if tide_mask.shape is None:
            return
    # find ice shelf points
    is_els=tide_mask.interp(D.x, D.y) > 0.5
    # need to assign the 'floating' field in case the SMB routines need it
    D.assign(floating=is_els)
    if verbose:
        print(f"\t\t{np.mean(is_els)*100}% shelf data")
        print(f"\t\ttide model: {tide_model}")
        print(f"\t\ttide directory: {tide_directory}")
    # extrapolate tide estimate beyond model bounds
    # extrapolation cutoff is in kilometers
    if np.any(is_els.ravel()):
        D.tide_ocean = pyTMD.compute_tide_corrections(
                D.x, D.y, D.delta_time,
                DIRECTORY=tide_directory, MODEL=tide_model,
                EPOCH=(2018,1,1,0,0,0), TYPE='drift', TIME='utc',
                EPSG=EPSG, EXTRAPOLATE=extrapolate, CUTOFF=cutoff)
    # use a flexure mask to adjust estimated tidal values
    if np.any(is_els.ravel()) and tide_adjustment:
        print(f"\t\t{tide_adjustment_file}") if verbose else None
        D.assign({'tide_adj_scale': np.ones_like(D.x)})
        # read adjustment grid and interpolate to values
        tide_adj_scale = pc.grid.data().from_file(tide_adjustment_file,
            file_format=tide_adjustment_format, field_mapping=dict(z='tide_adj_scale'),
            bounds=[np.array([-0.6, 0.6])*W+xy0[0], np.array([-0.6, 0.6])*W+xy0[1]])
        # interpolate tide adjustment to coordinate values
        D.tide_adj_scale[:] = tide_adj_scale.interp(D.x, D.y)
        # mask out scaling factors where grounded
        D.tide_adj_scale[is_els==0]=0.0
        # multiply tides and dynamic atmospheric correction by adjustments
        ii, = np.nonzero(D.tide_adj_scale != 0)
        D.tide_ocean[ii] *= D.tide_adj_scale[ii]
        D.dac[ii] *= D.tide_adj_scale[ii]
        print(f'mean tide adjustment scale={np.nanmean(D.tide_adj_scale)}')

    # replace invalid tide and dac values
    D.tide_ocean[is_els==0] = 0
    D.dac[is_els==0] = 0
    D.tide_ocean[~np.isfinite(D.tide_ocean)] = 0
    D.dac[~np.isfinite(D.dac)] = 0
    # apply tide and dac to heights
    D.z -= (D.tide_ocean + D.dac)
    return D

def read_bedmachine_greenland(mask_file, xy0, Wxy):
    from xarray import open_dataset
    x0=np.arange(xy0[0]-0.6*Wxy, xy0[0]+0.6*Wxy,100)
    y0=np.arange(xy0[1]-0.6*Wxy, xy0[1]+0.6*Wxy,100)
    with  open_dataset(mask_file,'r') as md:
        mask_data=pc.grid.data().from_dict({'x':x0,'y':y0,
            'z':np.array(md['mask'].interp(x0, y0))})
    mask_data.z[~np.isfinite(mask_data.z)]=0
    # ice-shelves in the mask are represented by 3
    tide_mask_data=pc.grid.data().from_dict({'x':x0,'y':y0,
            'z':np.abs(mask_data.z-3)<0.5})
    # land ice is either 2 (grounded ice) or 3 (shelf ice)
    mask_data.z = ((mask_data.z-3)<0.5) | ((mask_data.z-2)<0.5)
    return mask_data, tide_mask_data

def decimate_data(D, N_target, W_domain,  W_sub, x0, y0):
    # reduce the data density to a target value in small bins
    ij_bin=np.round((D.x - (x0 - W_domain/2 + W_sub/2))/W_sub)+ \
        1j*np.round((D.y - (y0 - W_domain/2 + W_sub/2))/W_sub)
    rho_target = N_target / W_domain**2
    ind_buffer=[]
    print(f'decimate_data: N_target:{N_target}, N={D.size}, R_target={D.size/N_target}')
    for bin0 in np.unique(ij_bin):
        ii = np.flatnonzero((ij_bin==bin0) & D.along_track)
        this_rho = len(ii) / (W_sub**2)
        #print(f'bin0={bin0}, this_rho={this_rho}, ratio={this_rho/rho_target}')
        if this_rho < rho_target:
            # if there aren't too many points in this bin, continue
            #print(f'bin:{bin0}, N: {len(ii)}, N_target:{len(ii)*rho_target/this_rho}, N_out:{len(ii)}, R_target={this_rho/rho_target}, R=1.0')
            ind_buffer += [ii]
            continue
        # make a global reference point number (equal to the number of ref pts in an orbit * rgt + ref_pt)
        global_ref_pt=D.ref_pt[ii]+40130000/20*D.rgt[ii]
        u_ref_pts = np.unique(global_ref_pt)
        # select a subset of the global reference points (this skips the right number of points)
        # note that the np.arange() call can sometimes return a value of len(u_ref_pts) (?), so
        # its output needs to be checked (!)
        sel_ind=np.arange(0, len(u_ref_pts), this_rho/rho_target).astype(int)
        sel_ref_pts = u_ref_pts[sel_ind[sel_ind < len(u_ref_pts)]]
        # keep the points matching the selected ref pt numbers
        isub = np.in1d(global_ref_pt, sel_ref_pts)
        #print(f'bin:{bin0}, N: {len(ii)}, N_target:{len(ii)*rho_target/this_rho}, N_out:{np.sum(isub)}, R_target={this_rho/rho_target}, R={len(ii)/np.sum(isub)}')
        ind_buffer.append(ii[isub])
    print(f"Decimate_data: N_AT={len(np.concatenate(ind_buffer))}, N_XO={np.sum(D.along_track==0)}")
    ind_buffer.append(np.flatnonzero(D.along_track==0))
    D.index(np.concatenate(ind_buffer))

def set_three_sigma_edit_with_DEM(data, xy0, Wxy, DEM_file, DEM_tol, W_med=None):
    '''
    Check data against DEM

    inputs:
        data (pc.data): input data
        xy0  (list of floats) : tile center
        Wxy (float) : tile width
        DEM_file (string) : filename for DEM (if None, data.dem_h is used)
        DEM_tol (float) : tolerance for test
        W_med (float) : distance over which the DEM is corrected to the median of the data

    outputs:  None (modifies data.three_sigma_edit)
    '''

    if DEM_tol is None:
        return

    if W_med is None:
        W_med=Wxy/8
    if 'three_sigma_edit' not in data.fields:
        data.assign({'three_sigma_edit':np.ones_like(data.z, dtype=bool)})

    if DEM_file is None:
        r_DEM = data.z-data.dem_h
    else:
         r_DEM = data.z - pc.grid.data()\
                .from_geotif(DEM_file, bounds=[xy+0.6*np.array([-Wxy, Wxy]) for xy in xy0])\
                .interp(data.x, data.y)

    good = np.ones_like(data.z, dtype=bool)
    for x0 in np.arange(xy0[0]-Wxy/2, xy0[0]+Wxy/2, W_med):
        for y0 in np.arange(xy0[1]-Wxy/2, xy0[1]+Wxy/2, W_med):
            ii = (data.x>=x0) & (data.x <= x0+W_med) &\
                (data.y>=y0) & (data.y <= y0+W_med)
            good[ii] &= (np.abs(r_DEM[ii] - np.nanmedian(r_DEM[ii])) < DEM_tol)
    data.three_sigma_edit &= good

def ATL11_to_ATL15(xy0, Wxy=4e4, ATL11_index=None, E_RMS={}, \
            t_span=[2019.25, 2020.5], \
            spacing={'z0':2.5e2, 'dz':5.e2, 'dt':0.25},  \
            sigma_geo=6.5,\
            sigma_radial=0.03,\
            dzdt_lags=[1, 4],\
            hemisphere=1,\
            reference_epoch=None,\
            region=None,\
            reread_dirs=None, \
            sigma_extra_bin_spacing=None,\
            sigma_extra_max=None,\
            prior_edge_args=None,\
            data_file=None, \
            restart_edit=False, \
            max_iterations=5, \
            N_subset=8,  \
            W_edit=None, \
            out_name=None, \
            compute_E=False, \
            DEM_file=None,\
            DEM_tol=None,\
            geoid_tol=None, \
            sigma_tol=None,\
            mask_file=None,\
            rock_mask_file=None,\
            rock_mask_reject_value=None,\
            geoid_file=None,\
            E_d3zdx2dt_scale_file=None,\
            tide_mask_file=None,\
            tide_directory=None,\
            tide_adjustment=False,\
            tide_adjustment_file=None,\
            tide_adjustment_format=None,\
            tide_model=None,\
            firn_correction=None, \
            firn_directory=None,\
            firn_version=None,\
            firn_grid_file=None,\
            avg_scales=None,\
            edge_pad=None,\
            error_res_scale=None,\
            calc_error_file=None,\
            bias_params=['rgt','cycle'],\
            verbose=False,\
            write_data_only=False):
    '''
    Function to generate DEMs and height-change maps based on ATL11 surface height data.

    Inputs:
        xy0: 2 elements specifying the center of the domain
        Wxy: (float) Width of the domain
        ATL11_index: (string) Index file (from pointCollection.geoIndex) for ATL11 data
        E_RMS: (dict) dictionary specifying constraints on derivatives of the ice-sheet surface.
        t_span: (2-element list of floats) starting and ending times for the output grids (in years CE)
        spacing: (dict) dictionary specifying the grid spacing for z0, dz, and dt
        dzdt_lags: (list) lags over which elevation change rates and errors will be calculated
        hemisphere: (int) the hemisphere in which the grids will be generated. 1 for northern hemisphere, -1 for southern
        reference_epoch: (int) The time slice (counting from zero, in steps of spacing['dt']) corresponding to the DEM
        reread_dirs: (string) directory containing output files from which data will be read (if None, data will be read from the index file)
        data_file: (string) output file from which to reread data (alternative to reread_dirs)
        max_iterations: (int) maximum number of iterations in three-sigma edit of the solution
        N_subset: (int) If specified, the domain is subdivided into this number of divisions in x and y, and a three-sigma edit is calculated for each
        sigma_extra_bin_spacing: (float) grid spacing that will be used to calculate extra geophysical errors
        W_edit: (float) Only data within this distance of the grid center can be editied in the iterative step
        out_name: (string) Name of the output file
        compute_E: (bool) If true, errors will be calculated
        DEM_file: (string) DEM against which the data will be checked using the DEM_tol parameter
        DEM_tol: (float) Points different from the DEM by more than this value will be edited
        sigma_tol: (float) Points with a sigma parameter greater than this value will be edited
        mask_file: (string) File specifying areas for which data should be used and strong constraints should be applied
        tide_mask_file: (string)  File specifying the areas for which the tide model should be calculated
        tide_directory: (string)  Directory containing the tide model data
        tide_adjustment: (bool)  Adjust tides for ice shelf flexure
        tide_adjustment_file: (string)  File for adjusting tide and dac values for ice shelf flexure
        tide_adjustment_format: (string) file format of the scaling factor grid
        tide_model: (string)  Name of the tide model to use for a given domain
        E_d3zdx2dt_scale_file: (string) Filename containing scaling to be applied to the d3zdx2dt parameter
        avg_scales: (list of floats) scales over which the output grids will be averaged and errors will be calculated
        error_res_scale: (float) If errors are calculated, the grid resolution will be coarsened by this factor
        calc_error_file: (string) Output file for which errors will be calculated.
        verbose: (bool) Print progress of processing run
    outputs:
        S: dict containing fit output
    '''
    SRS_proj4, EPSG=get_SRS_info(hemisphere)

    E_RMS0={'d2z0_dx2':200000./3000/3000, 'd3z_dx2dt':3000./3000/3000, 'd2z_dt2':5000}
    E_RMS0.update(E_RMS)

    W={'x':Wxy, 'y':Wxy,'t':np.diff(t_span)}
    ctr={'x':xy0[0], 'y':xy0[1], 't':np.mean(t_span)}
    bds={ dim: c_i+np.array([-0.5, 0.5])*W[dim]  for dim, c_i in ctr.items()}

    # work out which mask to use based on the region
    tide_mask_data=None
    mask_data=None

    read_mask_file=None
    if data_file is not None:
        read_mask_file=data_file
    elif calc_error_file is not None:
        read_mask_file=calc_error_file

    mask_update_function=None
    if geoid_file is not None:
        ancillary_data={'geoid_file':geoid_file}
        from ATL1415.update_masks_with_geoid import update_masks_with_geoid
        mask_update_function=update_masks_with_geoid

    if read_mask_file is not None:
        print("READING MASK DATA")
        mask_data={'z0':pc.grid.data().from_h5(read_mask_file, group='z0', fields=['mask']),
                   'dz':pc.grid.data().from_h5(read_mask_file, group='dz', fields=['mask'])}
        for key, mask in mask_data.items():
            mask.assign({'z':mask.mask})
    elif region is not None:
        if region in ['AA', 'GL']:
            pad=np.array([-1.e4, 1.e4])
            mask_data=pc.grid.data().from_file(mask_file,
                                             bounds=[bds['x']+pad, bds['y']+pad],
                                             t_range=bds['t']+np.array([-1, 1]))
            if rock_mask_file is not None:
                rock_mask=pc.grid.data().from_file(rock_mask_file, bounds=mask_data.bounds())
                rock_mask.reject=rock_mask.mask==rock_mask_reject_value
                rock=rock_mask.interp(mask_data.x, mask_data.y, gridded=True, field='reject') > 0.5
                for band in range(mask_data.z.shape[2]):
                    mask_data.z[:,:,band] *= (rock==0)
            while mask_data.t[-1] < ctr['t']+W['t']/2:
                # append a copy of the last field in the mask data to the end of the mask data
                mask_data.z = np.concatenate([mask_data.z,mask_data.z[:,:,-1:]], axis=2)
                mask_data.t = np.concatenate([mask_data.t,mask_data.t[-1:]+1], axis=0)
            mask_data.__update_size_and_shape__()
            mask_data.z[~np.isfinite(mask_data.z)]=0.
            mask_file=None
        elif mask_file.endswith('.shp') or mask_file.endswith('.db'):
            mask_data=make_mask_from_vector(mask_file, W, ctr, spacing['z0'], srs_proj4=SRS_proj4)

        # check if mask_data is 2D or 3D
        if len(mask_data.z.shape)==2:
            #repeat the mask data to make a 3d field
            mask_data.t=np.arange(t_span[0]-1, t_span[1]+1)
            mask_data.z=np.tile(mask_data.z[:,:,None], [1, 1, len(mask_data.t)])
            mask_data.__update_size_and_shape__()


    # initialize file_list to empty in case we're rereading the data
    file_list=[]

    constraint_scaling_maps=None
    if E_d3zdx2dt_scale_file is not None:
        constraint_scaling_maps={'d3z_dx2dt': pc.grid.data().from_file(
            E_d3zdx2dt_scale_file,
            bounds = [ii+Wxy*np.array([-0.6, 0.6]) for ii in xy0])}

    # figure out where to get the data
    if data_file is not None:
        data=pc.data().from_h5(data_file, group='data')
    elif calc_error_file is not None:
        data=pc.data().from_h5(calc_error_file, group='data')
        max_iterations=0
        compute_E=True
    elif reread_dirs is not None:
        data, tile_reread_list = reread_data_from_fits(xy0, Wxy, reread_dirs, template='E%d_N%d.h5')
    else:
        data, file_list = read_ATL11(xy0, Wxy, ATL11_index, SRS_proj4,
                                     sigma_geo=sigma_geo, sigma_radial=sigma_radial)
        if sigma_tol is not None and data is not None:
            data.index(data.sigma < sigma_tol)
        if data is not None:
            N0=data.size
            decimate_data(data, 1.2e6, Wxy, 5000, xy0[0], xy0[1] )
            print(f'decimated {N0} to {data.size}')

    if data is None:
        print("No data present for region, returning.")
        return None

    # if any manual edits are needed, make them here:
    manual_edits(data)

    if data is None or data.size < 10:
        print("Fewer than 10 data points, returning")
        return None

    if data.time.max()-data.time.min() < 80./365.:
        print("time span too short, returning.")
        return None

    if edge_pad is not None:
        ctr_dist = np.max(np.abs(data.x-xy0[0]), np.abs(data.y-xy0[1]))
        in_ctr = ctr_dist < Wxy/2 - edge_pad
        if np.sum(in_ctr) < 50:
            print("After editing by edge pad, <50 data present, returning")
            return None

    if W_edit is not None:
        # this is used for data that we are rereading from a set of other files.
        # data that are far from the center of this file cannot be eliminated by
        # the three-sigma edit
        data.assign({'editable':  (np.abs(data.x-xy0[0])<=W_edit/2) & (np.abs(data.y-xy0[1])<=W_edit/2)})

    if restart_edit:
        if 'three_sigma_edit' in data.fields:
            data.three_sigma_edit[:]=True

    if firn_correction is not None:
        if 'h_firn' in data.fields:
            # if there is a firn varaible already, undo the correction
            data.z += data.h_firn
        if firn_grid_file is not None:
            if firn_correction=='MERRA2_hybrid':
                # defaults work here:
                SMB_corr_from_grid(data,
                    model_file=os.path.join(firn_directory,firn_grid_file))
        else:
            assign_firn_variable(data, firn_correction, firn_directory, hemisphere,
                             model_version=firn_version, subset_valid=False)
        # make the correction
        data.z -= data.h_firn

    # to reject ATL06/11 blunders
    set_three_sigma_edit_with_DEM(data, xy0, Wxy, DEM_file, DEM_tol)

    # apply the tides if a directory has been provided
    # NEW 2/19/2021: apply the tides only if we have not read the data from first-round fits.
    if (tide_mask_file is not None or tide_mask_data is not None) and reread_dirs is None and calc_error_file is None and data_file is None:
        apply_tides(data, xy0, Wxy,
                    tide_mask_file=tide_mask_file,
                    tide_mask_data=tide_mask_data,
                    tide_directory=tide_directory,
                    tide_model=tide_model,
                    tide_adjustment=tide_adjustment,
                    tide_adjustment_file=tide_adjustment_file,
                    tide_adjustment_format=tide_adjustment_format,
                    EPSG=EPSG, verbose=verbose)
    elif 'floating' not in data.fields:
        # fix for firn runs where the floating variable did not get assigned
        if tide_mask_file is not None and tide_mask_data is None:
            try:
                tide_mask = pc.grid.data().from_geotif(tide_mask_file,
                                                       bounds=[np.array([-0.6, 0.6])*Wxy+xy0[0], np.array([-0.6, 0.6])*Wxy+xy0[1]])
                if tide_mask.shape is not None:
                    data.assign(floating=tide_mask.interp(data.x, data.y) > 0.5)
            except IndexError:
                data.assign(floating=np.zeros_like(data.x, dtype=bool))

    if geoid_tol is not None:
        data.index((data.z - data.geoid_h) > geoid_tol)

    if write_data_only:
        return {'data':data}

    # call smooth_xytb_fitting
    S=smooth_fit(data=data,
                      ctr=ctr, W=W,
                      spacing=spacing, E_RMS=E_RMS0,
                      reference_epoch=reference_epoch,
                      compute_E=compute_E,
                      constraint_scaling_maps=constraint_scaling_maps,
                      sigma_extra_bin_spacing=sigma_extra_bin_spacing,
                      sigma_extra_max=sigma_extra_max,
                      prior_edge_args=prior_edge_args,
                      bias_params=bias_params,
                      max_iterations=max_iterations,
                      srs_proj4=SRS_proj4, VERBOSE=True,
                      dzdt_lags=dzdt_lags,
                      mask_file=mask_file, mask_data=mask_data, mask_scale={0:10, 1:1},
                      converge_tol_frac_edit=0.001,
                      error_res_scale=error_res_scale,
                      ancillary_data=ancillary_data,
                      mask_update_function=mask_update_function,
                      avg_scales=avg_scales)
    S['file_list'] = file_list
    return S

def save_fit_to_file(S,  filename, dzdt_lags=None, reference_epoch=0):
    if os.path.isfile(filename):
        os.remove(filename)
    with h5py.File(filename,'w') as h5f:
        h5f.create_group('/data')
        for key in S['data'].fields:
            h5f.create_dataset('/data/'+key, data=getattr(S['data'], key))
        # metadata:
        h5f.create_group('/meta')
        h5f.create_group('/meta/timing')
        for key in S['timing']:
            h5f['/meta/timing/'].attrs[key]=S['timing'][key]
        # write out list of ATL11 files so that lineage can be populated in ATL14/15
        if 'file_list' in S:
            h5f['meta'].attrs['input_files'] = ','.join([os.path.basename(Si) for Si in S['file_list']]).encode('ascii')
        h5f['meta'].attrs['first_delta_time']=np.nanmin(S['data'].delta_time)
        h5f['meta'].attrs['last_delta_time']=np.nanmax(S['data'].delta_time)
        h5f.create_group('/RMS')
        for key in S['RMS']:
            h5f.create_dataset('/RMS/'+key, data=S['RMS'][key])
        h5f.create_group('E_RMS')
        for key in S['E_RMS']:
            h5f.create_dataset('E_RMS/'+key, data=S['E_RMS'][key])
        for key in S['m']['bias']:
            h5f.create_dataset('/bias/'+key, data=S['m']['bias'][key])

    # if we have a 3d mask on z0, use it to replace the 'mask' in S['m']['dz']
    if S['grids']['dz'].mask_3d is not None:
        S['m']['dz'].mask = S['grids']['dz'].mask_3d

    for key , ds in S['m'].items():
        if isinstance(ds, pc.grid.data):
            ds.to_h5(filename, group=key)
    with h5py.File(filename,'r+') as h5f:
        h5f.create_dataset('/z0/mask', data=S['grids']['z0'].mask.astype(int), \
                           chunks=True, compression="gzip", fillvalue=-1)
        h5f.create_dataset('/dz/mask', data=S['grids']['dz'].mask_3d.z.astype(int),\
                           chunks=True, compression="gzip", fillvalue=-1)
    return

def mask_components_by_time(dz):
    """
    identify the connected components in the data, mark unconstrained epochs as invalid

    inputs
    dz: pc.grid.data() instance containing fields dz, cell_area
    outputs:
    modified inputs
    """
    from scipy.ndimage import label

    components, n_components = label(dz.cell_area>0)
    first_epoch=np.zeros(n_components, dtype=int)+n_components
    last_epoch=np.zeros(n_components, dtype=int)

    for comp in range(1, n_components):
        these = components==comp
        for t_slice in range(dz.shape[2]):
            sampled=np.any(dz.count[:,:,t_slice][these]>1)
            if t_slice <= first_epoch[comp]:
                if sampled:
                    first_epoch[comp]=t_slice
            if t_slice >= last_epoch[comp]:
                if sampled:
                    last_epoch[comp]=t_slice

    last_epoch_map=np.zeros_like(dz.cell_area)+np.NaN
    first_epoch_map=np.zeros_like(dz.cell_area)+np.NaN

    for comp in range(1, n_components):
        last_epoch_map[components==comp]=last_epoch[comp]
        first_epoch_map[components==comp]=first_epoch[comp]

    for t_slice in range(dz.dz.shape[2]):
        dz.dz[:,:,t_slice][t_slice < first_epoch_map]=np.NaN
        dz.dz[:,:,t_slice][t_slice > last_epoch_map]=np.NaN


def interp_ds(ds, scale):

    delta_xy=[(ds.x[1]-ds.x[0])/scale, (ds.y[1]-ds.y[0])/scale]
    xi=np.arange(ds.x[0], ds.x[-1]+delta_xy[0], delta_xy[0])
    yi=np.arange(ds.y[0], ds.y[-1]+delta_xy[1], delta_xy[1])
    out=pc.grid.data().from_dict({'x':xi,'y':yi})
    if len(ds.shape)==2:
        out.time=ds.time
    for field in ds.fields:
        z0=getattr(ds, field)
        # for sigma fields, fill in the areas around the valid points with
        # a smoothed version of the error field
        if len(ds.shape)==2:
            if 'sigma' in field:
                z0=pc.grid.fill_edges( z0, w_smooth=1, ndv=0)
            out.assign({field:pc.grid.data().from_dict({'x':ds.x, 'y':ds.y, 'z':z0}).interp(xi, yi, gridded=True)})
        else:
            if 'sigma' in field:
                z0=pc.grid.fill_edges(z0, w_smooth=1, dim=2, ndv=0)
            zi=np.zeros([xi.size, yi.size, ds.time.size])
            for epoch in range(ds.time.size):
                temp=pc.grid.data().from_dict({'x':ds.x, 'y':ds.y, 'z':np.squeeze(z0[:,:,epoch])})
                zi[:,:,epoch] = temp.interp(xi, yi, gridded=True)
            out.assign({field:zi})
    return out
def save_errors_to_file( S, filename, dzdt_lags=None, reference_epoch=None, grid_datasets=None):

    for key, ds in S['E'].items():
        if isinstance(ds, pc.grid.data):
            print(key)
            ds.to_h5(filename, group=key.replace('sigma_',''))

    with h5py.File(filename,'r+') as h5f:
        for key in S['E']['sigma_bias']:
            if 'bias/sigma' in h5f and  key in h5f['/bias/sigma']:
                print(f'{key} already exists in sigma_bias')
                h5f['/bias/sigma/'+key][...]=S['E']['sigma_bias'][key]
            else:
                h5f.create_dataset('/bias/sigma/'+key, data=S['E']['sigma_bias'][key])

    return

def main(argv):
    # account for a bug in argparse that misinterprets negative agruents
    for i, arg in enumerate(argv):
        if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

    import argparse
    parser=argparse.ArgumentParser(\
        description="function to fit ICESat-2 data with a smooth elevation-change model", \
        fromfile_prefix_chars="@")
    parser.add_argument('--xy0', type=float, nargs=2, help="fit center location")
    parser.add_argument('--ATL11_index', type=lambda p: os.path.abspath(os.path.expanduser(p)), required=True, help="ATL11 index file")
    parser.add_argument('--Width','-W',  type=float, help="Width of grid")
    parser.add_argument('--time_span','-t', type=str, help="time span, first year,last year AD (comma separated, no spaces)")
    parser.add_argument('--grid_spacing','-g', type=str, help='grid spacing:DEM (meters),dh maps xy (meters),dh_maps time (years): comma-separated, no spaces', default='250.,4000.,1.')
    parser.add_argument('--Hemisphere','-H', type=int, default=1, help='hemisphere: -1=Antarctica, 1=Greenland')
    parser.add_argument('--base_directory','-b', type=lambda p: os.path.abspath(os.path.expanduser(p)), help='base directory')
    parser.add_argument('--out_name', '-o', type=lambda p: os.path.abspath(os.path.expanduser(p)), help="output file name")
    parser.add_argument('--prelim', action='store_true')
    parser.add_argument('--centers', action="store_true")
    parser.add_argument('--edges', action="store_true")
    parser.add_argument('--corners', action="store_true")
    parser.add_argument('--matched', action="store_true")
    parser.add_argument('--prior_edge_include', type=float, help='include prior constraints over this width at the edge of each tile')
    parser.add_argument('--prior_sigma_scale', type=float, default=1, help='scale prior error estimates by this value')
    parser.add_argument('--tile_spacing', type=float, default=60000.)
    parser.add_argument('--sigma_extra_bin_spacing', type=float, default=None, help='width over which data are collected to calculate the extra error estimates')
    parser.add_argument('--sigma_extra_max', type=float, default=None, help='maximum value for sigma_extra,')
    parser.add_argument('--W_edit', type=int)
    parser.add_argument('--E_d2zdt2', type=float, default=5000)
    parser.add_argument('--E_d2z0dx2', type=float, default=0.02)
    parser.add_argument('--E_d3zdx2dt', type=float, default=0.0003)
    parser.add_argument('--E_d2z0dx2_file', type=lambda p: os.path.abspath(os.path.expanduser(p)), help='file from which to read the expected d2z0dx2 values')
    parser.add_argument('--E_d3zdx2dt_scale_file', type=str, help='grid file containing scaling values for the E_d3zdx2dt parameter')
    parser.add_argument('--data_gap_scale', type=float,  default=2500)
    parser.add_argument('--sigma_geo', type=float,  default=6.5)
    parser.add_argument('--sigma_radial', type=float,  default=0.03)
    parser.add_argument('--dzdt_lags', type=str, default='1,4', help='lags for which to calculate dz/dt, comma-separated list, no spaces')
    parser.add_argument('--avg_scales', type=str, default='10000,40000', help='scales at which to report average errors, comma-separated list, no spaces')
    parser.add_argument('--N_subset', type=int, default=None, help="number of pieces into which to divide the domain for (cheap) editing iterations.")
    parser.add_argument('--max_iterations', type=int, default=6, help="maximum number of iterations used to edit the data.")
    parser.add_argument('--map_dir','-m', type=lambda p: os.path.abspath(os.path.expanduser(p)))
    parser.add_argument('--DEM_file', type=lambda p: os.path.abspath(os.path.expanduser(p)), help='DEM file to use with the DEM_tol parameter')
    parser.add_argument('--DEM_tol', type=float, default=50, help='points different from the DEM by more than this value will be edited in the first iteration')
    parser.add_argument('--geoid_tol', type=float, help='points closer than this to the geoid will be rejected')
    parser.add_argument('--sigma_tol', type=float, help='points with sigma greater than this value will be edited')
    parser.add_argument('--mask_file', type=lambda p: os.path.abspath(os.path.expanduser(p)))
    parser.add_argument('--rock_mask_file', type=lambda p: os.path.abspath(os.path.expanduser(p)), help='mask indicating exposed rock')
    parser.add_argument('--rock_mask_reject_value', type=float, default=1, help='value within the rock mask file that indicates rock')
    parser.add_argument('--geoid_file', type=lambda p: os.path.abspath(os.path.expanduser(p)), help="file containing geoid information")
    parser.add_argument('--tide_mask_file', type=lambda p: os.path.abspath(os.path.expanduser(p)))
    parser.add_argument('--tide_directory', type=lambda p: os.path.abspath(os.path.expanduser(p)))
    parser.add_argument('--tide_adjustment', action="store_true", help="Adjust tides for ice shelf flexure")
    parser.add_argument('--tide_adjustment_file', type=lambda p: os.path.abspath(os.path.expanduser(p)), help="File for adjusting tide and dac values for ice shelf flexure")
    parser.add_argument('--tide_adjustment_format', type=str, choices=('geotif','h5','nc'), default='h5', help="File format of the scaling factor grid")
    parser.add_argument('--tide_model', type=str)
    parser.add_argument('--firn_directory', type=str, help='directory containing firn model')
    parser.add_argument('--firn_model', type=str, help='firn model name')
    parser.add_argument('--firn_version', type=str, help='firn version')
    parser.add_argument('--firn_grid_file', type=str, help='gridded firn model file that can be interpolated directly.')
    parser.add_argument('--reference_epoch', type=int, default=0, help="Reference epoch number, for which dz=0")
    parser.add_argument('--data_file', type=lambda p: os.path.abspath(os.path.expanduser(p)), help='read data from this file alone')
    parser.add_argument('--restart_edit', action='store_true')
    parser.add_argument('--calc_error_file','-c', type=lambda p: os.path.abspath(os.path.expanduser(p)), help='file containing data for which errors will be calculated')
    parser.add_argument('--calc_error_for_xy', action='store_true', help='calculate the errors for the file specified by the x0, y0 arguments')
    parser.add_argument('--error_res_scale','-s', type=float, nargs=2, default=[4, 2], help='if the errors are being calculated (see calc_error_file), scale the grid resolution in x and y to be coarser')
    parser.add_argument('--bias_params', type=str, default="rgt,cycle", help='one bias parameter will be assigned for each unique combination of these ATL11 parameters (comma-separated list with no spaces)')
    parser.add_argument('--region', type=str, help='region for which calculation is being performed')
    parser.add_argument('--verbose','-v', action="store_true")
    parser.add_argument('--write_data_only', action='store_true', help='save data without processing')
    args, unknown=parser.parse_known_args()

    args.grid_spacing = [np.float64(temp) for temp in args.grid_spacing.split(',')]
    args.dzdt_lags = [np.int64(temp) for temp in args.dzdt_lags.split(',')]
    args.time_span = [np.float64(temp) for temp in args.time_span.split(',')]
    args.avg_scales = [np.int64(temp) for temp in args.avg_scales.split(',')]

    spacing={'z0':args.grid_spacing[0], 'dz':args.grid_spacing[1], 'dt':args.grid_spacing[2]}
    E_RMS={'d2z0_dx2':args.E_d2z0dx2, 'd3z_dx2dt':args.E_d3zdx2dt, 'd2z_dt2':args.E_d2zdt2}

    if args.data_gap_scale > 0:
        E_RMS[ 'd2z_dxdt'] = args.E_d3zdx2dt*args.data_gap_scale
    print("E_RMS="+str(E_RMS))

    reread_dirs=None
    dest_dir=args.base_directory
    if args.W_edit is None:
        W_edit=args.Width/2
    else:
        W_edit=args.W_edit

    if args.centers:
        dest_dir += '/centers'
        prior_dirs=None
        W_edit=None
    elif args.prelim:
        dest_dir += '/prelim'
        prior_dirs=None
        W_edit=None
    elif args.edges or args.corners:
        reread_dirs=[args.base_directory+'/centers']
        if args.edges:
            dest_dir += '/edges'
        prior_dirs=reread_dirs
    elif args.corners:
        reread_dirs += [args.base_directory+'/edges']
        dest_dir +='/corners'
        prior_dirs=reread_dirs
    elif args.matched:
        dest_dir += '/matched'
        args.max_iterations=1
        prior_dirs = [args.base_directory+'/'+ii for ii in ['prelim','centers','edges','corners']]

    prior_edge_args=None
    if args.prior_edge_include is not None:
        prior_edge_args={'prior_dir': prior_dirs,
                         'edge_include':args.prior_edge_include,
                         'sigma_scale':args.prior_sigma_scale,
                         'tile_spacing':args.tile_spacing}

    if args.xy0 is None and args.calc_error_file is not None or args.data_file is not None:
        # get xy0 from the filename
        if args.calc_error_file is not None:
            re_match=re.compile('E(.*)_N(.*).h5').search(os.path.basename(args.calc_error_file))
        elif args.data_file is not None:
            re_match=re.compile('E(.*)_N(.*).h5').search(os.path.basename(args.data_file))
        args.xy0=[float(re_match.group(ii))*1000 for ii in [1, 2]]

    if args.E_d2z0dx2_file is not None and args.calc_error_file is None:
        E_d2z0dx2=pc.grid.data().from_geotif(args.E_d2z0dx2_file)#, bounds=[args.xy0[0]+np.array([-1, 1])*args.Width, args.xy0[1]+np.array([-1, 1])*args.Width])
        col = np.argmin(np.abs(E_d2z0dx2.x-args.xy0[0]))
        row = np.argmin(np.abs(E_d2z0dx2.y-args.xy0[1]))
        args.E_d2z0dx2 = np.minimum(1.e-2, np.maximum(1.e-4, E_d2z0dx2.z[row,col]))
        if np.isnan(args.E_d2z0dx2):
            args.E_d2z0dx2=1.e-2

    if args.calc_error_file is not None:
        dest_dir=os.path.dirname(args.calc_error_file)
        args.out_name=args.calc_error_file
        with h5py.File(args.calc_error_file) as h5f:
            args.E_d2z0dx2 = h5f['E_RMS']['d2z0_dx2'][()]
            E_RMS['d2z0_dx2'] = args.E_d2z0dx2
        if not os.path.isfile(args.out_name):
            print(f"{args.out_name} not found, returning")
            return 1

    if args.out_name is None:
        args.out_name=dest_dir + '/E%d_N%d.h5' % (args.xy0[0]/1e3, args.xy0[1]/1e3)

    if args.calc_error_for_xy:
        args.calc_error_file=args.out_name
        if not os.path.isfile(args.out_name):
            print(f"{args.out_name} not found, returning")
            return 1

    if args.tide_adjustment_file is not None:
        args.tide_adjustment=True

    if args.error_res_scale is not None:
        if args.calc_error_file is not None:
            for ii, key in enumerate(['z0','dz']):
                spacing[key] *= args.error_res_scale[ii]

    args.bias_params=args.bias_params.split(',')

    if not os.path.isdir(args.base_directory):
        os.mkdir(args.base_directory)
    try:
        os.mkdir(dest_dir)
    except FileExistsError:
        pass

    print("ATL11_to_ATL15: working on "+args.out_name)

    S=ATL11_to_ATL15(args.xy0, ATL11_index=args.ATL11_index,
           Wxy=args.Width, E_RMS=E_RMS, t_span=args.time_span, spacing=spacing, \
           E_d3zdx2dt_scale_file=args.E_d3zdx2dt_scale_file,\
           bias_params=args.bias_params,\
           prior_edge_args=prior_edge_args, \
           sigma_geo=args.sigma_geo, \
           sigma_radial=args.sigma_radial, \
           hemisphere=args.Hemisphere, reread_dirs=reread_dirs, \
           data_file=args.data_file, \
           restart_edit=args.restart_edit, \
           out_name=args.out_name,
           dzdt_lags=args.dzdt_lags, \
           N_subset=args.N_subset,\
           mask_file=args.mask_file, \
           rock_mask_file=args.rock_mask_file, \
           rock_mask_reject_value=args.rock_mask_reject_value,\
           region=args.region, \
           geoid_file=args.geoid_file,\
           tide_mask_file=args.tide_mask_file, \
           tide_directory=args.tide_directory, \
           tide_adjustment=args.tide_adjustment, \
           tide_adjustment_file=args.tide_adjustment_file, \
           tide_adjustment_format=args.tide_adjustment_format, \
           tide_model=args.tide_model, \
           firn_directory=args.firn_directory,\
           firn_version=args.firn_version,\
           firn_correction=args.firn_model,\
           firn_grid_file=args.firn_grid_file,\
           max_iterations=args.max_iterations, \
           sigma_extra_bin_spacing=args.sigma_extra_bin_spacing,\
           sigma_extra_max=args.sigma_extra_max,\
           reference_epoch=args.reference_epoch, \
           W_edit=W_edit,\
           calc_error_file=args.calc_error_file, \
           error_res_scale=args.error_res_scale, \
           avg_scales=args.avg_scales, \
           verbose=args.verbose, \
           DEM_file=args.DEM_file, \
           DEM_tol=args.DEM_tol, \
           geoid_tol=args.geoid_tol,\
           sigma_tol=args.sigma_tol, \
           write_data_only=args.write_data_only)

    status=0
    if S is None:
        # indicates not enough data, but typically not a failure
        return status
    if args.write_data_only:
        S['data'].to_h5(args.out_name, group='data')
        return status
    status=1
    if args.calc_error_file is None and 'm' in S and len(S['m'].keys()) > 0:
        # if this isn't an error-calculation run, save the gridded fit data to the output file
        save_fit_to_file(S, args.out_name, dzdt_lags=args.dzdt_lags, reference_epoch=args.reference_epoch)
        status=0
    elif 'E' in S and len(S['E'].keys()) > 0:
        # If this is an error-calculation run, save the errors to the output file
        S['E']['sigma_z0']=interp_ds(S['E']['sigma_z0'], args.error_res_scale[0])
        for field in ['sigma_dz'] + [ f'sigma_dzdt_lag{lag}' for lag in args.dzdt_lags ]:
            S['E'][field] = interp_ds( S['E'][field], args.error_res_scale[1] )
        save_errors_to_file(S, args.out_name, dzdt_lags=args.dzdt_lags, reference_epoch=args.reference_epoch)
        status=0
    print(f"done with {args.out_name}")
    return status

if __name__=='__main__':
    status=main(sys.argv)
    if status is None:
        status=1
    sys.exit(status)

#-160000 -1800000 --centers @/home/ben/git_repos/surfaceChange/default_args/test.txt
#-160000 -1800000 --centers @/home/ben/git_repos/surfaceChange/default_args_z03xlooser_dt10xlooser_errors.txt -c /Volumes/ice2/ben/ATL14_test/IS2//U07/z03xlooser_dt10xlooser_40km/centers/E-160_N-1800.h5
