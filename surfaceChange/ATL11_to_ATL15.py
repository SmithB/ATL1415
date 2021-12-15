#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 17:15:08 2019

@author: ben
 """
import os
os.environ["MKL_NUM_THREADS"]="1"  # multiple threads don't help that much and tend to eat resources
os.environ["OPENBLAS_NUM_THREADS"]="1"  # multiple threads don't help that much and tend to eat resources

import numpy as np
from LSsurf.smooth_xytb_fit import smooth_xytb_fit
import pointCollection as pc

import re
import sys
import h5py
import traceback
import matplotlib.pyplot as plt
from surfaceChange.reread_data_from_fits import reread_data_from_fits
import pyTMD
import scipy.optimize
import pdb

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

def read_ATL11(xy0, Wxy, index_file, SRS_proj4, sigma_geo=6.5):
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
                        'ref_surf':['e_slope','n_slope', 'x_atc', 'fit_quality', 'dem_h']}
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
    file_list= [ Di.filename for Di in D11_list ]
    for D11 in D11_list:
        D11.get_xy(proj4_string=SRS_proj4)
        # select the subset of the data within the domain
        D11.index((np.abs(D11.x[:,0]-xy0[0]) <= Wxy/2) &\
                 (np.abs(D11.y[:,0]-xy0[1]) <= Wxy/2))
        if D11.size==0:
            continue
        sigma_corr=np.sqrt((sigma_geo*np.abs(np.median(D11.n_slope)))**2+\
                           (sigma_geo*np.abs(np.median(D11.e_slope)))**2+0.03**2)
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
           'delta_time': D11.delta_time,
           'time':D11.delta_time/24/3600/365.25+2018,
           'n_slope':D11.n_slope,
           'e_slope':D11.e_slope,
           'along_track':np.ones_like(D11.x, dtype=bool)})]

        if len(D11.ref_pt) == 0:
            continue
        # N.B.  D11 is getting indexed in this step, and it's leading to the warning in
        # line 76.  Can fix by making crossover_data.from_h5 copy D11 on input
        D_x = pc.ATL11.crossover_data().from_h5(D11.filename, pair=D11.pair, D_at=D11)
        if D_x is None:
            continue
        # fit_quality and dem_h D11 apply to all cycles, but are mapped only to some specific
        # cycle of the reference cycle
        temp={'fit_quality':np.nanmax(D_x.fit_quality[:,:,0], axis=1),
              'dem_h':np.nanmax(D_x.dem_h[:,:,0], axis=1)}       
        for cycle in range(D_x.shape[1]):
            D_x.fit_quality[:,cycle,1]=temp['fit_quality']
            D_x.dem_h[:,cycle,1]=temp['dem_h']
        
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
        D.index(np.isfinite(D.z))
        
    except ValueError:
        # catch empty data
        return None, file_list

    D.index(( D.fit_quality ==0 ) | ( D.fit_quality == 2 ))

    return D, file_list

def apply_tides(D, xy0, W,
                tide_mask_file=None,
                tide_mask_data=None,
                tide_directory=None,
                tide_model=None,
                tide_adjustment=False,
                tide_adjustment_file=None,
                sigma_tolerance=0.5,
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
        tide_adjustment: use bounded least-squares fit to adjust tides for extrapolated points
        tide_adjustment_file: file specifying the areas for which the tide model should be empirically adjusted
        sigma_tolerance: maximum height sigmas allowed in tidal adjustment fit
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

    is_els=tide_mask.interp(D.x, D.y) > 0.5
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
    # use a bounded least-squares fit to adjust tides
    # this is a purely empirical correction that
    # does not take into account ice shelf flexure physics
    if np.any(is_els.ravel()) and tide_adjustment:
        D.assign({'tide_adj_scale': np.ones_like(D.x)})
        D.tide_adj_scale[is_els==0]=0.0
        tide_adj_sigma = np.zeros_like(D.x) + np.inf
        # adjust indices that are extrapolated of within a defined mask
        if not tide_adjustment_file:
            # check if point is within model domain
            inmodel = pyTMD.check_tide_points(D.x, D.y,
                DIRECTORY=tide_directory, MODEL=tide_model,
                EPSG=EPSG)
            # find where points have an extrapolated tide value
            # only calculate adjustments for ice shelf values
            adjustment_indices, = np.nonzero(np.isfinite(D.tide_ocean) &
                np.logical_not(inmodel) & is_els.ravel() & D.along_track)
        else:
            # read adjustment mask and calculate indices
            try:
                adjustment_mask = pc.grid.data().from_geotif(tide_adjustment_file,
                    bounds=[np.array([-0.6, 0.6])*W+xy0[0], np.array([-0.6, 0.6])*W+xy0[1]])
                adjustment_indices, = np.nonzero((adjustment_mask.interp(D.x, D.y) > 0.5) &
                    np.isfinite(D.tide_ocean) & is_els.ravel() & D.along_track)
            except IndexError:
                adjustment_indices = []
        # make a global reference point number combining ref_pt, rgt and pair
        # convert both pair and rgt to zero-based indices
        global_ref_pt = 3*1387*D.ref_pt + 3*(D.rgt-1) + (D.pair-1)
        u_ref_pt = np.unique(global_ref_pt[adjustment_indices])
        print(f"\t\ttide adjustment: {len(u_ref_pt)} refs") if verbose else None
        print(f"\t\t{tide_adjustment_file}") if tide_adjustment_file else None
        # calculate temporal fit of tide points with model phases
        # only for reference points that are extrapolated
        for ref_pt in u_ref_pt:
            # indices for along-track coordinates for reference point
            iref, = np.nonzero((global_ref_pt == ref_pt) & D.along_track)
            # calculate distance from central point
            x = np.median(D.x[iref])
            y = np.median(D.y[iref])
            dist = np.sqrt((D.x - x)**2 + (D.y - y)**2)
            # indices of nearby points (include nearby crossover points)
            # reduce to ice shelf points with errors less than tolerance
            ii, = np.nonzero(((global_ref_pt == ref_pt) | (dist <= 100)) &
                np.isfinite(D.tide_ocean) & is_els.ravel() &
                (D.sigma < sigma_tolerance))
            # calculate differences in spatial coordinates
            dx = D.x[ii] - x
            dy = D.y[ii] - y
            # check if minimum number of values for fit
            if (len(ii) < 4):
                # continue to next global reference point
                continue
            elif np.any(dx**2 + dy**2) and (len(ii) < 6):
                # continue to next global reference point
                continue
            # reduce time and ocean amplitude to global reference point
            t = np.copy(D.time[ii])
            tide = np.copy(D.tide_ocean[ii])
            # correct heights for ocean variability
            dac = np.copy(D.dac[ii])
            dac[~np.isfinite(dac)] = 0.0
            # reduce DAC-corrected heights to global reference point
            h = D.z[ii] - dac
            # use linear least-squares with bounds on the variables
            # create design matrix
            p0 = np.ones_like(t)
            p1 = t - np.median(t)
            DMAT = np.c_[tide,p0,p1]
            # check if there are enough unique dates for fit
            u_days = np.unique(np.round(p1*365.25))
            if (len(u_days) <= 3):
                continue
            # tuple for parameter bounds (lower and upper)
            # output tidal amplitudes must be between 0 and 1 of model
            # average height must be between minimum and maximum of h
            # elevation change rate must be between range
            lb,ub = ([0.0,np.nanmin(h),-2.0],[1.0,np.nanmax(h),2.0])
            # check if there are coordinates away from central point
            if np.any(dx**2 + dy**2):
                # append horizontal coordinates to design matrix
                DMAT = np.c_[DMAT,dx,dy]
                # calculate min and max of surface slopes
                e_slope_min = np.nanmin(D.e_slope[ii])
                e_slope_max = np.nanmax(D.e_slope[ii])
                n_slope_min = np.nanmin(D.n_slope[ii])
                n_slope_max = np.nanmax(D.n_slope[ii])
                # append lists for parameter bounds (lower and upper)
                # fit slopes must be within the range of ATL11 slopes
                lb.extend([e_slope_min,n_slope_min])
                ub.extend([e_slope_max,n_slope_max])
            # calculate degrees of freedom
            n_max,n_terms = np.shape(DMAT)
            nu = np.float64(n_max - n_terms)
            # attempt bounded linear least squares
            try:
                # results from bounded least squares
                results = scipy.optimize.lsq_linear(DMAT, h,
                    bounds=(lb,ub))
                # model covariance matrix
                Hinv = np.linalg.inv(np.dot(DMAT.T,DMAT))
            except:
                # print exceptions
                if verbose:
                    traceback.print_exc()
                # continue to next global reference point
                # and use original tide values
                continue
            else:
                # extract tidal adjustment estimate
                adj,*_ = np.copy(results['x'])
                # calculate mean square error
                MSE = np.sum(results['fun']**2)/nu
                # standard error from covariance matrix
                adj_sigma,*_ = np.sqrt(MSE*np.diag(Hinv))
            # use best case fits for each point
            imin, = np.nonzero(tide_adj_sigma[ii] >= adj_sigma)
            # copy adjustments and estimated uncertainties
            D.tide_adj_scale[ii[imin]] = adj
            tide_adj_sigma[ii[imin]] = adj_sigma
        # multiply tides by adjustments
        ii, = np.nonzero((D.tide_adj_scale < 1) & (tide_adj_sigma < 1))
        D.tide_ocean[ii] *= D.tide_adj_scale[ii]
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
    for bin0 in np.unique(ij_bin):
        ii = np.flatnonzero((ij_bin==bin0) & D.along_track)
        this_rho = len(ii) / (W_sub**2)
        #print(f'bin0={bin0}, this_rho={this_rho}, ratio={this_rho/rho_target}')
        if this_rho < rho_target:
            # if there aren't too many points in this bin, continue
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
        ind_buffer.append(ii[isub])
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
            dzdt_lags=[1, 4],\
            hemisphere=1, reference_epoch=None,\
            region=None, reread_dirs=None, \
            data_file=None, \
            max_iterations=5, \
            N_subset=8,  \
            W_edit=None, \
            out_name=None, \
            compute_E=False, \
            DEM_file=None,\
            DEM_tol=None,\
            sigma_tol=None,\
            mask_file=None,\
            tide_mask_file=None,\
            tide_directory=None,\
            tide_adjustment=False,\
            tide_adjustment_file=None,\
            tide_model=None,\
            avg_scales=None,\
            edge_pad=None,\
            error_res_scale=None,\
            calc_error_file=None,\
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
        W_edit: (float) Only data within this distance of the grid center can be editied in the iterative step
        out_name: (string) Name of the output file
        compute_E: (bool) If true, errors will be calculated
        DEM_file: (string) DEM against which the data will be checked using the DEM_tol parameter
        DEM_tol: (float) Points different from the DEM by more than this value will be edited
        sigma_tol: (float) Points with a sigma parameter greater than this value will be edited
        mask_file: (string) File specifying areas for which data should be used and strong constraints should be applied
        tide_mask_file: (string)  File specifying the areas for which the tide model should be calculated
        tide_directory: (string)  Directory containing the tide model data
        tide_adjustment: (bool)  Use bounded least-squares fit to adjust tides for extrapolated points
        tide_adjustment_file: (string)  File specifying the areas for which the tide model should be empirically adjusted
        tide_model: (string)  Name of the tide model to use for a given domain
        avg_scales: (list of floats) scales over which the output grids will be averaged and errors will be calculated
        error_res_scale: (float) If errors are calculated, the grid resolution will be coarsened by this factor
        calc_error_file: (string) Output file for which errors will be calculated.
        verbose: (bool) Print progress of processing run
    outputs:
        S: dict containing fit output
    '''
    SRS_proj4, EPSG=get_SRS_info(hemisphere)

    E_RMS0={'d2z0_dx2':200000./3000/3000, 'd3z_dx2dt':3000./3000/3000, 'd2z_dxdt':3000/3000, 'd2z_dt2':5000}
    E_RMS0.update(E_RMS)

    W={'x':Wxy, 'y':Wxy,'t':np.diff(t_span)}
    ctr={'x':xy0[0], 'y':xy0[1], 't':np.mean(t_span)}

    # initialize file_list to empty in case we're rereading the data
    file_list=[]
    
    # figure out where to get the data
    if data_file is not None:
        data=pc.data().from_h5(data_file, group='data')
    elif calc_error_file is not None:
        data=pc.data().from_h5(calc_error_file, group='data')
        max_iterations=0
        compute_E=True
        N_subset=None
    elif reread_dirs is not None:
        data = reread_data_from_fits(xy0, Wxy, reread_dirs, template='E%d_N%d.h5')
    else:
        data, file_list = read_ATL11(xy0, Wxy, ATL11_index, SRS_proj4, sigma_geo=sigma_geo)
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
            return None

    if W_edit is not None:
        # this is used for data that we are rereading from a set of other files.
        # data that are far from the center of this file cannot be eliminated by
        # the three-sigma edit
        data.assign({'editable':  (np.abs(data.x-xy0[0])<=W_edit/2) & (np.abs(data.y-xy0[1])<=W_edit/2)})


    # work out which mask to use based on the region
    tide_mask_data=None
    mask_data=None
    if region is not None:
        if region=='AA':
            mask_data=pc.grid.data().from_geotif(mask_file, bounds=[xy0[0]+np.array([-1.2, 1.2])*Wxy/2, xy0[1]+np.array([-1.2, 1.2])*Wxy/2])
            import scipy.ndimage as snd
            mask_data.z=snd.binary_erosion(snd.binary_erosion(mask_data.z, np.ones([1,3])), np.ones([3,1]))
            mask_file=None
        if region=='GL' and mask_file.endswith('.nc'):
            mask_data, tide_mask_data = read_bedmachine_greenland(mask_file, xy0, Wxy)

    # new 8/11/2021: use a DEM (if provided) to reject ATL06/11 blunders
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
                    EPSG=EPSG, verbose=verbose)

    if write_data_only:
        return {'data':data}

    # call smooth_xytb_fitting
    S=smooth_xytb_fit(data=data, ctr=ctr, W=W, spacing=spacing, E_RMS=E_RMS0,
                     reference_epoch=reference_epoch, N_subset=N_subset, compute_E=compute_E,
                     bias_params=['rgt','cycle'],  max_iterations=max_iterations,
                     srs_proj4=SRS_proj4, VERBOSE=True, dzdt_lags=dzdt_lags,
                     mask_file=mask_file, mask_data=mask_data, mask_scale={0:10, 1:1},
                     converge_tol_frac_edit=0.001,
                     error_res_scale=error_res_scale,
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
    for key , ds in S['m'].items():
        if isinstance(ds, pc.grid.data):
            ds.to_h5(filename, group=key)
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
    for field in ds.fields:
        delta_xy=[(ds.x[1]-ds.x[0])/scale, (ds.y[1]-ds.y[0])/scale]
        xi=np.arange(ds.x[0], ds.x[-1]+delta_xy[0], delta_xy[0])
        yi=np.arange(ds.y[0], ds.y[-1]+delta_xy[1], delta_xy[1])
        z0=getattr(ds, field)
        if len(ds.shape)==2:
            zi=pc.grid.data().from_dict({'x':ds.x, 'y':ds.y, 'z':z0}).interp(xi, yi, gridded=True)
            return pc.grid.data().from_dict({'x':xi, 'y':yi, field:zi})
        else:
            zi=np.zeros([xi.size, yi.size, ds.time.size])
            for epoch in range(ds.time.size):
                temp=pc.grid.data().from_dict({'x':ds.x, 'y':ds.y, 'z':np.squeeze(z0[:,:,epoch])})
                zi[:,:,epoch] = temp.interp(xi, yi, gridded=True)
            return pc.grid.data().from_dict({'x':xi, 'y':yi, 'time':ds.time, field:zi})

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
    parser.add_argument('--ATL11_index', type=str, required=True, help="ATL11 index file")
    parser.add_argument('--Width','-W',  type=float, help="Width of grid")
    parser.add_argument('--time_span','-t', type=str, help="time span, first year,last year AD (comma separated, no spaces)")
    parser.add_argument('--grid_spacing','-g', type=str, help='grid spacing:DEM (meters),dh maps xy (meters),dh_maps time (years): comma-separated, no spaces', default='250.,4000.,1.')
    parser.add_argument('--Hemisphere','-H', type=int, default=1, help='hemisphere: -1=Antarctica, 1=Greenland')
    parser.add_argument('--base_directory','-b', type=str, help='base directory')
    parser.add_argument('--out_name', '-o', type=str, help="output file name")
    parser.add_argument('--centers', action="store_true")
    parser.add_argument('--edges', action="store_true")
    parser.add_argument('--corners', action="store_true")
    parser.add_argument('--W_edit', type=int)
    parser.add_argument('--E_d2zdt2', type=float, default=5000)
    parser.add_argument('--E_d2z0dx2', type=float, default=0.02)
    parser.add_argument('--E_d3zdx2dt', type=float, default=0.0003)
    parser.add_argument('--E_d2z0dx2_file', type=str, help='file from which to read the expected d2z0dx2 values')
    parser.add_argument('--data_gap_scale', type=float,  default=2500)
    parser.add_argument('--sigma_geo', type=float,  default=6.5)
    parser.add_argument('--dzdt_lags', type=str, default='1,4', help='lags for which to calculate dz/dt, comma-separated list, no spaces')
    parser.add_argument('--avg_scales', type=str, default='10000,40000', help='scales at which to report average errors, comma-separated list, no spaces')
    parser.add_argument('--N_subset', type=int, default=None, help="number of pieces into which to divide the domain for (cheap) editing iterations.")
    parser.add_argument('--max_iterations', type=int, default=6, help="maximum number of iterations used to edit the data.")
    parser.add_argument('--map_dir','-m', type=str)
    parser.add_argument('--DEM_file', type=str, help='DEM file to use with the DEM_tol parameter')
    parser.add_argument('--DEM_tol', type=float, default=50, help='points different from the DEM by more than this value will be edited in the first iteration')
    parser.add_argument('--sigma_tol', type=float, help='points with sigma greater than this value will be edited')
    parser.add_argument('--mask_file', type=str)
    parser.add_argument('--tide_mask_file', type=str)
    parser.add_argument('--tide_directory', type=str)
    parser.add_argument('--tide_adjustment', action="store_true", help="Use bounded least-squares fit to adjust tides")
    parser.add_argument('--tide_adjustment_file', type=str, help="File specifying the areas for which the tide model should be empirically adjusted")
    parser.add_argument('--tide_model', type=str)
    parser.add_argument('--reference_epoch', type=int, default=0, help="Reference epoch number, for which dz=0")
    parser.add_argument('--data_file', type=str, help='read data from this file alone')
    parser.add_argument('--calc_error_file','-c', type=str, help='file containing data for which errors will be calculated')
    parser.add_argument('--calc_error_for_xy', action='store_true', help='calculate the errors for the file specified by the x0, y0 arguments')
    parser.add_argument('--error_res_scale','-s', type=float, nargs=2, default=[4, 2], help='if the errors are being calculated (see calc_error_file), scale the grid resolution in x and y to be coarser')
    parser.add_argument('--region', type=str, help='region for which calculation is being performed')
    parser.add_argument('--verbose','-v', action="store_true")
    parser.add_argument('--write_data_only', action='store_true', help='save data without processing')
    args, unknown=parser.parse_known_args()

    args.grid_spacing = [np.float64(temp) for temp in args.grid_spacing.split(',')]
    args.dzdt_lags = [np.int64(temp) for temp in args.dzdt_lags.split(',')]
    args.time_span = [np.float64(temp) for temp in args.time_span.split(',')]
    args.avg_scales = [np.int64(temp) for temp in args.avg_scales.split(',')]

    if args.E_d2z0dx2_file is not None and args.calc_error_file is None:
        E_d2z0dx2=pc.grid.data().from_geotif(args.E_d2z0dx2_file)#, bounds=[args.xy0[0]+np.array([-1, 1])*args.Width, args.xy0[1]+np.array([-1, 1])*args.Width])
        col = np.argmin(np.abs(E_d2z0dx2.x-args.xy0[0]))
        row = np.argmin(np.abs(E_d2z0dx2.y-args.xy0[1]))
        args.E_d2z0dx2 = np.minimum(1.e-2, np.maximum(1.e-4, E_d2z0dx2.z[row,col]))
        if np.isnan(args.E_d2z0dx2):
            args.E_d2z0dx2=1.e-2

    spacing={'z0':args.grid_spacing[0], 'dz':args.grid_spacing[1], 'dt':args.grid_spacing[2]}
    E_RMS={'d2z0_dx2':args.E_d2z0dx2, 'd3z_dx2dt':args.E_d3zdx2dt, 'd2z_dxdt':args.E_d3zdx2dt*args.data_gap_scale,  'd2z_dt2':args.E_d2zdt2}

    reread_dirs=None
    dest_dir=args.base_directory
    if args.W_edit is None:
        W_edit=args.Width/2
    else:
        W_edit=args.W_edit

    if args.centers:
        dest_dir += '/centers'
        W_edit=None
    if args.edges or args.corners:
        reread_dirs=[args.base_directory+'/centers']
        if args.edges:
            dest_dir += '/edges'
    if args.corners:
        reread_dirs += [args.base_directory+'/edges']
        dest_dir +='/corners'

    if args.calc_error_file is not None:
        dest_dir=os.path.dirname(args.calc_error_file)
        # get xy0 from the filename
        re_match=re.compile('E(.*)_N(.*).h5').search(args.calc_error_file)
        args.xy0=[float(re_match.group(ii))*1000 for ii in [1, 2]]
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

    if args.error_res_scale is not None:
        if args.calc_error_file is not None:
            for ii, key in enumerate(['z0','dz']):
                spacing[key] *= args.error_res_scale[ii]

    if not os.path.isdir(args.base_directory):
        os.mkdir(args.base_directory)
    try:
        os.mkdir(dest_dir)
    except FileExistsError:
        pass
    S=ATL11_to_ATL15(args.xy0, ATL11_index=args.ATL11_index,
           Wxy=args.Width, E_RMS=E_RMS, t_span=args.time_span, spacing=spacing, \
           sigma_geo=args.sigma_geo, \
           hemisphere=args.Hemisphere, reread_dirs=reread_dirs, \
           data_file=args.data_file, \
           out_name=args.out_name,
           dzdt_lags=args.dzdt_lags, \
           N_subset=args.N_subset,\
           mask_file=args.mask_file, \
           region=args.region, \
           tide_mask_file=args.tide_mask_file, \
           tide_directory=args.tide_directory, \
           tide_adjustment=args.tide_adjustment, \
           tide_adjustment_file=args.tide_adjustment_file, \
           tide_model=args.tide_model, \
           max_iterations=args.max_iterations, \
           reference_epoch=args.reference_epoch, \
           W_edit=W_edit,\
           calc_error_file=args.calc_error_file, \
           error_res_scale=args.error_res_scale, \
           avg_scales=args.avg_scales, \
           verbose=args.verbose, \
           DEM_file=args.DEM_file, \
           DEM_tol=args.DEM_tol, \
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
