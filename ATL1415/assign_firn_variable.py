#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 14:08:42 2019.

@author: ben
"""

import numpy as np
import SMBcorr
import os


def assign_firn_variable(data, firn_correction, firn_dir, hemisphere,\
                         model_version='1.0', subset_valid=False, 
                         variables=None, infer_FAC=True, rho_water=1):

    EPSG={-1:'EPSG:3031', 1:'EPSG:3413'}[hemisphere]
    
    if firn_correction == 'MAR':
        if hemisphere==1:
            data.assign({'h_firn' : SMBcorr.interp_MAR_firn(data.x, data.y, data.time)})
    elif firn_correction == 'MERRA2_hybrid':
        KWARGS={}
        # get MERRA-2 version and major version
        MERRA2_VERSION = model_version
        # MERRA-2 hybrid directory
        DIRECTORY=os.path.join(firn_dir,'MERRA2_hybrid',MERRA2_VERSION)
        # MERRA-2 region name from ATL11 region
        MERRA2_REGION = {-1:'ais',1:'gris'}[hemisphere]
        # keyword arguments for MERRA-2 interpolation programs
        if MERRA2_VERSION in ('v0','v1','v1.0'):
            KWARGS['VERSION'] = MERRA2_VERSION
            DEFAULT_VARIABLES = {'FAC':'FAC', 'SMB_a':'cum_smb_anomaly', 'h_a':'height'}
        else:
            KWARGS['VERSION'] = MERRA2_VERSION.replace('.','_')
            DEFAULT_VARIABLES = {'FAC':'FAC', 'SMB_a':'SMB_a','h_a':'h_a'}
        # use compressed files
        if hemisphere==1:
            KWARGS['GZIP'] = False
        else:
            KWARGS['GZIP'] = False
        if variables is None:
            VARIABLES = DEFAULT_VARIABLES
        else:
            VARIABLES=variables
        # output variable keys for both direct and derived fields
        SMB_data={out_var: SMBcorr.interpolate_merra_hybrid(DIRECTORY, EPSG,
                            MERRA2_REGION, data.time, data.x, data.y,
                            VARIABLE=model_var, **KWARGS) for out_var, model_var in VARIABLES.items()}
        if infer_FAC:
            SMB_data['FAC']=SMB_data['h_a']-SMB_data['SMB_a']
            
        if 'floating' in data.fields:
            # assume SMB is in ice equivalent.  
            #  For grounded ice, dh = FAC + SMB
            #  For floating ice, dh = FAC + (rho_water-rho_i)/(rho_w) SMB
            float_scale = (data.floating==0) + (rho_water-.917)/rho_water*(data.floating==1)
            data.assign({'h_firn':SMB_data['FAC'] + float_scale*SMB_data['SMB_a']})
        else:
            data.assign({'h_firn':SMB_data['FAC'] + SMB_data['SMB_a']})
            
        # use the mask values to set the ouput to NaN if the model is invalid
        if infer_FAC:
            data.h_firn[SMB_data['FAC'].mask]=np.NaN
        data.h_firn[SMB_data['SMB_a'].mask]=np.NaN
    if subset_valid:
        data.index(np.isfinite(data.h_firn))
