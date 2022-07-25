#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 15:21:52 2022

@author: ben
"""
import numpy as np
import pointCollection as pc

def update_masks_with_geoid(grids, m, args):
    
    # read the geoid:
    geoid_file=args['ancillary_data']['geoid_file']
    geoid = pc.grid.data().from_nc(geoid_file,
        xname='lon', yname='lat', field_mapping=dict(z='geoid_h'))

    # interpolate the geoid to the z0 grid
    z0=m['z0'].copy()
    z0.get_latlon(srs_proj4=args['srs_proj4'])
    P2 = 0.5*(3.0*np.sin(z0.latitude*np.pi/180.0)**2 - 1.0)
    z0.assign({'geoid_z' : geoid.interp(z0.longitude % 360.0, z0.latitude) \
               + ( -0.198*P2*(1.0 + 0.3))})
    # mask the z0 nodes for which z0 is below the geoid
    below_geoid = z0.z0 < z0.geoid_z
    grids['z0'].mask[below_geoid] = 0

    # interpolate the geoid to the dz grid
    dz_geoid = z0.interp(grids['dz'].ctrs[1],
                         grids['dz'].ctrs[0], gridded=True, field='geoid_z')
    # interpolate z0 to the dz grid
    dz_z0 = z0.interp(grids['dz'].ctrs[1],
                         grids['dz'].ctrs[0], gridded=True, field='z0')

    # zero the mask and cell_area for nodes that are below the geoid
    if grids['dz'].mask_3d is not None:
        for ti in range(grids['dz'].mask_3d.shape[2]):
            above_geoid = (dz_z0 + m['dz'].dz[:,:,ti]) > dz_geoid
            grids['dz'].mask_3d.z[:,:,ti] &= above_geoid
            grids['dz'].cell_area[:,:,ti] *= above_geoid.astype(float)

