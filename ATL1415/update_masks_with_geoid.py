#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 15:21:52 2022

@author: ben
"""
import numpy as np
import pointCollection as pc

def _geoid_read_bounds(z0, pad_deg=1.0):
    """
    The [[lon0, lon1], [lat0, lat1]] window of the geoid this tile needs.

    EGM2008_geoid_h.nc is a global 2.5-arcminute grid: 4321 x 8641 float32,
    81.7 MiB, chunked (721, 1441).  A 60 km tile covers well under a degree of
    latitude, so reading the whole file to interpolate one tile fetches roughly
    twenty times more than the window does, on every tile of a run of
    thousands.

    Returns None only if the tile has no finite coordinates.  A tile that
    crosses the 0/360 seam, or that encircles a pole, has no contiguous window
    in this coordinate and reads every longitude -- still a large saving,
    because the latitude band is what dominates.

    The pad is generous on purpose: interpolating outside the window returns
    NaN, and a NaN geoid would silently mis-mask the tile rather than fail.
    """
    lon = np.asarray(z0.longitude) % 360.0
    lat = np.asarray(z0.latitude)
    good = np.isfinite(lon) & np.isfinite(lat)
    if not np.any(good):
        return None

    lat_range = [max(-90.0, float(np.min(lat[good])) - pad_deg),
                 min(90.0, float(np.max(lat[good])) + pad_deg)]

    lon_min = float(np.min(lon[good]))
    lon_max = float(np.max(lon[good]))
    if (lon_max - lon_min) > 180.0:
        return [[0.0, 360.0], lat_range]
    return [[max(0.0, lon_min - pad_deg),
             min(360.0, lon_max + pad_deg)], lat_range]


def update_masks_with_geoid(grids, m, args):

    # interpolate the geoid to the z0 grid
    z0=m['z0'].copy()
    z0.get_latlon(srs_proj4=args['srs_proj4'])

    # read the geoid, windowed to this tile.  The read has to come AFTER
    # get_latlon(), which is what says which window is needed.
    geoid_file=args['ancillary_data']['geoid_file']
    geoid = pc.grid.data().from_nc(geoid_file,
        xname='lon', yname='lat', field_mapping=dict(z='geoid_h'),
        bounds=_geoid_read_bounds(z0))

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

