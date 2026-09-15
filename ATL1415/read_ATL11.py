#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 09:14:51 2025

@author: ben
"""
import pointCollection as pc
import numpy as np
import os
import re


def _lonlat_bounding_box(bounds, SRS_proj4):
    '''
    Compute a geographic (lon/lat) bounding box for a projected [x,y] box,
    for use in an earthaccess.search_data(bounding_box=...) call.

    Samples the box's 4 corners, 4 edge midpoints, and center (9 points) --
    plain corners can under-estimate the true lon/lat extent under
    polar-stereographic curvature near the pole; a little slack costs
    nothing since this only feeds a CMR candidate search over full
    orbit-length granules.

    Handles the antimeridian correctly: near-pole tiles routinely have
    sample points on both sides of the +-180 deg discontinuity even when
    nowhere near the true pole singularity. In that case the returned
    lon_min > lon_max, which is how a CMR bounding_box search signals
    "crosses the antimeridian" -- this is intentional, not a bug. If the
    box appears to span (nearly) the full longitude range even after
    accounting for that, it's assumed to contain the pole itself, and the
    full (-180, 180) longitude range is returned.

    inputs:
        bounds: 2-element iterable of 2-element iterables, [[xmin,xmax],[ymin,ymax]]
        SRS_proj4: proj4 string for the coordinate system of bounds
    output:
        (lon_min, lat_min, lon_max, lat_max)
    '''
    xs = [bounds[0][0], bounds[0][1], np.mean(bounds[0])]
    ys = [bounds[1][0], bounds[1][1], np.mean(bounds[1])]
    xg, yg = np.meshgrid(xs, ys)
    D = pc.data().from_dict({'x': xg.ravel(), 'y': yg.ravel()}).get_latlon(proj4_string=SRS_proj4)
    lat_min, lat_max = D.latitude.min(), D.latitude.max()
    lon = D.longitude
    lon_min, lon_max = lon.min(), lon.max()
    if lon_max - lon_min > 180:
        # points are likely wrapped around +-180, not genuinely spanning
        # most of the globe (true for any tile far smaller than the earth,
        # unless it actually contains the pole) -- recompute in a shifted
        # [0, 360) frame and convert back
        lon_shifted = np.mod(lon, 360)
        lo, hi = lon_shifted.min(), lon_shifted.max()
        if hi - lo > 180:
            # genuinely spans (near-)all longitudes -- the tile contains
            # or nearly surrounds the pole; search the whole longitude range
            lon_min, lon_max = -180., 180.
        else:
            # convert back to signed longitude; lon_min > lon_max here is
            # intentional and correct, see docstring
            lon_min = ((lo + 180) % 360) - 180
            lon_max = ((hi + 180) % 360) - 180
    return (lon_min, lat_min, lon_max, lat_max)


def select_best_xover_index(D):
    """
    Select the best crossing-track data for each crossover location

    For each crossover (defined by a unique ref rgt, crossing rgt, ref pair, and crossing pair),
    select the measurement from each cycle with the smallest error in corrected height

    input:
        D: pointCollection.data structure containing crossing track data
    output:
        ii: numpy array of indexes of minimum-error points
    """

    _, i_pts = pc.unique_by_rows(np.c_[D.rgt, D.ref_rgt, D.pair_track, D.ref_pair, D.cycle_number], return_dict=True)
    ii = np.zeros(len(i_pts), dtype=int)
    for count, (pt, i_pt) in enumerate(i_pts.items()):
        ii[count]=i_pt[np.argsort(D.h_corr_sigma[i_pt])[0]]
    return ii


# The ATL11 crossover generation string, e.g. '007_cycle_03_30_v03'.  The two
# captured groups are the release (007) and the version (03), which is how they
# appear in an ATL11XO granule name; the cycle range in the middle is part of
# the generation's name but not of any granule's.
XOVER_VERSION_RE = re.compile(r'(\d\d\d)_cycle_\d\d_\d\d_v(\d\d)')


def parse_ATL11xo_version(ATL11xo_version):
    """
    split an ATL11 crossover generation into its release and version

    Shared with scripts/setup_ATL11_xover.py, which writes the local schema
    files: the two must agree character for character, because the strings
    this produces end up in granule names that are matched exactly against
    CMR (cloud) or against files on disk (discover).

    input:
        ATL11xo_version: str, e.g. '007_cycle_03_30_v03'
    output:
        (release, version), e.g. ('007', '03')
    """
    try:
        return XOVER_VERSION_RE.search(ATL11xo_version).groups()
    except AttributeError:
        raise AttributeError('ATL11xo version did not match pattern '
                             '(rrr)_cycle_(cc)_(cc)_v(vv)')


def xover_tiling_schema(x_cycle, hemi, xover_tile_dir=None, ATL11xo_version=None):
    """
    get the tiling schema for one crossover cycle

    Local mode (xover_tile_dir given) reads the schema file that
    setup_ATL11_xover.py wrote next to the tiles.  Cloud mode builds the
    identical schema in memory and gives it an EarthAccess source, so that
    resolve_files_for_box() turns the same tile names into ATL11XO granule
    URLs via CMR instead of looking for them on a filesystem.  Nothing is
    staged for the cloud case -- that is the point: there is no crossover
    tree on the MAAP bucket and no schema file to write one next to.

    inputs:
        x_cycle: int, crossover cycle number
        hemi: str, 'AA' or 'AR'
        xover_tile_dir: str, optional. Local tile directory holding
            cycle_xx/200km_tiling_<hemi>.json.  If None, a cloud schema is
            built and ATL11xo_version is required.
        ATL11xo_version: str, optional. Crossover generation, e.g.
            '007_cycle_03_30_v03'.  Required in cloud mode.
    output:
        pc.tilingSchema
    """
    if xover_tile_dir is not None:
        schema_file = os.path.join(xover_tile_dir,
                                   f'cycle_{x_cycle:02d}',
                                   f'200km_tiling_{hemi}.json')
        return pc.tilingSchema().from_file(schema_file)

    if ATL11xo_version is None:
        raise ValueError('xover_tiling_schema: ATL11xo_version is required when '
                         'no xover_tile_dir is given')
    release, version = parse_ATL11xo_version(ATL11xo_version)
    # These values mirror scripts/setup_ATL11_xover.py exactly.  Tiles are
    # labelled by their CENTERS, hence mapping_function_name='round', and the
    # labels are in km, hence scale=1000.
    return pc.tilingSchema(
        tile_spacing=200e3,
        mapping_function_name='round',
        format_str=f'ATL11XO_{hemi}_E%d_N%d_c{x_cycle:02d}_{release}_{version}',
        scale=1000,
        extension='.h5',
        directory=None,
        source={'type': 'EarthAccess', 'short_name': 'ATL11XO', 'daac': 'NSIDC'})


def read_ATL11(xy0, Wxy, index_file, SRS_proj4, xover_tile_root=None,
               sigma_geo=6.5, sigma_radial=0.03, xover_cycles=[1,2],
               verbose=False, hemisphere=None, fs=None, earthaccess=False,
               ATL11_release=None, ATL11xo_version=None):


    bounds = [xy0[0]+np.array([-Wxy/2, Wxy/2]), xy0[1]+np.array([-Wxy/2, Wxy/2])]

    D_at, ATL11_file_list = read_ATL11_at(bounds, index_file, SRS_proj4,
                  sigma_geo=sigma_geo,
                  sigma_radial=sigma_radial,
                  earthaccess=earthaccess, fs=fs, verbose=verbose,
                  ATL11_release=ATL11_release)

    # exit if no data returned
    if D_at is None:
        return None, []

    # Two ways to get crossovers, and each mode has exactly one switch:
    # locally, xover_tile_root names the tile tree setup_ATL11_xover.py built;
    # in the cloud there is no tree to name, so ATL11xo_version turns them on
    # and the tiles are resolved from CMR.  Neither given means no crossovers,
    # which is a legitimate configuration -- but it used to be the ONLY cloud
    # outcome, silently, because xover_tile_root cannot be set in that mode.
    if xover_tile_root is None and not (earthaccess and ATL11xo_version is not None):
        return D_at, ATL11_file_list

    # Otherwise, read the crossover tiles
    D_xo, xover_file_list = read_ATL11_xovers(bounds, SRS_proj4,
                                              xover_tile_dir = xover_tile_root,
                                              ATL11xo_version = ATL11xo_version,
                                              xover_cycles = xover_cycles,
                                              verbose=verbose, hemisphere=hemisphere,
                                              fs=fs)
    return pc.data().from_list([D_at, D_xo]), ATL11_file_list + xover_file_list


def read_ATL11_at(bounds, index_file, SRS_proj4,
               sigma_geo=6.5, sigma_radial=0.03,
               earthaccess=False, fs=None, verbose=False, ATL11_release=None):
    '''
    read ATL11 data from an index file, or from NASA Earthdata Cloud

    inputs:
        xy0 : 2-element iterable specifying the domain center
        Wxy : Width of the domain
        index_file : local mode (earthaccess=False, default): file made by
            pointCollection.geoindex pointing at ATL11 data, covering the
            whole archive.
            cloud mode (earthaccess=True): root directory of per-granule
            geoIndex files (one 'ATL11_index_<cycles>_<release>_<version>/
            <granule>.h5' subtree per release combo -- see
            pointCollection.scripts.query_ATL11_cloud.index_path_for_granule).
        SRS_proj4: projection information for the data
        earthaccess: if True, search NASA Earthdata Cloud (via earthaccess)
            for ATL11 granules intersecting the domain and read them
            directly from S3, using index_file as the per-granule geoIndex
            root. Uses short_name='ATL11' and version_mismatch='error'.
        fs: s3fs.S3FileSystem, optional (cloud mode only); reused across
            granules to avoid re-deriving S3 credentials. If None, one is
            obtained automatically.
        verbose: if True (cloud mode only), report the number of candidate
            granules found.

    output:
        D: data structure
        file_list: list of ATL11 files read
    '''

    field_dict_11={None:['latitude','longitude','delta_time',\
                        'h_corr','h_corr_sigma','h_corr_sigma_systematic', 'ref_pt'],\
                        '__calc_internal__' : ['rgt'],
                        'cycle_stats' : {'tide_ocean','dac'},
                        'ref_surf':['e_slope','n_slope', 'x_atc', 'fit_quality', 'dem_h', 'geoid_h']}

    if earthaccess:
        from pointCollection.scripts.query_ATL11_cloud import (
            find_ATL11_granules, index_path_for_granule, read_ATL11_granule_cloud_items)

        bbox = _lonlat_bounding_box(bounds, SRS_proj4)
        # Filter to the generation whose per-granule index is staged: index_file
        # is the ROOT, holding one ATL11_index_<cycles>_<release>_<version>/
        # subtree per generation, and an unfiltered search returns every
        # generation CMR holds -- so a second staged subtree would make this
        # tile read both and double-count its data.  ATL11_release is the ATL11
        # generation ONLY; the crossovers have their own (--ATL11xo_version).
        granules = find_ATL11_granules(bbox, granule_release=ATL11_release)
        # Do NOT hoist the filesystem out of the loop when we own it.  The
        # brokered NSIDC credentials expire (roughly four hours), and get_s3fs()
        # re-derives when a cached session is close to that; asking it per
        # granule is a dict lookup in the normal case and is what lets the
        # refresh happen at all.  A near-pole Antarctic tile, whose bounding box
        # spans all longitudes, reads enough granules to get there.  A caller
        # that passed its own fs keeps it -- that is the caller's to manage.
        caller_supplied_fs = fs is not None
        if verbose:
            print(f'read_ATL11_at: found {len(granules)} candidate granules')
        D11_list = []
        for granule in granules:
            if not caller_supplied_fs:
                fs = pc.io_utils.get_s3fs(daac='NSIDC')
            s3_url = granule.data_links(access='direct')[0]
            idx_file = index_path_for_granule(os.path.basename(s3_url), index_file)
            items = read_ATL11_granule_cloud_items(s3_url, idx_file, bounds[0], bounds[1],
                                                    fields=field_dict_11, fs=fs,
                                                    version_mismatch='error')
            if items:
                D11_list.extend(items)
    else:
        try:
            # catch empty data
            D11_list=pc.geoIndex().from_file(index_file).query_xy_box(
                *bounds, fields=field_dict_11)
        except ValueError:
            return None, []
    if D11_list is None:
        return None, []
    D_list=[]

    if len(D11_list) == 0:
        # NO DATA HERE, which is a normal outcome on the outer ring of the
        # dilated tile grid, not an error.  This has to return None the way the
        # local branch's `except ValueError` above does: pc.data().from_list([])
        # returns a LIVE object with fields==[] rather than None (it returns
        # early only for D_list is None), so the caller's `if data is not None`
        # test at ATL11_to_ATL15.py:639 passes and data.sigma then raises
        # AttributeError.  Returning None instead lets ATL11_to_ATL15 take the
        # insufficient-data path it already has, which is what run.sh's
        # "no fit written ... skipping error calculation" guard expects.
        return None, []

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

    if len(D_list) == 0:
        # Granules intersected the bounding box, but no point survived the
        # in-bounds filter above -- same conclusion, same reason.
        return None, D11_files

    return pc.data().from_list(D_list), D11_files

def read_ATL11_xovers(bounds, SRS_proj4, xover_tile_dir=None, ATL11xo_version=None,
                      xover_cycles=[1,2], hemisphere=None, verbose=True, fs=None):
    '''
    read crossover data from tiles

    Parameters
    ----------
    bounds : 2-iterable of 2-iterables
        iter.
    SRS_proj4 : str
        proj4 string for the spatial reference system to be used.
    xover_tile_dir : str, optional
        local tile directory to be read. The default is None, which selects
        cloud mode and requires ATL11xo_version.
    ATL11xo_version : str, optional
        crossover generation, e.g. '007_cycle_03_30_v03'. Used in cloud mode
        to build the ATL11XO granule names; ignored when xover_tile_dir is
        given, since the schema file there already carries them.
    xover_cycles : iterble of ints, optional
        crossover cycles to be read. The default is [1,2]. Reading only the
        first two cycles is a deliberate, longstanding design choice, not a
        truncation to be widened later.
    verbose : bool, optional
        if True, report status
    fs : s3fs.S3FileSystem, optional
        filesystem to reuse for remote tiles (only relevant if the tiling
        schema's 'source' specifies a remote source). If None and a remote
        source is used, one is obtained automatically and reused across
        every tile/cycle in this call.

    Returns
    -------
    D_xo : pointCollection.data
        data object containing crossover data.
    xover_files_used : list
        crossover files read.

    '''

    if hemisphere==-1:
        hemi='AA'
    else:
        hemi='AR'

    D_x=[]
    D_d=[]
    D_r=[]
    xover_files_used = []
    for x_cycle in xover_cycles:
        schema = xover_tiling_schema(x_cycle, hemi,
                                     xover_tile_dir=xover_tile_dir,
                                     ATL11xo_version=ATL11xo_version)
        resolved, fs = schema.resolve_files_for_box(bounds, fs=fs, verbose=verbose)
        for tile_name, xover_file in resolved.items():
            with pc.io_utils.open_h5(xover_file, fs=fs) as h5f:
                D_ri = pc.data().from_h5(xover_file, h5_f=h5f).get_xy(proj4_string=SRS_proj4)
                keep = (D_ri.x >= bounds[0][0]) & (D_ri.x <= bounds[0][1]) &\
                     (D_ri.y >= bounds[1][0]) & (D_ri.y <= bounds[1][1])
                if not np.any(keep):
                    continue
                D_ri.index(keep)
                D_xi = pc.data().from_h5(xover_file, group='crossing_track', h5_f=h5f)
                D_xi.index(keep)
                D_di = pc.data().from_h5(xover_file, group='datum_track',
                                         fields=['rgt','ref_pt','pair_track',
                                                 'dem_h','geoid_h','fit_quality',
                                                 'n_slope','e_slope'], h5_f=h5f)
                D_di.index(keep)
            D_xi.assign(cycle_number=np.zeros_like(D_ri.x)+x_cycle)
            D_x += [D_xi]
            D_d += [D_di]
            D_r += [D_ri]
            xover_files_used += [xover_file]
    if len(D_x)==0 or not hasattr(D_x[0],'rgt'):
        return None, []
    D_x = pc.data().from_list(D_x)
    D_d = pc.data().from_list(D_d)
    D_r = pc.data().from_list(D_r)
    D_x.assign(ref_rgt = D_d.rgt, ref_pair=D_d.pair_track)
    # choose the smallest_sigma xover for each rgt and pair
    ii = select_best_xover_index(D_x)
    D_x = D_x[ii]
    D_d = D_d[ii]
    D_r = D_r[ii]

    blank = np.zeros_like(D_x.h_corr) + np.nan
    D_xo = pc.data().from_dict({
        'z':D_x.h_corr,
        'sigma':D_x.h_corr_sigma,
        'sigma_corr': D_x.h_corr_sigma_systematic,
        'x':D_r.x,
        'y':D_r.y,
        'latitude':D_r.latitude,
        'longitude':D_r.longitude,
        'dem_h':D_r.dem_h,
        'geoid_h':D_r.geoid_h,
        'rgt':D_x.rgt,
        'pair':D_x.pair_track,
        'ref_pt':D_d.ref_pt,
        'cycle':D_x.cycle_number,
        'n_cycles':blank,
        'fit_quality':D_d.fit_quality,
        'tide_ocean':D_x.tide_ocean,
        'dac':D_x.dac,
        'delta_time':D_x.delta_time,
        'n_slope':D_r.n_slope,
        'e_slope':D_r.e_slope,
        'time': D_x.delta_time/24/3600/365.25+2018,
        'along_track':np.zeros_like(D_r.x, dtype=bool)})
    if verbose:
        print(f"read_ATL11_xovers: read {D_xo.size} crossing_track measurements from {len(xover_files_used)} files")
    return D_xo, xover_files_used
