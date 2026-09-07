#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find the ATL14/ATL15 files of a previous release, locally or in the cloud.

Split out of ATL11_to_ATL15.set_three_sigma_edit_from_previous_product so that
discovery (which differs between a discover run and a MAAP run) is separate
from the residual/sigma_extra arithmetic (which does not differ at all).

See docs/Transition_to_maap.md, Q27, work items W1 and W3.
"""
import glob
import os
import re

import pointCollection as pc

# ATL14_IS_0329_100m_005_02.nc  /  ATL15_A1_0329_01km_005_02.nc
# The two products share a field layout: <product>_<region>_<cycles>_
# <resolution>_<release>_<revision>.nc .  The cycles and release fields are
# what Q27 F2 requires filtering on: the published collection mixes cycle
# ranges, so ATL14_A1_0328_100m_005_01.nc sits beside the 0329 granules under
# the same short_name, and mosaicking the two together is wrong data with no
# error.
GRANULE_RE = re.compile(
    r'^(?P<product>ATL1[45])_(?P<region>[A-Za-z0-9]+)_(?P<cycles>\d{4})_'
    r'(?P<resolution>[0-9a-z]+)_(?P<release>\d{3})_(?P<revision>\d+)\.nc$')

# ATL15 is published at four averaging scales; only the finest is used here,
# matching the 'ATL15_*1km_*.nc' glob the local path has always used.
ATL15_RESOLUTION = '01km'


RELEASE_SPEC_RE = re.compile(r'^(?:rel)?(?P<release>\d{3})_(?P<cycles>\d{4})$')


def previous_product_arg(value):
    '''
    argparse type for --previous_product: a release spec, or a path/URI.

    --previous_product is normally a directory and is normalized by
    paths.path_or_uri(), which makes it absolute against the working directory.
    A cloud release spec must survive that untouched -- '005_0329' would
    otherwise become '<cwd>/005_0329' and no longer parse as a release at all.
    Shape decides: '<release>_<cycles>' is a spec, anything else is a path.
    '''
    from ATL1415.paths import path_or_uri
    if value is not None and RELEASE_SPEC_RE.match(str(value).strip()):
        return str(value).strip()
    return path_or_uri(value)


def parse_release_spec(spec):
    '''
    Split a cloud previous-product spec into (release, cycles).

    The spec replaces a directory when --previous_product_earthaccess is set,
    and is written '<release>_<cycles>' -- '005_0329' -- which is the same pair
    the discover directory name rel005_0329 carries.  Both fields are required:
    a search on release alone would mix cycle ranges (Q27 F2).

    inputs:
        spec (str): '<release>_<cycles>', e.g. '005_0329'
    output:
        (release, cycles) as strings, e.g. ('005', '0329')
    '''
    m = RELEASE_SPEC_RE.match(str(spec).strip())
    if m is None:
        raise ValueError(
            f"--previous_product={spec!r} is not a valid cloud previous-product spec.  "
            "With --previous_product_earthaccess set, --previous_product is the release "
            "and cycle range to search for, written '<release>_<cycles>' (e.g. '005_0329'), "
            "not a directory.")
    return m['release'], m['cycles']


def _latest_revisions(urls):
    '''
    Keep one granule per (product, region, resolution): the highest revision.

    CMR returns every revision of a granule, so ATL14_A1_0329_100m_005_01.nc and
    ..._005_02.nc can both come back; mosaicking both would let the older one
    fill holes in the newer.
    '''
    best = {}
    for url in urls:
        m = GRANULE_RE.match(os.path.basename(url))
        if m is None:
            continue
        key = (m['product'], m['region'], m['resolution'])
        rev = int(m['revision'])
        if key not in best or rev > best[key][0]:
            best[key] = (rev, url)
    return [url for _, url in sorted(best.values(), key=lambda kv: kv[1])]


def _search_cloud(short_name, release, cycles, bbox, resolution=None):
    '''
    Search CMR for one product's granules of a given release and cycle range.

    bbox=None searches the whole collection, which is how a misconfigured
    release is told apart from a tile that is simply outside coverage.
    '''
    # find_ATL11_granules is the same earthaccess.login(strategy='netrc') +
    # search_data() path read_ATL11_at uses; short_name is a parameter, so it
    # is not ATL11-specific despite the name.  It always passes bounding_box
    # through, so the unfiltered search calls earthaccess itself.
    if bbox is not None:
        from pointCollection.scripts.query_ATL11_cloud import find_ATL11_granules
        granules = find_ATL11_granules(bbox, short_name=short_name, version=release)
    else:
        import earthaccess
        earthaccess.login(strategy='netrc')
        granules = earthaccess.search_data(short_name=short_name, version=release)

    urls = []
    for granule in granules:
        try:
            url = granule.data_links(access='direct')[0]
        except (IndexError, KeyError):
            continue
        m = GRANULE_RE.match(os.path.basename(url))
        if m is None:
            continue
        if m['cycles'] != cycles or m['release'] != release:
            continue
        if resolution is not None and m['resolution'] != resolution:
            continue
        urls.append(url)
    return _latest_revisions(urls)


def find_previous_product_files(previous_product, bbox=None, earthaccess=False,
                                verbose=False):
    '''
    Resolve --previous_product into lists of ATL14 and ATL15 files to read.

    LOCAL MODE (earthaccess=False) is what it has always been: each entry of
    previous_product is a directory, globbed for ATL14_*.nc and ATL15_*1km_*.nc.
    A URI is rejected rather than globbed -- glob.glob() over an s3:// path
    returns [] with no error, which used to disable the three-sigma pre-edit
    silently on every tile (Q27 W1).

    CLOUD MODE (earthaccess=True) searches CMR for the release and cycle range
    named by previous_product, restricted to the tile's bounding box, so only
    the sectors that actually intersect the tile come back (Q27 W3).

    inputs:
        previous_product (list of str): directories, or '<release>_<cycles>'
                         specs in cloud mode
        bbox     (tuple): (lon_min, lat_min, lon_max, lat_max) for the tile;
                          cloud mode only
        earthaccess (bool): search NASA Earthdata Cloud instead of globbing
        verbose  (bool): report what was found
    outputs:
        (ATL14_files, ATL15_files): lists of paths (local) or s3:// URLs (cloud)

    raises RuntimeError if the search finds no files at all anywhere, which is a
    misconfiguration rather than an uncovered tile -- see the module docstring.
    '''
    if isinstance(previous_product, str):
        previous_product = [previous_product]

    if not earthaccess:
        remote = [d for d in previous_product if pc.io_utils.is_remote_path(d)]
        if remote:
            raise ValueError(
                f"--previous_product entries {remote} are URIs, but "
                "--previous_product_earthaccess is not set.  glob() cannot list a URI: "
                "it returns nothing, and the three-sigma pre-edit would be skipped "
                "silently on every tile.  Set --previous_product_earthaccess and give "
                "the release instead (e.g. --previous_product=005_0329).")
        ATL14_files, ATL15_files = [], []
        for directory in previous_product:
            ATL14_files += sorted(glob.glob(os.path.join(directory, 'ATL14_*.nc')))
            ATL15_files += sorted(glob.glob(os.path.join(directory, 'ATL15_*1km_*.nc')))
        if not ATL14_files and not ATL15_files:
            raise RuntimeError(
                "no ATL14_*.nc or ATL15_*1km_*.nc files under any of "
                f"{list(previous_product)}.  The three-sigma pre-edit would be skipped "
                "on every tile; drop --previous_product if that is what you want.")
        if verbose:
            print(f'\tfind_previous_product_files: {len(ATL14_files)} ATL14, '
                  f'{len(ATL15_files)} ATL15 files')
        return ATL14_files, ATL15_files

    releases = {parse_release_spec(spec) for spec in previous_product}
    if len(releases) > 1:
        raise ValueError('--previous_product may name only one release in cloud mode, '
                         f'got {sorted(releases)}')
    release, cycles = releases.pop()

    ATL14_files = _search_cloud('ATL14', release, cycles, bbox)
    ATL15_files = _search_cloud('ATL15', release, cycles, bbox,
                                resolution=ATL15_RESOLUTION)
    if verbose:
        print(f'\tfind_previous_product_files: ATL14 {[os.path.basename(f) for f in ATL14_files]}')
        print(f'\tfind_previous_product_files: ATL15 {[os.path.basename(f) for f in ATL15_files]}')

    if not ATL14_files and not ATL15_files and bbox is not None:
        # Nothing intersects this tile.  That is normal at the edge of the
        # previous product's domain, and a fatal misconfiguration if the
        # release does not exist at all -- one unfiltered search tells them
        # apart, and only runs on the tiles that found nothing.
        anywhere = (_search_cloud('ATL14', release, cycles, None)
                    or _search_cloud('ATL15', release, cycles, None,
                                     resolution=ATL15_RESOLUTION))
        if not anywhere:
            raise RuntimeError(
                f"no ATL14/ATL15 granules anywhere for release {release}, cycles {cycles} "
                f"(--previous_product={release}_{cycles}).  Check the release and cycle "
                "range: CMR mixes cycle ranges under one short_name, so a wrong cycles "
                "string matches nothing rather than falling back to another.")
    return ATL14_files, ATL15_files
