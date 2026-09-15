#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tide elevations from either a local tide model or one of the zarr stores pyTMD
publishes in the cloud.

pyTMD's compute.tide_elevations() reaches its model through
pyTMD.io.model.open_dataset(), which has no zarr branch, so it cannot be pointed
at a zarr store.  What it does *after* opening the model is short and generic --
coords_as -> interp -> predict + infer, all on an ordinary xarray Dataset -- so
the cloud path here substitutes the dataset and reuses that sequence rather than
reimplementing any tide computation.  See
https://pytmd.readthedocs.io/en/latest/user_guide/Cloud-Access.html

Why this is worth the extra code path, measured on a 50k-point 60 km Ross tile:
the published zarr stores are chunked, so a tile reads 1.6 MiB (CATS2008-v2023)
or 6.3 MiB (Gr1km-v2).  Reaching the same model as the hosted netCDF instead
costs 582 MiB, because every variable in CATS2008_v2023.nc is stored contiguously
and uncompressed and no spatial subset can avoid reading a variable whole.  The
two agree to max 7.2e-05 m over those points.
"""
import os

import numpy as np


# The stores pyTMD publishes are named '<model>.zarr' under the bucket, which is
# the convention its own documentation uses.
ZARR_SUFFIX = '.zarr'


def is_remote_tide_directory(tide_directory):
    """
    True if `tide_directory` names a cloud location holding zarr stores rather
    than a directory of local model files.
    """
    return isinstance(tide_directory, str) and tide_directory.startswith('s3://')


def _store_path(tide_directory, tide_model):
    """s3 path (no scheme) of the zarr store for `tide_model`"""
    return f"{tide_directory[len('s3://'):].rstrip('/')}/{tide_model}{ZARR_SUFFIX}"


def open_tide_dataset(tide_directory, tide_model, group='z', anon=None):
    """
    Open a cloud tide model as an xarray Dataset, plus the pyTMD model object
    that carries its nodal corrections.

    Parameters
    ----------
    tide_directory : str
        s3:// prefix holding '<model>.zarr' stores, e.g. 's3://pytmd'.
    tide_model : str
        pyTMD model name, e.g. 'CATS2008-v2023' or 'Gr1km-v2'.
    group : str, default 'z'
        Model group; 'z' is elevations.
    anon : bool or None
        Whether to read the store anonymously.  Defaults to True, overridable
        with ATL1415_TIDE_ANON=0.  This matters: pyTMD's bucket is publicly
        readable but REJECTS signed requests from another account, so a DPS
        worker -- which has its own AWS credentials -- must ask for the store
        anonymously or it gets AccessDenied.

    Returns
    -------
    ds : xarray.Dataset
    m : pyTMD.io.model
        Built from the database only; it reads no model files.
    """
    import s3fs
    import zarr
    import xarray as xr
    import pyTMD

    if anon is None:
        anon = os.environ.get('ATL1415_TIDE_ANON', '1') not in ('0', 'false', 'False')
    fs = s3fs.S3FileSystem(anon=anon)
    store = zarr.storage.FsspecStore(fs, path=_store_path(tide_directory, tide_model))
    ds = xr.open_zarr(store, group=group, zarr_format=3)
    # verify=False: nothing local to verify, and the model object is wanted only
    # for its metadata (corrections), not to locate files
    m = pyTMD.io.model(verify=False).from_database(tide_model)
    return ds, m


def tide_elevations(x, y, delta_time, tide_directory=None, tide_model=None,
                    crs=4326, epoch=(2018, 1, 1, 0, 0, 0), standard='GPS',
                    method='linear', extrapolate=False, cutoff=10.0,
                    type='drift', infer_minor=True, chunks=None, **kwargs):
    """
    Tide elevations at (x, y, delta_time), from a local model directory or a
    cloud zarr store.

    A local `tide_directory` is passed straight to pyTMD.compute.tide_elevations
    and behaves exactly as before, so local and SLURM runs are unaffected.  An
    s3:// `tide_directory` takes the zarr path.

    Parameters mirror pyTMD.compute.tide_elevations.

    Returns
    -------
    np.ndarray
        Tide elevations, in meters.
    """
    import pyTMD

    if not is_remote_tide_directory(tide_directory):
        return np.array(pyTMD.compute.tide_elevations(
            x, y, delta_time, directory=tide_directory, model=tide_model,
            infer_minor=infer_minor, chunks=chunks, crs=crs, type=type,
            epoch=epoch, standard=standard, method=method,
            extrapolate=extrapolate, cutoff=cutoff, **kwargs))

    ds, m = open_tide_dataset(tide_directory, tide_model)

    # the same sequence compute.tide_elevations runs once its model is open
    import timescale
    ts = timescale.from_deltatime(delta_time, epoch=epoch, standard=standard)
    X, Y = ds.tmd.coords_as(x, y, type=type, time=delta_time, crs=crs)
    local = ds.tmd.interp(X, Y, method=method,
                          extrapolate=extrapolate, cutoff=cutoff)
    corrections = kwargs.get('corrections') or m.corrections
    tide = np.array(local.tmd.predict(ts.tide, deltat=ts.tt_ut1,
                                      corrections=corrections))
    if infer_minor:
        tide += np.array(local.tmd.infer(ts.tide, deltat=ts.tt_ut1,
                                         corrections=corrections))
    return tide
