#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Path handling for arguments that may name either a local file or a cloud object.

A DPS worker has no workspace mount, so the masks and ancillary grids a tile
solve reads are named by s3:// URI rather than by path.  The helpers here keep
those URIs intact where the usual local-path normalization would corrupt them.
"""
import os

import pointCollection as pc


def path_or_uri(path):
    """
    Normalize a path argument that may be a local path or a URI.

    Arguments naming input files are normally passed through
    os.path.abspath(os.path.expanduser(...)), so that a run does not depend on
    the working directory it was started from.  That is actively wrong for a
    URI: abspath() would rewrite 's3://bucket/key' as '<cwd>/s3:/bucket/key',
    collapsing the '//' and prefixing the cwd.  A URI is already absolute, so
    it is returned unchanged.

    Parameters
    ----------
    path : str or None

    Returns
    -------
    str or None
        `path` unchanged if it is a URI (or None), otherwise the absolute,
        user-expanded local path.
    """
    if path is None:
        return None
    if pc.io_utils.is_remote_path(path):
        return path
    return os.path.abspath(os.path.expanduser(path))


def join_path_or_uri(root, path):
    """
    Join `path` onto `root` unless `path` already names an absolute location.

    os.path.join() alone is not enough here: it treats 's3://bucket/key' as a
    relative path (it does not start with '/') and would produce
    '<root>/s3://bucket/key'.  A URI is absolute and must win outright, the
    same way an absolute local path does.

    Parameters
    ----------
    root : str
        Directory or prefix to join onto.  May itself be a local path or a URI.
    path : str
        Relative path, absolute path, or URI.

    Returns
    -------
    str
    """
    if pc.io_utils.is_remote_path(path) or os.path.isabs(path):
        return path
    if pc.io_utils.is_remote_path(root):
        # posix join: a URI's separator is '/' whatever the local os.sep is
        return root.rstrip('/') + '/' + path
    return os.path.join(root, path)


def exists(path):
    """
    True if `path` names an existing local file or directory, or a URI whose
    object exists.

    A local path is checked with os.path.exists(); a URI is checked against the
    bucket, which costs a network round trip, so call this at setup time rather
    than per tile.
    """
    if pc.io_utils.is_remote_path(path):
        fs = pc.io_utils.get_s3fs(daac=None)
        return fs.exists(path)
    return os.path.exists(path)


def is_tile_manifest(tiles_dir):
    """True if tiles_dir names a manifest (a .txt list of tiles), not a directory"""
    return tiles_dir.rstrip('/').endswith('.txt')


def read_tile_manifest(manifest):
    """
    The tiles a manifest lists, SORTED: one local path or URI per line.

    An Antarctic sector's tiles live under two prefixes on S3 (the 60 km and
    44 km halves) and S3 has no symlinks to gather them in one directory as
    setup_AA_sectors.py does on discover, so a sector names its tiles in a
    manifest instead (plan_AA_mosaic_dps.sh M0).  The writers key tiles by
    file name, so a name listed twice is an error, not a silent choice.
    """
    if pc.io_utils.is_remote_path(manifest):
        text = pc.io_utils.get_s3fs(daac=None).cat(manifest).decode('utf-8')
    else:
        with open(manifest) as fh:
            text = fh.read()
    tiles = sorted(line.strip() for line in text.splitlines() if line.strip())
    seen, dup = set(), set()
    for name in (os.path.basename(t) for t in tiles):
        (dup if name in seen else seen).add(name)
    if dup:
        raise ValueError(f'{manifest} lists {len(dup)} tile name(s) more than once: {sorted(dup)[:10]}')
    return tiles


def list_tiles(tiles_dir, pattern='*.h5'):
    """
    The tile files in a local directory, under a URI prefix, or listed in a
    manifest (read_tile_manifest), SORTED.

    On MAAP the netCDF writers read the prelim tiles in place on the bucket
    (docs/plan_dps_mosaic.sh D2): nothing is copied to /home, whose quota is
    150 GB.  Sorted in all cases, so what is read first -- and so which tile
    set_lineage names for a granule -- does not depend on the filesystem.
    """
    if is_tile_manifest(tiles_dir):
        import fnmatch
        return [t for t in read_tile_manifest(tiles_dir)
                if fnmatch.fnmatch(os.path.basename(t), pattern)]
    if pc.io_utils.is_remote_path(tiles_dir):
        return pc.io_utils.glob_remote(tiles_dir.rstrip('/') + '/' + pattern,
                                       fs=pc.io_utils.get_s3fs(daac=None))
    import glob
    return sorted(glob.glob(os.path.join(tiles_dir, pattern)))


def open_tile(path):
    """
    Open a tile for reading, from a local path or a URI; use as
    `with open_tile(path) as h5f:`.

    A remote tile is read with pointCollection's small block size and block
    cache: these readers take a few datasets from each tile, and s3fs's
    default block (50 MiB) would fetch most of the file for them.  The
    session is the default AWS credential chain, since the tiles are in our
    bucket, not a DAAC's.  Leaving the `with` closes the remote file object as
    well as the HDF5 file (h5py does not close a file object it was given).
    """
    import contextlib
    import h5py
    if not pc.io_utils.is_remote_path(path):
        return h5py.File(path, 'r')

    @contextlib.contextmanager
    def remote():
        with pc.io_utils.open_remote(path, fs=pc.io_utils.get_s3fs(daac=None),
                                     block_size=pc.io_utils.DEFAULT_REMOTE_BLOCK_SIZE,
                                     daac=None) as fileobj, h5py.File(fileobj, 'r') as h5f:
            yield h5f
    return remote()
