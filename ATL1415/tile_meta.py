#!/usr/bin/env python3
"""
Per-tile metadata for the netCDF writers, read from each prelim tile ONCE.

ATL14_write2nc and ATL15_write2nc need two things from every prelim tile: the
tile_stats values (make_tile_stats_group) and the ATL11 lineage inputs
(ATL1415_attrs_meta.set_lineage).  Read in place on S3 that is ~1.3 s a tile,
and ATL15 used to do it twice for each of its four files: >2 h for Greenland.

So the 200 km mosaic jobs, which already work through the tiles a block at a
time, save the records of the tiles they own (`write`, below) as one small
JSON file each, and the writers read those (--tile_meta_dir) instead of the
tiles.  Without --tile_meta_dir the writers read the tiles, as on discover --
once per run either way (tile_records caches).

A tile belongs to the 200 km tile whose square holds its center
(mosaic_groups.centers_200km), so the files partition the tiles.  The writers
check that against a listing of tiles_dir and stop, naming the tiles, on any
gap, overlap, or tile rewritten since its record was made.

usage:
    python -m ATL1415.tile_meta write <tiles_dir> <out.json> [--center X Y] [--workers N]
"""
import argparse
import json
import os
import re
import sys

import numpy as np

META_VERSION = 1


def tile_xy_km(name):
    """
    (x, y) in km from a tile name, parsed as make_tile_stats_group always has
    (loosely: E<x>_..., N<y>.<ext>); None where the x part does not parse, the
    files the stats reader skips.
    """
    try:
        x = int(re.match(r'^.*E(.*)\_.*$', name).group(1))
    except Exception:
        return None
    return x, int(re.match(r'^.*N(.*)\..*$', name).group(1))


def _fs_for(path):
    import pointCollection as pc
    return pc.io_utils.get_s3fs(daac=None) if pc.io_utils.is_remote_path(path) else None


def parse_input_files(raw):
    """meta/input_files as set_lineage has always parsed it: a list of names."""
    inputs = str(raw)
    if inputs[:1] == 'b':
        inputs = inputs[1:]
    inputs = inputs.replace("'", '')
    # a tile that read no ATL11 (a matched tile) has input_files == ''
    return list(filter(None, inputs.split(',')))


def read_tile_stats(h5):
    """The tile_stats values of one open tile (make_tile_stats_group's fields)."""
    return {
        'N_data': int(np.sum(h5['data']['three_sigma_edit'][:])),
        'RMS_data': float(h5['RMS']['data'][()]),
        'RMS_bias': float(np.sqrt(np.mean((h5['bias']['val'][:] / h5['bias']['expected'][:])**2))),
        'N_bias': int(len(h5['bias']['val'][:])),
        'RMS_d2z0dx2': float(h5['RMS']['grad2_z0'][()]),
        'RMS_d2zdt2': float(h5['RMS']['d2z_dt2'][()]),
        'RMS_d2zdx2dt': float(h5['RMS']['grad2_dzdt'][()]),
        'sigma_xx0': float(h5['E_RMS']['d2z0_dx2'][()]),
        'sigma_tt': float(h5['E_RMS']['d2z_dt2'][()]),
        'sigma_xxt': float(h5['E_RMS']['d3z_dx2dt'][()]),
    }


def read_tile_lineage(h5):
    """
    {'input_files': [...], 'stored': {granule: {attr: text}}} for one open tile.
    Raises KeyError for a tile without meta/input_files, as set_lineage expects.
    """
    from ATL1415.ATL1415_attrs_meta import lineage_from_tile
    files = parse_input_files(h5['/meta/'].attrs['input_files'])
    stored = {}
    for granule in files:
        if granule not in stored:
            this = lineage_from_tile(h5, granule)
            if this:
                stored[granule] = this
    return {'input_files': files, 'stored': stored}


def tile_record(path):
    """
    Everything the writers need from one tile, from one open.

    'stats' is None for a file whose name is not E<x>_N<y>.h5 (the stats reader
    skips those); a stats read that fails raises, as it always has.  'lineage'
    is None, with 'lineage_error', for a tile set_lineage would skip
    (unreadable, or without meta/input_files).
    """
    from ATL1415.paths import open_tile
    rec = {'name': os.path.basename(path), 'stats': None, 'lineage': None, 'lineage_error': None}
    want_stats = tile_xy_km(rec['name']) is not None
    try:
        with open_tile(path) as h5:
            if want_stats:
                rec['stats'] = read_tile_stats(h5)
            try:
                rec['lineage'] = read_tile_lineage(h5)
            except (OSError, KeyError) as e:
                rec['lineage_error'] = f'{type(e).__name__}: {e}'
    except OSError as e:
        if want_stats:
            raise
        rec['lineage_error'] = f'{type(e).__name__}: {e}'
    return rec


def tile_listing(tiles_dir, pattern='*.h5'):
    """
    {name: stamp} for the tiles in tiles_dir; a stamp is 'size|modified',
    which changes when the tile is rewritten.  Local or remote.
    """
    from ATL1415.paths import list_tiles
    paths = list_tiles(tiles_dir, pattern)
    fs = _fs_for(tiles_dir)
    out = {}
    if fs is None:
        for p in paths:
            st = os.stat(p)
            out[os.path.basename(p)] = (p, f'{st.st_size}|{st.st_mtime_ns}')
        return out
    details = {os.path.basename(d['name']): d
               for d in fs.ls(tiles_dir.rstrip('/'), detail=True)}
    for p in paths:
        d = details.get(os.path.basename(p), {})
        mod = d.get('LastModified', d.get('mtime', d.get('created', '')))
        out[os.path.basename(p)] = (p, f"{d.get('size', '')}|{mod}")
    return out


def owned_by(name, center, tile_W=200e3):
    """True if tile `name` belongs to the 200 km tile centered at `center` (m)."""
    from ATL1415.mosaic_groups import centers_200km
    c = centers_200km([name], tile_W=tile_W)
    return bool(c) and np.allclose(c[0], center)


def write_meta(tiles_dir, out_path, center=None, workers=8):
    """
    Read the records of the tiles in tiles_dir (those owned by the 200 km tile
    at `center`, or all) and write them, with their listing stamps, to out_path
    (local or s3://).  Returns the number of tiles.
    """
    listing = tile_listing(tiles_dir)
    names = sorted(n for n in listing if center is None or owned_by(n, center))
    paths = [listing[n][0] for n in names]
    if workers > 1 and len(paths) > 1:
        # PROCESSES, not threads: h5py holds one lock around every HDF5 call,
        # including the S3 reads it makes through the file object, so threads
        # read one tile at a time (measured: 8 threads 0.9 s/tile, 8 processes
        # 0.26 s/tile).  'spawn': an s3fs session cannot cross a fork.
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(min(workers, len(paths)), mp_context=mp.get_context('spawn')) as ex:
            records = list(ex.map(tile_record, paths, chunksize=2))
    else:
        records = [tile_record(p) for p in paths]
    for rec in records:
        rec['stamp'] = listing[rec['name']][1]
    doc = {'version': META_VERSION, 'tiles_dir': tiles_dir,
           'center': None if center is None else [float(c) for c in center],
           'tiles': records}
    text = json.dumps(doc, indent=0)
    fs = _fs_for(out_path)
    if fs is None:
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        with open(out_path, 'w') as fh:
            fh.write(text)
    else:
        fs.pipe(out_path, text.encode('utf-8'))
    return len(records)


def load_meta(meta_dir, tiles_dir):
    """
    [(tile path, record)] for the tiles in tiles_dir, in list_tiles order, from
    the JSON files in meta_dir.  Stops on any tile without exactly one current
    record.
    """
    listing = tile_listing(tiles_dir)
    fs = _fs_for(meta_dir)
    if fs is None:
        import glob
        files = sorted(glob.glob(os.path.join(meta_dir, '*.json')))
        read = lambda f: open(f).read()
    else:
        import pointCollection as pc
        files = pc.io_utils.glob_remote(meta_dir.rstrip('/') + '/*.json', fs=fs)
        read = lambda f: fs.cat(f).decode('utf-8')
    if not files:
        raise RuntimeError(f'tile_meta: no *.json in {meta_dir}: run the 200 km jobs '
                           '(they write it), or omit --tile_meta_dir to read the tiles')
    records, where = {}, {}
    for f in files:
        doc = json.loads(read(f))
        if doc.get('version') != META_VERSION:
            raise RuntimeError(f'tile_meta: {f} is version {doc.get("version")}, need {META_VERSION}')
        for rec in doc['tiles']:
            if rec['name'] in records:
                raise RuntimeError(f'tile_meta: {rec["name"]} is in both {where[rec["name"]]} and {f}')
            records[rec['name']], where[rec['name']] = rec, f
    problems = []
    missing = sorted(set(listing) - set(records))
    extra = sorted(set(records) - set(listing))
    stale = sorted(n for n in set(listing) & set(records)
                   if records[n].get('stamp') != listing[n][1])
    for label, names in (('no record', missing), ('record but no tile', extra),
                         ('rewritten since its record was made', stale)):
        if names:
            problems.append(f'{len(names)} with {label}: {", ".join(names[:10])}'
                            f'{" ..." if len(names) > 10 else ""}')
    if problems:
        raise RuntimeError(f'tile_meta: {meta_dir} does not match the tiles in {tiles_dir}: '
                           + '; '.join(problems)
                           + '.  Rerun the 200 km jobs, or omit --tile_meta_dir.')
    return [(listing[n][0], records[n]) for n in sorted(listing, key=lambda n: listing[n][0])]


_CACHE = {}


def tile_records(args):
    """
    [(tile path, record)] for args.tiles_dir, in list_tiles order: from
    args.tile_meta_dir if it is set, else from the tiles themselves.  Read once
    per process -- ATL15's four files share it.
    """
    meta_dir = getattr(args, 'tile_meta_dir', None)
    key = (args.tiles_dir, meta_dir)
    if key not in _CACHE:
        if meta_dir:
            _CACHE[key] = load_meta(meta_dir, args.tiles_dir)
        else:
            from ATL1415.paths import list_tiles
            _CACHE[key] = [(p, tile_record(p)) for p in list_tiles(args.tiles_dir)]
    return _CACHE[key]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = parser.add_subparsers(dest='cmd', required=True)
    w = sub.add_parser('write', help="save the records of a 200 km tile's prelim tiles (or all)")
    w.add_argument('tiles_dir')
    w.add_argument('out')
    w.add_argument('--center', type=float, nargs=2, metavar=('X', 'Y'))
    w.add_argument('--workers', type=int, default=4, help='processes reading tiles (default 4)')
    args = parser.parse_args(argv)
    n = write_meta(args.tiles_dir, args.out, center=args.center, workers=args.workers)
    print(f'tile_meta: {n} tile records -> {args.out}')
    if n == 0:
        print(f'tile_meta: ERROR: no tiles in {args.tiles_dir}'
              + (f' for center {args.center}' if args.center else ''), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
