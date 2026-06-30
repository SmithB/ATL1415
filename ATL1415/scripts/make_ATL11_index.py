#!/usr/bin/env python3

import os
import sys
import glob
import argparse

HEMI_CONFIG = {
    'north': {'regions': ['03', '04', '05'], 'hemisphere':  1, 'prefix': 'Arctic'},
    'south': {'regions': ['10', '11', '12'], 'hemisphere': -1, 'prefix': 'Antarctic'},
}

def parse_atl11_release(atl11_release):
    """
    Parse ATL11 release string into (release, dir_suffix).

    'ATL11_007_0329_v02' -> ('007', '007_cycle_03_29_v02')
    '005_0314'           -> ('005', '005_cycle_03_14')
    """
    ver = atl11_release[len('ATL11_'):] if atl11_release.startswith('ATL11_') else atl11_release
    parts = ver.split('_')
    release = parts[0]
    cycles = parts[1]
    cycle_start, cycle_end = cycles[:2], cycles[2:]
    version_tag = parts[2] if len(parts) > 2 else ''
    suffix = f"{release}_cycle_{cycle_start}_{cycle_end}"
    if version_tag:
        suffix += f"_{version_tag}"
    return release, suffix


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        description='Build ATL11 index directory tree, symlinks, and index queue files.',
        fromfile_prefix_chars='@')
    parser.add_argument('--ATL14_root', type=str, required=True,
                        help='ATL14 processing root (from discover.txt)')
    parser.add_argument('--ATL11_root', type=str, required=True,
                        help='ATL11 source data root (from discover.txt)')
    parser.add_argument('--ATL11_release', type=str, required=True,
                        help='ATL11 release string, e.g. ATL11_007_0329_v02 (from latest_release.txt)')
    parser.add_argument('--ATL11_north', type=str, default=None,
                        help='explicit north source dir override; "None" to skip hemisphere')
    parser.add_argument('--ATL11_south', type=str, default=None,
                        help='explicit south source dir override; "None" to skip hemisphere')
    args = parser.parse_args(argv)

    release, suffix = parse_atl11_release(args.ATL11_release)

    atl11_release_dir = args.ATL11_release if args.ATL11_release.startswith('ATL11_') \
        else f'ATL11_{args.ATL11_release}'

    # resolve source directories, applying any overrides
    source_dirs = {}
    for hemi, cfg in HEMI_CONFIG.items():
        override = getattr(args, f'ATL11_{hemi}')
        if override == 'None':
            source_dirs[hemi] = None
        elif override is not None:
            source_dirs[hemi] = override
        else:
            source_dirs[hemi] = os.path.join(
                args.ATL11_root, f"{cfg['prefix']}_{suffix}", release)

    with open('index_queue.txt', 'w') as index_fh, \
         open('post_queue.txt', 'w') as post_fh:

        for hemi, cfg in HEMI_CONFIG.items():
            source_dir = source_dirs[hemi]
            if source_dir is None:
                print(f"Skipping {hemi}")
                continue

            print(f"{hemi}: {source_dir}")

            hemi_dir  = os.path.join(args.ATL14_root, atl11_release_dir, hemi)
            index_dir = os.path.join(hemi_dir, 'index')
            os.makedirs(index_dir, exist_ok=True)

            for region in cfg['regions']:
                pattern = os.path.join(source_dir, f'ATL11_????{region}_????_???_??.h5')
                files = sorted(glob.glob(pattern))
                if not files:
                    print(f"  WARNING: no files found matching {pattern}")
                for filepath in files:
                    base = os.path.basename(filepath)
                    link_path = os.path.join(hemi_dir, base)
                    if not os.path.islink(link_path):
                        os.symlink(filepath, link_path)
                    if os.path.isfile(os.path.join(index_dir, base)):
                        continue
                    index_fh.write(
                        f"pushd {index_dir}; index_glob.py"
                        f" --hemisphere {cfg['hemisphere']} --type ATL11"
                        f" -g ../{base} --index_file ./{base} --Relative; popd\n"
                    )

            post_fh.write(
                f"pushd {index_dir}; index_glob.py"
                f" --hemisphere {cfg['hemisphere']} --type h5_geoindex"
                f" -g 'ATL11*.h5' --index_file GeoIndex.h5 --Relative; popd\n"
            )

    print("wrote index_queue.txt and post_queue.txt")


if __name__ == '__main__':
    main()
