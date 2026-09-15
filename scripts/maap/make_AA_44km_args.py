#!/usr/bin/env python3
"""
Derive the Antarctic SOUTH-half args file from the north-half one.

The south half solves a smaller footprint -- -W 44000 against 60000, same
40 km spacing -- and writes to its own region directory, which is how the two
halves stay separable on a bucket where a tile is named for its center alone.

WHY THIS IS NOT A default_args OVERRIDES FILE, which was tried first and
failed: setup_ATL1415_region.py derives BOTH the region directory and the args
file name from --region, so getting input_args_AA_44km.txt out of it means
passing --region=AA_44km -- and --region is not a label, it is a behavioural
switch.  ATL11_to_ATL15.py:595-596 loads the gridded mask only for
region in ['AA', 'GL']; anything else falls through both branches with
mask_data left None, and the solve dies at line 619 with

    AttributeError: 'NoneType' object has no attribute 'z'

All four 44 km jobs of the 2026-09-08 cost queue failed exactly that way, while
their 60 km counterparts -- including the SAME tile center, E420_N20, in the
overlap band -- succeeded.  So the south half keeps --region=AA and differs
only in the two lines that describe its geometry and its output location.

Usage:  make_AA_44km_args.py <input_args_AA.txt> <output_args_AA_44km.txt>
"""
import os
import sys

WIDTH_44KM = '44000'


def derive(src_lines, out_base_directory):
    """Rewrite -W and -b; leave every other line, --region included, alone."""
    out, seen_W, seen_b = [], False, False
    for line in src_lines:
        stripped = line.strip()
        if stripped.startswith('-W='):
            out.append(f'-W={WIDTH_44KM}\n'); seen_W = True
        elif stripped.startswith('-b='):
            out.append(f'-b={out_base_directory}\n'); seen_b = True
        else:
            out.append(line)
    if not seen_W:
        out.append(f'-W={WIDTH_44KM}\n')
    if not seen_b:
        out.append(f'-b={out_base_directory}\n')
    return out


def main(argv):
    if len(argv) != 3:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    src, dst = argv[1], argv[2]
    base_directory = os.path.dirname(os.path.abspath(dst))
    os.makedirs(base_directory, exist_ok=True)
    with open(src) as fh:
        lines = fh.readlines()

    if not any(l.strip() == '--region=AA' for l in lines):
        print(f'ERROR: {src} does not set --region=AA.  The south half must keep '
              'the north half\'s region, or ATL11_to_ATL15 will not load the '
              'gridded mask (see this file\'s docstring).', file=sys.stderr)
        return 1

    with open(dst, 'w') as fh:
        fh.writelines(derive(lines, base_directory))
    print(f'wrote {dst}')
    print(f'  -W={WIDTH_44KM}, -b={base_directory}, --region=AA unchanged')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
