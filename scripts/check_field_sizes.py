#!/usr/bin/env python3
"""
Check the per-tile field-size reports of one step directory.

Each tile job writes <step_dir>/field_sizes/<tile>_report.json beside its tile
(ATL11_to_ATL15.save_field_size_report), holding the shapes of dz/dz and
dz/sigma_dz, with null for a field that was not in the file.  This reads those
reports -- never the tiles -- and checks, per docs/plan_check_field_sizes.sh:

  1 SHAPE    dz/dz equals the shape derived from the run's args file:
               nx = ny = W / dz_spacing + 1    (-W, second value of -g)
               nt = (t1 - t0) / dt + 1         (-t, third value of -g)
             -t, NOT --t_crop: the tiles span -t (IS: 32 epochs, not 30).
  2 PRELIM   dz/sigma_dz is present and equals dz/dz.
  3 MATCHED  dz/sigma_dz is null -- matched tiles carry no sigma by design;
             the uncertainties come from the prelim step.  Checks 2 and 3 are
             inverted between the steps, so the step is never guessed.
  4 PAIRING  every tile has a report and every report has a tile, matched by
             basename.  The report's "file" value is a path inside the DPS
             worker and is never used.

-W and -t are REQUIRED: with no derived shape there is nothing to check
check 1 against, and whether tiles should instead be compared with each other
is an open question for Ben (plan C1).  AA's two halves are separate step
directories with separate args files, so each gets its own derived shape.

Usage:
  check_field_sizes.py <step_dir> @<input_args.txt> [--step prelim|matched]
  check_field_sizes.py <step_dir> -W=60000 -g=100,1000,0.25 -t=2018.75,2026.5

@file is read the way the solver reads it (argparse fromfile, '@' includes
followed), with the solver's names, aliases and -g default.  Every other line
of the file is ignored; an unknown option typed on the command line itself is
an error, so a mistyped --step cannot vanish into the ignored lines.

--step defaults to the step directory's name when that is literally 'prelim'
or 'matched'; otherwise it is required.

Exit status:
  0  every check passed
  1  at least one check failed (one line per problem, then a summary)
  2  the check did not happen: no such directory, no reports, arguments
     that do not give a shape, or a step that cannot be determined
"""
import argparse
import glob
import json
import os
import sys

REPORT_SUFFIX = '_report.json'
REPORT_KEYS = {'file', 'dz/dz', 'dz/sigma_dz'}
STEPS = ('prelim', 'matched')


class CannotCheck(Exception):
    """The check could not run at all (exit 2)."""


def _count(span, spacing, what):
    """span/spacing + 1, which must come out a whole number."""
    n = span / spacing
    if spacing <= 0 or abs(n - round(n)) > 1e-6:
        raise CannotCheck(f'{what}: {span} / {spacing} is not a whole number of cells')
    return int(round(n)) + 1


def _spacing(text):
    """One -g entry, which may be a fraction: monthly runs pass dt as 1/12.

    make_mosaic_jobs.py and ATL15_write2nc.py both read -g this way.
    """
    if '/' in text:
        num, den = [float(v) for v in text.split('/')]
        return num / den
    return float(text)


def expected_shape(Width, time_span, grid_spacing):
    """
    Derive the dz shape [nx, ny, nt] from the solver's -W, -t and -g values.

    Returns (shape, derivation text).
    """
    if Width is None or time_span is None:
        raise CannotCheck('-W and -t are both required (pass @<input_args file>)')
    try:
        t0, t1 = [float(t) for t in time_span.split(',')]
        _, dz_spacing, dt = [_spacing(g) for g in grid_spacing.split(',')]
    except (ValueError, ZeroDivisionError):
        raise CannotCheck(f'cannot read -t={time_span} and -g={grid_spacing}')
    nxy = _count(Width, dz_spacing, '-W / dz spacing')
    nt = _count(t1 - t0, dt, '-t span / dt')
    dt_text = grid_spacing.split(',')[2]
    derivation = (f'-W={Width:g} / {dz_spacing:g} + 1 = {nxy};  '
                  f'-t={t0:g},{t1:g}: ({t1:g} - {t0:g}) / {dt_text} + 1 = {nt}')
    return [nxy, nxy, nt], derivation


def infer_step(step_dir, step):
    if step is not None:
        return step
    name = os.path.basename(os.path.normpath(step_dir))
    if name in STEPS:
        return name
    raise CannotCheck(f'cannot tell the step from the directory name {name!r}; '
                      'pass --step prelim or --step matched')


def check_report(report, step, shape):
    """Problems with one parsed report, as a list of strings."""
    problems = []
    missing = REPORT_KEYS - set(report)
    extra = set(report) - REPORT_KEYS
    if missing:
        problems.append(f'missing keys {sorted(missing)}')
    if extra:
        problems.append(f'unexpected keys {sorted(extra)}')
    dz = report.get('dz/dz')
    sigma = report.get('dz/sigma_dz')
    if 'dz/dz' in report and dz != shape:
        problems.append(f'dz/dz is {dz}, expected {shape}')
    if 'dz/sigma_dz' in report:
        if step == 'prelim' and sigma != dz:
            problems.append(f'prelim dz/sigma_dz is {sigma}, expected {dz} (same as dz/dz)')
        if step == 'matched' and sigma is not None:
            problems.append(f'matched dz/sigma_dz is {sigma}, expected null (no sigma in matched)')
    return problems


def check_step_dir(step_dir, shape, derivation, step=None):
    """
    Run all four checks against an expected dz shape.
    Returns (problems, summary lines); raises CannotCheck.
    """
    if not os.path.isdir(step_dir):
        raise CannotCheck(f'step directory {step_dir} not found')
    step = infer_step(step_dir, step)

    tiles = {os.path.basename(f)[:-len('.h5')]
             for f in glob.glob(os.path.join(step_dir, '*.h5'))}
    report_files = sorted(glob.glob(os.path.join(step_dir, 'field_sizes', '*' + REPORT_SUFFIX)))
    if not report_files:
        raise CannotCheck(f'no reports in {os.path.join(step_dir, "field_sizes")} '
                          f'({len(tiles)} tiles in {step_dir})')

    problems = []
    reports = set()
    bad_tiles = set()
    for report_file in report_files:
        name = os.path.basename(report_file)[:-len(REPORT_SUFFIX)]
        reports.add(name)
        try:
            with open(report_file) as fh:
                report = json.load(fh)
            if not isinstance(report, dict):
                raise ValueError('not a JSON object')
        except (OSError, ValueError) as e:
            problems.append(f'{name}: cannot read {report_file}: {e}')
            bad_tiles.add(name)
            continue
        for problem in check_report(report, step, shape):
            problems.append(f'{name}: {problem}')
            bad_tiles.add(name)
    for name in sorted(reports - tiles):
        problems.append(f'{name}: report has no tile')
        bad_tiles.add(name)
    for name in sorted(tiles - reports):
        problems.append(f'{name}: tile has no report')
        bad_tiles.add(name)

    names = tiles | reports
    summary = [f'step {step}: {step_dir}',
               f'expected dz/dz {shape}  ({derivation})',
               f'expected dz/sigma_dz ' + ('== dz/dz' if step == 'prelim' else 'null'),
               f'{len(reports)} reports, {len(tiles)} tiles, '
               f'{len(names) - len(bad_tiles)} of {len(names)} passed, {len(problems)} problems']
    return problems, summary


def parse_args(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser(
        description='Check the per-tile field-size reports of a prelim or matched directory.',
        fromfile_prefix_chars='@', allow_abbrev=False)
    parser.add_argument('step_dir')
    parser.add_argument('--step', choices=STEPS,
                        help="default: the step directory's name, if prelim or matched")
    # the solver's names, aliases and -g default (ATL11_to_ATL15.parse_args)
    parser.add_argument('--Width', '-W', type=float, help='tile width (m)')
    parser.add_argument('--time_span', '-t', type=str, help='first year,last year')
    parser.add_argument('--grid_spacing', '-g', type=str, default='250.,4000.,1.',
                        help='z0 spacing,dz spacing,dt')
    args, unknown = parser.parse_known_args(argv)
    # the rest of an @file is expected; anything typed directly is a mistake
    typed = [arg for arg in unknown if arg in argv]
    if typed:
        parser.error(f'unrecognized arguments: {" ".join(typed)}')
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        shape, derivation = expected_shape(args.Width, args.time_span, args.grid_spacing)
        problems, summary = check_step_dir(args.step_dir, shape, derivation, args.step)
    except CannotCheck as e:
        print(f'check_field_sizes.py: CHECK NOT DONE: {e}', file=sys.stderr)
        return 2
    for problem in problems:
        print(f'PROBLEM {problem}')
    for line in summary:
        print(line)
    print('FAILED' if problems else 'OK')
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
