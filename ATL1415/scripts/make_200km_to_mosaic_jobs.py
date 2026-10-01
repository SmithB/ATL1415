import os
import sys
import argparse
import ATL1415
import subprocess


def mosaic_commands(base, lags, in_base=None, crop='', skip_z0=False, workers=1):
    """
    The make_mosaic.py commands that join a region's 200 km tiles into its
    mosaics, by group, in task order.

    Parameters
    ----------
    base : str
        region directory the mosaics are written into (local).
    lags : iterable of int
        dzdt lags.
    in_base : str, optional
        where the 200 km tiles are read from, <in_base>/200km_tiles/<group>/,
        if not base -- an s3:// prefix on DPS (docs/plan_dps_mosaic.sh D3a-2).
    crop : str, optional
        '-c <bounds>' for every command, or ''.
    skip_z0 : bool, optional
        leave out z0 even if its 200 km tiles exist.
    workers : int, optional
        make_mosaic.py -j: read each command's tiles in this many processes.

    Returns
    -------
    dict
        group name (z0, dz, dzdt_lag<N>, avg_dz_<S>m, avg_dzdt_<S>m_lag<N>)
        -> list of commands, run in order (the first creates the file, -R).
    """
    from ATL1415.paths import exists, join_path_or_uri
    in_base = base if in_base is None else in_base
    j = f' -j {workers}' if workers and workers > 1 else ''

    def lines(group, out_file, fields):
        d = join_path_or_uri(in_base, f'200km_tiles/{group}')
        return [f"make_mosaic.py {crop} {'-R' if ii == 0 else ''} -d {d} -g '*.h5' "
                f"-O {base}/{out_file} --in_group {group}/ -F {field}{j}"
                for ii, field in enumerate(fields)]

    commands = {}
    commands['dz'] = lines('dz', 'dz.h5',
                           ["dz", "sigma_dz", "count", "misfit_rms", "misfit_scaled_rms", "mask", "cell_area"])
    for lag_num in lags:
        field = f"dzdt_lag{lag_num}"
        commands[field] = lines(field, f'{field}.h5', [field, f"sigma_{field}", "cell_area"])
    for group in ["avg_dz_40000m", "avg_dz_20000m", "avg_dz_10000m"]:
        out = group.replace("000m", "km").replace("avg_", "")
        commands[group] = lines(group, f'{out}.h5', [group, f"sigma_{group}", "cell_area"])
        group_dt = group.replace("dz", "dzdt")
        out_dt = out.replace("dz", "dzdt")
        for lag_num in lags:
            field = f"{group_dt}_lag{lag_num}"
            commands[field] = lines(field, f'{out_dt}_lag{lag_num}.h5', [field, f"sigma_{field}", "cell_area"])
    # z0 only where its 200 km tiles were made (a local check, or a listing
    # of the prefix when the tiles are on the bucket)
    if not skip_z0 and exists(join_path_or_uri(in_base, '200km_tiles/z0')):
        commands['z0'] = lines('z0', 'z0.h5',
                               ["z0", "misfit_rms", "misfit_scaled_rms", "mask", "cell_area", "count", "sigma_z0"])
    return commands


def make_mosaic_jobs(base, region, lags, skip_z0=False, run=False, in_base=None,
                     group=None, environment='IS2', workers=1):

    mosaic_run = f"mosaic_run_{region}"

    # Check if bounds.txt exists and set crop variable
    bounds_file = os.path.join(base, "bounds.txt")
    if os.path.isfile(bounds_file):
        with open(bounds_file, 'r') as f:
            crop = f"-c {f.readline().strip()}"
    else:
        crop = ""

    # Create mosaic_run directory and subdirectories
    os.makedirs(mosaic_run, exist_ok=True)
    for thedir in ['queue', 'running', 'done', 'logs', 'active_logs', 'error_logs']:
        os.makedirs(os.path.join(mosaic_run, thedir), exist_ok=True)

    commands = mosaic_commands(base, lags, in_base=in_base, crop=crop, skip_z0=skip_z0, workers=workers)
    if group is not None:
        # one group: one DPS job's worth
        if group not in commands:
            raise SystemExit(f"make_200km_to_mosaic_jobs.py: no group {group!r}; "
                             f"the groups are {', '.join(commands)}")
        commands = {group: commands[group]}

    for task, (this_group, these_commands) in enumerate(commands.items(), start=1):
        with open(f"{mosaic_run}/queue/task_{task}", 'w') as f:
            if environment:
                f.write(f"source activate {environment}\n")
            for command in these_commands:
                f.write(command + "\n")

    ATL1415.make_slurm_file(os.path.join(mosaic_run, 'slurm_run.sh'),
                    subs={'JOB_NAME': f'mosaic_{region}',
                          'TIME': "04:00:00",
                          'NUM_TASKS': 6,
                          'JOB_NUMBERS':f'{1}-{task}'})
    if run:
        os.chdir(mosaic_run)
        subprocess.run(["sbatch", "slurm_run.sh"])

def main():
    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,  fromfile_prefix_chars='@')
    parser.add_argument('-b','--base_dir', type=str, default=os.getcwd(), help='directory in which to look for mosaicked .h5 files')
    parser.add_argument('-rr','--region', type=str, help='2-letter region indicator \n'
                                                         '\t A(1-4): Antarctica, by quadrant \n'
                                                         '\t AK: Alaska \n'
                                                         '\t CN: Arctic Canada North \n'
                                                         '\t CS: Arctic Canada South \n'
                                                         '\t GL: Greeland and peripheral ice caps \n'
                                                         '\t IS: Iceland \n'
                                                         '\t SV: Svalbard \n'
                                                         '\t RA: Russian Arctic')
    parser.add_argument('--grid_spacing','-g', type=str, help='grid spacing:DEM (meters),dh maps xy (meters),dh_maps time (years): comma-separated, no spaces', default='100.,1000.,1/4')
    parser.add_argument('--time_span','-t', type=str, help='time span, first year,last year AD (comma separated, no spaces); used to infer --dzdt_lags if not given explicitly')
    parser.add_argument('--dzdt_lags', type=str, default=None, help='comma-separated list of dzdt lags to process; inferred from --time_span and --grid_spacing if omitted')
    parser.add_argument('--run', action='store_true', help="run the script")
    parser.add_argument('--in_base', type=str, help='where the 200 km tiles are read from (<in_base>/200km_tiles/<group>/), if not --base_dir; may be s3://')
    parser.add_argument('--group', type=str, help='write the task for this one group only (z0, dz, dzdt_lag4, avg_dz_40000m, ...)')
    parser.add_argument('--environment','-e', type=str, default='IS2', help="environment each task activates; '' for none")
    parser.add_argument('--workers', '-j', type=int, default=1, help='make_mosaic.py -j: processes reading each command\'s tiles')
    args, unknown = parser.parse_known_args()

    # get the time interval:
    spacing={}
    for dim, this_sp in zip(['z0','dz','dt'], args.grid_spacing.split(',')):
        if '/' in this_sp:
            # this is a fractional spacing (e.g. 1/12 year)
            this_sp = this_sp.split('/')
            this_sp = float(this_sp[0])/float(this_sp[1])
        else:
            this_sp = float(this_sp)
        spacing[dim] = this_sp
    args.grid_spacing = [spacing['z0'], spacing['dz'], spacing['dt']]

    if args.grid_spacing[0] > 1000:
        skip_z0 = True
    else:
        skip_z0 = False

    if args.dzdt_lags is not None:
        lags = [*map(int, args.dzdt_lags.split(','))]
    else:
        time_span = [*map(float, args.time_span.split(','))]
        lags = ATL1415.infer_dzdt_lags(args.grid_spacing[2], time_span)

    make_mosaic_jobs(args.base_dir, args.region, lags,
                     run = args.run, skip_z0=skip_z0, in_base=args.in_base,
                     group=args.group, environment=args.environment, workers=args.workers)

if __name__=="__main__":
    main()
