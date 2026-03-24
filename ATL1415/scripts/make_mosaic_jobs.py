#!/usr/bin/env python3
import os
import argparse
import subprocess
import ATL1415

def write_task(task_id, content, mosaic_run, environment, append=False):
    """
    Write a task to a script that can be run in parallel

    Parameters
    ----------
    task_id : int
        number of task to write.
    content : str
        content to write to script.
    mosaic_run : str
        directory into which to write the output.
    environment : str
        environment to be activated.
    append : bool, optional
        if true, append to task file. The default is False.

    Returns
    -------
    None.

    """

    mode = "a" if append else "w"
    task_file = os.path.join(mosaic_run, "queue" , f"task_{task_id}")
    with open(task_file, mode) as f:
        if mode == 'w':
            f.write(f"source activate {environment};\n")
        f.write(content + "\n")

def make_mosaic_jobs(base, region, lags,
                     step='matched',
                     t_res=0.25,
                     skip_z0=False,
                     tasks=4,
                     environment='IS2'):
    """
    make a set of jobs for mosaicking a region

    Parameters
    ----------
    base: str
        directory containing regions.
    region : str
        region name.
    lags : iterable
        lags for which to calculate dz/dt.
    step : str, optional
        'prelim' or 'matched. The default is 'matched'.
    t_res : float, optional
        time interval between steps. The default is 0.25.
    skip_z0 : bool, optional
        if true, do not mosaic z0. The default is False.
    tasks : int, optional
        tasks per slurm job. The default is 4.
    environment : str, optional
        environment to be activated. The default is 'IS2'.

    Returns
    -------
    mosaic_run: mosaic directory.
    tasks: tasks created

    """

    mosaic_run = f"mosaic_run_{region}"

    # crop handling
    bounds_file = os.path.join(base, "bounds.txt")
    if os.path.exists(bounds_file):
        with open(bounds_file) as f:
            first_line = f.readline().strip()
        crop = f"-c {first_line}"
    else:
        crop = ""

    # initial parameters
    pad = 5000
    feather = 10000
    compute_sigma = True
    compute_SMB = False

    glob_str = f"'{step}/*.h5'"

    # make directories
    os.makedirs(mosaic_run, exist_ok=True)
    for sub in ["queue", "running", "done", "logs"]:
        os.makedirs(os.path.join(mosaic_run, sub), exist_ok=True)

    task = 0
    if not skip_z0:
        field='z0'
        group='z0'
        field_list = ["z0", "misfit_rms", "misfit_scaled_rms", "mask",
                      "cell_area", "count","sigma_z0'"]
        task += 1
        this_replace = '-R'
        append=False
        for field in field_list:
            cmd = (
                f"make_mosaic.py {crop} {this_replace} "
                f"-d {base} -g {glob_str} -p {pad} -f {feather} "
                f"-O {base}/z0.h5 --in_group {group}/ -F {field}"
            )
            this_replace=""
            write_task(task, cmd, mosaic_run, environment, append=append)
            append = True

    # ---- first loop over groups ----
    for group in ["avg_dz_40000m", "avg_dz_20000m", "avg_dz_10000m"]:
        field = group

        this_pad = pad
        this_feather = feather
        this_S = ""
        this_w = "-w"

        if group == "avg_dz_40000m":
            this_pad = 0
            this_feather = 0
            this_S = "-S 40000 40000"
            this_w = ""
        elif group == "avg_dz_20000m":
            this_pad = 0
            this_feather = 0
            this_S = ""
            this_w = ""

        out = group.replace("000m", "km").replace("avg_", "")
        TR = "--t_range 2019 2050"

        task += 1
        cmd = (
            f"make_mosaic.py {crop} {this_w} {TR} -R "
            f"-d {base} -g {glob_str} -p {this_pad} -f {this_feather} {this_S} "
            f"-O {base}/{out}.h5 --in_group {group}/ -F {field} cell_area"
        )
        write_task(task, cmd, mosaic_run, environment)

        if compute_sigma:
            cmd2 = (
                f"make_mosaic.py {crop} {this_w} {TR} -d {base} "
                f"-g 'prelim/*.h5' -p {this_pad} -f {this_feather} {this_S} "
                f"-O {base}/{out}.h5 --in_group {group}/ -F sigma_{group}"
            )
            write_task(task, cmd2,  mosaic_run, environment, append=True)

        group_dt = group.replace("dz", "dzdt")
        out_dt = out.replace("dz", "dzdt")

        for lag_num in lags:
            lag = f'_lag{lag_num}'
            field = f"{group_dt}{lag}"
            field_list = f"{field} cell_area"

            task += 1

            start_year = 2019 + t_res * lag_num / 2
            TR = f"--t_range {start_year} 2050"

            cmd = (
                f"make_mosaic.py {crop} {TR} -R {this_w} "
                f"-d {base} -g {glob_str} -p {this_pad} -f {this_feather} {this_S} "
                f"-O {base}/{out_dt}{lag}.h5 --in_group {field}/ -F {field_list}"
            )
            write_task(task, cmd,  mosaic_run, environment)

            if compute_sigma:
                sigma_field = f"sigma_{group_dt}{lag}"
                cmd2 = (
                    f"make_mosaic.py {crop} {TR} {this_w} -d {base} "
                    f"-g 'prelim/*.h5' -p {this_pad} -f {this_feather} {this_S} "
                    f"-O {base}/{out_dt}{lag}.h5 --in_group {field}/ -F {sigma_field}"
                )
                write_task(task, cmd2, mosaic_run, environment, append=True)

    # ---- dz base task ----
    field = "dz"
    TR = "--t_range 2019 2050"

    task += 1
    cmd = (
        f"make_mosaic.py {crop} {TR} -R -w "
        f"-d {base} -g {glob_str} -p {pad} -f {feather} "
        f"-O {base}/dz.h5 --in_group dz/ "
        f"-F count misfit_rms misfit_scaled_rms mask cell_area {field}"
    )
    write_task(task, cmd, mosaic_run, environment)

    if compute_sigma:
        cmd2 = (
            f"make_mosaic.py {crop} {TR} -w -d {base} "
            f"-g 'prelim/*.h5' -p {pad} -f {feather} "
            f"-O {base}/dz.h5 --in_group dz/ -F sigma_dz"
        )
        write_task(task, cmd2, mosaic_run, environment, append=True)

    if compute_SMB:
        cmd3 = (
            f"make_mosaic.py {crop} -w {TR} -d {base} "
            f"-g 'prelim/*.h5' -p {pad} -f {feather} "
            f"-O {base}/dz.h5 --in_group dz/ -F SMB_a FAC"
        )
        write_task(task, cmd3, mosaic_run, environment, append=True)

    # ---- lag loop for dzdt ----
    for lag_num in lags:
        lag = f'_lag{lag_num}'
        start_year = 2019 + t_res * lag_num / 2
        TR = f"--t_range {start_year} 2050"

        task += 1
        field = f"dzdt{lag}"

        cmd = (
            f"make_mosaic.py {crop} {TR} -R -w "
            f"-d {base} -g {glob_str} -p {pad} -f {feather} "
            f"-O {base}/dzdt{lag}.h5 --in_group dzdt{lag}/ "
            f"-F {field} cell_area"
        )
        write_task(task, cmd, mosaic_run, environment)

        if compute_sigma:
            cmd2 = (
                f"make_mosaic.py {crop} {TR} -w -d {base} "
                f"-g 'prelim/*.h5' -p {pad} -f {feather} "
                f"-O {base}/dzdt{lag}.h5 --in_group dzdt{lag}/ "
                f"-F sigma_{field}"
            )
            write_task(task, cmd2, mosaic_run, environment, append=True)

    return mosaic_run, task


def main():

    parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,  fromfile_prefix_chars='@')
    parser.add_argument('-b','--base_dir', type=str, default=os.getcwd(), help='region in which to mosaic')
    parser.add_argument('-rr','--region', type=str,
                        help='2-letter region indicator \n'
                                 '\t CN: Arctic Canada North \n'
                                 '\t CS: Arctic Canada South \n'
                                 '\t GL: Greeland and peripheral ice caps \n'
                                 '\t IS: Iceland \n'
                                 '\t SV: Svalbard \n'
                                 '\t RA: Russian Arctic')
    parser.add_argument('--grid_spacing','-g', type=str, help='grid spacing:DEM (meters),dh maps xy (meters),dh_maps time (years): comma-separated, no spaces', default='100.,1000.,1/4')
    parser.add_argument('--dzdt_lags', type=str, help='comma-separated list of dzdt lags to process')
    parser.add_argument('--run', action='store_true', help="run the script")
    parser.add_argument('--num_tasks', type=int, default=4, help='number of slurm tasks to assign')
    parser.add_argument('--environment','-e', type=str, default='IS2', help='environment to activate for each job')
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

    lags = [*map(int, args.dzdt_lags.split(','))]
    mosaic_run, tasks = make_mosaic_jobs(args.base_dir, args.region,
                                         lags,
                                         t_res = args.grid_spacing[2],
                                         skip_z0=skip_z0,
                                         environment=args.environment)
    ATL1415.make_slurm_file(os.path.join(mosaic_run, 'slurm_run.sh'),
                subs={'JOB_NAME': f'mosaic_{args.region}',
                      'TIME': "04:00:00",
                      'NUM_TASKS': tasks,
                      'JOB_NUMBERS':f'{1}-{tasks}'})

    if args.run:
        os.chdir(mosaic_run)
        subprocess.run(["sbatch", "slurm_run.sh"])

if __name__ == "__main__":
    main()
