import os
import sys
import argparse
import ATL1415
import subprocess

def make_mosaic_jobs(base, region, lags, t_res=0.25, skip_z0=False, run=False):

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
    for thedir in ['queue', 'running', 'done', 'logs']:
        os.makedirs(os.path.join(mosaic_run, thedir), exist_ok=True)

    compute_sigma = 'True'
    field = 'dz'
    glob_str = "'matched/*.h5'"
    task = 0

    task += 1
    this_replace = '-R'
    field_list = ["dz", "sigma_dz", "count", "misfit_rms", "misfit_scaled_rms", "mask", "cell_area"]

    with open(f"{mosaic_run}/queue/task_{task}", 'w') as f:
        f.write("source activate IS2\n")
        for field in field_list:
            print(field)
            f.write(f"make_mosaic.py {crop} {this_replace} -d {base}/200km_tiles/dz -g '*.h5' -O {base}/dz.h5 --in_group dz/ -F {field}\n")
            this_replace = ''

    for lag_num in lags:
        lag = f'_lag{lag_num}'
        task += 1
        TR = f"--t_range {2019 + t_res * int(lag_num) / 2} 2050"
        field = f"dzdt{lag}"
        with open(f"{mosaic_run}/queue/task_{task}", 'w') as f:
            f.write("source activate IS2\n")
            this_replace = '-R'
            for ff in [field, f"sigma_{field}", "cell_area"]:
                f.write(f"make_mosaic.py {crop} {TR} {this_replace} -d {base}/200km_tiles/{field} -g '*.h5' -O {base}/{field}.h5 --in_group {field}/ -F {ff}\n")
                this_replace = ''

    for group in ["avg_dz_40000m", "avg_dz_20000m", "avg_dz_10000m"]:
        field = group
        print(group, field)
        this_replace = '-R'
        out = group.replace("000m", "km").replace("avg_", "")

        task += 1
        with open(f"{mosaic_run}/queue/task_{task}", 'w') as f:
            f.write("source activate IS2\n")
            this_replace = '-R'
            for ff in [field, f"sigma_{field}", "cell_area"]:
                f.write(f"make_mosaic.py {crop} {this_replace} -d {base}/200km_tiles/{group} -g '*.h5' -O {base}/{out}.h5 --in_group {group}/ -F {ff}\n")
                this_replace = ''

        group = group.replace("dz", "dzdt")
        out = out.replace("dz", "dzdt")
        for lag_num in lags:
            lag = f'_lag{lag_num}'
            field = f"{group}{lag}"
            TR = f"--t_range {2019 + t_res * lag_num / 2} 2050"
            task += 1
            with open(f"{mosaic_run}/queue/task_{task}", 'w') as f:
                f.write("source activate IS2\n")
                this_replace = '-R'

                print(f"lag={lag}, group={group}, task={task}")
                for ff in [field, f"sigma_{field}", "cell_area"]:
                    f.write(f"make_mosaic.py {crop} {TR} {this_replace} -d {base}/200km_tiles/{field} -g '*.h5' -O {base}/{out}{lag}.h5 --in_group {field}/ -F {ff}\n")
                    this_replace = ''

    if os.path.isdir(os.path.join(base, "200km_tiles", "z0")) and not skip_z0:
        field_list = ["z0", "misfit_rms", "misfit_scaled_rms", "mask", "cell_area", "count"]
        task += 1
        this_replace = '-R'

        field_list.append("sigma_z0")
        for field in field_list:
            print(field)
            with open(f"{mosaic_run}/queue/task_{task}", 'w') as f:
                f.write(f"make_mosaic.py {crop} {this_replace} -d {base}/200km_tiles/z0 -g '*.h5' -O {base}/z0.h5 --in_group z0/ -F {field}\n")
                this_replace = ''

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
    parser.add_argument('--dzdt_lags', type=str, help='comma-separated list of dzdt lags to process')
    parser.add_argument('--run', action='store_true', help="run the script")
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
    
    make_mosaic_jobs(args.base_dir, args.region, [ *map(int, args.dzdt_lags.split(',')) ], 
                     t_res = args.grid_spacing[2],
                     run = args.run, skip_z0=skip_z0)

if __name__=="__main__":
    main()
