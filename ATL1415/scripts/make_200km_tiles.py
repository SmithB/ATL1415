#! /usr/bin/env python


import glob
import numpy as np
import re
import pointCollection as pc
import os
import stat
import ATL1415
from ATL1415.scripts.setup_slurm_run import setup_directories


def make_fields(dzdt_lags, t_res=0.25, skip_z0=False ):
    """
    Build the field lists and time ranges to mosaic for each output group.

    Parameters
    ----------
    dzdt_lags : iterable
        dzdt lags (in grid-spacing units) to generate dzdt/avg_dzdt groups for.
    t_res : float, optional
        dz/dt grid time resolution, used to convert lags to years. The default is 0.25.
    skip_z0 : bool, optional
        if true, omit the z0 group. The default is False.

    Returns
    -------
    fields : dict
        mapping of group name to list of fields to mosaic.
    time_ranges : dict
        mapping of group name to [start, end] year range for the mosaic.

    """

    fields={}
    if not skip_z0:
        fields['z0']="z0 sigma_z0 misfit_rms misfit_scaled_rms mask cell_area count".split(' ')

    fields['dz']="dz sigma_dz count misfit_rms misfit_scaled_rms mask cell_area".split(' ')

    time_ranges={}
    time_ranges['dz']=[2019, 2050]

    lags = [ f'_lag{lag}' for lag in dzdt_lags ]

    for lag in lags:
        field_str='dzdt'+lag
        fields[field_str] = ["dzdt"+lag, "sigma_dzdt"+lag, "cell_area"]
        time_ranges[field_str] = [2019 + t_res * int(lag.replace('_lag',''))/2, 2050]
    for res in ["_40000m", "_20000m", "_10000m"]:
        fields['avg_dz'+res] = ["avg_dz"+res, "sigma_avg_dz"+res,'cell_area']
        time_ranges['avg_dz'+res] = [2019, 2050]
        for lag in lags:
            field_str='avg_dzdt'+res+lag
            fields[field_str]=[field_str, 'sigma_'+field_str, 'cell_area']
            time_ranges[field_str] = [2019 + t_res * int(lag.replace('_lag',''))/2, 2050]
    #for key, item in fields.items():
    #print(key+" : "+str(item))
    #print(fields)
    return fields, time_ranges

def make_200km_tiles(region_dir, tile_W=200e3):
    """
    Find or build the list of 200km-tile centers for a region.

    Parameters
    ----------
    region_dir : str
        directory containing the region's prelim tile output.
    tile_W : float, optional
        width of the tiles into which the small tiles are grouped, in meters.
        The default is 200e3.

    Returns
    -------
    xyc : list
        list of [x, y] tile-center coordinates.

    """
    print("looking for tiles for "+region_dir)
    tile_ctr_file=os.path.join(region_dir,'200km_tile_list.txt')

    tile_re = re.compile('E(.*)_N(.*).h5')

    if os.path.isfile(tile_ctr_file):
        with open(tile_ctr_file) as fh:
            xyc=[ [*map(float, line.rstrip().split(' '))] for line in fh]
        return xyc

    tile_files=[]
    for sub in ['prelim']:
        tile_files += glob.glob(os.path.join(region_dir, sub, 'E*N*.h5'))

    tile_list=[]
    for tile_name in tile_files:
        try:
            tile_list += [np.array([*map(int, tile_re.search(tile_name).groups())])]
        except Exception as e:
            print(tile_name)
            print(e)

    xy0=np.c_[tile_list]*1000
    xyc=pc.unique_by_rows(np.floor(xy0/tile_W)*tile_W+tile_W/2)
    with open(tile_ctr_file,'w') as fh:
        for line in xyc:
            fh.write(str(line[0])+' '+str(line[1])+'\n')
    return xyc

def main():
    tile_W=2.e5
    avg_re = re.compile('_(\d+)m')
    import argparse
    parser = argparse.ArgumentParser(fromfile_prefix_chars="@")
    parser.add_argument('region_dir', type=str)
    parser.add_argument('region', type=str)
    parser.add_argument('--grid_spacing','-g', type=str, help='grid spacing:DEM (meters),dh maps xy (meters),dh_maps time (years): comma-separated, no spaces', default='100.,1000.,1/4')
    parser.add_argument('--step', type=str, default='matched')
    parser.add_argument('--pad', type=float)
    parser.add_argument('--feather', type=float)
    parser.add_argument('--W', type=int, default=60000)
    parser.add_argument('--spacing', type=int, default=40000)
    parser.add_argument('--skip_sigma', action='store_true')
    parser.add_argument('--skip_z0', action='store_true')
    parser.add_argument('--time_span','-t', type=str, help='time span, first year,last year AD (comma separated, no spaces); used to infer --dzdt_lags if not given explicitly')
    parser.add_argument('--dzdt_lags', type=str, default=None, help='comma-separated list of dzdt lags to process; inferred from --time_span and --grid_spacing if omitted')
    parser.add_argument('--name', type=str)
    parser.add_argument('--lags_only', action='store_true')
    parser.add_argument('--environment','-e', type=str, default='ATL14', help="environment that each job will activate")
    args, _ =parser.parse_known_args()

    region_dir=args.region_dir
    region=args.region

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
        args.skip_z0 = True

    step=args.step

    # make the pad and feather work for 44 km tiles:
    overlap=args.W-args.spacing
    if args.pad is None:
        args.pad = overlap/4
    if args.feather is None:
        args.feather = overlap/2
    print(f"***pad={args.pad}, feather={args.feather}, overlap={overlap}")

    print(f"Skip sigma is {args.skip_sigma}, step is {step}")
    print("region_dir is " +region_dir)

    if args.dzdt_lags is not None:
        dzdt_lags = [*map(int, args.dzdt_lags.split(','))]
    else:
        time_span = [*map(float, args.time_span.split(','))]
        dzdt_lags = ATL1415.infer_dzdt_lags(args.grid_spacing[2], time_span)

    fields, time_ranges=make_fields(dzdt_lags,
                                    t_res = args.grid_spacing[2],
                                    skip_z0 = args.skip_z0)

    xyc=make_200km_tiles(region_dir)

    tile_dir_200km=os.path.join(region_dir,'200km_tiles')
    if not os.path.isdir(tile_dir_200km):
        os.mkdir(tile_dir_200km)

    if args.name is None:
        args.name=region
    run_dir=f'tile_run_{args.name}'
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)

    if os.path.isdir(run_dir+'/logs'):
        N=len(glob.glob(run_dir+'/logs_round_*'))
        for sub in ['logs','done','active_logs','error_logs']:
            os.rename(run_dir+'/'+sub, run_dir+f'/{sub}_round_{N+1}')

    setup_directories(run_dir)

    non_sigma_fields={}
    sigma_fields={}
    for group, field_list in fields.items():
        non_sigma_fields[group]=" ".join([field for field in field_list if 'sigma' not in field])
        sigma_fields[group]=" ".join([field for field in field_list if 'sigma' in field])

    for count, xy in enumerate(xyc):
        search_bounds =[xy[0]-tile_W/2-1.e4, xy[0]+tile_W/2+1.e4, xy[1]-tile_W/2-1.e4, xy[1]+tile_W/2+1.e4]
        search_bounds_str = " ".join([str(ii) for ii in search_bounds])
        tile_bounds = [xy[0]-tile_W/2, xy[0]+tile_W/2, xy[1]-tile_W/2, xy[1]+tile_W/2]
        tile_bounds_1km = "_".join([str(int(ii/1000)) for ii in tile_bounds])
        tile_bounds_str = " ".join([str(ii) for ii in tile_bounds])

        task_file=f'{run_dir}/queue/task_{count+1}'
        with open(task_file,'w') as fh:
            for group in fields.keys():
                pad=args.pad
                feather=args.feather
                spacing_str=""
                avg_scale=avg_re.search(group)
                if avg_scale is not None:
                    avg_scale=float(avg_scale.groups()[0])
                    # NOTE - this is to deal with the truncated 20-km averages in
                    # release 003.  May need to be fixed in the future (avg_scale > overlap makes more sense)
                    if avg_scale >= overlap:
                        pad=0
                        feather=0
                    if "40000m" in group:
                        spacing_str="-S 40000 40000"
                out_dir = os.path.join(tile_dir_200km, group)
                if not os.path.isdir(out_dir):
                    os.mkdir(out_dir)
                out_file = os.path.join(out_dir, f"{group}{tile_bounds_1km}.h5")
                fh.write("#\n")
                if group=='z0':
                    time_str=''
                else:
                    time_str = f"--t_range {time_ranges[group][0]} {time_ranges[group][1]}"
                if args.lags_only and 'lag' not in group:
                    continue
                # bundled lines run one after another regardless of each
                # other's exit status; self-report failure via a grep-able
                # marker so an early line's crash isn't masked by a later
                # line's clean exit (see packable_job.txt's error_logs check)
                rc_suffix = '; rc=$?; [ $rc -ne 0 ] && echo "##TASK_LINE_FAILED## rc=$rc"; (exit $rc)\n'
                fh.write(f"make_mosaic.py -w -R -d {region_dir} -g '{step}/E*.h5' -r {search_bounds_str} -f {feather} -p {pad} -c {tile_bounds_str} -G {group} -F {non_sigma_fields[group]} -O {out_file} {spacing_str} {time_str}"+rc_suffix)
                if not args.skip_sigma:
                    fh.write(f"make_mosaic.py -w  -d {region_dir} -g 'prelim/E*.h5' -r {search_bounds_str} -f {feather} -p {pad} -c {tile_bounds_str} -G {group} -F {sigma_fields[group]} -O {out_file} {spacing_str} {time_str}"+rc_suffix)
        st=os.stat(task_file)
        os.chmod(task_file, st.st_mode | stat.S_IEXEC)

    ATL1415.make_slurm_file(os.path.join(run_dir, 'slurm_mos_run'),
                subs={'JOB_NAME': f'1415_{args.name}_mos',
                      'TIME': "04:00:00",
                      'NUM_TASKS': "4",
                      'ENVIRONMENT': args.environment,
                      'JOB_NUMBERS': f'1-{count+1}'})

if __name__=='__main__':
    main()
