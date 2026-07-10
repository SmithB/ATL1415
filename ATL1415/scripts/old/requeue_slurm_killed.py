#! /usr/bin/env python

import os
import re
import glob
import sys

def make_killed_list(list_file, done_dir, run_name):
    job_re=re.compile('\d+_(\d+):')
    job_list=[]
    if list_file is not None:
        with open(list_file,'r') as fh:
            for line in fh:
                try:
                    job_num=job_re.search(line).group(1)
                    job_list += [job_num]
                except Exception:
                    continue
    else:
        all_dones=glob.glob(done_dir+'/*_*')
        job_list=[]
        for done_file in all_dones:
            with open('logs/'+os.path.basename(done_file)+'.log','r') as fh:
                for line in fh:
                    if 'illed' in line:
                        job_num=done_file.split('_')[-1]
                        job_list += [job_num]
                        break
    if run_name is not None:
        slurm_log_files=glob.glob(run_name+'*_*')
        print(f"found {len(slurm_log_files)} log files")
        for file in slurm_log_files:
            with open(file,'r') as fh:
                for line in fh:
                    line=line.lower()
                    if 'kill' in line or 'cancelled' in line:
                        job_num=file.split('_')[-1]
                        job_list += [job_num]
                        print(job_num)
        print(f"found {len(job_list)} killed or cancelled files")
    done_list=[]
    for job_num in set(job_list):
        done_file=glob.glob(f'{done_dir}/*_{job_num}')
        if len(done_file)==0:
            print(f"missing done file for {job_num} in {done_dir}")
        done_list += done_file
    return done_list

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--done_dir', type=str, default='done')
    parser.add_argument('--log_dir', type=str, default='logs')
    parser.add_argument('--run_name', type=str)
    parser.add_argument('--list_file', type=str)
    parser.add_argument('--new_run_name', type=str)
    parser.add_argument('--dry_run', action='store_true')
    args=parser.parse_args()

    done_dirs=glob.glob('done*')
    N_done_dirs=len(done_dirs)

    while os.path.isdir(f'done_round_{N_done_dirs}'):
        N_done_dirs += 1
    done_old_dir=f'done_round_{N_done_dirs}'

    if args.new_run_name is None:
        args.new_run_name=f'taskr{N_done_dirs+1}'

    killed_list=make_killed_list(args.list_file, args.done_dir, args.run_name)


    count=len(glob.glob('queue/*'))

    for file in killed_list:
        status='not run'

        log_file='logs/'+os.path.basename(file)+'.log'

        with open(log_file,'r') as fh:
            for line in fh:
                if 'done with' in line:
                    if status=='not run':
                        status='run but no errors'
                    else:
                        print(f'file: {file} appears to have completed')
                        break
        count += 1
        out_file=f'queue/{args.new_run_name}_{count}'
        if args.dry_run:
            print(f'\t---reprocessing {file} as {out_file}')
            continue
        with open(out_file,'w') as fh_out:
            fh_out.write('#! /usr/bin/env bash\n')
            fh_out.write('# reprocess ' +file+'\n')

            with open(file) as fh_in:
                for line in fh_in:
                    if 'ATL11_to_ATL15' in line:
                        if status=='not run':
                            fh_out.write(line)
                        elif status=='run but no errors' and 'calc_error' in line:
                            fh_out.write(line)
                    else:
                        fh_out.write(line)


    print(f'moving {args.done_dir}/ to ', done_old_dir)
    if not args.dry_run:
        os.rename(args.done_dir,done_old_dir)
        os.mkdir('done')

if __name__=='__main__':
    main()
