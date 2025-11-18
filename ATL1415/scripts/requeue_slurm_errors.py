#! /usr/bin/env python

import os
import re
import glob
import sys

def make_error_list(list_file, done_dir):
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
                to_rerun=False
                for line in fh:
                    line=line.lower()
                    if 'killed' in line:
                        to_rerun=False
                        break
                    if not to_rerun and 'error' in line:
                        to_rerun=True
                if to_rerun:
                    job_num=done_file.split('_')[-1]
                    job_list += [job_num]
    done_list=[]
    for job_num in job_list:
        done_file=glob.glob(f'{done_dir}/*_{job_num}')
        done_list += done_file
    return done_list

done_dirs=glob.glob('done*')
N=len(done_dirs)

while os.path.isdir(f'done_round_{N}'):
    N += 1
done_old_dir=f'done_round_{N}'


list_file=None
run_name=sys.argv[1]
if len(sys.argv) > 2:
    list_file=sys.argv[2]

error_list=make_error_list(list_file, 'done')
print(error_list)

if len(error_list) ==0:
    print("no errors found, returning")
    sys.exit()


count=len(glob.glob('queue/*'))

for file in error_list:
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
    out_file=f'queue/{run_name}_{count}'
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

os.rename('done',done_old_dir)
os.mkdir('done')

