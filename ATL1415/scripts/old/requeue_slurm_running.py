#! /usr/bin/env python

import os
import re
import glob
import sys

if not os.path.isdir('running_old'):
    os.rename('running','running_old')
    os.mkdir('running')

run_name=sys.argv[1]

count=len(glob.glob('queue/*'))

for file in glob.glob('running_old/*'):
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

            
 
