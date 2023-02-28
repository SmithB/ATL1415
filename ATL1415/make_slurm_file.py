#! /usr/bin/env python


from importlib import resources
import os
import re

def make_slurm_file(dst_file, subs=None, source_file='packable_job.txt', css=False):
    # get the template file from [the current package source directory]/templates
    in_file=os.path.join(resources.files('ATL1415'),'resources','slurm_templates',source_file)

    print(css)
    if subs is None:
        subs={}
    print(subs)
    in_re=re.compile('(\[\[(.*)=(.*)\]\])')
    with open(in_file,'r') as fh_in:
        with open(dst_file,'w') as fh_out:
            for line in fh_in:
                  if css:
                      if "##SBATCH --constraint" in line and "cssro" in line:
                          line=line.replace('##','#')

                  m=in_re.search(line)
                  if m is None:
                      fh_out.write(line)
                      continue
                  pattern=m.group(1)
                  key=m.group(2)
                  print(key)
                  val=m.group(3)
                  if key in subs:
                      val=subs[key]
                  line=line.replace(pattern, val)

                  fh_out.write(line)
