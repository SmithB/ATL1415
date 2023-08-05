


#This notebook contains code to add another year of dh/dt values onto the ATL15_output_attrs.csv file. It duplicates the fields for a specified lag from each group in the file, replaces the name of the lag with a larger lag, and writes the results out to a new csv file.

#The original CSV files did not have a newline character on the last line. Without this, this script will omit the last, largest average, so it is important to manually add the newline if it is not present.

#After running this, move the original ATL15_output_attrs.csv into a cache directory, and move the new version to ATL15_output_attrs.csv

import re


import re

in_file='ATL15_output_attrs.csv'
out_file='ATL15_output_attrs_rel003.csv'

old_lag='lag12'
new_lag='lag16'
old_name='Triennial'
new_name='Quadrennial'

with open(in_file,'r') as fh_in:
    with open(out_file,'w') as fh_out:
        for line in fh_in:
            #print(line)
            if old_lag in line:
                # check if lagxx is in the line, if so create an updated version with the new lag
                temp=line.replace(old_lag, new_lag).replace(old_name.lower(), new_name.lower()).replace(old_name, new_name)
                field=line.split(',')[1]
                if 'time' in field or 'ice_area' in field:
                    # time_lagxx or ice_area_lagxx fields: write immediately
                    fh_out.write(line)
                    fh_out.write(temp)
                elif 'sigma' not in line:
                    # delta_h_lagxx or dhdt_lagxx fields: write out the field, cache the
                    # updated field
                    fh_out.write(line)
                    not_sigma_line = temp
                else:
                    # sigma line: write the field, followed by the cached field and the updated field.
                    fh_out.write(line)
                    fh_out.write(not_sigma_line)
                    fh_out.write(temp)
            else:
                # Everything else gets written immediately
                fh_out.write(line)
                