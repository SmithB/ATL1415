#! /usr/bin/env python3


import pointCollection as pc
import matplotlib.pyplot as plt
import numpy as np
import sys
import json
import os

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
           return int(obj)
        if isinstance(obj, np.floating):
           return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

pad=1.e3
tols=[0, 25, 50, 100]
tols.sort()
tols=tols[::-1]
print(tols)

tile_file=sys.argv[1]
DEM_file=sys.argv[2]
dest_dir=sys.argv[3]

if not os.path.isdir(dest_dir):
	os.mkdir(dest_dir)
for tol in tols:
	tol_dir=os.path.join(dest_dir, str(tol))
	if not os.path.isdir(tol_dir):
		os.mkdir(tol_dir)

D=pc.data().from_h5(tile_file, group='data')

bounds=[[np.min(D.x)-pad, np.max(D.x)+pad], [np.min(D.y)-pad, np.max(D.y)+pad]]

D.assign({'DEM':pc.grid.data().from_geotif(DEM_file, bounds=bounds).interp(D.x, D.y)})
D.assign({'r':D.z-D.DEM})
D_ed=D[(D.three_sigma_edit==1) & (np.isfinite(D.sigma))]



report_file=os.path.basename(tile_file).replace('.h5','.DEM_deltas.json')
for tol in tols:
	if np.any(np.abs(D_ed.r)>tol):
		this_tol=tol
		break
report_file=os.path.join(dest_dir, str(this_tol), report_file)
sigma_bins=[0] + list(2**np.arange(0, 7))
report={'tile_file':tile_file, 
				'DEM_file':DEM_file, 
				'tols':[tol for tol in tols],
				'N':[np.sum(np.abs(D_ed.r)>tol) for tol in tols],
                                'N_raw':[np.sum(np.abs(D.r)>tol) for tol in tols], 
				'sigma':np.nanstd((D.z-D.DEM)[D.three_sigma_edit==1]),
	                        'std_r_data':np.nanstd((D.z-D.z_est)[D.three_sigma_edit==1]),
                'sigma_bins':sigma_bins,
                'N_sigma_gt':[np.sum(D_ed.sigma>sigma_i) for sigma_i in sigma_bins],
                'N_sigma_raw_gt':[np.sum(D.sigma>sigma_i) for sigma_i in sigma_bins],
                'rss_sigma_ed':[np.sqrt(np.mean(D_ed.sigma[np.abs(D_ed.r)>tol]**2)) for tol in tols]}
print(report)

with open(report_file, 'w', encoding='utf-8') as f:
    json.dump(report, f, ensure_ascii=False, indent=4, cls=NpEncoder)


