import csv, numpy as np, rasterio, pyproj, os
from rasterio.windows import from_bounds
os.environ.setdefault('AWS_REGION','us-west-2')
B='s3://maap-ops-workspace/ben_smith/ATL1415/masks/Antarctic/'
ice  = rasterio.open(B+'AntarcticIceMask_2018.00_2026.25_240m_v4.1.tif')
tide = rasterio.open(B+'BedMachineAntarcticaOceanv2.tif')
to_ll = pyproj.Transformer.from_crs(ice.crs,'EPSG:4326',always_xy=True)
W=60000.0
def frac(src,x0,y0):
    a=src.read(1,window=from_bounds(x0-W/2,y0-W/2,x0+W/2,y0+W/2,src.transform),
               boundless=True,fill_value=0)
    return float(np.mean(a>0)) if a.size else 0.0

rows=[]
for fn, group in (('AA_transect_xy.txt','transect'), ('AA_tide_xy.txt','tide')):
    for ln in open(fn):
        if not ln.strip(): continue
        x,y = (int(float(v)) for v in ln.split())
        i,t = frac(ice,x,y), frac(tide,x,y)
        lon,lat = to_ll.transform(x,y)
        r_km = int(round(np.hypot(x,y)/1000))
        if   group=='tide' and t>0.99: why='floating ice: full tide correction'
        elif group=='tide':            why='grounding line: partial tide correction'
        elif i==0:                     why='inside ICESat-2 pole hole: expect a clean empty skip'
        elif i<0.95:                   why='pole-hole edge: peak track convergence'
        elif r_km>=2000:               why='coastal margin: lowest track density'
        else:                          why=f'grounded interior, r={r_km} km'
        rows.append(dict(x0=x, y0=y, lat=round(lat,2), r_km=r_km,
                         ice_frac=round(i,3), tide_frac=round(t,3),
                         group=group, purpose=why))

rows.sort(key=lambda r: (r['group'], r['r_km']))
with open('AA_queue_manifest.csv','w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
open('AA_queue_xy.txt','w').write('\n'.join(f"{r['x0']} {r['y0']}" for r in rows)+'\n')

print(f"{'x0':>9} {'y0':>9} {'lat':>7} {'ice':>5} {'tide':>5} {'group':>8}  purpose")
for r in rows:
    print(f"{r['x0']:9d} {r['y0']:9d} {r['lat']:7.2f} {r['ice_frac']:5.2f} "
          f"{r['tide_frac']:5.2f} {r['group']:>8}  {r['purpose']}")
print(f'\n{len(rows)} tiles -> AA_queue_xy.txt, AA_queue_manifest.csv')
