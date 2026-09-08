import numpy as np, rasterio, pyproj, os
from rasterio.windows import from_bounds
os.environ.setdefault('AWS_REGION','us-west-2')
B='s3://maap-ops-workspace/ben_smith/ATL1415/masks/Antarctic/'
ice  = rasterio.open(B+'AntarcticIceMask_2018.00_2026.25_240m_v4.1.tif')
tide = rasterio.open(B+'BedMachineAntarcticaOceanv2.tif')
to_ll = pyproj.Transformer.from_crs(ice.crs,'EPSG:4326',always_xy=True)
HALF, W = 20000.0, 60000.0
def frac(src,x0,y0):
    a=src.read(1,window=from_bounds(x0-W/2,y0-W/2,x0+W/2,y0+W/2,src.transform),
               boundless=True,fill_value=0)
    return float(np.mean(a>0)) if a.size else 0.0
def on_grid(v):
    n=np.round(v/HALF)
    return (n+1 if n%2==0 else n)*HALF

# named shelf search boxes (EPSG:3031 km), kept coarse on purpose
BOXES = {'Ross':   (-600, 400, -1300, -500),
         'Ronne':  (-1500, -500, 0, 900),
         'Amery':  (1500, 2300, 500, 1100)}
shelf, grounding = [], []
for name,(x0k,x1k,y0k,y1k) in BOXES.items():
    for xk in range(x0k, x1k+1, 80):
        for yk in range(y0k, y1k+1, 80):
            x,y = on_grid(xk*1000), on_grid(yk*1000)
            i = frac(ice,x,y)
            if i < 0.99: continue
            t = frac(tide,x,y)
            if t > 0.99:            shelf.append((name,x,y,i,t))
            elif 0.25 < t < 0.75:   grounding.append((name,x,y,i,t))

def pick(rows, n):
    """spread the picks across shelves rather than clustering in one"""
    out, seen = [], {}
    for r in rows:
        if seen.get(r[0],0) < n:
            out.append(r); seen[r[0]] = seen.get(r[0],0)+1
    return out

sel = pick(shelf,1) + pick(grounding,1)
print(f'{"shelf":8} {"x0":>9} {"y0":>9} {"lat":>7} {"ice":>5} {"tide":>5}  class')
lines=[]
for name,x,y,i,t in sel:
    lon,lat = to_ll.transform(x,y)
    cls = 'floating (full tide correction)' if t>0.99 else 'grounding line (mixed)'
    print(f'{name:8} {x:9.0f} {y:9.0f} {lat:7.2f} {i:5.2f} {t:5.2f}  {cls}')
    lines.append(f'{x:.0f} {y:.0f}')
open('AA_tide_xy.txt','w').write('\n'.join(lines)+'\n')
print(f'\nwrote AA_tide_xy.txt ({len(lines)} tiles); '
      f'candidates: {len(shelf)} floating, {len(grounding)} grounding-line')
