#!/usr/bin/env python3
"""
Count ice-bearing AA tiles.

Spacing is 40 km for BOTH halves (--tile_spacing); -W (60/44 km) is the tile
WIDTH, so tiles overlap.  Read the mask ONCE, decimated, instead of issuing a
windowed S3 read per tile -- the per-tile version was OOM-killed.
"""
import os, numpy as np, rasterio
os.environ.setdefault('AWS_REGION','us-west-2')
B='s3://maap-ops-workspace/ben_smith/ATL1415/masks/Antarctic/'
STEP=40000; MIN_XY, MAX_XY = 360000, 440000
DEC=16   # 240 m -> 3.84 km pixels

with rasterio.open(B+'AntarcticIceMask_2018.00_2026.25_240m_v4.1.tif') as ice:
    b=ice.bounds
    h=int(np.ceil(ice.height/DEC)); w=int(np.ceil(ice.width/DEC))
    a=ice.read(1, out_shape=(h,w))>0
    px=(b.right-b.left)/w; py=(b.top-b.bottom)/h
print(f'mask {a.shape}, {px/1000:.2f} km px, {a.sum()} ice px', flush=True)

def has_ice(x0,y0,W):
    c0=int(np.floor((x0-W/2-b.left)/px)); c1=int(np.ceil((x0+W/2-b.left)/px))
    r0=int(np.floor((b.top-(y0+W/2))/py)); r1=int(np.ceil((b.top-(y0-W/2))/py))
    c0=max(c0,0); r0=max(r0,0); c1=min(c1,w); r1=min(r1,h)
    if c1<=c0 or r1<=r0: return False
    return bool(a[r0:r1, c0:c1].any())

xs=np.arange(np.round(b.left/STEP)*STEP, b.right+STEP, STEP)
ys=np.arange(np.round(b.bottom/STEP)*STEP, b.top+STEP, STEP)
tot=0
for W,label,sel in ((60000,'60 km','outer'),(44000,'44 km','inner')):
    n=cand=0
    for x in xs:
        for y in ys:
            mx=max(abs(x),abs(y))
            if sel=='outer' and mx < MIN_XY: continue
            if sel=='inner' and mx > MAX_XY: continue
            cand+=1
            if has_ice(x,y,W): n+=1
    print(f'{label}: {n} ice-bearing of {cand} candidate centres', flush=True)
    tot+=n
print(f'TOTAL: {tot} solve jobs')
