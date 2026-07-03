HEMI_EPSG = {'north':3413,'south':3031}
HEMI_ABBREV = {'north':'AR','south':'AA'}
HEMI_NAME = {'north':'Arctic','south':'Antarctic'}
import os
import glob
import pointCollection as pc
import re
def setup_ATL11_xover(dst_dir, ATL11xo_top=None, ATL11xo_version=None, hemi=None, cycles=['01','02']):
    print(f"\tsetup_ATL11xo: \n\t\tdst_dir={dst_dir},\n\t\themi={hemi}")
    xover_src = os.path.join(ATL11xo_top, HEMI_NAME[hemi]+'_'+ATL11xo_version)
    xover_dst = os.path.join(dst_dir, 'xover_tiles')


    try:
        release, version = re.compile('(\d\d\d)_cycle_\d\d_\d\d_v(\d\d)')\
                        .search(ATL11xo_version)\
                        .groups()
    except AttributeError:
        raise AttributeError('ATL11xo version did not match pattern (rrr)_cycle_(cc)_(cc)_v(vv)')

    for cycle in cycles:
        cycle_src = os.path.join(xover_src, 'xover_tiles', 'cycle_'+cycle)
        if not os.path.isdir(cycle_src):
            raise PathError(cycle_src+' not found')
        cycle_dst = os.path.join(xover_dst, os.path.basename(cycle_src))
        os.makedirs(cycle_dst, exist_ok=True)
        schema_file = os.path.join(cycle_dst,f'200km_tiling_{HEMI_ABBREV[hemi]}.json')
        this_format_str = f"ATL11XO_{HEMI_ABBREV[hemi]}_E%d_N%d_c{cycle}_{release}_{version}"

        pc.tilingSchema(tile_spacing=200e3,
                        mapping_function_name='round',
                        format_str=this_format_str,
                        scale=1000,
                        directory=cycle_src,
                        extension='.h5').to_json(schema_file)
