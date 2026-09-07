from .ATL11_to_ATL15 import *
from .ATL1415_attrs_meta import *
from .make_slurm_file import *
from .lags import infer_dzdt_lags
from .assign_firn_variable import *
from .SMB_corr_from_grid import SMB_corr_from_grid
from .read_ATL11 import *
from  .tile_names import tile_centers_from_files
from  .tile_names import tile_centers_from_scripts
from .make_nc_projection_variable import *
from .make_tile_stats_group import *

# pyTMD sets AWS_NO_SIGN_REQUEST=YES in the environment when pyTMD.io is first
# imported (v3.0.9, pyTMD/io/__init__.py), so that its own reads of the public
# s3://pytmd stores are anonymous.  It is a PROCESS-WIDE variable, and GDAL
# honours it for every /vsis3 read -- including the masks, geoid and tide masks
# this code reads from s3://maap-ops-workspace, which then fail with HTTP 403
# and surface as ogr.Open() returning None.  The first line of this file
# imports ATL11_to_ATL15, which imports pyTMD at module level, so by the time
# execution reaches here pyTMD.io has run and this pop is the last word.
# Our tide reads do not need the variable: ATL1415/tides.py asks for anonymity
# explicitly with s3fs.S3FileSystem(anon=...), which ignores it.
import os as _os
_os.environ.pop('AWS_NO_SIGN_REQUEST', None)
