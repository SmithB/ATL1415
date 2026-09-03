"""
Tests for ATL1415.read_ATL11._lonlat_bounding_box(), which turns a
projected [x,y] box into a geographic bounding box for an
earthaccess.search_data(bounding_box=...) call.

No network access is required. Imports read_ATL11.py directly via
importlib, bypassing ATL1415/__init__.py (which pulls in LSsurf, not
installed in every environment this runs in, and unrelated to this file).
"""
import os
import importlib.util
import numpy as np
import pytest

_HERE = os.path.dirname(__file__)
_MODULE_PATH = os.path.join(_HERE, '..', 'ATL1415', 'read_ATL11.py')
_spec = importlib.util.spec_from_file_location('read_ATL11', _MODULE_PATH)
read_ATL11_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(read_ATL11_mod)
_lonlat_bounding_box = read_ATL11_mod._lonlat_bounding_box

SRS_PROJ4_AA = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'


def test_mid_latitude_box_no_wraparound():
    # a modest box well away from the pole and the antimeridian (which, for
    # this south-polar-stereographic lon_0=0 projection, runs along the
    # y-axis -- pick an off-axis point instead)
    xy0 = (1.5e6, -1.5e6)
    W = 6.e4
    bounds = [xy0[0]+np.array([-W/2, W/2]), xy0[1]+np.array([-W/2, W/2])]
    lon_min, lat_min, lon_max, lat_max = _lonlat_bounding_box(bounds, SRS_PROJ4_AA)
    assert lon_min < lon_max
    assert lat_min < lat_max
    assert -90 <= lat_min <= lat_max <= -60


def test_antimeridian_crossing_box_reports_west_gt_east():
    # SRS_PROJ4_AA has lon_0=0 -- the negative-y axis in this south-polar
    # stereographic projection sits on the +-180 deg meridian, so a box
    # centered there straddles the antimeridian without containing the pole
    xy0 = (0., -2.2e5)
    W = 4.4e4
    bounds = [xy0[0]+np.array([-W/2, W/2]), xy0[1]+np.array([-W/2, W/2])]
    lon_min, lat_min, lon_max, lat_max = _lonlat_bounding_box(bounds, SRS_PROJ4_AA)
    # crossing the antimeridian is signalled by west > east
    assert lon_min > lon_max
    assert lon_min > 90
    assert lon_max < -90
    assert lat_min < lat_max
    assert lat_max < -85


def test_box_containing_pole_returns_full_longitude_range():
    xy0 = (0., 0.)
    W = 4.4e4
    bounds = [xy0[0]+np.array([-W/2, W/2]), xy0[1]+np.array([-W/2, W/2])]
    lon_min, lat_min, lon_max, lat_max = _lonlat_bounding_box(bounds, SRS_PROJ4_AA)
    assert lon_min == pytest.approx(-180.)
    assert lon_max == pytest.approx(180.)
    assert lat_max <= -89
