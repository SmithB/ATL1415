"""
ATL11_to_ATL15 --no_data_group (docs/plan_dps_mosaic.sh D2b): a matched tile
can leave out its per-point /data group (80-92% of the file), which nothing
downstream reads.  A prelim tile must keep it -- the matched and error runs
reread it -- so the flag is refused, at parse time, without --matched.
"""
from types import SimpleNamespace

import h5py
import numpy as np
import pytest
import pointCollection as pc

from ATL1415.ATL11_to_ATL15 import parse_args, save_fit_to_file

BASE = ['ATL11_to_ATL15.py', '--ATL11_index', '/tmp/index.h5', '--xy0', '0', '0']


def test_off_by_default():
    assert parse_args(BASE + ['--matched']).no_data_group is False


def test_matched_may_drop_data():
    assert parse_args(BASE + ['--matched', '--no_data_group']).no_data_group is True


@pytest.mark.parametrize('step', [['--prelim'], []], ids=['prelim', 'no step'])
def test_refused_without_matched(step):
    with pytest.raises(SystemExit) as exit_info:
        parse_args(BASE + step + ['--no_data_group'])
    assert exit_info.value.code == 2


def fit_result():
    """the smallest S that save_fit_to_file writes"""
    x = np.arange(0., 5.)
    t = np.arange(3.)
    data = pc.data().from_dict({'x': x, 'y': x, 'z': x, 'delta_time': 10. + x})
    z0 = pc.grid.data().from_dict({'x': x, 'y': x, 'z0': np.ones((5, 5))})
    dz = pc.grid.data().from_dict({'x': x, 'y': x, 't': t, 'dz': np.ones((5, 5, 3))})
    mask_3d = pc.grid.data().from_dict({'x': x, 'y': x, 't': t, 'z': np.ones((5, 5, 3))})
    return {'data': data, 'timing': {'fit': 1.}, 'RMS': {'data': 0.5}, 'E_RMS': {'d2z0_dx2': 1.},
            'm': {'bias': {}, 'z0': z0, 'dz': dz},
            'grids': {'z0': SimpleNamespace(mask=np.ones((5, 5))), 'dz': SimpleNamespace(mask_3d=mask_3d)}}


@pytest.mark.parametrize('write_data', [True, False])
def test_save_fit_to_file(tmp_path, write_data):
    filename = str(tmp_path / 'E0_N0.h5')
    save_fit_to_file(fit_result(), filename, write_data=write_data)
    with h5py.File(filename, 'r') as h5f:
        assert ('data' in h5f) is write_data
        # everything else is written either way, including the times that
        # come from the data
        assert {'meta', 'RMS', 'E_RMS', 'z0', 'dz'} <= set(h5f)
        assert h5f['meta'].attrs['first_delta_time'] == 10.
        assert h5f['meta'].attrs['last_delta_time'] == 14.
        if write_data:
            assert np.array_equal(h5f['data/z'][:], np.arange(0., 5.))
