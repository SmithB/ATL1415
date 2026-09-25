"""
ATL11_to_ATL15 --solver (docs/plan_cholmod_fit.sh, QC2: opt-in).  SPQR is the
default; --solver=cholmod is checked at parse time, before any data are read,
because LSsurf's smooth_fit would silently ignore a solver it does not know.
"""
import sys

import pytest

from ATL1415.ATL11_to_ATL15 import parse_args

BASE = ['ATL11_to_ATL15.py', '--ATL11_index', '/tmp/index.h5', '--xy0', '0', '0']


def test_spqr_is_the_default():
    assert parse_args(list(BASE)).solver == 'spqr'


def test_cholmod_is_opt_in():
    assert parse_args(BASE + ['--solver=cholmod']).solver == 'cholmod'


def test_unknown_solver_is_refused():
    with pytest.raises(SystemExit):
        parse_args(BASE + ['--solver=lsqr'])


@pytest.mark.parametrize('module, message', [('sksparse.cholmod', 'scikit-sparse'),
                                             ('LSsurf.ls_solvers', 'ls_solvers')])
def test_cholmod_without_its_parts_stops_loudly(monkeypatch, module, message):
    monkeypatch.setitem(sys.modules, module, None)
    with pytest.raises(ImportError, match=message):
        parse_args(BASE + ['--solver=cholmod'])
    # spqr needs neither
    assert parse_args(list(BASE)).solver == 'spqr'
