"""
Tests for the lazy attribute resolution in ATL1415/__init__.py.

Seven submodules define a function of their own name, and the eager __init__
that preceded 0c6ea35 bound ATL1415.<name> to that FUNCTION.  The lazy
__getattr__ returned the module instead, which broke every queue builder
(ATL1415.make_slurm_file(...) -> "'module' object is not callable") and the
write2nc scripts, and nothing noticed until the mosaic step ran.

Each case runs in a fresh interpreter: lazy resolution depends on what has
already been imported, so one test's imports must not satisfy another's.
A case is skipped only when a third-party dependency is missing from the
environment running the tests (the solver chain is not installed everywhere).
"""
import subprocess
import sys

import pytest

FUNCTION_OVER_MODULE = [
    'make_slurm_file',
    'make_nc_projection_variable',
    'make_tile_stats_group',
    'read_ATL11',
    'assign_firn_variable',
    'SMB_corr_from_grid',
    'ATL11_to_ATL15',
]


def _run(code):
    result = subprocess.run([sys.executable, '-c', code],
                            capture_output=True, text=True)
    if result.returncode != 0 and 'ModuleNotFoundError' in result.stderr \
            and "No module named 'ATL1415" not in result.stderr:
        pytest.skip('dependency missing: ' + result.stderr.strip().splitlines()[-1])
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize('name', FUNCTION_OVER_MODULE)
def test_attribute_is_the_function(name):
    _run(f"""
import inspect, ATL1415
value = ATL1415.{name}
assert inspect.isfunction(value), f'ATL1415.{name} is {{value!r}}'
assert value.__name__ == '{name}'
""")


@pytest.mark.parametrize('name', FUNCTION_OVER_MODULE)
def test_from_import_is_the_function(name):
    _run(f"""
import inspect
from ATL1415 import {name} as value
assert inspect.isfunction(value), f'from ATL1415 import {name} gave {{value!r}}'
""")


def test_write2nc_import_line():
    """The exact import in ATL14_write2nc.py and ATL15_write2nc.py."""
    _run("""
import inspect
from ATL1415 import ATL1415_attrs_meta, make_nc_projection_variable, make_tile_stats_group
assert inspect.ismodule(ATL1415_attrs_meta)
assert inspect.isfunction(make_nc_projection_variable)
assert inspect.isfunction(make_tile_stats_group)
""")


def test_submodule_still_importable_by_path():
    _run("""
import inspect
import ATL1415.make_slurm_file as module
assert inspect.ismodule(module)
""")
