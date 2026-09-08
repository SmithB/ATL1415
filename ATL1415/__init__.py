"""
ATL1415: ICESat-2 ATL14/ATL15 gridded land-ice height.

SUBMODULES ARE IMPORTED LAZILY.  `import ATL1415` used to pull in the whole
package, and the first line of it reached ATL11_to_ATL15 -> LSsurf -> sparseqr,
a CFFI extension built against SuiteSparse.  That made the solver's compiled
toolchain a hard requirement for every entry point, including the ones that
never solve anything: setup_ATL1415_region.py and make_ATL1415_queue.py only
need ATL1415.paths, but could not run at all without a full solver environment.
On MAAP that is the difference between composing an args file in the ADE's
stock environment and needing a bespoke conda env for it.

Attribute access resolves on demand (PEP 562), so `ATL1415.read_ATL11`,
`from ATL1415 import ATL1415_attrs_meta` and a bare `ATL1415.<name>`
re-exported from any submodule all still work; they simply import what they
need when they are first touched.  `from ATL1415.paths import ...` never
triggers any of it.
"""

import importlib as _importlib

# Searched in this order for a re-exported name, so the cheap modules resolve
# first and the solver chain is reached only by something that genuinely needs
# it.  ATL11_to_ATL15 is last on purpose: it is the expensive one.
_LAZY_SUBMODULES = (
    'paths',
    'tile_names',
    'lags',
    'make_nc_projection_variable',
    'make_tile_stats_group',
    'ATL1415_attrs_meta',
    'make_slurm_file',
    'read_ATL11',
    'assign_firn_variable',
    'SMB_corr_from_grid',
    'ATL11_to_ATL15',
)


def _suppress_pytmd_anonymous_reads():
    """
    Stop pyTMD's AWS_NO_SIGN_REQUEST from making every GDAL read anonymous.

    pyTMD (v3.0.9, pyTMD/io/__init__.py) sets AWS_NO_SIGN_REQUEST=YES in
    os.environ when pyTMD.io is imported, so its own reads of the public
    s3://pytmd stores are unsigned.  The variable is PROCESS-WIDE and GDAL
    honours it for every /vsis3 read -- including the masks, geoid and tide
    masks read from s3://maap-ops-workspace, which then come back HTTP 403 and
    surface as ogr.Open() returning None or from_geotif() yielding an object
    with no .z.

    This used to be a bare os.environ.pop() at the bottom of this file, which
    worked only because the eager `from .ATL11_to_ATL15 import *` above it had
    already triggered pyTMD -- so the pop was guaranteed to run last.  Under
    lazy imports that ordering is gone: pyTMD may be imported at any later
    point, after the pop, and the bug would come back silently.

    gdal.SetConfigOption() does not depend on ordering at all -- a config
    option set programmatically takes precedence over the environment
    variable, verified 2026-09-08 by reading a private object with
    AWS_NO_SIGN_REQUEST=YES still in the environment.  The pop is kept as well,
    for any non-GDAL consumer that reads the variable directly.

    Our own tide reads are unaffected: ATL1415/tides.py passes anon= to s3fs
    explicitly, and s3fs ignores both the GDAL option and the variable.
    """
    import os
    os.environ.pop('AWS_NO_SIGN_REQUEST', None)
    try:
        from osgeo import gdal
    except ImportError:
        # A GDAL-less environment can still use the path helpers; there is no
        # /vsis3 read to protect.
        return
    gdal.SetConfigOption('AWS_NO_SIGN_REQUEST', 'NO')


def __getattr__(name):
    """Resolve a submodule, or a name re-exported from one, on first access."""
    if name in _LAZY_SUBMODULES:
        # Asked for by name, so an ImportError is the honest answer and
        # propagates: the caller wanted this module specifically.
        module = _importlib.import_module('.' + name, __name__)
        _suppress_pytmd_anonymous_reads()
        globals()[name] = module
        return module

    # A submodule that will not import must not take the whole lookup down
    # with it.  In a light environment -- one without the solver toolchain, or
    # without an optional dependency like timescale -- some of these raise
    # ImportError, and a name that lives in a module which DID import should
    # still resolve.  Failures are collected rather than swallowed, and
    # reported if the name is not found anywhere, so "no such attribute" is
    # never confused with "the module holding it could not be imported".
    unavailable = []
    for submodule in _LAZY_SUBMODULES:
        try:
            module = _importlib.import_module('.' + submodule, __name__)
        except ImportError as exc:
            unavailable.append(f'{submodule} ({exc})')
            continue
        _suppress_pytmd_anonymous_reads()
        if hasattr(module, name):
            value = getattr(module, name)
            globals()[name] = value
            return value

    message = f'module {__name__!r} has no attribute {name!r}'
    if unavailable:
        message += ('.  These submodules could not be imported and were not '
                    'searched: ' + '; '.join(unavailable))
    raise AttributeError(message)


def __dir__():
    return sorted(set(globals()) | set(_LAZY_SUBMODULES))


# Applied eagerly too, so a process that imports ATL1415 and goes straight to
# GDAL is covered even before any submodule is touched.
_suppress_pytmd_anonymous_reads()
