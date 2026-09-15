#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 14:08:42 2019.

@author: ben
"""

# This module used to carry its own copy of assign_firn_variable(), identical to
# the one in SMBcorr apart from whitespace and a stale np.NaN (removed in numpy
# 2.0).  It now delegates, so the two cannot drift.  The import is deferred
# because SMBcorr is an optional dependency -- see the [firn] extra.


def assign_firn_variable(*args, **kwargs):
    """
    Assign a firn-height variable to a data structure.

    Thin wrapper around SMBcorr.assign_firn_variable; see that function for the
    signature.  SMBcorr is optional, so it is imported only when called.
    """
    try:
        from SMBcorr.assign_firn_variable import assign_firn_variable as _assign_firn_variable
    except ImportError as exc:
        raise ImportError(
            'assign_firn_variable requires SMBcorr, which is not installed.  '
            'Install it with "pip install ATL1415[firn]".') from exc
    return _assign_firn_variable(*args, **kwargs)
