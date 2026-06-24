#! /usr/bin/env python

def infer_dzdt_lags(t_res, time_span):
    """
    Infer the set of dz/dt lags (in units of t_res) to mosaic/tile.

    Includes the native time-resolution step (lag=1), the quarter- and
    half-year steps if t_res is finer than quarterly, then one step per
    whole year out to the time span.

    Parameters
    ----------
    t_res : float
        time resolution, in years, of the dz grids (the 'dt' component of
        --grid_spacing).
    time_span : iterable of two floats
        [first_year, last_year] (the --time_span / -t argument).

    Returns
    -------
    list of int
        sorted, deduplicated lags.
    """

    max_years = int(time_span[1] - time_span[0])

    lags = {1}
    if t_res < 0.25:
        # finer than quarterly: also include the quarter- and half-year steps
        lags.add(round(0.25 / t_res))
        lags.add(round(0.5 / t_res))
    for year in range(1, max_years + 1):
        lags.add(round(year / t_res))

    return sorted(lags)
