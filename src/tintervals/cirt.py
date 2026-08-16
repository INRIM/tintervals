from datetime import datetime, timedelta, timezone

import numpy as np

import tintervals as ti


def _next_month(year, month):
    if month == 12:
        return year + 1, 1
    return year, month + 1


def cirt2mjd(year, month):
    """Convert the BIPM Circular T year-month in MJD start-stop.

    Parameters
    ----------
    year : int
            Circular-T year
    month : int
            Circular-T month

    Returns
    -------
    1d array
            Circular-T start and stop as MJD
    """
    # year + month -> circular T
    next_year, next_month = _next_month(year, month)

    first = datetime(year, month, 1, tzinfo=timezone.utc)
    next_first = datetime(next_year, next_month, 1, tzinfo=timezone.utc)

    s_mjd = ti.datetime2mjd(first - timedelta(days=1))
    e_mjd = ti.datetime2mjd(next_first - timedelta(days=1))

    s_mod = (s_mjd + 1) % 5
    e_mod = (e_mjd + 1) % 5

    return np.array([s_mjd - s_mod, e_mjd - e_mod])


def mjd2cirt(mjd):
    """Return the Circular T (year/month) corresponding to the given MJD.

    Parameters
    ----------
    mjd : float
            MJD date

    Returns
    -------
    year, month
            corresponding Circular-T containing the MJD.
    """
    mod = (mjd + 1) % 5
    dt = ti.mjd2datetime(mjd - mod + 5).replace(day=1)
    return dt.year, dt.month


def cirt_vals(start, stop=None):
    """Calculate MJD intervals corresponding to the BIPM Circular Ts between start and stop.


    Parameters
    ----------
    start : float or datetime
            starting MJD or date
    stop : float or datetime, optional
            stop MJD or date, by default (present day).

    Returns
    -------
    2d array
            start, stop intervals in MJD for each Circular T
    """
    if isinstance(start, datetime):
        start = ti.datetime2mjd(start)

    if stop is None:
        stop = datetime.now(timezone.utc)

    if isinstance(stop, datetime):
        stop = ti.datetime2mjd(stop)

    year, month = mjd2cirt(start)
    stop_year, stop_month = mjd2cirt(stop)

    intervals = []

    while (year, month) <= (stop_year, stop_month):
        intervals.append(cirt2mjd(year, month))
        year, month = _next_month(year, month)

    return np.array(intervals)


def yymm2cirt_number(yymm):
    """Return the Circular T number corresponding to YYMM."""
    year = 2000 + int(yymm[:2])
    month = int(yymm[2:])

    return 312 + 12 * (year - 2014) + month
