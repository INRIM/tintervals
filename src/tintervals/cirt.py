from datetime import datetime, timezone
from dateutil import relativedelta, rrule
import tintervals as ti
import numpy as np


def cirt2mjd(year, month):
	""" Convert the BIPM Circular T year-month in MJD start-stop.

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
	x = datetime(year, month, 1, tzinfo=timezone.utc) 
	s = x + relativedelta.relativedelta(days=-1)
	e = x + relativedelta.relativedelta(months=1) + relativedelta.relativedelta(days=-1)
	


	s_mjd = ti.epoch2mjd(ti.datetime2epoch(s))
	e_mjd = ti.epoch2mjd(ti.datetime2epoch(e))

	s_mod = (s_mjd + 1)%5
	e_mod = (e_mjd + 1)%5

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
	mod = (mjd + 1)%5
	dt = ti.mjd2datetime(mjd -mod + 5).replace(day=1)
	return dt.year, dt.month

def cirtvals(start, stop=None):
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
	# if isinstance(start, datetime):
	# 	start_mjd = ti.datetime2mjd(start)
	# else:
	# 	start_mjd = start
	# 	start = ti.mjd2datetime(start)
		

	# if stop is None:
	# 	stop = datetime.utcnow()
		
	# if isinstance(stop, datetime):
	# 	stop_mjd = ti.datetime2mjd(stop)
	# else:
	# 	stop_mjd = stop
	# 	stop = ti.mjd2datetime(stop)

	# # for rrule start and stop should be extened to get all circularT intervals, inclusive
	# dtstart = start.replace(day=1).replace(tzinfo=None) #rrule does not like tzinfo
	# until = (stop + relativedelta.relativedelta(months=2)).replace(tzinfo=None)

	# months = list(rrule.rrule(rrule.MONTHLY, bymonthday=1, dtstart=dtstart, until=until))

	# # convert each 1st of month to start of circularT
	# # timezone.utc assure proper conversion at integer mjds
	# s_mjd = np.array([ti.epoch2mjd(ti.datetime2epoch(s.replace(tzinfo=timezone.utc))) for s in months])
	# s_mod = (s_mjd + 1)%5	

	# breaks = s_mjd-s_mod

	# vals = np.column_stack((breaks[:-1], breaks[1:]))

	# # mask result as the extension done before may have been too extreme
	# mask = (vals[:,1] > start_mjd) & (vals[:,0] < stop_mjd)
	# return vals[mask]

	if isinstance(start, datetime):
		start = ti.datetime2mjd(start)
	
	if stop is None:
		stop = datetime.utcnow()
		
	if isinstance(stop, datetime):
		stop = ti.datetime2mjd(stop)
	
	starty, startm = mjd2cirt(start)
	stopy, stopm = mjd2cirt(stop)

	dtstart = datetime(starty,startm,1)
	until = datetime(stopy,stopm,1) + relativedelta.relativedelta(months=1)
	months = list(rrule.rrule(rrule.MONTHLY, bymonthday=1, dtstart=dtstart, until=until))

	# # convert each 1st of month to start of circularT
	# # timezone.utc assure proper conversion at integer mjds
	s_mjd = np.array([ti.epoch2mjd(ti.datetime2epoch(s.replace(tzinfo=timezone.utc))) for s in months])  # first of every month
	s_mjd -= 1  # last of every month
	s_mod = (s_mjd + 1)%5	

	breaks = s_mjd-s_mod

	vals = np.column_stack((breaks[:-1], breaks[1:]))
	return vals
	
