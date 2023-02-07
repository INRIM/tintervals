
from datetime import datetime, timezone
import numpy as np
# import julian
# fast data import from iso format
import ciso8601 as ciso




# K+K -> ISO 	   : string manipulation
# ISO -> datetime  : ciso.parse_datetime(), naive=system
# datetime -> ISO  : datetime.isoformat + some manipulations
# datetime -> epoch: datetime.timestamp(), naive=system
# epoch -> datetime: fromtimestamp + timezone.utc to avoid naive date
# epoch -> MJD     : astropy.time.Time
# MJD -> epoch     : astropy.time.Time


# MJD and unix time zero in the other scale
# note that both are aligned with UTC and basically ignore leap seconds 
mjd_epoch_0 = 40587.0 #julian.to_jd(epoch2datetime(0), fmt='mjd')
epoch_mjd_0 = -3506716800.0 #datetime2epoch(julian.from_jd(0, fmt='mjd').replace(tzinfo=timezone.utc))



def kk2iso(s, year_digits='20'):
	"""Convert a dateime string from K+K to ISO format.

	Parameters
	----------
	s : str	
		datetime in K+K format.
	year_digits : str, optional
		digits of year to be prepended to K+K format, by default '20'.

	Returns
	-------
	str
		datetime in ISO format.
	"""
	return year_digits + s.replace("*", " ")


def iso2datetime(s, tzinfo=None):
	"""Convert a datetime string in ISO format to a datime object.

	Parameters
	----------
	s : str
		datetime in ISO format.
	tzinfo : timezone object, optional
		if given, interprets the datetime as naive and replace the tzinfo of the result, by default None.

	Returns
	-------
	datetime
		datetime object.
	"""
	
	if tzinfo:
		return ciso.parse_datetime_as_naive(s).replace(tzinfo=tzinfo)
	else:
		return ciso.parse_datetime(s)



# inspiration
# https://stackoverflow.com/questions/2150739/iso-time-iso-8601-in-python

def datetime2iso(d):
	"""Convert a datetime object in a datetime string in ISO format.

	Parameters
	----------
	d : datetime
		datetime object

	Returns
	-------
	str
		datetime in ISO format

	Notes
	-----
	Result is rounded to the nearest integer second.
	The ouput uses 'Z' for UTC datetimes.

	"""
	return d.replace(microsecond=0).isoformat().replace('+00:00','Z')


def datetime2epoch(d):
	"""Convert a datetime object to seconds from the epoch.

	Parameters
	----------
	d : datetime
		datetime object.

	Returns
	-------
	float
		seconds from the epoch.
	"""
	return d.timestamp()


def epoch2datetime(t):
	"""Convert seconds from the epoch to a datetime object.

	Parameters
	----------
	t : float
		seconds from the epoch.

	Returns
	-------
	datetime
		datetime object.
	"""
	return datetime.fromtimestamp(t, tz=timezone.utc)








# def mjd_from_epoch(epoch):
# 	return Time(epoch, format='unix').to_value('mjd') 

# def epoch_from_mjd(mjd):
# 	return Time(mjd, format='mjd', scale='utc').to_value('unix') 
	


def kk2epoch(s, year_digits='20'):
	"""Convert a datetime string in K+K format to seconds from the epoch.

	Parameters
	----------
	s : str
		datetime in K+K format.
	year_digits : str, optional
		digits of year to be prepended to K+K format, by default '20'

	Returns
	-------
	float
		seconds from the epoch.
	"""
	return datetime2epoch(ciso.parse_datetime_as_naive(kk2iso(s, year_digits)))
	
def iso2epoch(s):
	"""Convert a datetime string in ISO format to seconds from the epoch.

	Parameters
	----------
	s : str
		datetime in ISO format.
	
	Returns
	-------
	float
		seconds from the epoch.
	"""
	return datetime2epoch(iso2datetime(s))


def epoch2iso(t):
	"""Convert seconds from the epoch to a datetime string in ISO format.

	Parameters
	----------
	t : float
		seconds from the epoch.
	
	Returns
	-------
	str
		datetime in ISO format.
	"""
	return datetime2iso(epoch2datetime(round(t)))



def epoch2mjd(t):
	"""Convert seconds from the epoch to MJD.

	Parameters
	----------
	t : float or ndarray
		seconds from the epoch.

	Returns
	-------
	float or ndarray
		MJD.
	"""
	return t/86400. + mjd_epoch_0
	#return julian.to_jd(epoch2datetime(t), fmt='mjd')



def mjd2epoch(d):
	"""Convert MJD to seconds from the epoch.
	
	Parameters
	----------
	d : float or ndarray
		MJD.

	Returns
	-------
	float or ndarray
		seconds from the epoch.
	"""
	return d*86400 + epoch_mjd_0
	#return datetime2epoch(julian.from_jd(d, fmt='mjd').replace(tzinfo=timezone.utc))




def mjd2iso(mjd):
	"""Convert MJD to a datetime string in ISO format.

	Parameters
	----------
	mjd : float
		MJD.
	
	Returns
	-------
	str
		datetime in ISO format.
	"""
	return epoch2iso(epoch_from_mjd(mjd))
	

def iso2mjd(s):
	"""Convert a datetime string in ISO format to MJD.

	Parameters
	----------
	s : str
		datetime in ISO format.
	
	Returns
	-------
	float
		MJD.
	"""
	return epoch2mjd(iso2epoch(s))



def datetime2mjd(d):
	"""Convert a datetime object to seconds from the epoch.

	Parameters
	----------
	d : datetime
		datetime object.

	Returns
	-------
	float
		MJD.
	"""
	return epoch2mjd(datetime2epoch(d))


def mjd2datetime(d):
	"""Convert seconds from the epoch to a datetime object.

	Parameters
	----------
	d : float
		MJD.

	Returns
	-------
	datetime
		datetime object.
	"""
	return epoch2datetime(mjd2epoch(d))

def kk2mjd(s, year_digits='20'):
	"""Convert a datetime string in K+K format to MJD.

	Parameters
	----------
	s : str
		datetime in K+K format.
	year_digits : str, optional
		digits of year to be prepended to K+K format, by default '20'

	Returns
	-------
	float
		MJD.
	"""
	return datetime2mjd(ciso.parse_datetime_as_naive(kk2iso(s, year_digits)))



def myvectorize(f):
# https://stackoverflow.com/questions/32766210/making-a-vectorized-numpy-function-behave-like-a-ufunc?rq=1
	vf = np.vectorize(f)

	def newfunc(*args, **kwargs):
		return vf(*args, **kwargs)[()]

	newfunc.__doc__ = """Vectorized version of {}.
	""".format(f.__name__)
	newfunc.__doc__ +=  f.__doc__.replace('str\n', 'str or ndarray\n').replace('float', 'float or ndarray')

	return newfunc


def mjd_from_epoch(t):
	"""Same as epoch2mjd (already vectorized)."""
	return epoch2mjd(t)

def epoch_from_mjd(mjd):
	"""Same as mjd2epoch (already vectorized)."""
	return mjd2epoch(mjd)

iso_from_epoch = myvectorize(epoch2iso)

iso_from_mjd = myvectorize(mjd2iso)