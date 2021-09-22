
from datetime import datetime, timezone
from astropy.time import Time

import numpy as np

# fast data import from iso format
import ciso8601 as ciso




# K+K -> ISO 	   : string manipulation
# ISO -> datetime  : ciso.parse_datetime(), naive=system
# datetime -> ISO  : datetime.isoformat + some manipulations
# datetime -> epoch: datetime.timestamp(), naive=system
# epoch -> datetime: fromtimestamp + timezone.utc to avoid naive date
# epoch -> MJD     : astropy.time.Time
# MJD -> epoch     : astropy.time.Time



def myvectorize(f):
# https://stackoverflow.com/questions/32766210/making-a-vectorized-numpy-function-behave-like-a-ufunc?rq=1
    vf = np.vectorize(f)

    def newfunc(*args, **kwargs):
        return vf(*args, **kwargs)[()]
    return newfunc



def kk2iso(s, year_digits='20'):
	return year_digits + s.replace("*", " ")


def iso2datetime(s, tzinfo=None):
	
	if tzinfo:
		return ciso.parse_datetime_as_naive(s).replace(tzinfo=tzinfo)
	else:
		return ciso.parse_datetime(s)



# inspiration
# https://stackoverflow.com/questions/2150739/iso-time-iso-8601-in-python

def datetime2iso(d):
	return d.replace(microsecond=0).isoformat().replace('+00:00','Z')


def datetime2epoch(d):
	return d.timestamp()


def epoch2datetime(t):
	return datetime.fromtimestamp(t, tz=timezone.utc)


def mjd_from_epoch(epoch):
	return Time(epoch, format='unix').to_value('mjd') 

def epoch_from_mjd(mjd):
	return Time(mjd, format='mjd', scale='utc').to_value('unix') 
	

def kk2epoch(s, year_digits='20'):
	return datetime2epoch(ciso.parse_datetime_as_naive(kk2iso(s, year_digits)))
	
def iso2epoch(s):
	return datetime2epoch(iso2datetime(s))

@myvectorize
def iso_from_epoch(t):
	return datetime2iso(epoch2datetime(t))
	
@myvectorize
def iso_from_mjd(mjd):
	return iso_from_epoch(epoch_from_mjd(mjd))
	





