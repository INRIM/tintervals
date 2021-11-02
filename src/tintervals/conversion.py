
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


# MJD and unix time zero in the other scale
# note that both are aligned with UTC and basically ignore leap seconds 
mjd_epoch_0 = 40587.0 #julian.to_jd(epoch2datetime(0), fmt='mjd')
epoch_mjd_0 = -3506716800.0 #datetime2epoch(julian.from_jd(0, fmt='mjd').replace(tzinfo=timezone.utc))

def epoch2mjd(t):
	return t/86400. + mjd_epoch_0
	#return julian.to_jd(epoch2datetime(t), fmt='mjd')

mjd_from_epoch = epoch2mjd

def mjd2epoch(d):
	return d*86400 + epoch_mjd_0
	#return datetime2epoch(julian.from_jd(d, fmt='mjd').replace(tzinfo=timezone.utc))

epoch_from_mjd = mjd2epoch



# def mjd_from_epoch(epoch):
# 	return Time(epoch, format='unix').to_value('mjd') 

# def epoch_from_mjd(mjd):
# 	return Time(mjd, format='mjd', scale='utc').to_value('unix') 
	


def kk2epoch(s, year_digits='20'):
	return datetime2epoch(ciso.parse_datetime_as_naive(kk2iso(s, year_digits)))
	
def iso2epoch(s):
	return datetime2epoch(iso2datetime(s))


def epoch2iso(t):
	return datetime2iso(epoch2datetime(round(t)))

iso_from_epoch = myvectorize(epoch2iso)

	

def mjd2iso(mjd):
	return epoch2iso(epoch_from_mjd(mjd))
	
iso_from_mjd = myvectorize(mjd2iso)






