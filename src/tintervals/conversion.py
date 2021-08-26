import time
from datetime import datetime
from astropy.time import Time
import pytz


# fast data import from iso format
import ciso8601 as ciso



## convert the notation of the K+K counter file to seconds since the epoch
#def kk2epoch(s, summertime=True):
#	"""convert the K+K data format (as string) in Unix timestamp. Take as input the K+K timetag and if summerime is in effect.	
#	"""
#	date = datetime.strptime(s, "%y%m%d*%H%M%S.%f")
#	epoch = datetime.fromtimestamp(0) # in localtime, as date
#	return (date - epoch).total_seconds() -3600.*summertime #the calculation above is not automatically aware of summertime


# set local timezone and local epoch!
cet = pytz.timezone('CET')
epoch = cet.localize(datetime.fromtimestamp(0))

# convert the notation of the K+K counter file to seconds since the epoch
# using astropy
#def kk2epoch(s):
#	# if spaces replace them with the *
#	s = s.replace(" ", "*")
#	date = datetime.strptime(s, "%y%m%d*%H%M%S.%f")
#	
#	# lets try to use ciso8601 to get faster parsing of huge comb file
##	# i neet to add the first 2 digits of the year and replace the * with a space (or a T)
##	s = "20" + s.replace("*", " ")
##	date = ciso.parse_datetime(s)
#	
#	# set local timezone
#	date_tzaware = cet.localize(date)
#	date_utc = date_tzaware.astimezone(pytz.utc)
#	
#	out = Time(date_utc, scale='utc').unix
#	return out

# version about 10 time faster
def kk2epoch(s):
	# lets try to use ciso8601 to get faster parsing of huge comb file
#	# i neet to add the first 2 digits of the year and replace the * with a space (or a T)
	s = "20" + s.replace("*", " ")
	date = ciso.parse_datetime(s)
	
	# set local timezone
	date_tzaware = cet.localize(date)
	#date_utc = date_tzaware.astimezone(pytz.utc)
	
	# Time is slow, so lets do the math by hand
	out = (date_tzaware - epoch).total_seconds()
	return out



def epoch2mjd(epoch):
	return Time(epoch, format='unix').mjd

def mjd2epoch(mjd):
	return Time(mjd, format='mjd', scale='utc').unix
