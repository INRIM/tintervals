
"""
Functions to work with time intervals in the form of start-stop ranges.

"""


import numpy as np
import scipy.interpolate
import warnings


def array2intervals(t, tgap=1., tblock=0.):
	"""
	Calculate from an array of timetags t a 2-d array in the form (start,stop),
	including gaps > tgap and removing intervals < tblock
	
	Parameters
	----------
	t       : 1d array of timetags (assumed sorted)
	tgap    : resulting intervals will mark gaps > tgap
	tblock  : intervals shorter than tblock will be removed
	
	Returns
	-------	
	out     : 2-d array of start, stop intervals
	
	Example
	-------
	
	>>> array2intervals(np.array([1,2,3,7,8]))                                     
	array([[1, 3],[7, 8]])
	
	"""
	t = np.atleast_1d(t)
	submask = (t[1:]-t[:-1] > tgap)
	
	mask = np.ones_like(t, dtype=bool)
	mask[1:] = submask
	tstarts = t[mask]

	mask =  np.ones_like(t, dtype=bool)
	mask[:-1] = submask
	tstops = t[mask]

	mask = (tstops - tstarts) >= tblock
	tstarts = tstarts[mask]
	tstops = tstops[mask]

	out = np.column_stack((tstarts, tstops))
	return out
	

def intervals2weights(a, step=1, min=None, max=None, norm=False):
	"""
	Given some interval start-stop return an arange of timetags and an array of weights.
	
	Parameters
	----------
	a       : 2-d array of intervals in the form (start,stop)
	step    : float, default = 1
			  spacing between timetags
	min     : float, default = min of a
	          fix the minimum of the timetags
	max     : float, default = max of a
	          fix the maximum of the timetags 
	norm    : bool, default = False
	          if True normalize the weights
	
	Returns
	-------	
	t	: array of timetags
	w	: weights (positive only between each start,stop)
	        
	Example
	-------
		
	"""
	if min == None:
		min = np.amin(a)

	if max == None:
		max = np.amax(a)	

	start = np.ceil(min/step)*step
	stop = np.floor(max/step)*step
	
	t = np.arange(start,stop+step,step)
	

	res = np.zeros_like(t)
	for x in a:
		res +=(t>x[0]) & (t<x[1])

	if norm:
		res = res/sum(res)

	return t, res




def mix(a, b):
	"""
	Intersection of two  2-d arrays in the form (start,stop),
	
	Parameters
	----------
	a       : 2-d array of intervals in the form (start,stop)
	b       : 2-d array of intervals in the form (start,stop)
	
	Returns
	-------	
	out     : intersection of a and b

	"""
	res = []
	for x in a:
		for y in b:
			start = max(x[0], y[0])
			stop = min(x[1], y[1])
			if stop > start:
				res += [[start, stop]]
				
	return np.array(res)			


# 
def split(a, base = 10.):
	"""
	Given some interval start-stop return all the included intervals with start and stop multiple of base
	
	Parameters
	----------
	a       : 2-d array of intervals in the form (start,stop)
	base    : float, 
	
	Returns
	-------	
	out     : 2-d array of intervals in the form (start,stop),
	        with start and stop every multiple of base

	Example
	-------
	
	>>> split(np.array([[1,31]]), 10)                                              
	array([[10., 20.],[20., 30.]])

	"""
	res = []
	for x in a:
		start = np.ceil(x[0] / base) * base
		stop = np.floor(x[1] / base) * base
		tags = np.arange(start, stop+base, base)
		res += list(np.column_stack((tags[:-1], tags[1:])))
	
	return np.array(res)



def csaverage(f, tistart, tistop, tostart, tostop):
	"""
	Perform the average of data given in (tistart,tistop) ranges 
	in the intervals (tostart, tostop).
	It uses linear interpolation of the data and can work with gaps in the data.
	
	The average is made by a "cumsum" algorithm
	
	Parameters
	----------
	f       : input data to be averaged
	tistart : start time for each measurement
	tistop  : stop time for each measuremnt
	tostart : start time for the output average
	tostop  : stop time for the ouput average
	
	Returns
	-------	
	res     : averaged array
	num     : number of points averaged for each result
	"""
	
	# lets look for every tstart that is not a tstop
	tiadd = np.array([t for t in tistart if t not in tistop])
	
	#print len(tiadd)
	
	# here are all the timestamps I care for
	t2 = np.append(tistop, tiadd)
	t3 = np.sort(t2)
	
	
	# in the cumsum nothing should be added until we get to a timetag in tiadd
	f2 = np.append(f, np.zeros(len(tiadd)))
	
	# now I sort f2 temporarly and perform the cumsum
	f3 = f2[t2.argsort()]
	cs = np.cumsum(f3)
	
	# lets do the same for the number of averaged points
	datapoints = np.ones(len(tistop))
	datapoints2 = np.append(datapoints, np.zeros(len(tiadd)))
	datapoints3 = datapoints2[t2.argsort()]
	ct = np.cumsum(datapoints3)
	
	
	
	# interpolate linearly
	csf = scipy.interpolate.interp1d(t3, cs, kind='linear')
	ctf = scipy.interpolate.interp1d(t3, ct, kind='linear')



	# now for the averaging	
	num = csf(tostop)-csf(tostart)
	den = ctf(tostop)-ctf(tostart)

	return num/den, den



