
import numpy as np
import scipy.interpolate
#import warnings


def array2intervals(t, tgap=1., tblock=0.):
	"""
	Convert from timetags to time intervals.

	Calculate from an array of timetags t a 2-d array in the form (start,stop),
	including gaps > tgap and removing intervals < tblock.
	
	Parameters
	----------
	t       : 1-d array
		timetags (assumed sorted).
	tgap    : float
		resulting intervals will mark gaps > tgap, by default 1.
	tblock  : float
		intervals shorter than tblock will be removed, by default 0.
	
	Returns
	-------	
	2-d array
		start, stop intervals
	
	Examples
	--------
	
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
	Convert from time intervals to timetags and weights.
	
	Given some interval start-stop return an arange of timetags and an array of weights.
	
	Parameters
	----------
	a       : 2-d array
		 	  intervals in the form (start,stop).
	step    : float
			  spacing between timetags.
	min     : float
	          fix the minimum of the timetags, by default min of a.
	max     : float
	          fix the maximum of the timetags, by default max of a.
	norm    : bool
	          if True normalize the weights, by default False
	
	Returns
	-------	
	1-d array	
		timetags.
	1-d array
		weights (positive only between each start,stop).
	        
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




def intersect(a, b):
	"""
	Calculate the intersection of two arrays of time intervals.

	Parameters
	----------
	a : 2-d array
		intervals in the form (start,stop).
	b : 2-d array
		intervals in the form (start,stop).
	
	Returns
	-------	
	2-d array
		intersection of a and b

	Examples
	--------
	>>> intersect([[1,3],[5,7]],[[2,6]])                                           
	array([[2, 3], [5, 6]])

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
	Split time intervals in regular steps.

	Given some interval start-stop return all the included intervals with start and stop multiple of base
	
	Parameters
	----------
	a       : 2-d array
			intervals in the form (start,stop).
	base    : float 
			base interval (generated intervals will be mutiple of base).
	Returns
	-------	
	2-d array
		intervals in the form (start,stop),
	    with start and stop every multiple of base

	Examples
	--------
	
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


def regvals(tstart, tstop, base=1., offset=0., extend=True):
	"""Generate regular time intervals.

	Parameters
	----------
	tstart : float
		starting time.
	tstop : float
		stopping time.
	base : float, optional
		base interval (generated intervals will be mutiple of base), by default 1.
	offset : float, optional
		time offset in the generated intervals, by default 0.
	extend : bool, optional
		if True, extend the intervals to include the extremes, by default True.

	Returns
	-------
	2-d array
		 regular start/stop intervals
	"""
	if extend:
		tstart = np.floor((tstart-offset)/base)*base + offset
		tstop = np.ceil((tstop-offset)/base)*base + offset
	else:
		tstart = np.ceil((tstart-offset)/base)*base + offset
		tstop = np.floor((tstop-offset)/base)*base + offset
	
	tags = np.arange(tstart, tstop+base, base)
	return np.column_stack((tags[:-1], tags[1:]))




def raverage(data, t, base=1., offset=0., step=None, gap_ext=10):
	"""Average data in regular intervals.

	The algorithm perform the average by reshaping data.
	
	Parameters
	----------
	data : ndarray
		data to be averaged (along the first axis).
	t : 1-d array
		timetags of the data to be averaged. The length should match the first dimension of data.
		Timetags are assumed sorted and it should be a subset of an arange with step = base
	base : float, optional
		base intervals (average will be at the mutiple of base, if offset=0), by default 1.
	offset : float, optional
		time offset for the averages (averages will be at mutiple of base + offset), by default 0.
	step : float, optional
		spacing between timetags, by default the median of the differences of t.
	gap_ext : int, optional
		algorithm will automatically skip gaps > gap_ext*base, by default 10.
		Changing this value may change the efficiency of the algorithm dependning on the input data.
		(Lower values will split the data more, but will have less gaps to fill.)

	Returns
	-------
	2-d array
		start/stop intervals of the averages.
	ndarray
		averaged data.
	1-d array
		number of averaged point for each intervals.

	Notes
	-----
	Timetags are included in the start intervals but excluded in the stop intervals.


	"""
	step = step if step else np.median(np.diff(t))

	if not (base/step).is_integer():
		raise ValueError(f'Not integer ratio between base {base} and step {step}.')
	
	N = int(base/step)

	invals = array2intervals(t, tgap=gap_ext*base)

	vals = []
	res = []
	count = []

	for a, b in invals:
		# both >=, as there is at least a big gap and I want to keep both the first and last timetag
		mask = (t>=a) & (t<=b)


		# b+step ccount the convention about timetag inclusion (see below)
		outvals = regvals(a, b+step, base=base, offset=offset, extend=True)

		ext_t = np.arange(outvals[0,0], outvals[-1,-1], step)
		isin = np.isin(ext_t, t)

		shape = list(data.shape)
		shape[0] = len(ext_t)

		ext_data = np.zeros(shape)
		ext_data[isin] = data[mask]

		#ext_data[:,0] = ext_t

		new_shape = [N, -1] + shape[1:]

		# note that this reshaping corresponds to including the start timetag but excluding the stop timetag
		num = np.sum(ext_data.reshape(new_shape, order='F'), axis=0)
		den = np.sum(isin.reshape((N,-1), order='F'), axis=0)

		vals += [outvals[den>0]]
		res += [(num[den>0].T/den[den>0]).T]
		count += [den[den>0]]
	
	return np.concatenate(vals, axis=0), np.concatenate(res), np.concatenate(count)



def maverage(data, t, intervals):
	"""Average data in given intervals.

	The average is done masking the data.


	Parameters
	----------
	data : ndarray
		data to be averaged (along the first axis).
	t : 1-d array
		timetags of the data to be averaged. The length should match the first dimension of data.
	intervals : 2-d array
		start/stop intervals where average the data.

	Returns
	-------
	2-d array 
		start/stop intervals of the averages.
	ndarray 
		averaged data.
	1-d array
		number of averaged point for each intervals.

	Notes
	-----
	Timetags are included in the start intervals but excluded in the stop intervals.

	"""
	res = []
	count = []
	
	intervals = np.atleast_2d(intervals)


	for a, b in intervals:
		mask = (t >= a) & (t < b)

		num = np.sum(data[mask], axis=0)
		den = np.sum(mask)

		if den > 0:
			res += [num/den]
		else:
			res += [num*0.] # garbage data, will be masked away later but I avoid a warning

		count += [den]

	res = np.stack(res, axis=0)
	count= np.atleast_1d(count)

	return intervals[count>0], res[count>0], count[count>0]




def csaverage(f, ti, to, axis=0):
	"""
	Average data in different time intervals.

	Perform the average of data given in (tistart,tistop) ranges 
	in the intervals (tostart, tostop).
	It uses linear interpolation of the data and can work with gaps in the data.
	
	The average is made by a "cumsum" algorithm
	
	Parameters
	----------
	f  : ndarray
		input data to be averaged.
	ti : 2-d array
		 start/stop time for each measurement.
	to : 2-d array 
		start/stop time for the output average.
	axis : int, optional
		axis of interpolation, by default 0.
	
	
	Returns
	-------	
	ndarray
		averaged array.
	1-d array
		number of points averaged for each result.
	"""
	
	tistart, tistop = np.atleast_2d(ti).T
	tostart, tostop = np.atleast_2d(to).T

	# lets look for every tstart that is not a tstop
	tiadd = np.array([t for t in tistart if t not in tistop])
	
	#print len(tiadd)
	
	# here are all the timestamps I care for
	t2 = np.append(tistop, tiadd)
	t3 = np.sort(t2)
	
	
	# in the cumsum nothing should be added until we get to a timetag in tiadd
	# f2 = np.append(f, np.zeros(len(tiadd)))
	# ndim version
	f = np.atleast_1d(f)
	npad = [(0, 0)] * f.ndim # https://stackoverflow.com/questions/19349410/how-to-pad-with-zeros-a-tensor-along-some-axis-python
	pad_size = len(tiadd)
	npad[axis] = (0, pad_size)
	f2 = np.pad(f, pad_width=npad)



	
	# now I sort f2 temporarly and perform the cumsum
	f3 = f2[t2.argsort()]
	cs = np.cumsum(f3, axis=axis)
	
	# lets do the same for the number of averaged points
	datapoints = np.ones(len(tistop))
	datapoints2 = np.append(datapoints, np.zeros(len(tiadd)))
	datapoints3 = datapoints2[t2.argsort()]
	ct = np.cumsum(datapoints3)
	
	
	
	# interpolate linearly
	csf = scipy.interpolate.interp1d(t3, cs, kind='linear', axis=axis, bounds_error=False, fill_value=(cs[0],cs[-1]))
	ctf = scipy.interpolate.interp1d(t3, ct, kind='linear', bounds_error=False, fill_value=(ct[0],ct[-1]))



	# now for the averaging	
	num = csf(tostop)-csf(tostart)
	den = ctf(tostop)-ctf(tostart)

	# proper (?) broadcasting
	res = np.moveaxis(np.moveaxis(num, axis, -1)/den, -1, axis)

	return res, den



def add_break(vals, brk):
	"""Add a break to given start-stop intervals.

	Parameters
	----------
	vals : 2d array
		Start/stop intervals
	brk : float
		Break to be inserted.

	Returns
	-------
	2d array
		New intervals with added break.

	Examples
	--------

	>>> vals = np.array([[0,1],[1,2],[2,3]])                                    
	>>> add_break(vals, 1.5)                                                    
	array([[0. , 1. ],
		   [1. , 1.5],
		   [1.5, 2. ],
		   [2. , 3. ]])


	"""
	out = []
	for start,stop in vals:
		if (brk > start) and (brk < stop):
			out += [[start, brk],
			[brk, stop]]
		else:
			out += [[start, stop]]
	return np.atleast_2d(out)

def remove_break(vals, brk):
	"""Remove a break from given start/stop intervals.

	Parameters
	----------
	vals : 2d array
		Start/stop intervals
	brk : float
		Break to be removed

	Returns
	-------
	2d array
		New intervals with break removed.

	Returns
	-------
	2d array
		New intervals with added break.

	Examples
	--------

	>>> vals = np.array([[0,1],[1,2],[2,3]])                                    
	>>> remove_break(vals, 2)                                                    
	array([[0. , 1. ],
		   [1. , 3. ]])

	"""
	out = []
	begin = vals[0,0]
	end = vals[-1,-1]
	gaps = np.column_stack((vals[:-1,1], vals[1:,0]))
	
	for start,stop in gaps:
		if (brk >= start) and (brk <= stop):
			pass
		else:
			out += [[start,stop]]
	
	out = np.atleast_2d(out)
	stops = np.append(out[:,0], end)
	starts = np.append(begin, out[:,1])

	return np.atleast_2d(np.column_stack((starts,stops)))