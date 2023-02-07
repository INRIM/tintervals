
import tintervals as ti


from datetime import datetime
import glob
import os
import numpy as np
from functools import reduce
import decimal
import sys
import copy
from pandas import read_csv 
#from ruamel.yaml import YAML
import yaml

# class Oscillator():
# 	name
# 	description
# 	lab
# 	code
# 	nominal_frequency
# 	grs_correction
# 	systematic_uncertainty

# class Link():
# 	oscA
# 	oscB
# 	sB
# 	r0

# 	data = t, delta, flag, extra

class Oscillator():
	""" A simple class to store oscillator informations.
	
	Attributes
	----------
	name : str
		Name of the oscillator
	v0 : decimal.Decimal
		Nominal frequency (could be zero)
	description : str
		Description of the oscillator (if any)
	grs_correction : float
		If the oscillator is a clock, gravitational redshift still to be applied to the clock (if any)
	systematic_uncertainty : float
		If the oscillator is a clock, systematic uncertainty of the clock (if any)
	"""
	def __init__(self, name='', v0=0., description='', grs_correction=0., systematic_uncertainty=0.):
		self.name = name
		
		self.v0 = decimal.Decimal(v0) if v0 else decimal.Decimal(0)
		
		self.description = description
		self.grs_correction = float(grs_correction)
		self.systematic_uncertainty = float(systematic_uncertainty)


	def __repr__(self):
		return str(self.__dict__)
	


class Link():
	""" The main class for handling link data.

	Parameters
	----------
	data : ndarray, optional
		data of the link (with columns timetags, comparator, flag), by default np.empty(0)
	r0 : decimal.Decimal, optional
		Nominal ratio, by default ratio of oscillator nominal frequencies (if any) or zero
	oscA : Oscillator, optional
		Oscillator A (denominator), by default None
	oscB : Oscillator, optional
		Oscillator B (numerator), by default None
	sB : float, optional
		Scaling factor, by default oscB nominal frequency (if any) else 1
	name : string, optional
		Name of the link, by default None
	step : float, optional
		Time step of the data, by default 1.
	
	"""
	def __init__(self, data=np.empty(0), r0=None,  oscA=None, oscB=None, sB=None, name=None, step=1.):
		if oscA:
			self.oscA = Oscillator(oscA) if isinstance(oscA, str) else oscA
		else:
			self.oscA = Oscillator()

		if oscB:
			self.oscB = Oscillator(oscB) if isinstance(oscB, str) else oscB
		else:
			self.oscB = Oscillator()


		implicit_r0 = self.oscB.v0/self.oscA.v0 if self.oscA.v0 else 0

		self.r0 = decimal.Decimal(r0) if r0 else implicit_r0 if implicit_r0 else decimal.Decimal(1.)


		# scaling as proposed by Lodewyck (= oscB.v0 for fractional frequencies, +/- 1 for transfer beats)
		self.sB = float(sB) if sB else float(self.oscB.v0) if self.oscB.v0 else 1

		

		data = np.atleast_2d(data)

		# get at least 3 columns
		if data.shape[1] < 3:
			pad_size = 3 - data.shape[1]
			npad = [(0, 0)] * data.ndim # https://stackoverflow.com/questions/19349410/how-to-pad-with-zeros-a-tensor-along-some-axis-python
			npad[1] = (0, pad_size)
			data = np.pad(data, npad, constant_values=1)

		self.data = data

		self._set_view()
		self.step = step if step else np.median(np.diff(self.t))
 

		self.name = name if name else self.oscB.name + '-' + self.oscA.name


	def normalize(self, new_r0=None):
		"""Renormalize link data to a new r0

		Parameters
		----------
		new_r0 : deciaml.Decimal, optional
			explicit new ratio for the normalization, by default None.
			If none, the link is renormalized to the ratio of oscillators nominal frequencies (if any). 


		Notes
		-----
		It requires that the nominal frequency of OscA is set, and either new_r0 or the nominal frequeny of oscB set.

		"""
		if not self.oscA.v0:
			raise ValueError('Cannot change r0 if the first oscillator has no nominal value (not accurate oscillator)')
		
		if new_r0:
			new_r0 = decimal.Decimal(new_r0)
			self.oscB.v0 = new_r0*self.oscA.v0
		elif self.oscB.v0:
			new_r0 = self.oscB.v0/self.oscA.v0
		else:
			raise ValueError('Cannot normalize link without explicit or implicit nominal ratio r0.')

		new_sB = float(self.oscB.v0)


		# Lodewyck2020 appendix C
		fT = self.delta*self.sB
		new_fT = fT + float((self.r0 - new_r0)*self.oscA.v0)
		self.data[:,1]  =new_fT/new_sB
		self.r0 = new_r0
		self.sB = new_sB

	def _set_view(self):
		# note that here delta is always the comparator output (notation sligly different than Lodewyck2020 )
		self.t = self.data[:,0].view()
		self.delta = self.data[:,1].view()
		self.flag = self.data[:,2].view()
		
	def drop_invalid(self, drop_flag=0):
		"""Remove invalid points

		Parameters
		----------
		drop_flag : int or float, optional
			drop point whose flag is < drop_flag, by default 0
		"""
		mask = self.flag > drop_flag
		self.data = self.data[mask]
		self._set_view()

	def crop(self, start=None, stop=None):
		if start is None:
			start = np.amin(self.t)
		if stop is None:
			stop = np.amax(self.t) + self.step

		mask = (self.t>=start) & (self.t<=stop)
		self.data = self.data[mask]
		self._set_view()

	def extend(self, start=None, stop=None):
		"""Fill in missing timetags with invalid data.

		Parameters
		----------
		start : float, optional
			start point of the extension, by default the minimum timetag
		stop : float, optional
			stop point of the extension, by default the maximum timetag

		Timetags are filled in the np.arange(start,stop,link.step).
		Inserted data has delta and flag zero.

		"""
		if start is None:
			start = np.amin(self.t)
		if stop is None:
			stop = np.amax(self.t) + self.step
		
		if start > np.amin(self.t) or stop < (np.amax(self.t) + self.step):
			raise ValueError('Start/stop of extend have to be larger than the available data.')

		ext_t = np.arange(start, stop, self.step)

		isin = np.isin(ext_t, self.t)

		if sum(isin) != len(self.t):
			raise ValueError('Extended time scale does not include existing timetags!')
		else:
			axis=0
			npad = [(0, 0)] * self.data.ndim # https://stackoverflow.com/questions/19349410/how-to-pad-with-zeros-a-tensor-along-some-axis-python
			pad_size = len(ext_t) - len(self.t)
			npad[axis] = (pad_size, 0)
			ext_data = np.pad(self.data, pad_width=npad)

			idx = np.argsort(isin) # sort first the data in ext_t but not self.t
			idx2 = np.argsort(ext_t[idx])
			self.data = ext_data[idx2]
			self.data[:,0] = ext_t
			self._set_view()

	def copy(self):
		"""Return a copy of the object.
		"""
		cp =  copy.deepcopy(self)
		cp._set_view()
		return cp


	def average(self, base=1., offset=0., drop_flag=0, timetags_as_start=True, **kwargs):
		"""Reduce the link averaging data in place.

		Parameters
		----------
		base : float, optional
			Time base of average (average data every base seconds), by default 1.
		offset : float, optional
			offset in the time base (if zero, multiple of base are used), by default 0.
		drop_flag : int or float, optional
			drop point whose flag is < drop_flag, by default 0
		timetags_as_start : bool, optional
			If true, new timetags are the start of the averaging intervals; if false, the center, by default True
		**kwargs : 
			other keyword arguments passed to link_average

		Notes
		-----
		Method base on the function link_average.

		"""
		vals, average, count = link_average(self, base=base, offset=offset, drop_flag=drop_flag, timetags_as_start=timetags_as_start, **kwargs)


		self.data = np.atleast_2d(average.data)
		self._set_view()
		self.step = base



	def __repr__(self):
		return repr(self.data)


	def __add__(self, other):
		link, mask1, mask2 = chain2(self, other)
		return link

	def __neg__(self):
		# Lodewyck2020 eq 34
		
		new_r0 = self.r0**-1
		new_sB = float(self.oscA.v0) if self.oscA.v0 else 1.

		fT = self.delta*self.sB
		new_fT = -float(new_r0)*fT

		new_data = self.data.copy()
		new_data[:,1]=new_fT/new_sB

		return Link(oscA=self.oscB, oscB=self.oscA, data=new_data, r0=new_r0, sB = new_sB, step=self.step)
		

	def __sub__(self, other):
		link, mask1, mask2 = chain2(self, -other)
		return link


		


def chain2(link1, link2):
	"""Chain 2 links in sequence.

	Parameters
	----------
	link1 : link B/A
		link object to be chained.
	link2 : link C/B
		link object to be chained.

	Returns
	-------
	res 
		chained link object C/A.
	mask1
		mask to be applied to link1 to select timetags common to link2
	mask2
		mask to be applied to link2 to select timetags common to link1

	"""
	if link1.step != link2.step:
	 	raise ValueError('Cannot chain links with inconsitent time steps.')

	if not link1.oscA.v0:
		raise ValueError('Cannot chain links in case the first oscillator in the chain has no nominal value (not accurate oscillator)')


	common_t = np.intersect1d(link1.t, link2.t)

	mask1 = np.isin(link1.t, common_t)
	mask2 = np.isin(link2.t, common_t)

	common_flag = np.minimum(link1.flag[mask1], link2.flag[mask2])

	# from Lodewyck2020
	R1 = link1.sB/float(link1.oscA.v0*link1.r0)*link1.delta[mask1]
	R2 = link2.sB/float(link1.oscA.v0*(link1.r0*link2.r0))*link2.delta[mask2]

	new_data = np.column_stack((common_t, R1+R2, common_flag))
	
	

	res = Link(oscB=link2.oscB, oscA=link1.oscA, data=new_data, r0=link1.r0*link2.r0, sB=link1.r0*link2.r0*link1.oscA.v0, step=link1.step)
	if res.oscB.v0:
		res.normalize()
	
			

	return res, mask1, mask2


def chain(*links):
	"""Chain multiple links in sequence.

	Parameters
	----------
	*links :
		link objects to be chained.
	
	Returns
	-------
	res 
		chained link object.
	masks
		list of masks to select timetags in the common uptime (one for each link).

	"""
	step = links[0].step
	for l in links:
		if l.step != step:
			raise ValueError('Cannot chain links with inconsitent time steps.')

	if not links[0].oscA.v0:
		raise ValueError('Cannot chain links in case the first oscillator in the chain has no nominal value (not accurate oscillator)')


	common_t = reduce(np.intersect1d, [link.t for link in links]) # see https://numpy.org/doc/stable/reference/generated/numpy.intersect1d.html

	masks = [np.isin(link.t, common_t) for link in links]
	
	common_flag = np.min([link.flag[mask] for link, mask in zip(links, masks)], axis=0)

	den = np.cumprod([x.r0 for x in links])
	v0A = links[0].oscA.v0
	
	common_R = [link.sB*link.delta[mask]/(float(v0A*d)) for link,mask,d in zip(links, masks, den)]
	res = np.sum(common_R, axis=0)



	new_data = np.column_stack((common_t, res, common_flag))
	
	

	res = Link(oscB=links[-1].oscB, oscA=links[0].oscA, data=new_data, r0=den[-1], sB=den[-1]*links[0].oscA.v0, step=step)
	if res.oscB.v0:
		res.normalize()
	
			

	return res, masks



def average(link, intervals):
	"""Average a link in given intervals

	Parameters
	----------
	link : Link
		link object to be averaged.
	intervals : str or 2d arrays
		One of 'hour','day' or 'bipm' to average every hour, day or BIPM 5-day intervals.
		Or a 2d-array of start-stop intervals.

	Returns
	-------
	vals
		start/stop intervals of the averages.
	res
		averaged link data (all columns, so average timetag, average delta and average flag).
	count
		number of averaged point for each intervals.


	Notes
	-----
	Based on tintervals.maverage.

	"""
	if isinstance(intervals, str):
		base, offset = _string_to_base(intervals)
		intervals = ti.regvals(np.amin(link.t), np.amax(link.t), base, offset)
	else:
		intervals = np.atleast_2d(intervals)

	vals, ave, count = ti.maverage(link.data, link.t, intervals)
	return vals, ave, count

def link_average(link, base=1., offset=0., drop_flag=0, timetags_as_start=True,  **kwargs):
	"""Average a link in regular intervals.

	Parameters
	----------
	link : Link
		link object to be averaged
	base : float, optional
		Time base of average (average data every base seconds), by default 1.
	offset : float, optional
		offset in the time base (if zero, multiple of base are used), by default 0.
	drop_flag : int or float, optional
		drop point whose flag is < drop_flag, by default 0
	timetags_as_start : bool, optional
		If true, new timetags are the start of the averaging intervals; if false, the center, by default True

	Returns
	-------
	vals
		start/stop intervals of the averages.
	res
		averaged link object.
	count
		number of averaged point for each intervals.

	"""
	link.drop_invalid(drop_flag)


	if isinstance(base, str):
		base, offset = _string_to_base(base, offset)

	vals, average, count = ti.raverage(link.data, link.t, base=base, offset=offset, step=link.step, **kwargs)

	# new timatags are starts
	if timetags_as_start:
		average[:,0] = vals[:,0]

	# new flags are percentage
	average[:,2] = count/(base/link.step)
	
	return vals, Link(average, r0=link.r0, oscA=link.oscA, oscB=link.oscB, sB=link.sB, step=base), count

	

def _string_to_base(short, offset=0):
	if short == 'hour': base=3600.
	elif short == 'day': base=86400.
	elif short == 'bipm':
		base=86400.*5
		offset=+2*86400.
	else:
		raise NotImplementedError('{} is unsupported. Use hour/day/bipm.'.format(short))

	return base, offset




HEADER_STD = """# Data for {linkname} 
# File generated on: {now}
# With the script: {command} 
"""

HEADER_MESSAGE = """# 
# {}
"""

def _decimal2string(d):
	"""Format a decimal to a  human readable string."""
	# https://stackoverflow.com/questions/11227620/drop-trailing-zeros-from-decimal/18769210#18769210
	normalized = d.normalize()
	sign, digits, exponent = normalized.as_tuple()
	if exponent > 0:
		normalized =  decimal.Decimal((sign, digits + (0,) * exponent, 0))
	
	out = "{}".format(normalized)
	return out #.replace(',', '_')



def save_link_to_dir(dir, link, extra_names=[], message='', time_format='mjd', yfmt='{:.10e}'):
	"""Save Link data to a directory.

	Parameters
	----------
	dir : string
		ouput directory
	link : Link
		link to be saved
	extra_names : list, optional
		if Link data has more than 3 columns, this can be used to specify their name, by default []
	message : str, optional
		message to be written in the header, by default ''
	time_format : ['mjd', 'iso', 'unix'], optional
		Specify ouput time format, by default 'mjd'
	digits : int, optional
		Number of digits for the Link delta, by default 10


	Notes
	-----
	Data is saved in separate files, one file for each day.
	Link metadata is saved in a yaml file in the directory.


	"""
	sub = os.path.join(dir,link.name)
	if not os.path.exists(sub):
		os.makedirs(sub)

	if time_format == 'iso':
		time_converter = ti.epoch2iso
		time_format = '{}'
	elif time_format == 'mjd':
		time_converter = ti.epoch2mjd
		time_format = '{:.6f}'
	elif time_format == 'unix':
		time_converter = lambda x: x
		time_format = '{}'
	else:
		raise ValueError("Unrecognized time_format. Valid formats are 'iso', 'mjd' and 'unix'.")

	if link.step <= 1:
		ffmt = '{:.0f}'
	else:
		ffmt = '{:.6f}'

	# save metadata to yaml
	metadata = {'name': link.name}
	if link.oscA.v0:
		metadata['numrhoBA'] = _decimal2string(link.oscA.v0*link.r0)
		metadata['denrhoBA'] = _decimal2string(link.oscA.v0)
	elif link.oscB.v0:
		metadata['numrhoBA'] = _decimal2string(link.oscB.v0)
		metadata['denrhoBA'] = _decimal2string(link.oscb.v0/link.r0)
	else:	
		metadata['numrhoBA'] = _decimal2string(link.r0)
		metadata['denrhoBA'] = _decimal2string(decimal.Decimal(1))

	metadata['sB'] = link.sB

	if link.oscA.v0:
		metadata['nu0A'] = _decimal2string(link.oscA.v0)
	if link.oscA.grs_correction or link.oscA.systematic_uncertainty:
		metadata['grsA'] = link.oscA.grs_correction
		metadata['uA_sys'] = link.oscA.systematic_uncertainty
		
	if link.oscB.v0:
		metadata['nu0B'] = _decimal2string(link.oscB.v0)
	if link.oscB.grs_correction or link.oscB.systematic_uncertainty:
		metadata['grsB'] = link.oscB.grs_correction
		metadata['uB_sys'] = link.oscB.systematic_uncertainty

	if link.step != 1:
		metadata['interval'] = float(link.step)

	with open(os.path.join(sub, link.name+'.yml'), 'w') as metafile:
		# print(metadata)
		#YAML().dump([metadata], metafile)
		yaml.dump([metadata], metafile, sort_keys=False)


	mjd = ti.epoch2mjd(link.t)
	
	mjd_with_data = np.unique(np.floor(mjd))
	now = datetime.now().astimezone().replace(microsecond=0).isoformat()

	for day in mjd_with_data:
		date = ti.iso_from_mjd(day)[:10]
		mask = (mjd >= day) & (mjd < day + 1)
	
		

		filename = date + '_' + link.name + '.dat'

		header = HEADER_STD.format(linkname=link.name, now=now, command=' '.join(sys.argv))

		if message:
			header += HEADER_MESSAGE.format(message)


		names = ['t', 'ΔA→B', 'flag'] + extra_names 
		#yfmt = '{{:.{}e}}'.format(digits)
		fmtlst = [time_format, yfmt, ffmt] + ['{}']*(link.data.shape[1]-3)
		fmt = '\t'.join(fmtlst) + '\n'

		col_len = [len(fmtlst[0].format(time_converter(link.data[0][0])))] + [len(f.format(d)) for f, d in zip(fmtlst[1:], link.data[0,1:])]
		col_tit = [x.ljust(y) for x, y in zip(names, col_len)]
		header += '# \n# ' + '\t'.join(col_tit) + '\n'

		with open(os.path.join(sub, filename), 'w') as filo:
			filo.write(header)
			for d in link.data[mask]:
				filo.write(fmt.format(time_converter(d[0]), *d[1:]))	


def load_link_from_dir(dir, meta=None, start=None, stop=None, discard_invalid=True, time_format='mjd', ext=['.dat','.txt'], round=True, remove_not_unique=True, verbose=False, **kwargs):
	"""Load link data from a directory.

	Parameters
	----------
	dir : str
		Directory to be read
	meta : str
		yaml file with metadata, by default search for a yaml file with the same name of dir.
	start : float
		if given, crop the data from this point, by default None
	stop : float
		if given, crop the data to this point, by default None
	discard_invalid : bool, optional
		if True, discard data with flag=0, by default True
	time_format : str, optional
		Specify input time format, 'mjd', 'iso' or 'unix', by default 'mjd'
	ext : list or string, optional
		extension of the file to read, by default ['.dat','.txt']
	round : bool, optional
		if True, round timetags to the nearest second, by default True
	remove_not_unique : bool, optional
		if True, remove data corresponding to not-unique timetags, by default True
	verbose : bool, optional
		if True, verbose ouput
	**kwargs :
		other keyword argument to pass to np.genfromtxt

	Returns
	-------
	Link
		Link object with data from the directory.
	"""
	if time_format == 'iso':
		time_converter = ti.iso2epoch
	elif time_format == 'mjd':
		time_converter = lambda x: ti.epoch_from_mjd(float(x))
	elif time_format == 'unix':
		time_converter = lambda x: x
	else:
		raise ValueError("Unrecognized time_format. Valid formats are 'iso', 'mjd' and 'unix'.")
	

	if isinstance(ext, str):
		ext = [ext]

	# normpath avoid problems with trailing slashes
	name = os.path.basename(os.path.normpath(dir))
	parts = name.split('-')

	
	# search yaml file if not specified
	if not meta:
		metafiles = glob.glob(os.path.join(dir,name + '.y*ml'))
		if len(metafiles) == 0:
			metafiles = glob.glob(os.path.join('**',name + '.y*ml'), recursive=True)			

		if len(metafiles) > 0:
			meta = metafiles[0]
			
	if meta:
		if verbose:
			print(f'YAML metadata file: {meta}')

		with open(meta, 'r') as mf:
			#metadata = YAML(typ='base').load(mf)	
			# baseloader : everything is a string
			metadata = yaml.load(mf, Loader=yaml.BaseLoader)
			metadata = {x['name']: x for x in metadata}
			metadata = metadata.get(name, {})
	else:
			if verbose:
				print('No metadata file found.')
			metadata = {}

	


	if len(parts) >=2:
		nameA, nameB = parts[-1], parts[0]
	else:
		nameA, nameB = name + 'A', name+'B'

	# conversions handled by the Oscilator class
	oscA = Oscillator(nameA, v0=metadata.get('nu0A', 0), grs_correction=metadata.get('grsA', 0), systematic_uncertainty=metadata.get('uA_sys', 0))
	oscB = Oscillator(nameB, v0=metadata.get('nu0B', 0), grs_correction=metadata.get('grsB', 0), systematic_uncertainty=metadata.get('uB_sys', 0))
	
	r0 = decimal.Decimal(metadata.get('numrhoBA', 1))/decimal.Decimal(metadata.get('denrhoBA', 1))
	sB = float(metadata.get('sB', 1))
	step = float(metadata.get('interval', '0').strip('s'))

	datafiles = []
	for e in ext:
		datafiles += glob.glob(os.path.join(dir, '*' + e)) 

	datafiles = np.sort(datafiles) # sort also tranform in array


	if len(datafiles) == 0:
		raise IOError('No datafiles with extension {} found in {}'.format(ext, dir))
	else:	
		load = []

		for fili in datafiles:
			#d = np.genfromtxt(fili, encoding='UTF-8', converters={0:time_converter}, **kwargs)

			# faster!
			d = read_csv(fili, sep='\s+', header=None, converters={0:time_converter}, comment='#', **kwargs)
			d = np.atleast_2d(d)

			if verbose:
				print(f'Loading {fili}')

			if np.size(d) >0:
				# TODO: add flag if missing
				
				if discard_invalid:
					d = d[d[:,2] > 0]

			
				load += [d]
			
		load = np.concatenate(load)
		if round:
			load[:,0] = np.around(load[:,0])

		if remove_not_unique:
			uniq, idx, count = np.unique(load[:,0], return_index=True, return_counts=True)
			if (count>1).any():
				load = load[idx]

				if verbose:
					print('Removing {} not unique timetags'.format(sum(count>1)))




		out = Link(data=load, r0=r0, sB=sB, oscA=oscA, oscB=oscB, name=name, step=step)
		out.crop(start, stop)
		return out



def load_links_from_osc_names(link_dos, dir='.',  **kwargs):
	"""Loads a chain of links from oscillators names.

	Parameters
	----------
	link_dos : list of strings
		names of the oscillators to be chained
	dir : str, optional
		directory to be read, by default '.'
	**kwargs : 
		other keyword arguments as in load_link_from_dir


	Returns
	-------
	list of Link objects
		loaded data

	Notes
	-----
	Loads a chain of links whose names are in the form OscA-OscB.
	For example:
	
	>>> load_links_from_osc_names(['A', 'B', 'C', 'D'])


	would load directories 'B-A', 'C-B' and 'D-C' 
	(or load and invert the links from directories 'A-B', 'B-C' and 'C-D').

	"""
	link_names = [b + '-' + a for b, a in zip(link_dos[::-1][:-1], link_dos[::-1][1:])]

	links = []

	for name in link_names[::-1]:
		
		if os.path.exists(os.path.join(dir, name)):
			links += [load_link_from_dir(os.path.join(dir, name),  **kwargs)]
		else:
			parts = name.split('-')
			invname = parts[1] + '-' + parts[0]
			links += [-load_link_from_dir(os.path.join(dir, invname),  **kwargs)]
	
	return links






if __name__ == '__main__':
	A = Oscillator('A','100.')
	B = Oscillator('B', '100.', grs_correction=1e-16)
	C = Oscillator('C', None)
	D = Oscillator('D','100.3')


	BA = Link(np.array([0, 0.1, 1]),sB=1, r0=1, oscA=A, oscB=B)
	CB = Link(np.array([0, 0.2, 1]), sB=1, r0=1, oscA=B, oscB=C)
	DC = Link(np.array([0, 0., 1]), sB=1, r0=1, oscA=C, oscB=D)

	save_link_to_dir('TEST',BA)
	BARET = load_link_from_dir('./TEST/B-A')



# 	print((BA +CB + DC).data)
# 	print(chain(BA,CB,DC)[0].data)

# 	A = Oscillator('test_A','10.')
# 	B = Oscillator('test_B', '20.')
# 	C = Oscillator('test_C', '20.')
# 	D = Oscillator('test_D','10.03')


# 	BA = Link(np.array([0, 1e-12, 1]), oscA=A, oscB=B)
# 	CB = Link(np.array([0, 1e-12, 1]), oscA=B, oscB=C)
# 	DC = Link(np.array([0, 2e-12, 1]), oscA=C, oscB=D)

# 	print((BA +CB + DC).data)
# 	print(chain(BA,CB,DC)[0].data)

# 	FE = Link(np.array([[0, 2e-12, 1],[1, 2.2e-12,1]]))



# 	N = 100
# 	data = np.column_stack((np.arange(N), np.arange(N), np.ones(N)))

# 	ame = Link(data)

# 	vals, ave, count = average(ame, [[0,20],[20,100],[1000,2000]])
# 	print(vals)
# 	print(ave)
# 	print(count)

# #	print(ame)
# 	ame.average(base=10)
# 	print(ame)
