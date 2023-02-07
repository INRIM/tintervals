# I am running a venv now!
# # https://stackoverflow.com/questions/4761041/python-import-src-modules-when-running-tests
# import sys
# sys.path.insert(0, '../src') 


import tintervals as ti
from datetime import datetime, timezone


assert ti.kk2iso('210922*161141.810') == '20210922 161141.810'
assert ti.kk2epoch('210922*161141.810') ==  1632319901.810

# check that this works also on solar hour
test = datetime(2021,12,12,12,12,12, tzinfo=timezone.utc)
assert ti.kk2epoch('211212*131212') == test.timestamp()



assert ti.epoch2datetime(ti.epoch_from_mjd(60_000))  == datetime(2023, 2, 25, 0, 0, tzinfo=timezone.utc)

assert ti.datetime2epoch(datetime(2021, 9, 22, 14, 25, 7, tzinfo=timezone.utc)) == 1632320707

assert ti.mjd_from_epoch(ti.datetime2epoch(datetime(2021, 9, 22, 0, 0, 0, tzinfo=timezone.utc))) == 59479.


# some timing

# from IPython import get_ipython
# ipython = get_ipython()

# test circular T functions
print(ti.cirtt2mjd(2023,2))
print(ti.mjd2cirt(59974))
print(ti.cirtvals(59823,59849) )



import timeit
import_module = """
import tintervals as ti
from datetime import datetime, timezone
"""

cmds = [
"ti.kk2iso('210922*161141.810')",
"ti.kk2epoch('210922*161141.810')",
"ti.iso2datetime('20210922T161141.810Z')",
"ti.iso2epoch('20210922T161141.810Z')",
"ti.datetime2iso(datetime(2021,9,22,11,12,tzinfo=timezone.utc))",
"ti.datetime2epoch(datetime(2021,9,22,11,12,tzinfo=timezone.utc))",
"ti.epoch2datetime(1632319901)",
"ti.mjd_from_epoch(1632319901)",
"ti.epoch_from_mjd(59479.356)",
"ti.iso_from_epoch(1632319901)",
"ti.iso_from_mjd(59479.356)",
"ti.kk2mjd('210922*161141.810')",
"ti.iso2mjd('20210922T161141.810Z')",
"ti.datetime2mjd(datetime(2021,9,22,11,12,tzinfo=timezone.utc))",
"ti.mjd2datetime(59479.356)",
]

for cmd in cmds:
	print(cmd)	
	print(f'{timeit.timeit(cmd, setup=import_module)*100:.2f} ns')




print(ti.deadtime.unc_fft([0,5], [0,10], wfm=1))