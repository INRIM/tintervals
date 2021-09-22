# https://stackoverflow.com/questions/4761041/python-import-src-modules-when-running-tests
import sys
sys.path.insert(0, '../src') 


import tintervals as ti
from datetime import datetime, timezone


assert ti.kk2iso('210922*161141.810') == '20210922 161141.810'
assert ti.kk2epoch('210922*161141.810') ==  1632319901.810

# check that this works also on solar hour
test = datetime(2021,12,12,12,12,12, tzinfo=timezone.utc)
assert ti.kk2epoch('211212*131212') == test.timestamp()



assert ti.epoch2datetime(ti.mjd2epoch(60_000))  == datetime(2023, 2, 25, 0, 0, tzinfo=timezone.utc)

assert ti.datetime2epoch(datetime(2021, 9, 22, 14, 25, 7, tzinfo=timezone.utc)) == 1632320707

assert ti.epoch2mjd(ti.datetime2epoch(datetime(2021, 9, 22, 0, 0, 0, tzinfo=timezone.utc))) == 59479.


# some timing

from IPython import get_ipython
ipython = get_ipython()


cmds = [
"ti.kk2iso('210922*161141.810')",
"ti.kk2epoch('210922*161141.810')",
"ti.iso2datetime('20210922T161141.810Z')",
"ti.iso2epoch('20210922T161141.810Z')",
"ti.datetime2iso(datetime(2021,9,22,11,12,tzinfo=timezone.utc))",
"ti.datetime2epoch(datetime(2021,9,22,11,12,tzinfo=timezone.utc))",
"ti.epoch2datetime(1632319901)",
"ti.epoch2mjd(1632319901)",
"ti.mjd2epoch(59479.356)",
]

for cmd in cmds:
	print(cmd)	
	ipython.magic("timeit " + cmd)



