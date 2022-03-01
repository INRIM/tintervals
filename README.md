# tintervals

A python package collecting  functions and tools to work with time intervals and to convert time data.

The package also provides utilities for handling optical link data developed for the [EMPIR project ROCIT](http://empir.npl.co.uk/rocit/).

Development is available at 
- https://github.com/INRIM/tintervals
- https://gitlab.ininrim.it/m.pizzocaro/tintervals (INRIM only)

Documentation is available at https://tintervals.readthedocs.io

Package is available at https://pypi.org/project/tintervals/


## Installation

The package can be installed using pip:

`pip install tintervals`

or directly from github:

`pip install git+https://github.com/INRIM/tintervals.git`

## Requirements

* `numpy`, `scipy`
* `ciso8601` (for fast ISO format reading)
* `pandas` (used for fast file loading)
* `pyyaml`

## Basic usage

`import tintervals as ti`

Functions to convert timetag formats.
Fast functions used as converters when importing files:

| Function        | From               | To    | 
| --------------- | ------------------ | ----- |
| `kk2epoch`      | K+K counter format | Epoch |
| `iso2epoch`     | ISO format*        | Epoch |
| `kk2iso`        | K+K counter format | ISO format (naive) |
| `iso2datetime`  | ISO format*        | `datetime` |
| `datetime2iso`  | `datetime`         | ISO format ('Z' notation for UTC) |
| `datetime2epoch`| `datetime`*        | Epoch |
| `epoch2datetime`| Epoch              | Datetime |
| `epoch2iso`     | Epoch              | ISO format ('Z' notation for UTC) |
| `mjd2iso`       | MJD                | ISO format ('Z' notation for UTC) |
| `epoch2mjd`     | Epoch              | MJD   |
| `mjd2epoch`     | MJD                | Epoch |


Starred  inputs (*), if naive are considered as system/local time.
K+K format is always naive and interpreted as sytem/local time.
ISO format can be read with microseconds but it is printed without.
Conversion from Epoch time (Unix) and MJD is done by simple affine function, 
as both are aligned with UTC and basically ignore leap seconds.

Vectorized functions:

| Function        | From               | To    | 
| --------------- | ------------------ | ----- |
| `mjd_from_epoch`| Epoch              | MJD   |
| `epoch_from_mjd`| MJD                | Epoch |
| `iso_from_epoch`| Epoch              | ISO format ('Z' notation for UTC) |
| `iso_from_mjd`  | MJD                | ISO format ('Z' notation for UTC) |


Functions to manipulate array of timetags or array of start/stop intervals:

| Function | Description | 
| ------ | ------ |
| `array2intervals` | convert from an array of timetags to an array of start/stop intervals |
| `intervals2weights` | convert from start/stop intervals to timetags and weights  |
| `intersect` | take the intersection of two arrays of start/stop intervals |
| `split` | Subdivide an array of start/stop intervals to a fixed scale (e.g., every 10 s)|
| `regvals` | retrun regular intervals between a start and stop |
| `raverage` | Average data with timetags in regular intervals (reshape algorithm)|
| `maverage` | Average data with timetags in given intervals (mask algorithm)|
| `csaverage` | Average data  start/stop intervals in different start/stop intervals (cumsum algorithm)|

Functions to calculate deadtime uncertainty:

| Function | Description | 
| ------ | ------ |
| `deadtime.unc_fft` | calculate deadtime uncertainty from given maser noise (FFT algorithm) |


## Advanced usage
For handling optical links:

`import ti.rocitlinks as rl`

## License

[MIT](https://opensource.org/licenses/MIT)


## Acknowledgments
This work is partially funded by the European Metrology Program for Innovation and Research (EMPIR) project 18SIB05 ROCIT.
The EMPIR initiative is cofunded by the European Union’s Horizon 2020 research and innovation programme and the EMPIR Participating States.

## Authors

(c) 2021, 2022 Marco Pizzocaro - Istituto Nazionale di Ricerca Metrologica (INRIM)

## Speed test
```
ti.kk2iso('210922*161141.810')
166 ns ± 2.09 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)
ti.kk2epoch('210922*161141.810')
658 ns ± 1.24 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
ti.iso2datetime('20210922T161141.810Z')
222 ns ± 18 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
ti.iso2epoch('20210922T161141.810Z')
473 ns ± 5.11 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
ti.datetime2iso(datetime(2021,9,22,11,12,tzinfo=timezone.utc))
2.51 µs ± 9.92 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
ti.datetime2epoch(datetime(2021,9,22,11,12,tzinfo=timezone.utc))
613 ns ± 3.99 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
ti.epoch2datetime(1632319901)
406 ns ± 2.87 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
ti.mjd_from_epoch(1632319901)
131 ns ± 0.334 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)
ti.epoch_from_mjd(59479.356)
99.4 ns ± 0.228 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)
ti.iso_from_epoch(1632319901)
19.8 µs ± 86.8 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
ti.iso_from_mjd(59479.356)
23 µs ± 266 ns per loop (mean ± std. dev. of 7 runs, 10000 loops each)
```

