# tintervals

A python package collecting  functions and tools to work with time intervals and to convert time data.

## Installation

`pip3 install git+https://gitlab.ininrim.it/m.pizzocaro/tintervals.git --upgrade`

## Requirements

* `numpy`, `scipy`
* `datetime`
* `astropy`
* `ciso8601`

## Basic usage

`import tintervals as ti`

Functions to convert timetag formats (based on `ciso8601` and `astropy.time`).
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

Starred inputs, if naive are considered as system/local time.
K+K format is always naive and interpreted as sytem/local time.

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
| array2intervals | convert from an array of timetags to an array of start/stop intervals |
| mix | take the intersection of two arrays of start/stop intervals |
| split | Subdivide an array of start/stop intervals to a fixed scale (e.g., every 10 s)|
| csaverage | Average data given for some start/stop intervals in different start/stop intervals|


## License

[MIT](https://opensource.org/licenses/MIT)

## Authors

(c) 2021 Marco Pizzocaro - Istituto Nazionale di Ricerca Metrologica (INRIM)

## Speed test
```
ti.kk2iso('210922*161141.810')
165 ns ± 0.687 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)
ti.kk2epoch('210922*161141.810')
691 ns ± 0.775 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
ti.iso2datetime('20210922T161141.810Z')
216 ns ± 2.85 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
ti.iso2epoch('20210922T161141.810Z')
477 ns ± 3.17 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
ti.datetime2iso(datetime(2021,9,22,11,12,tzinfo=timezone.utc))
2.6 µs ± 114 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
ti.datetime2epoch(datetime(2021,9,22,11,12,tzinfo=timezone.utc))
621 ns ± 4.96 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
ti.epoch2datetime(1632319901)
407 ns ± 12 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
ti.mjd_from_epoch(1632319901)
260 µs ± 733 ns per loop (mean ± std. dev. of 7 runs, 1000 loops each)
ti.epoch_from_mjd(59479.356)
194 µs ± 402 ns per loop (mean ± std. dev. of 7 runs, 10000 loops each)
ti.iso_from_epoch(1632319901)
16.7 µs ± 95.8 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)
ti.iso_from_mjd(59479.356)
502 µs ± 2.34 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
```



