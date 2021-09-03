# tintervals

A python package collecting  functions and tools to work with time intervals and to convert time data.

## Installation

`pip3 install git+https://gitlab.ininrim.it/m.pizzocaro/tintervals.git`

## Requirements

* `numpy`, `scipy`
* `datetime`
* `astropy`
* `tzlocal`
* `ciso8601`

## Basic usage

`import tintervals as ti`

Functions to convert timetag formats (based on `ciso8601` and `astropy.time`):
| Function | Description | 
| ------ | ------ |
| epoch2mjd | convert timetag from seconds from the epoch to MJD  |
| mjd2epoch | convert timetag from MJD to seconds from the epoch  |
| kk2epoch  | convert the timetag ouput of K+K counter (local timezone) to epoch, cannot accept arrays |

Functions to manipulate array of timetags or array of start/stop intervals
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
