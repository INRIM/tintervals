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

Functions to convert timetag formats (based on `ciso8601` and `astropy.time`):
| Function | Description | 
| ------ | ------ |
| epoch2mjd | convert timetag from seconds from the epoch to MJD  |
| mjd2epoch | convert timetag from MJD to seconds from the epoch  |
| kk2epoch  | convert the timetag ouput of K+K counter (local timezone) to epoch, cannot accept arrays |

## License

[MIT](https://opensource.org/licenses/MIT)

## Authors

(c) 2021 Marco Pizzocaro - Istituto Nazionale di Ricerca Metrologica (INRIM)
