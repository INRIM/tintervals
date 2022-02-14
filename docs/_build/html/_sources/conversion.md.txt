# Time conversion functions

`tintervals` provides conversion function from different time formats.

Supported formats are:

| Format       | Description               | Example    | 
| --------------- | ------------------ | ----- |
| Epoch (Unix) | Seconds from the Unix epoch | `time.time()` |
| MJD          | Modified Julian Date        | `59604.3`     |
| ISO          | ISO format                  | `20210922T161141Z` |
| K*K          | Format used by [K+K counters](http://www.kplusk-messtechnik.de/products/fxe_19.htm) | `210922*161141.810`|
|              | (either with `*` or space)  |  |
| datetime     | Python datetime object      | |


If naive, ISO and datetime inputs are considered as system/local time.
K+K format is always naive and interpreted as sytem/local time.
ISO format can be read with microseconds but it is printed without.
Conversion from Epoch time (Unix) and MJD is done by simple affine function, 
as both are aligned with UTC and basically ignore leap seconds.

Functions with a `2` in their name are optimized for speed and are meant to be used as converters, e.g.:
`np.genfromtxt('somedata.dat', converters={0: iso2epoch})`.
Function with a `from` in their name are numpy universal (vectorized) functions that work either on single values or on numpy arrays.
For a more general handling of time formats you should consider [`astropy.time`](https://docs.astropy.org/en/stable/time/index.html).


```{eval-rst}
.. automodule:: tintervals.conversion
  :members:
```