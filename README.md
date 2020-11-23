# darkspot

This module provides the script `find_darkspot`,
which gives you the darkest spot in the sky close
to zenith for a given date and time.

## installation:

```
pip install git+https://github.com/fact-project/darkspot 
```

Optionally, for plotting support, also install `cartopy`:
```
conda install geos proj4
pip install cartopy
```

## Usage:

```
$ find_darkspot --help
Usage: find_darkspot [OPTIONS]

Options:
  -t, --time TIME       ISO8601 datestring for the observation time. If not
                        given, now is used.

  --site TEXT           Site name (must be known to astropy)
  --lat FLOAT           Latitude of observatory in decimal degrees, used to
                        override --site

  --lon FLOAT           Longitude of observatory in decimal degrees, used to
                        override --site

  --min-altitude FLOAT  Minimum altitude to consider in degrees
  --fov FLOAT           Diameter of field of view in degrees
  --plot                Show plot
  --band [V|BT]         Which optical band to use.
```

![img](showflux.png)


## Background:

This module uses the Hipparcos Star catalogue and
takes all stars into account that are brighter than 10 mag.

Planets are added using the `ephem` library.
