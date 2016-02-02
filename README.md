# darkspot

This module provides the script `find_darkspot`,
which gives you the darkest spot in the sky close
to zenith for a given date and time.

Darkest spot means area with the size of the FACT FoV
(4.5 degrees) that has the lowest light flux.

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

Just start the script
```
find_darkspot
```
and you will be prompted for date and time.


## Background:

This module uses the Hipparcos Star catalogue and
takes all stars into account that are brighter than 
10 mag.

Planets are added using the `ephem` library.
