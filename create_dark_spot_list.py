#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
Creates a csv with all spots in the sky that have no brighter star
than 6.5 visual magnitude

Usage:
    create_darkspot_list.py [options]

Options:
    --plot                    show the selected position, does not work on gui
'''
from __future__ import division, print_function
from docopt import docopt

args = docopt(__doc__)
print(args)

import pandas as pd
import numpy as np
from numpy import sin, cos, pi, arccos
from progressbar import ProgressBar


def angular_distance(ra1, dec1, ra2, dec2):
    '''
    Calculates angular distance between two objects
    Formula taken from http://www.gyes.eu/calculator/calculator_page1.htm

    Parameters
    ----------
    ra1 : number or array-like
        right ascension of first object in radians
    dec1 : number or array-like
        declination of first object in radians
    ra2 : number or array-like
        right ascension of second object in radians
    dec2 : number or array-like
        declination of second object in radians

    Returns
    -------
    angular_distance : number or numpy.ndarray
        the angular distance in radians between the objects
    '''
    term1 = sin(dec1) * sin(dec2)
    term2 = cos(dec1) * cos(dec2) * cos(ra1 - ra2)
    return arccos(term1 + term2)


def get_stars_in_fov(ra, dec, stars, fov=4.6):
    '''
    returns all the stars which are in FACTs field of view pointing
    for given coords ra, dec
    '''
    dist = angular_distance(ra, dec, stars.ra, stars.dec)
    mask = dist <= np.deg2rad(fov/2)
    return stars[mask]


def light_in_fov(ra, dec, stars, fov=4.6):
    '''
    returns the total star light in the fov

    Arguments
    ---------

    ra : float
        right ascension of pointing position in radians
    dec : float
        declination of pointing position in radians
    stars : pd.DataFrame
        the dataframe of the star catalogue containing flux, vmag, ra, dec
    '''
    infov = get_stars_in_fov(ra, dec, stars, fov)
    light = infov.flux.sum()
    brightest = infov.vmag.min()
    return light, brightest


def create_dataframe():
    stars = pd.read_csv(
        'hipparcos_vmag_less_10.tsv',
        delimiter='|',
        skiprows=40,  # the first 40 rows contain comments
        names=['ra', 'dec', 'vmag'],
    )

    # transform degrees to radians
    stars['ra'] = np.deg2rad(stars.ra)
    stars['dec'] = np.deg2rad(stars.dec)
    stars['flux'] = stars.vmag.apply(lambda x: 10**(-x/2.5)).sum()

    return stars


def create_darkspot_list(n_samples,
                         stars,
                         ):

    darkspots = pd.DataFrame(index=np.arange(n_samples))

    sin_dec = np.random.uniform(0, 1, n_samples//2)
    dec = np.arcsin(sin_dec)
    darkspots['dec'] = np.append(dec, -dec)
    darkspots['ra'] = np.random.uniform(0, 2*pi, n_samples)

    min_vmags = []
    fluxes = []
    prog = ProgressBar(maxval=len(darkspots.index))
    for idx, row in darkspots.iterrows():
        flux, min_vmag = light_in_fov(row.ra, row.dec, stars)
        min_vmags.append(min_vmag)
        fluxes.append(flux)
        prog.update(idx+1)

    darkspots['flux'] = fluxes
    darkspots['min_vmag'] = min_vmags

    return darkspots


def plot(stars, darkspots):
    import matplotlib.pyplot as plt
    from cartopy import crs

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=crs.Mollweide())
    ax.background_patch.set_facecolor('black')

    ax.scatter(
        np.rad2deg(darkspots.ra) - 180,
        np.rad2deg(darkspots.dec),
        transform=crs.PlateCarree(),
        linewidth=0,
        s=10,
        c='blue',
    )
    ax.scatter(
        np.rad2deg(darkspots.query('min_vmag > 6.5').ra) - 180,
        np.rad2deg(darkspots.query('min_vmag > 6.5').dec),
        transform=crs.PlateCarree(),
        linewidth=0,
        s=10,
        c='red',
    )

    ax.scatter(
        np.rad2deg(stars.ra) - 180,
        np.rad2deg(stars.dec),
        c=stars.vmag,
        vmax=15,
        s=0.3*(-stars.vmag + stars.vmag.max())**2,
        transform=crs.PlateCarree(),
        cmap='gray_r',
        linewidth=0,
        vmin=-1,
    )

    fig.tight_layout()
    plt.show()


def main():
    stars = create_dataframe()
    darkspots = create_darkspot_list(20000, stars)

    if args['--plot']:
        plot(stars, darkspots)

    darkspots = darkspots.query('min_vmag > 6.5')
    darkspots.to_csv('darkspots.csv', index=False)


if __name__ == '__main__':
    try:
        main()
    except (KeyboardInterrupt, SystemExit):
        pass
