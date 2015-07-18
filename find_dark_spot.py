#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
Calculates the darkest spot at <date> for a zenith angle below 10 degrees
(altitude of 80 degrees).
Use YYYY-MM-DD hh:mm format for date and time.


Usage:
    find_dark_spot.py [<date> <timeutc>] [options]

Options:
    --max-zenith=<degrees>    maximal zenith for the dark spot [default: 10]
    --plot                    show the selected position, does not work on gui
'''
from __future__ import division, print_function
from docopt import docopt

args = docopt(__doc__)
min_altitude = 90 - int(args['--max-zenith'])
max_zenith = int(args['--max-zenith'])

import pandas as pd
import numpy as np
from numpy import sin, cos, tan, arctan2, arcsin, pi, arccos
import ephem
from progressbar import ProgressBar
from datetime import datetime
from blessings import Terminal
term = Terminal()


def fact_setup(date):
    ''' creates an ephem.Observer for FACT at given date '''
    fact = ephem.Observer()
    fact.lon = '-17:53:05'
    fact.lat = '28:45:15'
    fact.elevation = 2200
    fact.date = date
    fact.epoch = ephem.J2000

    return fact


def enter_datetime():
    ''' a small cli utility to get date and time from the user '''
    print('\nPlease enter date and time for the ratescan')
    print(term.red('This is the real date, be aware for times after 0:00'))
    date = raw_input('Date (YYYY-MM-DD): ')
    time = raw_input('Time UTC: (hh:mm): ')
    return date, time


def equatorial2horizontal(ra, dec, observer):
    '''
    Transforms from right ascension, declination to azimuth, altitude for
    the given observer.
    Formulas are taken from https://goo.gl/1wMU4u

    Parameters
    ----------

    ra : number or array-like
        right ascension in radians of the object of interest

    dec : number or array-like
        declination in radians of the object of interest

    observer : ephem.Observer
        the oberserver for which azimuth and altitude are calculated

    Returns
    -------
    az : number or numpy.ndarray
        azimuth in radians for the given ra, dec
    alt : number or numpy.ndarray
        altitude in radians for the given ra, dec
    '''

    obs_lat = float(observer.lat)

    h = observer.sidereal_time() - ra
    alt = arcsin(sin(obs_lat) * sin(dec) + cos(obs_lat) * cos(dec) * cos(h))
    az = arctan2(sin(h), cos(h) * sin(obs_lat) - tan(dec)*cos(obs_lat))

    # add 180 degrees to azimuth to have the same definition as stellarium
    az = np.mod(az + pi, 2*pi)

    return az, alt


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


def get_stars_in_fov(az, alt, stars, observer, fov=4.6):
    '''
    returns all the stars which are in FACTs field of view pointing
    for given coords az, alt
    '''
    ra, dec = observer.radec_of(az, alt)
    dist = angular_distance(ra, dec, stars.ra, stars.dec)
    mask = dist <= np.deg2rad(fov/2)
    return stars[mask]


def light_in_fov(az, alt, stars, observer):
    ''' returns the total star light in the fov at azimuth az, altitude alt '''
    infov = get_stars_in_fov(az, alt, stars, observer)
    light = infov.vmag.apply(lambda x: 10**(-x/2.5)).sum()
    return light


def create_dataframe(observer):
    stars = pd.read_csv(
        'hipparcos_vmag_less_10.tsv',
        delimiter='|',
        skiprows=40,  # the first 40 rows contain comments
        names=['ra', 'dec', 'vmag'],
    )

    # transform degrees to radians
    stars.ra = np.deg2rad(stars.ra)
    stars.dec = np.deg2rad(stars.dec)

    # add the planets
    sol_objects = [
        ephem.Mercury(),
        ephem.Venus(),
        ephem.Moon(),
        ephem.Mars(),
        ephem.Jupiter(),
        ephem.Saturn(),
        ephem.Uranus(),
        ephem.Neptune(),
    ]

    for sol_object in sol_objects:
        sol_object.compute(observer)
        data = {
            'ra': float(sol_object.a_ra),
            'dec': float(sol_object.a_dec),
            'vmag': float(sol_object.mag),
        }
        stars = stars.append(data, ignore_index=True)


    stars['azimuth'], stars['altitude'] = equatorial2horizontal(
        stars.ra, stars.dec, observer
    )
    stars = stars.query('altitude > {}'.format(np.deg2rad(min_altitude - 5)))

    return stars


def dark_spot_gridsearch(observer,
                         stars,
                         min_altitude,
                         alt_resolution=0.5,
                         az_resolution=4,
                         ):

    azs = np.deg2rad(np.arange(0, 360, az_resolution))
    alts = np.deg2rad(np.arange(min_altitude, 91, 0.5))
    light = []
    coords = []

    prog = ProgressBar(maxval=len(azs) * len(alts)).start()
    for i, az in enumerate(azs):
        for j, alt in enumerate(alts):
            coords.append((az, alt))
            light.append(light_in_fov(az, alt, stars, observer))
            prog.update(i*len(alts) + j)

    light = np.array(light)
    coords = np.array(coords)
    azs = coords[:, 0]
    alts = coords[:, 1]

    min_index = np.argmin(light)
    az = azs[min_index]
    alt = alts[min_index]

    ra, dec = observer.radec_of(az, alt)

    return az, alt, ra, dec


def plot_dark_spot(stars, az, alt, min_altitude):
    import matplotlib.pyplot as plt
    from cartopy import crs

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=crs.NorthPolarStereo())
    ax.background_patch.set_facecolor('black')

    ax.set_extent([-180, 180, min_altitude - 5, 90], crs.PlateCarree())

    ax.plot(np.rad2deg(az), np.rad2deg(alt), 'ro', transform=crs.PlateCarree())

    plot = ax.scatter(
        np.rad2deg(stars.azimuth),
        np.rad2deg(stars.altitude),
        c=stars.vmag,
        vmax=15,
        s=0.5*(-stars.vmag + stars.vmag.max())**2,
        transform=crs.PlateCarree(),
        cmap='gray_r',
    )

    ax.gridlines(color='blue', ylocs=[min_altitude,])
    fig.colorbar(plot, ax=ax, label='visual magnitude')

    fig.tight_layout()
    plt.show()


def main():

    if args['<date>']:
        date = args['<date>']
        time = args['<timeutc>']
    else:
        date, time = enter_datetime()

    valid = False
    while not valid:
        try:
            date = datetime.strptime(date + ' ' + time, '%Y-%m-%d %H:%M')
            valid = True
        except ValueError:
            print('Could not parse date/time, please use the given notation\n')
            date, time = enter_datetime()

    fact = fact_setup(date)
    stars = create_dataframe(fact)
    az, alt, ra, dec = dark_spot_gridsearch(fact, stars, min_altitude)


    print(u'best ratescan position:')
    print(u'RA: {:1.3f} h'.format(np.rad2deg(ra) * 24/360))
    print(u'DEC: {:1.3f}°'.format(np.rad2deg(dec)))
    print(u'Az: {:1.3f}°'.format(np.rad2deg(az)))
    print(u'Alt: {:1.3f}°'.format(np.rad2deg(alt)))

    if args['--plot']:
        plot_dark_spot(stars, az, alt, min_altitude)



if __name__ == '__main__':
    try:
        main()
    except (KeyboardInterrupt, SystemExit):
        pass
