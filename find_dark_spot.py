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
    --show-flux               show the flux as overlay
'''
from __future__ import division, print_function
from docopt import docopt

args = docopt(__doc__)
min_altitude = 90 - float(args['--max-zenith'])
max_zenith = float(args['--max-zenith'])

import pandas as pd
import numpy as np
from numpy import sin, cos, tan, arctan2, arcsin, pi, arccos
import ephem
from six.moves import input
from datetime import datetime
from blessings import Terminal
from tqdm import tqdm
term = Terminal()


def fact_setup(date):
    ''' creates an ephem.Observer for FACT at given date '''
    fact = ephem.Observer()
    fact.lon = '-17:53:28'
    fact.lat = '28:45:42'
    fact.elevation = 2200
    fact.date = date
    fact.epoch = ephem.J2000

    return fact


def enter_datetime():
    ''' a small cli utility to get date and time from the user '''
    print('\nPlease enter date and time for the ratescan')
    print(term.red('This is the real date, be aware for times after 0:00'))
    date = input('Date (YYYY-MM-DD): ')
    time = input('Time UTC: (hh:mm): ')
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

    azs = np.deg2rad(np.arange(0, 360+0.1*az_resolution, az_resolution))
    alts = np.arange(90, min_altitude, -alt_resolution)
    alts = np.deg2rad(alts)
    light = []
    coords = []

    with tqdm(total=len(azs) * len(alts)) as pbar:
        for i, az in enumerate(azs):
            for j, alt in enumerate(alts):
                coords.append((az, alt))
                light.append(light_in_fov(az, alt, stars, observer))
                pbar.update(1)

    light = np.array(light)
    coords = np.array(coords)
    azs_flat = coords[:, 0]
    alts_flat = coords[:, 1]

    min_index = np.argmin(light)
    az = azs_flat[min_index]
    alt = alts_flat[min_index]

    ra, dec = observer.radec_of(az, alt)
    darkspot = {'ra': ra, 'dec': dec, 'az': az, 'alt': alt}
    darkspot_data = {'az': azs,  'alt': alts, 'flux': light}

    return darkspot, darkspot_data


def plot_dark_spot(stars, darkspot, darkspot_data, min_altitude):
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from cartopy import crs
    except ImportError:
        print('You need matplotlib and cartopy to plot')
        return

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=crs.NorthPolarStereo())
    ax.background_patch.set_facecolor('black')

    ax.set_extent([-180, 180, min_altitude - 5, 90], crs.PlateCarree())
    ax.plot(
        np.rad2deg(darkspot['az']),
        np.rad2deg(darkspot['alt']),
        'ro',
        transform=crs.PlateCarree(),
    )

    if args['--show-flux']:
        palt, paz = np.meshgrid(
            np.rad2deg(darkspot_data['alt']),
            np.rad2deg(darkspot_data['az']),
        )
        plot = plt.pcolormesh(
            paz,
            palt,
            np.reshape(darkspot_data['flux'], (len(darkspot_data['az']), -1)),
            cmap='gnuplot',
            norm=LogNorm(),
            transform=crs.PlateCarree()
        )
        fig.colorbar(plot, ax=ax, label='light flux / a.u.')

    plot = ax.scatter(
        np.rad2deg(stars.azimuth),
        np.rad2deg(stars.altitude),
        c=stars.vmag,
        vmax=15,
        s=0.5*(-stars.vmag + stars.vmag.max())**2,
        transform=crs.PlateCarree(),
        cmap='gray_r',
        linewidth=0,
    )

    ax.text(
        0,
        min_altitude - 0.2,
        u'{:2.1f}˚ zenith distance'.format(max_zenith),
        va='top',
        ha='center',
        color='blue',
        transform=crs.PlateCarree(),
    )

    # draw fov, ugly
    fov_az, fov_alt = np.meshgrid(
        np.linspace(0, 2*pi, 400),
        np.linspace(np.deg2rad(min_altitude - 10), pi/2, 200)
    )
    dist = angular_distance(fov_az, fov_alt, darkspot['az'], darkspot['alt'])
    ax.contour(np.rad2deg(fov_az),
               np.rad2deg(fov_alt),
               dist,
               [np.deg2rad(2.25), ],
               colors='red',
               transform=crs.PlateCarree(),
               )

    limit_az = np.linspace(0, 360, 100)
    plt.plot(
        limit_az, np.ones_like(limit_az) * min_altitude,
        'b-', transform=crs.PlateCarree(),
    )

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
    darkspot, data = dark_spot_gridsearch(fact, stars, min_altitude, 0.25, 2)
    stars_fov = get_stars_in_fov(darkspot['az'], darkspot['alt'], stars, fact)

    print(u'best ratescan position:')
    print(u'RA: {:2.2f} h'.format(np.rad2deg(darkspot['ra']) * 24/360))
    print(u'DEC: {:2.2f}°'.format(np.rad2deg(darkspot['dec'])))
    print(u'Az: {:2.1f}°'.format(np.rad2deg(darkspot['az'])))
    print(u'Alt: {:2.1f}°'.format(np.rad2deg(darkspot['alt'])))
    print(u'Brightest star in FOV: {:1.2f} mag'.format(stars_fov.vmag.min()))

    print('\nOutput for FACT schedule:')
    print('"ra":{:.3f}, "dec": {:.3f}'.format(
        np.rad2deg(darkspot['ra']) * 24/360,
        np.rad2deg(darkspot['dec']),
    ))

    if args['--plot']:
        plot_dark_spot(stars, darkspot, data, min_altitude)


if __name__ == '__main__':
    try:
        main()
    except (KeyboardInterrupt, SystemExit):
        pass
