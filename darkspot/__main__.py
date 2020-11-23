#!/usr/bin/env python
import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import EarthLocation
import click

from . import (
    create_catalogue,
    dark_spot_gridsearch,
    get_stars_in_fov,
    plot_dark_spot,
)


@click.command()
@click.option('-t', '--time', type=Time, help='ISO8601 datestring for the observation time. If not given, now is used.')
@click.option('--site', help='Site name (must be known to astropy)', default='Roque de los Muchachos')
@click.option('--lat', type=float, help='Latitude of observatory in decimal degrees, used to override --site')
@click.option('--lon', type=float, help='Longitude of observatory in decimal degrees, used to override --site')
@click.option('--min-altitude', type=float, help='Minimum altitude to consider in degrees', default=70)
@click.option('--fov', type=float, help='Diameter of field of view in degrees', default=4.5)
@click.option('--plot', type=float, help='Show plot', is_flag=True)
@click.option('--band', type=click.Choice(['V', 'BT']), default='BT', help='Which optical band to use.')
def main(time, site, lat, lon, min_altitude, fov, plot, band):
    min_altitude *= u.deg
    fov *= u.deg

    if time is None:
        time = Time.now()

    if lat is not None and lon is not None:
        lon *= u.deg
        lat *= u.deg
        location = EarthLocation(lon=lon, lat=lat)
    else:
        location = EarthLocation.of_site(site)

    stars = create_catalogue(time, location, min_altitude=min_altitude - fov)
    darkspot_table = dark_spot_gridsearch(stars, min_altitude=min_altitude, fov=fov, band=band)

    best = np.argmin(darkspot_table['total_light'])
    darkspot = darkspot_table[best]

    stars_fov = get_stars_in_fov(darkspot['altaz'], stars)

    print('Darkest position:')
    print('RA: {:3.2f}°'.format(darkspot['icrs'].ra.deg))
    print('DEC: {:3.2f}°'.format(darkspot['icrs'].dec.deg))
    print('Alt: {:3.2f}°'.format(darkspot['altaz'].alt.deg))
    print('Az: {:3.2f}°'.format(darkspot['altaz'].az.deg))
    print('Brightest star in FOV: {:1.2f}'.format(stars_fov['BTmag'].min()))

    if plot:
        plot_dark_spot(stars, darkspot, darkspot_table, min_altitude, fov=fov)

        import matplotlib.pyplot as plt
        plt.show()


if __name__ == '__main__':
    main()
