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
from blessings import Terminal
from datetime import datetime
import numpy as np

from . import (
    fact_setup,
    create_dataframe,
    dark_spot_gridsearch,
    get_stars_in_fov,
    plot_dark_spot
)

term = Terminal()


def enter_datetime():
    ''' a small cli utility to get date and time from the user '''
    print('\nPlease enter date and time for the ratescan')
    print(term.red('This is the real date, be aware for times after 0:00'))
    date = input('Date (YYYY-MM-DD): ')
    time = input('Time UTC: (hh:mm): ')
    return date, time


def main():
    args = docopt(__doc__)
    try:
        min_altitude = 90 - float(args['--max-zenith'])
        max_zenith = float(args['--max-zenith'])

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
        stars = create_dataframe(fact, min_altitude)
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
            np.rad2deg(darkspot['ra']) * 24 / 360,
            np.rad2deg(darkspot['dec']),
        ))

        if args['--plot']:
            plot_dark_spot(stars, darkspot, data, min_altitude, args['--show-flux'])

    except (KeyboardInterrupt, SystemExit):
        pass


if __name__ == '__main__':
        main()
