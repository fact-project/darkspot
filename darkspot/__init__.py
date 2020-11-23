import numpy as np
from pkg_resources import resource_filename
from tqdm import tqdm

import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz
from astropy.table import QTable
import ephem
from astropy.coordinates.angle_utilities import angular_separation


# astropy cannot compute magnitude, so we use ephem
SOLAR_SYSTEM_OBJECTS = [
    ephem.Mercury(),
    ephem.Venus(),
    ephem.Moon(),
    ephem.Mars(),
    ephem.Jupiter(),
    ephem.Saturn(),
    ephem.Uranus(),
    ephem.Neptune(),
]


def get_stars_in_fov(pointing, stars, fov=4.5 * u.deg):
    '''
    returns all the stars which are in FACTs field of view pointing
    for given coords az, alt
    '''
    return stars[stars['altaz'].separation(pointing) <= (fov / 2)]


def get_total_light(stars_in_fov, band='BT'):
    ''' returns the total star light in the fov at azimuth az, altitude alt '''
    return np.sum(10**(- stars_in_fov[band + 'mag'].to_value(u.mag) / 2.5))


def create_catalogue(obstime, location, min_altitude=0 * u.deg):
    stars = QTable.read(resource_filename('darkspot', 'hipparcos_bmag_10.fits.gz'))

    # prevent value cannot be set without precision loss error
    for col in stars.colnames:
        if not stars[col].dtype.isnative:
            stars[col] = stars[col].byteswap().newbyteorder()

        if stars[col].dtype == np.float32:
            stars[col] = stars[col].astype(np.float64)

    stars['icrs'] = SkyCoord(ra=stars['RAICRS'], dec=stars['DEICRS'])

    for sol_object in SOLAR_SYSTEM_OBJECTS:
        sol_object.compute(obstime.datetime)
        coord = SkyCoord(
            ra=float(sol_object.ra),
            dec=float(sol_object.dec),
            unit=u.rad,
        )
        mag = float(sol_object.mag) << u.mag

        stars.add_row({
            'icrs': coord,
            'RAICRS': coord.ra, 'DEICRS': coord.dec,
            'Vmag': mag, 'BTmag': mag,
            'HIP': -1, 'recno': -1,
        })

    stars['altaz'] = stars['icrs'].transform_to(
        AltAz(obstime=obstime, location=location)
    )

    return stars[stars['altaz'].alt >= min_altitude]


def dark_spot_gridsearch(
    stars,
    min_altitude,
    fov=4.5 * u.deg,
    alt_resolution=0.5 * u.deg,
    az_resolution=4 * u.deg,
    band='BT',
):

    az = np.arange(0, 360, az_resolution.to_value(u.deg)) << u.deg
    alt = np.arange(90, min_altitude.to_value(u.deg), -alt_resolution.to_value(u.deg)) << u.deg
    alt, az = np.meshgrid(alt, az)
    pointings = SkyCoord(alt=alt.ravel(), az=az.ravel(), frame=stars['altaz'].frame)

    light = []

    for pointing in tqdm(pointings):
        stars_in_fov = get_stars_in_fov(pointing, stars, fov)
        light.append(get_total_light(stars_in_fov, band=band))

    light = np.array(light)

    table = QTable({
        'icrs': pointings.transform_to('icrs'),
        'altaz': pointings,
        'total_light': light,
    })

    return table


def plot_dark_spot(stars, darkspot, darkspot_data, min_altitude, show_flux=True, fov=4.5 * u.deg):
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from cartopy import crs
    except ImportError:
        raise ImportError('You need matplotlib and cartopy to plot') from None

    min_altitude = min_altitude.to_value(u.deg)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=crs.NorthPolarStereo())
    ax.background_patch.set_facecolor('black')

    ax.set_extent([-180, 180, min_altitude - fov.to_value(u.deg), 90], crs.PlateCarree())
    ax.plot(
        darkspot['altaz'].az.deg,
        darkspot['altaz'].alt.deg,
        'ro',
        transform=crs.PlateCarree(),
    )

    if show_flux:
        alt = darkspot_data['altaz'].alt.deg
        az = darkspot_data['altaz'].az.deg
        n_alt = len(np.unique(alt))
        n_az = len(np.unique(az))

        shape = (n_az, n_alt)
        palt = alt.reshape(shape)
        paz = az.reshape(shape)
        light = darkspot_data['total_light'].reshape(shape)

        plot = plt.pcolormesh(
            paz,
            palt,
            light,
            cmap='inferno',
            norm=LogNorm(),
            transform=crs.PlateCarree()
        )
        fig.colorbar(plot, ax=ax, label='light flux / a.u.')

    # draw fov, ugly
    fov_az, fov_alt = np.meshgrid(
        np.linspace(0, 360, 400),
        np.linspace(min_altitude - 10, 90, 200)
    ) << u.deg
    dist = angular_separation(
        fov_az, fov_alt, darkspot['altaz'].az, darkspot['altaz'].alt
    )
    ax.contour(
        fov_az.to_value(u.deg),
        fov_alt.to_value(u.deg),
        dist.to_value(u.deg),
        [fov.to_value(u.deg) / 2],
        colors='red',
        transform=crs.PlateCarree(),
    )

    mag = stars['BTmag'].to_value(u.mag)
    ax.scatter(
        stars['altaz'].az.deg,
        stars['altaz'].alt.deg,
        c=mag,
        vmax=15,
        s=0.5*(-mag + mag.max())**2,
        transform=crs.PlateCarree(),
        cmap='gray_r',
        linewidth=0,
    )

    ax.text(
        0,
        min_altitude - 0.2,
        u'{:2.1f}˚ zenith distance'.format(90 - min_altitude),
        va='top',
        ha='center',
        color='blue',
        transform=crs.PlateCarree(),
    )

    limit_az = np.linspace(0, 360, 100)
    ax.plot(
        limit_az, np.full_like(limit_az, min_altitude),
        'b-', transform=crs.PlateCarree(),
    )
