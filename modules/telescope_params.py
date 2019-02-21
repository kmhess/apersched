from astropy.coordinates import EarthLocation
from astropy import units as u


def westerbork():
    westerbork = EarthLocation(lat=52.91460037 * u.deg, lon=6.60449982 * u.deg, height=82.2786 * u.m)
    return westerbork
