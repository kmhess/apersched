# make_imaging_sched: Make a schedule for Apertif imaging
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)

__author__ = "Kelley M. Hess"
__date__ = "$19-feb-2019 16:00:00$"
__version__ = "0.3"

import datetime

from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astropy import units as u
import numpy as np

from .calc_slewtime import calc_slewtime

westerbork = EarthLocation(lat=52.91460037*u.deg, lon=6.60449982*u.deg, height=82.2786*u.m)
names = ['3C138','3C147','CTD93']
calibrators = [SkyCoord.from_name(name) for name in names]

dowait = 5         # number of minutes to wait before checking source availability

_int,priority,lo,sub1,_type,weight,beam,sub2,freq1,freq2,freqcent,intent,person,switch_type = '30', \
     'A','4800','64','T','compound','0','320','1250.000','1450.000','1350.000','compound','KH','-'

###################################################################
# Required functions for observing

def observe_calibrator(obstimeUTC, obstime=15):
    return obstimeUTC + datetime.timedelta(minutes=obstime)

def wait_for_rise(obstimeUTC, waittime=5):
    return obstimeUTC + datetime.timedelta(minutes=waittime)
