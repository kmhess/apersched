# functions: useful things that I pulled out of the make_*_sched.py to make things less messy.
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)

__author__ = "Kelley M. Hess"
__date__ = "$19-feb-2019 16:00:00$"
__version__ = "0.3"

import datetime

from .calibrators import *

_int, priority, lo, sub1, _type, weight, beam, sub2, freq1, freq2, freqcent, intent, person, switch_type, freqmode = \
     '30', 'A', '4800', '64', 'T', 'compound', '0', '320', '1250.000', '1450.000', '1350.000', 'compound', 'KH', '-', '300'
# if args.calib_all_beams: _int = 10

###################################################################
# Required functions for observing and writing observations to csv file.


def write_to_csv(writer, source_name, source_pos, start_datetime, end_datetime, pulsar=False):
    # scan='{}{:03d}'.format(start_datetime.strftime('%Y%m%d'),i)
    # source='{}_{}'.format(source_name.split('_')[0],start_datetime.strftime('%Y%m%d'))
    source = '{}'.format(source_name.split('_')[0])
    ra = str(source_pos.to_string('hmsdms').split(' ')[0]).replace('h', ':').replace('m', ':').replace('s', '')
    dec = str(source_pos.to_string('hmsdms').split(' ')[1]).replace('d', ':').replace('m', ':').replace('s', '')
    date1, time1 = start_datetime.strftime('%Y-%m-%d'), start_datetime.strftime('%H:%M:%S')
    date2, time2 = end_datetime.strftime('%Y-%m-%d'), end_datetime.strftime('%H:%M:%S')
    all_cols=[source, ra, '', dec, date1, time1, date2, time2, '1400', 'square_39p1', '0', '39', str(pulsar), '0']
    writer.writerow(all_cols)


def observe_calibrator(obstimeUTC, obstime=20):
    return obstimeUTC + datetime.timedelta(minutes=obstime)


def observe_target(fields, obstimeUTC, name, obstime=3):
    # Use name to subtract 1 value from the weights column!
    fields['weights'][fields['name'] == name] -= 1
    return obstimeUTC + datetime.timedelta(hours=obstime)


def wait_for_rise(obstimeUTC, waittime=5):
    return obstimeUTC + datetime.timedelta(minutes=waittime)
