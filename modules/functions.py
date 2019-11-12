# functions: useful things that I pulled out of the make_*_sched.py to make things less messy.
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)

__author__ = "Kelley M. Hess"
__date__ = "$11-jul-2019 16:00:00$"
__version__ = "0.5"

import datetime

from .calibrators import *

_int, lo, sub1, _type, weight, beam, sub2, centfreq, intent, switch_type, freqmode = \
     '30', '4800', '64', 'T', 'compound', '0', '320', '1280', 'compound', '-', '300'
template = '/opt/apertif/share/parsets/parset_start_observation_atdb_SubbandPhaseCorrection.template'

###################################################################
# Required functions for observing and writing observations to csv file.


def write_to_csv(writer, source_name, source_pos, start_datetime, end_datetime, pulsar=False):
    # scan='{}{:03d}'.format(start_datetime.strftime('%Y%m%d'),i)
    # source='{}_{}'.format(source_name.split('_')[0],start_datetime.strftime('%Y%m%d'))
    source = source_name #'{}'.format(source_name.split('_')[0])
    ra = str(source_pos.to_string('hmsdms').split(' ')[0]).replace('h', ':').replace('m', ':').replace('s', '')
    dec = str(source_pos.to_string('hmsdms').split(' ')[1]).replace('d', ':').replace('m', ':').replace('s', '')
    date1, time1 = start_datetime.strftime('%Y-%m-%d'), start_datetime.strftime('%H:%M:%S')
    date2, time2 = end_datetime.strftime('%Y-%m-%d'), end_datetime.strftime('%H:%M:%S')
    if (source_name in flux_names) or (source_name in pol_names):
        all_cols = [source, ra, '', dec, date1, time1, date2, time2, '10', 'S*', weight, beam, 'system', freqmode, centfreq, template]
    elif ('S' in source_name) or ('M' in source_name):
        all_cols = [source, ra, '', dec, date1, time1, date2, time2, _int, _type, weight, beam, switch_type, freqmode, centfreq, template]
    elif (source_name == 'imaging_start') or (source_name == 'imaging_end'):
        all_cols = [source, ra, '', dec, date1, time1, date2, time2, '10', _type, weight, beam, switch_type, freqmode,
                    centfreq, template]
    else:
        all_cols = [source, ra, '', dec, date1, time1, date2, time2, '1280', 'square_39p1', '0', '39', str(pulsar), '0']
    writer.writerow(all_cols)


def observe_calibrator(obstimeUTC, obstime=20):
    return obstimeUTC + datetime.timedelta(minutes=obstime)


def observe_target(fields, obstimeUTC, name, obstime=3):
    # Use name to subtract 1 value from the weights column!
    fields['weights'][fields['name'] == name] -= 1
    return obstimeUTC + datetime.timedelta(hours=obstime)


def wait_for_rise(obstimeUTC, waittime=5):
    return obstimeUTC + datetime.timedelta(minutes=waittime)
