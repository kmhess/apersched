# make_imaging_sched: Make a schedule for Apertif imaging
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)
__author__ = "Kelley M. Hess"
__date__ = "$18-jul-2019 16:00:00$"
__version__ = "1.1.1"

import csv
import datetime
import os

from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.coordinates import Longitude, SkyCoord, get_sun, get_moon
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time
from astropy import units as u
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
import numpy as np

#import atdbquery
from modules.calc_slewtime import calc_slewtime  # Wants [ra,dec] start/end positions in radians; outputs seconds.
from modules.calibrators import *
from modules.functions import *
from modules.telescope_params import westerbork


# flux_names = ['3C147']
# flux_cal = [SkyCoord.from_name(name) for name in flux_names]
# pol_names = ['3C138']
# pol_cal = [SkyCoord.from_name(name) for name in pol_names]

###################################################################
# Survey specific functions for doing observations and calibration

# def do_calibration_40b(i, obstime_utc, telescope_position, csvfile, total_wait, next_cal, mins_per_beam):
# 
#     calib_sun_dist = 0.
# 
#     current_lst = Time(obstime_utc).sidereal_time('apparent', westerbork().lon)
#     # Consider HA limits for shadowing:
#     #     https://old.astron.nl/radio-observatory/astronomers/wsrt-guide-observations/3-telescope-parameters-and-array-configuration
#     # Note ha_limit[1:2] depend on length of calibration!
#     syswait = 2.0  # minutes
#     obstime = (mins_per_beam + syswait) * 40. - syswait  # minutes
#     sun_position = get_sun(Time(obstime_utc, scale='utc'))
#     if (i == 4):
#         next_cal = 'flux'
#         print(next_cal)
#     if next_cal == 'flux':
#         calibrators = flux_cal
#         names = flux_names
#         type_cal = 'Flux'
#         ha_limit = [-5.0, 5.0 - obstime / 60., 0.4]  # Entry 2 is hardcoded for 5 min per beam
#     if next_cal == 'pol':
#         calibrators = pol_cal
#         names = pol_names
#         type_cal = 'Polarization'
#         ha_limit = [-3.2, 3.3 - obstime / 60., - 1.0]  # Entry 2 is hardcoded for 5 min per beam
#     is_cal_up = np.array([(current_lst.hour - calibrators[0].ra.hour > ha_limit[0]) and (current_lst.hour - calibrators[0].ra.hour < ha_limit[1])])#,
#                  # (current_lst.hour - calibrators[1].ra.hour > ha_limit[0]) and (current_lst.hour - calibrators[1].ra.hour < ha_limit[2])])
#     is_sundist_okay = np.array([sun_position.separation(calibrators[0]).value > calib_sun_dist])#,
#                                 # sun_position.separation(calibrators[1]).value > calib_sun_dist])
#     calib_wait = 0
#     new_obstime_utc = obstime_utc
#     print(names)
# 
#     # Wait for calibrator to be at least an hour above the observing horizon, or not shadowed.
#     while not np.any(is_cal_up * is_sundist_okay):  # and (calib_wait < 6. * 60.): # and (not is3c286):
#         calib_wait += dowait
#         new_obstime_utc = wait_for_rise(new_obstime_utc, waittime=dowait)
#         new_lst = Time(new_obstime_utc).sidereal_time('apparent', westerbork().lon)
#         is_cal_up = [(new_lst.hour - calibrators[0].ra.hour > ha_limit[0]) and (new_lst.hour - calibrators[0].ra.hour < ha_limit[1])]#,
#                      # (new_lst.hour - calibrators[1].ra.hour > ha_limit[0]) and (new_lst.hour - calibrators[1].ra.hour < ha_limit[2])]
#     n = np.where(is_cal_up)[0][0]
#     ##### EDITABLE: Can change the number of hours the program will wait for a calibrator #####
#     if calib_wait != 0 and calib_wait < 4.0 * 60:
#         total_wait += calib_wait
#         n = np.where(is_cal_up)[0][0]
#         print("\tCalibrator not up, waiting {} minutes until LST: {}.".format(calib_wait, str(new_lst)))
#     # The commented part is hopefully obsolete with the new calibrator.py and observing strategy, but there's still a bad starting point in the sky
#     # elif calib_wait >= 4.0 * 60:
#     elif calib_wait >= 10.0 * 60:
#         after_cal = obstime_utc - datetime.timedelta(minutes=syswait)
#         new_telescope_position = telescope_position
#         i -= 1
#         # If can't do any pol at beginning, first must be a flux cal
#         # (if statement untested; trying to fix skipping pol cal when this happens in the middle!):
#         if i == 1:
#             next_cal = 'flux'
#         print("Must wait {} hours for calibrator to rise.  Instead, go directly to target.".format(calib_wait/60.))
#         print("\tIf this appears anywhere other than beginning of scheduling block, probably need to expanding "
#               "options in pointing file.")
#         # break
#         return i, after_cal, new_telescope_position, total_wait, next_cal
# 
#     slew_seconds = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian],
#                                      [calibrators[n].ra.radian, calibrators[n].dec.radian])
#     if slew_seconds < calib_wait * 60.:
#         new_obstime_utc = obstime_utc + datetime.timedelta(minutes=calib_wait)
#     else:
#         new_obstime_utc = obstime_utc + datetime.timedelta(seconds=slew_seconds)
# 
#     # Calculate appropriate observe time for the calibrator and observe it.
#     obstime = (mins_per_beam + syswait) * 40. - syswait    # <mins_per_beam> minutes per beam, 2 min wait
#     if n == 1:
#         obstime = (5.0 + syswait) * 40. - syswait  # force 5 minutes per beam, 2 min wait on calibs with natural gap before target
#     after_cal = observe_calibrator(new_obstime_utc, obstime=obstime)
#     if i == 1:
#         write_to_csv(csvfile, 'imaging_start', calibrators[n], new_obstime_utc - datetime.timedelta(minutes=3.0),
#                      new_obstime_utc - datetime.timedelta(minutes=2.0))
#     write_to_csv(csvfile, names[n], calibrators[n], new_obstime_utc, after_cal)
# 
#     print("Scan {} observed {} calibrator {}.".format(i, type_cal, names[n]))
#     check_sun = sun_position.separation(calibrators[n])
#     if check_sun.value < 30.0: #calib_sun_dist:
#         print("\tWARNING: {} is THIS close to Sun: {:5.2f}".format(names[n], check_sun))
#     new_telescope_position = calibrators[n]
# 
#     # Set up for next calibrator
#     if next_cal == 'flux':
#         next_cal = 'pol'
#     else:
#         next_cal = 'flux'
# 
#     return i, after_cal, new_telescope_position, total_wait, next_cal

def do_calibration_40b(i, obstime_utc, telescope_position, csvfile, total_wait, next_cal, mins_per_beam):

    calib_sun_dist = 30.

    current_lst = Time(obstime_utc).sidereal_time('apparent', westerbork().lon)
    # Consider HA limits for shadowing:
    #     https://old.astron.nl/radio-observatory/astronomers/wsrt-guide-observations/3-telescope-parameters-and-array-configuration
    # Note ha_limit[1:2] depend on length of calibration!
    syswait = 2.0  # minutes
    obstime = (mins_per_beam + syswait) * 40. - syswait  # minutes
    sun_position = get_sun(Time(obstime_utc, scale='utc'))
    if next_cal == 'flux':
        calibrators = flux_cal
        names = flux_names
        type_cal = 'Flux'
        ha_limit = [-5.0, 5.0 - obstime / 60., 0.4]  # Entry 2 is hardcoded for 5 min per beam
    if next_cal == 'pol':
        calibrators = pol_cal
        names = pol_names
        type_cal = 'Polarization'
        ha_limit = [-3.3, 3.3 - obstime / 60., - 1.0]  # Entry 2 is hardcoded for 5 min per beam
    is_cal_up = np.array([(current_lst.hour - calibrators[0].ra.hour > ha_limit[0]) and (current_lst.hour - calibrators[0].ra.hour < ha_limit[1]),
                 (current_lst.hour - calibrators[1].ra.hour > ha_limit[0]) and (current_lst.hour - calibrators[1].ra.hour < ha_limit[2])])
    is_sundist_okay = np.array([sun_position.separation(calibrators[0]).value > calib_sun_dist,
                                sun_position.separation(calibrators[1]).value > calib_sun_dist])
    calib_wait = 0
    new_obstime_utc = obstime_utc

    # Wait for calibrator to be at least an hour above the observing horizon, or not shadowed.
    while not np.any(is_cal_up * is_sundist_okay):  # and (calib_wait < 6. * 60.): # and (not is3c286):
        calib_wait += dowait
        new_obstime_utc = wait_for_rise(new_obstime_utc, waittime=dowait)
        new_lst = Time(new_obstime_utc).sidereal_time('apparent', westerbork().lon)
        is_cal_up = [(new_lst.hour - calibrators[0].ra.hour > ha_limit[0]) and (new_lst.hour - calibrators[0].ra.hour < ha_limit[1]),
                     (new_lst.hour - calibrators[1].ra.hour > ha_limit[0]) and (new_lst.hour - calibrators[1].ra.hour < ha_limit[2])]
    n = np.where(is_cal_up * is_sundist_okay)[0][0]
    ##### EDITABLE: Can change the number of hours the program will wait for a calibrator #####
    # if calib_wait != 0 and calib_wait < 4.0 * 60:
    if calib_wait != 0 and calib_wait < 9.0 * 60:
        total_wait += calib_wait
        # n = np.where(is_cal_up * is_sundist_okay)[0][0]
        print("\tCalibrator not up, waiting {} minutes until LST: {}.".format(calib_wait, str(new_lst)))
    # The commented part is hopefully obsolete with the new calibrator.py and observing strategy, but there's still a bad starting point in the sky
    # elif calib_wait >= 4.0 * 60:
    elif calib_wait >= 9.0 * 60:
        after_cal = obstime_utc - datetime.timedelta(minutes=syswait)
        new_telescope_position = telescope_position
        i -= 1
        # If can't do any pol at beginning, first must be a flux cal
        # (if statement untested; trying to fix skipping pol cal when this happens in the middle!):
        if i == 1:
            next_cal = 'flux'
        print("Must wait {} hours for calibrator to rise.  Instead, go directly to target.".format(calib_wait/60.))
        print("\tIf this appears anywhere other than beginning of scheduling block, probably need to expanding "
              "options in pointing file.")
        # break
        return i, after_cal, new_telescope_position, total_wait, next_cal

    slew_seconds = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian],
                                     [calibrators[n].ra.radian, calibrators[n].dec.radian])
    if slew_seconds < calib_wait * 60.:
        new_obstime_utc = obstime_utc + datetime.timedelta(minutes=calib_wait)
    else:
        new_obstime_utc = obstime_utc + datetime.timedelta(seconds=slew_seconds)

    # Calculate appropriate observe time for the calibrator and observe it.
    obstime = (mins_per_beam + syswait) * 40. - syswait    # <mins_per_beam> minutes per beam, 2 min wait
    # if (n == 1) & (next_cal == 'pol'):
    #     obstime = (5.0 + syswait) * 40. - syswait  # force 5 minutes per beam, 2 min wait on calibs with natural gap before target
    after_cal = observe_calibrator(new_obstime_utc, obstime=obstime)
    if i == 1:
        write_to_csv(csvfile, 'imaging_start', calibrators[n], new_obstime_utc - datetime.timedelta(minutes=3.0),
                     new_obstime_utc - datetime.timedelta(minutes=2.0))
    write_to_csv(csvfile, names[n], calibrators[n], new_obstime_utc, after_cal)

    print("Scan {} observed {} calibrator {}.".format(i, type_cal, names[n]))
    check_sun = sun_position.separation(calibrators[n])
    if check_sun.value < calib_sun_dist:
        print("\tWARNING: {} is THIS close to Sun: {:5.2f}".format(names[n], check_sun))
    new_telescope_position = calibrators[n]

    # Set up for next calibrator
    if next_cal == 'flux':
        next_cal = 'pol'
    else:
        next_cal = 'flux'

    return i, after_cal, new_telescope_position, total_wait, next_cal
    
def do_target_observation(i, obstime_utc, telescope_position, csvfile, total_wait): #, closest_field):
    # Get objects that are close to current observing horizon:
    current_lst = Time(obstime_utc).sidereal_time('apparent', westerbork().lon)

    proposed_ra = (current_lst + Longitude('6h')).wrap_at(360 * u.deg)
    test_slew_seconds = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian],
                                      [proposed_ra.radian, telescope_position.dec.radian])

    avail_fields = apertif_fields[apertif_fields['weights'] > 0]
    availability = SkyCoord(np.array(avail_fields['hmsdms'])).ra.hour - proposed_ra.hour - test_slew_seconds / 3600.
    availability[availability < -12.] += 24.

    sun_position = get_sun(Time(obstime_utc, scale='utc'))
    sun_okay = np.array(sun_position.separation(SkyCoord(avail_fields['hmsdms'])).value)

    targ_wait = 0.     # minutes
    ##### EDITABLE: Will control the maximum wait time for a target before it gives up and looks for a calibrator #####
    # wait_limit = 5.0  # hours
    #wait_limit = 10.0  # hours
    wait_limit = 23.0  # hours

    new_obstime_utc = obstime_utc
    # First check what is *already* up.  If nothing, then wait for something to rise.
    ##### If target is passing it's HA limit, adjust availability numbers (especially the second one) #####
    while not np.any((availability < -0.05) & (availability > -0.40) & (sun_okay > args.sun_distance)):
        targ_wait += dowait
        new_obstime_utc = wait_for_rise(new_obstime_utc, waittime=dowait)
        new_lst = Time(new_obstime_utc).sidereal_time('apparent', westerbork().lon)
        proposed_ra = (new_lst + Longitude('6h')).wrap_at(360 * u.deg)
        availability = SkyCoord(np.array(avail_fields['hmsdms'])).ra.hour - proposed_ra.hour
        availability[availability < -12.] += 24.
        # print(availability, sun_okay > args.sun_distance, sun_okay)
    if (targ_wait != 0.) and (targ_wait <= wait_limit * 60.):
        total_wait += targ_wait
        print("\tTarget not up (or sun issues), waiting {} minutes until LST: {}".format(targ_wait, str(new_lst)))

    if targ_wait <= wait_limit * 60.:
        # Choose M101 field first if available or within a 1 hour wait AND (before this function) if users had requested it in args.repeat_m101
        if 'M1403+5324' in avail_fields[(availability < -0.02) & (availability > -0.40)]['name']:
            first_field = avail_fields[avail_fields['name'] == 'M1403+5324'][0]
            print("*** M1403+5324 OBSERVED!  WAIT A MONTH TO SCHEDULE AGAIN! ***")
        elif 'M1403+5324' in avail_fields[(availability < 1.00) & (availability > -0.02)]['name']:
            m101_field = avail_fields[avail_fields['name'] == 'M1403+5324'][0]
            m101_availability = SkyCoord(m101_field['hmsdms']).ra.hour - proposed_ra.hour
            while not (m101_availability < -0.02) & (m101_availability > -0.40):
                targ_wait += dowait
                new_obstime_utc = wait_for_rise(new_obstime_utc, waittime=dowait)
                new_lst = Time(new_obstime_utc).sidereal_time('apparent', westerbork().lon)
                proposed_ra = (new_lst + Longitude('6h')).wrap_at(360 * u.deg)
                m101_availability = SkyCoord(m101_field['hmsdms']).ra.hour - proposed_ra.hour
                # m101_availability[m101_availability < -12] += 24
            first_field = m101_field
            print("*** M1403+5324 OBSERVED!  WAIT A MONTH TO SCHEDULE AGAIN! ***")
        else:
            first_field = avail_fields[(availability < -0.05) & (availability > -0.40) & (sun_okay > args.sun_distance)][0]
            check_sun = sun_position.separation(SkyCoord(first_field['hmsdms']))
            if check_sun.value < 50.:
                print("\tField is THIS close to Sun: {}".format(check_sun))

        # NOTE SLEW TIME IS CALCULATED TO THE *OBSERVING* HORIZON, NOT TO THE NEW RA!
        # TELESCOPE SHOULD MOVE TO HORIZON AND WAIT!!!
        slew_seconds = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian],
                                     [SkyCoord(first_field['hmsdms']).ra.radian,
                                      SkyCoord(first_field['hmsdms']).dec.radian])
        if slew_seconds < targ_wait * 60.:
            new_obstime_utc = obstime_utc + datetime.timedelta(minutes=targ_wait)
        else:
            new_obstime_utc = obstime_utc + datetime.timedelta(seconds=slew_seconds)
        after_target = observe_target(apertif_fields, new_obstime_utc, first_field['name'], obstime=11.5)
        write_to_csv(csvfile, first_field['name'], SkyCoord(first_field['hmsdms']), new_obstime_utc, after_target)
        print("Scan {} observed {}.".format(i, first_field['hmsdms']))
        return i, after_target, SkyCoord(first_field['hmsdms']), total_wait
    else:
        print("\tNo target for {} hours. Go to a calibrator instead.".format(targ_wait/60.))
        i -= 1
        return i, obstime_utc, telescope_position, total_wait

###################################################################

parser = ArgumentParser(description="Make observing schedule for the Apertif imaging surveys. Saves schedule in CSV file to be parsed by atdbspec. "
                                    "Outputs a png of the completed and scheduled pointings.",
                        formatter_class=RawTextHelpFormatter)

parser.add_argument('-f', '--filename', default='./ancillary_data/apertif_6mos2.03oct19.txt', #all_pointings.v5.29apr19.txt',
                    help='Specify the input file of pointings to choose from (default: %(default)s).')
parser.add_argument('-p', '--previous_obs', default='',
                    help='Specify a file of previously scheduled pointings. (No default.)')
parser.add_argument('-c', '--copy_previous',
                    help='Copy file of previously scheduled pointings and append to it. (Default: False.)',
                    action='store_true')
parser.add_argument('-o', '--output', default='temp',
                    help='Specify the root of output csv and png files. If file exists, append to it (default: imaging_%(default)s.csv).')
parser.add_argument('-b', "--all_beam_calib",
                    help="Default behavior is 15 minutes on a calibrator in the central beam. If option is included, run 40 beam calibration.",
                    action='store_true')
parser.add_argument('-m', "--mins_per_beam", default=3.0,
                    help="Number of minutes for calibrator in 40b scan (default: %(default)s).",
                    type=float)
parser.add_argument('-s', "--starttime_utc", default="2019-07-01 08:00:00",
                    help="The start time in ** UTC ** ! - format 'YYYY-MM-DD HH:MM:SS' (default: '%(default)s').",
                    type=datetime.datetime.fromisoformat)
parser.add_argument('-l', "--schedule_length", default=7.0,
                    help="Number of days to schedule (can be float; default: %(default)s).",
                    type=float)
parser.add_argument('-d', "--sun_distance", default=45.0,
                    help="Minimum allowed distance in decimal degrees to Sun (default: %(default)s).",
                    type=float)
parser.add_argument('-x', '--pol_cal_true',
                    help='Force program to start with a polarization calibrator rather than flux calibrator.',
                    action='store_true')
parser.add_argument('-r', "--repeat_m101",
                    help="If option is included, Try to schedule the M101 field once this time. Works until it's been observed to MDS depth.",
                    action='store_true')
parser.add_argument('-a', "--check_atdb",
                    help="If option is included, *DO NOT* check ATDB for previous observations.",
                    action='store_false')
parser.add_argument('-v', "--verbose",
                    help="If option is included, print updated UTC times after each scan.",
                    action='store_true')

# Parse the arguments above
args = parser.parse_args()

# Filename for the csv file of observed fields:
csv_filename = 'imaging_{}.csv'.format(args.output)

# Filename for the map of observed fields:
filename = 'imaging_map_{}.png'.format(args.output)

# Number of minutes to wait before checking source availability
dowait = 2

# Load all-sky pointing file and select the pointings with the label for the appropriate survey:

# labels: l=lofar; m=medium-deep; s=shallow; t=timing; g=Milky Way +/-5 in galactic latitude
#         h=NCP that will be covered with hexagonal compound beam arrangement
fields = Table(ascii.read(args.filename, format='fixed_width'))
apertif_fields = fields[((fields['label'] == 'm') | (fields['label'] == 'o')) ]#& (fields['ra']>11.*15) & (fields['ra']<18.0*15)]
weights = np.zeros(len(apertif_fields))
##### EDITABLE: Use different labels to control how areas of the Medium-deep are built up in H-ATLAS #####
weights[apertif_fields['label'] == 'm'] = 1
weights[apertif_fields['label'] == 'o'] = 0
weights[apertif_fields['label'] == 'l'] = 0

# Add "weights" column to table.
apertif_fields['weights'] = weights

print("\n##################################################################")
print("Number of all-sky fields are: {}".format(len(fields)))
print("Number of Apertif imaging fields are: {}".format(len(apertif_fields)))

# Retrieve names of observations from ATDB (excludes calibrator scans)
if args.check_atdb:
    observations = atdbquery.atdbquery(obs_mode='imaging', failures=False, transient=False)
    imaging_obs = [dict(observations[i])['name'] for i in range(len(observations))
                   if (dict(observations[i])['name'][0:2] != '3c') and (dict(observations[i])['name'][0:2] != '3C')
                   and (dict(observations[i])['name'][0:2] != 'CT')]
    # Adjust 'weights' field for objects that have been previously observed:
    for obs in imaging_obs:
        if obs in apertif_fields['name']:
            i = np.where(apertif_fields['name'] == obs)
            apertif_fields['weights'][i] -= 1
else:
    print("Not querying ATDB for previous observations.")

# Mostly for testing purposes, read in a list of previously scheduled pointings, remove those as possible fields,
# save info for plotting later.
try:
    scheduled = Table.read(args.previous_obs, format='csv')
    mask = ['3C' not in entry for entry in scheduled['source']]
    scheduled = scheduled[mask]
    scheduled_coords = SkyCoord(ra=scheduled['ra'], dec=scheduled['dec'], unit=(u.hourangle, u.deg))
    print("Reading table of previously scheduled pointings.")
    for sch in scheduled:
        if sch['source'] in apertif_fields['name']:
            i = np.where(apertif_fields['name'] == sch['source'])
            apertif_fields['weights'][i] -= 1
except IOError:
    print("If file of previously scheduled observations was requested, it does not exist. Continuing.")
    scheduled_coords = []

# Add back fields that were deem failed in order to schedule again:
# try:
#     failed = Table.read('./ancillary_data/failed_obs.csv')
#     for fail in failed:
#         if fail['name'] in apertif_fields['name']:
#             i = np.where(apertif_fields['name'] == fail['name'])
#             apertif_fields['weights'][i] += 1
# except IOError:
#     print("No list of failed observations exists. Continuing")

# Try to repeat M101 once at users request by appropriately modifying the weights after they are read in from the pointing file.
if args.repeat_m101:
    if apertif_fields['weights'][apertif_fields['name'] == 'M1403+5324'] > 0:
        apertif_fields['weights'][apertif_fields['name'] == 'M1403+5324'] = 1
    else:
        print("M1403+5324 with M101 & scintillating source already observed to full depth. Will not force to schedule.")
else:
    apertif_fields['weights'][apertif_fields['name'] == 'M1403+5324'] = 0

# Append schedule to end of previously existing file and tell the user what's happening.
if os.path.isfile(args.previous_obs) & args.copy_previous:
    print("-c selected: prepending previous file to new file {}".format(csv_filename))
    os.system("cp {} {}".format(args.previous_obs, csv_filename))
    header = False
elif os.path.isfile(csv_filename) & os.path.isfile(args.previous_obs):
    print("Specified output file exists; appending schedule to previous file {}".format(csv_filename))
    header = False
elif os.path.isfile(csv_filename):
    print("Output file {} exists but not specified with '-p', so overwriting.".format(csv_filename))
    os.remove(csv_filename)
    header = ['source', 'ra', 'ha', 'dec', 'date1', 'time1', 'date2', 'time2', 'int', 'type', 'weight', 'beam',
              'switch_type', 'freqmode', 'centfreq', 'template']
else:
    header = ['source', 'ra', 'ha', 'dec', 'date1', 'time1', 'date2', 'time2', 'int', 'type', 'weight', 'beam',
              'switch_type', 'freqmode', 'centfreq', 'template']

print("Will shift pointings if Sun is within {} degrees.".format(args.sun_distance))

# Estimate the telescope starting position as on the meridian (approximately parked)
telescope_position = SkyCoord(ra=Time(args.starttime_utc).sidereal_time('apparent', westerbork().lon), dec='50d00m00s')
current_lst = Time(args.starttime_utc).sidereal_time('apparent', westerbork().lon)
next_cal = 'flux'
if args.pol_cal_true:
    next_cal = 'pol'
# if (current_lst.hour - pol_cal[1].ra.hour > -5.0) and (current_lst.hour - pol_cal[1].ra.hour < 0.1):
#     next_cal = 'pol'

# Create a record of the positions planned to be observed so we can plot later.
observed_pointings = []

# Keep track of idle time.
total_wait = 0

# Do the observations: select calibrators and target fields, and write the output to a CSV file.
# Also, writes out a record of what is observed, when, and if the telescope has to wait for objects to rise.
print("\nStarting observations! UTC: " + str(args.starttime_utc))
print("                       LST: " + str(current_lst))
print("Calculating schedule for the following " + str(args.schedule_length) + " days.\n")

# Open & prepare CSV file to write parset parameters to, in format given by V.M. Moss.
# (This could probably be done better because write_to_parset is in modules/function.py)
with open(csv_filename, 'a') as csvfile:
    writer = csv.writer(csvfile)
    if header:
        writer.writerow(header)

    # Always start on a calibrator (unless the wait time is too long; modified in do_calibration_40b)
    i = 1
    if (telescope_position.ra.deg < 3.0 * 15) | (telescope_position.ra.deg > 17.0 * 15):  # Need enough time to do calibrator before first field.
        if args.all_beam_calib:
            i, new_obstime_utc, new_position, total_wait, next_cal = \
                do_calibration_40b(i, args.starttime_utc, telescope_position, writer, total_wait, next_cal, args.mins_per_beam)
        else:
            i, new_obstime_utc, new_position, total_wait = \
                do_calibration(i, args.starttime_utc, telescope_position, writer, total_wait)
        obstime_utc = new_obstime_utc + datetime.timedelta(minutes=2.0)
        total_wait += 2
        telescope_position = new_position
        if args.verbose:
            print("\tUTC: " + str(obstime_utc) + ",  LST: " + str(
                Time(obstime_utc).sidereal_time('apparent', westerbork().lon)) + " at end of scan.")
    else:
        obstime_utc = args.starttime_utc

    # Iterate between (HATLAS) target and calibrator for the specified amount of time & write to CSV file:
    while obstime_utc < args.starttime_utc + datetime.timedelta(days=args.schedule_length) - datetime.timedelta(hours=10.):
        # i += 1
        # i, new_obstime_utc, new_position, total_wait, next_cal = \
        #             do_calibration_40b(i, obstime_utc, telescope_position, writer, total_wait, next_cal, args.mins_per_beam)
        # if args.verbose:
        #     print("\tUTC: " + str(new_obstime_utc) + ",  LST: " + str(
        #         Time(new_obstime_utc).sidereal_time('apparent', westerbork().lon)) + " at end of scan.", end="")
        #     print("\tTotal time between end of scans: {:0.4} hours".format(
        #         (new_obstime_utc - obstime_utc).seconds / 3600.))
        # # Add datawriter wait time before considering next observation
        # obstime_utc = new_obstime_utc + datetime.timedelta(minutes=2.0)
        # total_wait += 2
        # telescope_position = new_position
        # observed_pointings.append(telescope_position)

        i += 1
        i, new_obstime_utc, new_position, total_wait = do_target_observation(i, obstime_utc, telescope_position, writer,
                                                                             total_wait)
        if args.verbose:
            print("\tUTC: " + str(new_obstime_utc) + ",  LST: " + str(
                Time(new_obstime_utc).sidereal_time('apparent', westerbork().lon)) + " at end of scan.", end="")
            print("\tTotal time between end of scans: {:0.4} hours".format(
                (new_obstime_utc - obstime_utc).seconds / 3600.))
        # Add datawriter wait time before considering next observation
        obstime_utc = new_obstime_utc + datetime.timedelta(minutes=2.0)
        total_wait += 2
        telescope_position = new_position
        observed_pointings.append(telescope_position)

        if (i==2) & (telescope_position.ra.deg > 4.75 * 15):  # If start w 22h and "past" 3C147, go to target in 12h.
            i += 1
            i, new_obstime_utc, new_position, total_wait = do_target_observation(i, obstime_utc, telescope_position,
                                                                                 writer,
                                                                                 total_wait)
            if args.verbose:
                print("\tUTC: " + str(new_obstime_utc) + ",  LST: " + str(
                    Time(new_obstime_utc).sidereal_time('apparent', westerbork().lon)) + " at end of scan.", end="")
                print("\tTotal time between end of scans: {:0.4} hours".format(
                    (new_obstime_utc - obstime_utc).seconds / 3600.))
            # Add datawriter wait time before considering next observation
            obstime_utc = new_obstime_utc + datetime.timedelta(minutes=2.0)
            total_wait += 2
            telescope_position = new_position
            observed_pointings.append(telescope_position)

        i += 1
        if (next_cal == 'flux'):
            next_cal == 'pol'
        i, new_obstime_utc, new_position, total_wait, next_cal = \
                    do_calibration_40b(i, obstime_utc, telescope_position, writer, total_wait, next_cal, args.mins_per_beam)
        if args.verbose:
            print("\tUTC: " + str(new_obstime_utc) + ",  LST: " + str(
                Time(new_obstime_utc).sidereal_time('apparent', westerbork().lon)) + " at end of scan.", end="")
            print("\tTotal time between end of scans: {:0.4} hours".format(
                (new_obstime_utc - obstime_utc).seconds / 3600.))
        # Add datawriter wait time before considering next observation
        obstime_utc = new_obstime_utc + datetime.timedelta(minutes=2.0)
        total_wait += 2
        telescope_position = new_position
        # observed_pointings.append(telescope_position)

# Write the last line of the csv file
ha = telescope_position.ra.hour - Time(obstime_utc).sidereal_time('apparent', westerbork().lon).hour
with open(csv_filename, 'a') as csvfile:
    writer = csv.writer(csvfile)
    write_to_csv(writer, 'imaging_end', telescope_position, obstime_utc, obstime_utc + datetime.timedelta(minutes=1.0))

print("\nEnding observations! UTC: " + str(obstime_utc))
print("Total number of survey fields observed is: {}".format(len(observed_pointings)))
print("Total time on survey fields is {:3.1f}% of total.".format(len(observed_pointings) * 11.5 * 3600 /
                                                            (obstime_utc - args.starttime_utc).total_seconds() * 100))
print("Total wait time: {} mins is {:3.1f}% of total.".format(total_wait,
                                                              total_wait * 60. / (obstime_utc - args.starttime_utc).total_seconds() * 100))
print("\nThe schedule has been written to " + csv_filename)
print("A map of the observed fields has been written to " + filename)
print("IF THIS IS NOT A REAL SURVEY OBSERVATION BUT WILL BE OBSERVED, EDIT THE TARGET NAMES IN THE csv FILE!")
print("##################################################################\n")

# Calculate ecliptic for plotting
year_arr = Time(args.starttime_utc) + np.linspace(0, 364, 365) * u.day
obs_arr = Time(args.starttime_utc) + np.linspace(0, args.schedule_length, np.int(np.ceil(args.schedule_length * 24))) * u.day
sun_year = get_sun(year_arr)
sun_obs = get_sun(obs_arr)
moon_obs = get_moon(obs_arr)

# Create and save a figure of all pointings selected for this survey, and which have been observed.
# plt.figure(figsize=[8, 8])
# m = Basemap(projection='nplaea', boundinglat=20, lon_0=310, resolution='l', celestial=True)
# m.drawparallels(np.arange(30, 90, 15), labels=[False, False, False, False], color='darkgray')
# m.drawmeridians(np.arange(0, 360, 15), labels=[True, True, False, True], color='darkgray', latmax=90)
# xsun_moll, ysun_moll = m(sun_year.ra.deg, sun_year.dec.deg)
# m.plot(xsun_moll, ysun_moll, 'o-', markersize=2, label='Ecliptic', color='orange')
# xsunobs_moll, ysunobs_moll = m(sun_obs.ra.deg, sun_obs.dec.deg)
# xmoonobs_moll, ymoonobs_moll = m(moon_obs.ra.deg, moon_obs.dec.deg)
# m.plot(xsunobs_moll, ysunobs_moll, 'o', markersize=6, label='Sun', color='orange')
# m.plot(xmoonobs_moll, ymoonobs_moll, 'o', markersize=5, label='Moon', color='gray')
# xpt_ncp, ypt_ncp = m(SkyCoord(np.array(apertif_fields['hmsdms'])).ra.deg,
#                      SkyCoord(np.array(apertif_fields['hmsdms'])).dec.deg)
# m.plot(xpt_ncp, ypt_ncp, 'o', markersize=7, label='SNS', mfc='none', color='0.1')
# # for i, f in enumerate(apertif_fields):
# #     if (f['label'] == 'm') & (f['weights'] != 1) & (f['name'] != "M1403+5324"):
# #         m.plot(xpt_ncp[i], ypt_ncp[i], 'o', markersize=7, mfc='red', color='0')
# #     elif (f['label'] == 's') & (f['weights'] == 0):
# #         m.plot(xpt_ncp[i], ypt_ncp[i], 'o', markersize=7, mfc='red', color='0')
# m.plot(0, 0, 'o', markersize=7, label='Already in ATDB', mfc='red', color='0')
# for p in scheduled_coords:
#     xpt_ncp, ypt_ncp = m(p.ra.deg, p.dec.deg)
#     m.plot(xpt_ncp, ypt_ncp, 'o', markersize=7, mfc='orange', color='0')
# if scheduled_coords:
#     m.plot(xpt_ncp, ypt_ncp, 'o', markersize=7, label='Previously scheduled', mfc='orange', color='0')
# for o in observed_pointings:
#     xpt_ncp, ypt_ncp = m(o.ra.deg, o.dec.deg)
#     m.plot(xpt_ncp, ypt_ncp, 'o', markersize=7, mfc='blue', color='0')
# m.plot(xpt_ncp, ypt_ncp, 'o', markersize=7, label='To be observed', mfc='blue', color='0')
# for f in flux_cal:
#     xcal_ncp, ycal_ncp = m(f.ra.deg, f.dec.deg)
#     m.plot(xcal_ncp, ycal_ncp, 'o', markersize=6, color='green')
# m.plot(xcal_ncp, ycal_ncp, 'o', label='Calibrators', markersize=6, color='green')
# for p in pol_cal:
#     xcal_ncp, ycal_ncp = m(p.ra.deg, p.dec.deg)
#     m.plot(xcal_ncp, ycal_ncp, 'o', markersize=6, color='green')
# plt.legend(loc=1)
# plt.title("Imaging survey fields {}".format(args.starttime_utc))
# plt.savefig(filename,bbox_inches='tight')
