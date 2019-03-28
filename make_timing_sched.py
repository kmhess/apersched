# make_timing_sched: Make a schedule for Apertif timing surveys (SC4)
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)

__author__ = "Kelley M. Hess"
__date__ = "$21-feb-2019 16:00:00$"
__version__ = "0.4"

import csv
import datetime

from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.coordinates import Longitude, SkyCoord
from astropy.io import ascii
from astropy.table import Table
from astropy.time import Time
from astropy import units as u
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import random

from atdbquery import atdbquery
from modules.calc_slewtime import calc_slewtime  # Wants [ra,dec] start/end positions in radians; outputs seconds.
from modules.calibrators import *
from modules.functions import *
from modules.telescope_params import westerbork

###################################################################
# Survey specific functions for doing observations and calibration

# Need to add a check to do_calib at the beginning of an observation (I think)
def do_calibration(i, obstime_utc, telescope_position, csvfile, total_wait):
    current_lst = Time(obstime_utc).sidereal_time('apparent', westerbork().lon)
    isB1933 = np.abs(current_lst.hour - psr_cal[0].ra.hour) < 5.0
    isB0531 = np.abs(current_lst.hour - psr_cal[1].ra.hour) < 5.0
    isB0329 = np.abs(current_lst.hour - psr_cal[2].ra.hour) < 5.0
    isB0950 = np.abs(current_lst.hour - psr_cal[3].ra.hour) < 5.0

    calib_wait = 180  # to make sure ATDB has enough time to start a new observation
    # Wait for calibrator to be at least an hour above the horizon.
    while (not isB1933) and (not isB0531) and (not isB0329) and (not isB0950):
        calib_wait += dowait
        total_wait += dowait
        new_obstime_utc = wait_for_rise(obstime_utc, waittime=dowait)
        obstime_utc = new_obstime_utc
        current_lst = Time(obstime_utc).sidereal_time('apparent', westerbork().lon)
        isB1933 = np.abs(current_lst.hour - psr_cal[0].ra.hour) < 5.0
        isB0531 = np.abs(current_lst.hour - psr_cal[1].ra.hour) < 5.0
        isB0329 = np.abs(current_lst.hour - psr_cal[2].ra.hour) < 5.0
        isB0950 = np.abs(current_lst.hour - psr_cal[3].ra.hour) < 5.0
    if calib_wait != 180:
        print("Calibrator not up, waiting {} minutes until LST: {}.".format(calib_wait, str(current_lst)))

    # Observe the calibrator(s) that is (are) up:
    if isB1933:
        name = psr_names[0]
        psr = psr_cal[0]
    elif isB0531:
        name = psr_names[1]
        psr = psr_cal[1]
    elif isB0329:
        name = psr_names[2]
        psr = psr_cal[2]
    elif isB0950:
        name = psr_names[3]
        psr = psr_cal[3]
    else:
        print("Error: No known test pulsar visible. This should not be possible")
        exit()

    slew_seconds = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian],
                                 [psr.ra.radian, psr.dec.radian])
    new_obstime_utc = obstime_utc + datetime.timedelta(seconds=slew_seconds)
    after_psr = observe_calibrator(new_obstime_utc, obstime=5)
    write_to_csv(csvfile, name, psr, new_obstime_utc, after_psr)
    print("Scan {} observed {}.".format(i, name))
    return i, after_psr, psr, total_wait

def do_target_observation(i, obstime_utc, telescope_position, csvfile, total_wait, closest_field):
    # Get first position or objects that are close to current observing horizon:
    current_lst = Time(obstime_utc).sidereal_time('apparent', westerbork().lon)
    proposed_ra = (current_lst + Longitude('1.5h')).wrap_at(360 * u.deg)

    avail_fields = timing_fields[timing_fields['weights'] > 0]
    if closest_field:
        availability = SkyCoord(timing_fields['hmsdms'][closest_field]).ra.hour - proposed_ra.hour
    else:
        availability = SkyCoord(np.array(avail_fields['hmsdms'])).ra.hour - proposed_ra.hour
    availability[availability < -12] += 24

    targ_wait = 180 # to make sure ATDB has enough time to start a new observation
    while not np.any((availability < 0.5) & (availability > -0.5)):
        targ_wait += dowait
        total_wait += dowait
        new_obstimeUTC = wait_for_rise(obstime_utc, waittime=dowait)
        obstime_utc = new_obstimeUTC
        current_lst = Time(obstime_utc).sidereal_time('apparent', westerbork().lon)
        proposed_ra = (current_lst + Longitude('1.5h')).wrap_at(360 * u.deg)
        if closest_field:
            availability = SkyCoord(timing_fields['hmsdms'][closest_field]).ra.hour - proposed_ra.hour
        else:
            availability = SkyCoord(np.array(avail_fields['hmsdms'])).ra.hour - proposed_ra.hour
        availability[availability < -12] += 24
    if targ_wait != 180:
        print("Target not up, waiting {} minutes until LST: {}".format(targ_wait, str(current_lst)))
    if closest_field:
        first_field = timing_fields[closest_field][0]
    else:
        first_field = random.choice(avail_fields[(availability < 0.5) & (availability > -0.5)])

    # NOTE SLEW TIME IS CALCULATED TO THE *OBSERVING* HORIZON, NOT TO THE NEW RA!
    # TELESCOPE SHOULD MOVE TO HORIZON AND WAIT!!!
    slew_seconds = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian],
                                 [SkyCoord(first_field['hmsdms']).ra.radian,
                                  SkyCoord(first_field['hmsdms']).dec.radian])
    new_obstimeUTC = obstime_utc + datetime.timedelta(seconds=slew_seconds) + datetime.timedelta(seconds=targ_wait)
    after_target = observe_target(timing_fields, new_obstimeUTC, first_field['name'], obstime=3)
    write_to_csv(csvfile, first_field['name'], SkyCoord(first_field['hmsdms']), new_obstimeUTC, after_target)
    print("Scan {} observed {}.".format(i, first_field['hmsdms']))
    return i, after_target, SkyCoord(first_field['hmsdms']), total_wait

###################################################################

parser = ArgumentParser(description="Make observing schedule for the Apertif timing SC4. Saves schedule in CSV file to be parsed by atdbspec. "
                                    "Outputs a png of the completed and scheduled pointings.",
                        formatter_class=RawTextHelpFormatter)

parser.add_argument('-f', '--filename', default='./ancillary_data/all_pointings.v4.13dec18.txt',
                    help='Specify the input file of pointings to choose from (default: %(default)s).')
parser.add_argument('-o', '--output', default='temp',
                    help='Specify the suffix of output csv and png files (default: timing_sched_%(default)s).csv.')
parser.add_argument('-s', "--starttime_utc", default="2019-03-11 08:00:00",
                    help="The start time in ** UTC ** ! - format 'YYYY-MM-DD HH:MM:SS' (default: '%(default)s').",
                    type=datetime.datetime.fromisoformat)
parser.add_argument('-l', "--schedule_length", default=2.0,
                    help="Number of days to schedule (can be float; default: %(default)s).",
                    type=float)
parser.add_argument('-a', "--check_atdb",
                    help="If option is included, *DO NOT* check ATDB for previous observations.",
                    action='store_false')

# Parse the arguments above
args = parser.parse_args()

# Filename for the csv file of observed fields:
csv_filename = 'timing_sched_{}.csv'.format(args.output)

# Filename for the map of observed fields:
filename = 'timing_map_{}.png'.format(args.output)

# Number of minutes to wait before checking source availability
dowait = 2

# Load all-sky pointing file and select the pointings with the label for the appropriate survey:

# labels: l=lofar; m=medium-deep; s=shallow; t=timing; g=Milky Way +/-5 in galactic latitude
#         h=NCP that will be covered with hexagonal compound beam arrangement
fields = Table(ascii.read(args.filename, format='fixed_width'))
timing_fields = fields[(fields['label'] == 't')]
weights = np.ones(len(timing_fields))  # Each timing field should be observed once.
# weights[apertif_fields['label']=='m'] = 10     # How to modify other fields.
timing_fields['weights'] = weights  # Add "weights" column to table.

print("\n##################################################################")
print("Number of all-sky fields are: {}".format(len(fields)))
print("Number of Timing fields are: {}".format(len(timing_fields)))

# Retrieve names of observations from ATDB (excludes calibrator scans)
if args.check_atdb:
    observations = atdbquery.atdbquery(obs_mode='SC4')
    timing_obs = [dict(observations[i])['name'] for i in range(len(observations))
                  if (dict(observations[i])['name'][0:2] != '3c') and (dict(observations[i])['name'][0:2] != '3C')
                  and (dict(observations[i])['name'][0:2] != 'CT')]
    # Adjust 'weights' field for objects that have been previously observed:
    for obs in timing_obs:
        if obs in fields['name']:
            i = np.where(timing_fields['name'] == obs)
            timing_fields['weights'][i] -= 1
else:
    print("Not querying ATDB for previous observations.")

# Estimate the telescope starting position as on the meridian (approximately parked)
telescope_position = SkyCoord(ra=Time(args.starttime_utc).sidereal_time('apparent', westerbork().lon), dec='50d00m00s')

# Create a record of the positions planned to be observed so we can plot later.
observed_pointings = []

print("\nStarting observations! UTC: " + str(args.starttime_utc))
print("Calculating schedule for the following " + str(args.schedule_length) + " days.\n")

# Keep track of observing efficiency
total_wait = 0
closest_field = None

# Open & prepare CSV file to write parset parameters to, in format given by V.M. Moss.
# (This could probably be done better because write_to_parset is in modules/function.py)
header = ['source', 'ra', 'dec', 'date1', 'time1', 'date2', 'time2', 'int', 'type', 'weight', 'beam', 'switch_type']

# Do the observations: select calibrators and target fields, and write the output to a CSV file.
# Also, writes out a record of what is observed, when, and if the telescope has to wait for objects to rise.
with open(csv_filename, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)

    # Always start on a calibrator
    i = 1
    i, new_obstime_utc, new_position, total_wait = do_calibration(i, args.starttime_utc, telescope_position, writer, total_wait)
    obstime_utc = new_obstime_utc
    telescope_position = new_position
    print("\tUTC: " + str(obstime_utc) + ",  LST: " + str(Time(obstime_utc).sidereal_time('apparent', westerbork().lon)) + " at end of scan.")

    # Iterate between target and calibrators for the specified amount of time & write to CSV file:
    while obstime_utc < args.starttime_utc + datetime.timedelta(days=args.schedule_length):
        i += 1
        i, new_obstime_utc, new_position, total_wait = do_target_observation(i, obstime_utc, telescope_position, writer, total_wait, closest_field)
        closest_field = None
        obstime_utc = new_obstime_utc
        telescope_position = new_position
        observed_pointings.append(new_position)
        print("\tUTC: " + str(obstime_utc) + ",  LST: " + str(Time(obstime_utc).sidereal_time('apparent', westerbork().lon)) + " at end of scan.")
        ask_calib = ((obstime_utc - args.starttime_utc).total_seconds() / 3600.) % 12   # Has it been a unit of ~12 hours since last calib?
        if (ask_calib > 10.5) | (ask_calib < 1.5):
            i += 1
            i, new_obstime_utc, new_position, total_wait = do_calibration(i, obstime_utc, telescope_position, writer, total_wait)
            obstime_utc = new_obstime_utc
            telescope_position = new_position
            print("\tUTC: " + str(obstime_utc) + ",  LST: " + str(Time(obstime_utc).sidereal_time('apparent', westerbork().lon)) + " at end of scan(s).")
        closest_field = None

print("Ending observations! UTC: " + str(obstime_utc))
print("Total wait time: {} mins is {:3.1f}% of total.".format(total_wait, total_wait * 60. / (obstime_utc - args.starttime_utc).total_seconds() * 100))
print("\nThe schedule has been written to " + csv_filename)
print("A map of the observed fields has been written to " + filename)
print("##################################################################\n")

# Create and save a figure of all pointings selected for this survey, and which have been observed.
plt.figure(figsize=[14, 10])
m = Basemap(projection='moll', lon_0=90, resolution='l', celestial=True)
m.drawparallels(np.arange(30, 90, 15), labels=[False, False, False, False], color='darkgray')
xpt_ncp, ypt_ncp = m(SkyCoord(np.array(timing_fields['hmsdms'])).ra.deg,
                     SkyCoord(np.array(timing_fields['hmsdms'])).dec.deg)
m.plot(xpt_ncp, ypt_ncp, 'o', markersize=6, label='SNS', mfc='none', color='0.3')
for i, f in enumerate(timing_fields):
    if (f['label'] == 't') & (f['weights'] == 0):
        m.plot(xpt_ncp[i], ypt_ncp[i], 'o', markersize=7, mfc='red', color='0.8')
m.plot(0, 0, 'o', markersize=7, label='Already in ATDB', mfc='red', color='0.8')
for o in observed_pointings:
    xpt_ncp, ypt_ncp = m(o.ra.deg, o.dec.deg)
    m.plot(xpt_ncp, ypt_ncp, 'o', markersize=6, mfc='blue', color='0.3')
m.plot(xpt_ncp, ypt_ncp, 'o', markersize=7, label='To be observed', mfc='blue', color='0.8')
plt.legend(loc=3)

plt.savefig(filename)
