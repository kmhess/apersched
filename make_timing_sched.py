# make_timing_sched: Make a schedule for Apertif timing surveys (SC4)
# K.M.Hess 19/02/2019 (hess@astro.rug.nl)

__author__ = "Kelley M. Hess"
__date__ = "$20-feb-2019 16:00:00$"
__version__ = "0.3"

import csv
import datetime

from astropy.coordinates import EarthLocation, Longitude, SkyCoord
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
#from modules.functions import *

###################################################################
# Required functions for writing to csv file

def write_to_parset(csvfile,i,source_name,source_pos,start_datetime,end_datetime):
    scan='{}{:03d}'.format(start_datetime.strftime('%Y%m%d'),i)
    #source='{}_{}'.format(source_name.split('_')[0],start_datetime.strftime('%Y%m%d'))
    source='{}'.format(source_name.split('_')[0])
    ra = str(source_pos.to_string('hmsdms').split(' ')[0]).replace('h',':').replace('m',':').replace('s','')
    dec = str(source_pos.to_string('hmsdms').split(' ')[1]).replace('d',':').replace('m',':').replace('s','')
    date1,time1 = start_datetime.strftime('%Y-%m-%d'),start_datetime.strftime('%H:%M:%S')
    date2,time2 = end_datetime.strftime('%Y-%m-%d'),end_datetime.strftime('%H:%M:%S')
    #all_cols=[scan,source,ra,dec,date1,time1,date2,time2,_int,priority,lo, \
    #          sub1,_type,weight,beam,sub2,freq1,freq2,freqcent,intent,person,switch_type]
    all_cols=[source,ra,dec,date1,time1,date2,time2,_int,_type,weight,beam,switch_type]
    writer.writerow(all_cols)

def observe_calibrator(obstimeUTC, obstime=20):
    return obstimeUTC + datetime.timedelta(minutes=obstime)

def observe_target(obstimeUTC, name, obstime=3):
    # Use name to subtract 1 value from the weights column!
    timing_fields['weights'][timing_fields['name'] == name] -= 1
    return obstimeUTC + datetime.timedelta(hours=obstime)

def wait_for_rise(obstimeUTC, waittime=5):
    return obstimeUTC + datetime.timedelta(minutes=waittime)


# Need to add a check to do_calib at the beginning of an observation (I think)
def do_calibration(i, obstimeUTC, telescope_position, csvfile, total_wait):
    currentLST = Time(obstimeUTC).sidereal_time('apparent', westerbork.lon)
    is3c147 = np.abs(currentLST.hour - calibrators[0].ra.hour < 5.0)
    is_ctd93 = np.abs(currentLST.hour - calibrators[2].ra.hour < 5.0)

    calib_wait = 0
    # Wait for calibrator to be at least an hour above the horizon.
    while (not is3c147) and (not is_ctd93):
        calib_wait += dowait
        total_wait += dowait
        new_obstimeUTC = wait_for_rise(obstimeUTC)
        obstimeUTC = new_obstimeUTC
        currentLST = Time(obstimeUTC).sidereal_time('apparent', westerbork.lon)
        is3c147 = np.abs(currentLST.hour - calibrators[0].ra.hour < 5.0)
        is_ctd93 = np.abs(currentLST.hour - calibrators[1].ra.hour < 5.0)
    if calib_wait != 0:
        print("Calibrator not up, waiting {} minutes until LST: {}.".format(calib_wait, str(currentLST)))

    # Observe the calibrator(s) that is (are) up:
    if is3c147:
        slew_seconds = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian],
                                     [calibrators[0].ra.radian, calibrators[0].dec.radian])
        new_obstimeUTC = obstimeUTC + datetime.timedelta(seconds=slew_seconds)
        after_3c147 = observe_calibrator(new_obstimeUTC)
        write_to_parset(csvfile, i, names[0], calibrators[0], new_obstimeUTC, after_3c147)
        print("Scan {} observed {}.".format(i, names[0]))

        i += 1
        slew_seconds = calc_slewtime([calibrators[0].ra.radian, calibrators[0].dec.radian],
                                     [calibrators[1].ra.radian, calibrators[1].dec.radian])
        new_obstimeUTC = after_3c147 + datetime.timedelta(seconds=slew_seconds)
        after_3c138 = observe_calibrator(new_obstimeUTC)
        write_to_parset(csvfile, i, names[1], calibrators[1], new_obstimeUTC, after_3c138)
        print("Scan {} observed {}.".format(i, names[1]))
        return i, after_3c138, calibrators[1], total_wait

    else:
        slew_seconds = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian], \
                                     [calibrators[2].ra.radian, calibrators[2].dec.radian])
        new_obstimeUTC = obstimeUTC + datetime.timedelta(seconds=slew_seconds)
        after_ctd93 = observe_calibrator(new_obstimeUTC)
        write_to_parset(csvfile, i, names[2], calibrators[2], new_obstimeUTC, after_ctd93)
        print("Scan {} observed {}.".format(i, names[2]))
        return i, after_ctd93, calibrators[2], total_wait

def do_target_observation(i, obstimeUTC, telescope_position, csvfile, total_wait, closest_field):
    # Get first position or objects that are close to current observing horizon:
    currentLST = Time(obstimeUTC).sidereal_time('apparent', westerbork.lon)
    proposed_ra = (currentLST + Longitude('1.5h')).wrap_at(360 * u.deg)
    slew_test = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian],  # seconds of time
                              [proposed_ra.radian, telescope_position.dec.radian])

    avail_fields = timing_fields[timing_fields['weights'] > 0]
    if closest_field:
        availability = SkyCoord(timing_fields['hmsdms'][closest_field]).ra.hour - proposed_ra.hour
    else:
        availability = SkyCoord(np.array(avail_fields['hmsdms'])).ra.hour - proposed_ra.hour
    availability[availability < -12] += 24

    targ_wait = 0
    while not np.any((availability < 0.5) & (availability > -0.5)):
        targ_wait += dowait
        total_wait += dowait
        new_obstimeUTC = wait_for_rise(obstimeUTC)
        obstimeUTC = new_obstimeUTC
        currentLST = Time(obstimeUTC).sidereal_time('apparent', westerbork.lon)
        proposed_ra = (currentLST + Longitude('1.5h')).wrap_at(360 * u.deg)
        if closest_field:
            availability = SkyCoord(timing_fields['hmsdms'][closest_field]).ra.hour - proposed_ra.hour
        else:
            availability = SkyCoord(np.array(avail_fields['hmsdms'])).ra.hour - proposed_ra.hour
        availability[availability < -12] += 24
    if targ_wait != 0:
        print
        "Target not up, waiting {} minutes until LST: {}".format(targ_wait, str(currentLST))
    if closest_field:
        first_field = timing_fields[closest_field][0]
    else:
        first_field = random.choice(avail_fields[(availability < 0.5) & (availability > -0.5)])

    # NOTE SLEW TIME IS CALCULATED TO THE *OBSERVING* HORIZON, NOT TO THE NEW RA!
    # TELESCOPE SHOULD MOVE TO HORIZON AND WAIT!!!
    slew_seconds = calc_slewtime([telescope_position.ra.radian, telescope_position.dec.radian], \
                                 [SkyCoord(first_field['hmsdms']).ra.radian,
                                  SkyCoord(first_field['hmsdms']).dec.radian])
    new_obstimeUTC = obstimeUTC + datetime.timedelta(seconds=slew_seconds) + datetime.timedelta(seconds=targ_wait)
    after_target = observe_target(new_obstimeUTC, first_field['name'])
    write_to_parset(csvfile, i, first_field['name'], SkyCoord(first_field['hmsdms']), new_obstimeUTC, after_target)
    print("Scan {} observed {}.".format(i, first_field['hmsdms']))
    return i, after_target, SkyCoord(first_field['hmsdms']), total_wait

###################################################################
westerbork = EarthLocation(lat=52.91460037*u.deg, lon=6.60449982*u.deg, height=82.2786*u.m)
names = ['3C138','3C147','CTD93']
calibrators = [SkyCoord.from_name(name) for name in names]

global dowait,obs_length
dowait = 5         # number of minutes to wait before checking source availability
obs_length = 0.5     # number of days to calculate observations for (can be a float)

global _int,priority,lo,sub1,_type,weight,beam,sub2,freq1,freq2,freqcent,intent,person,switch_type
_int,priority,lo,sub1,_type,weight,beam,sub2,freq1,freq2,freqcent,intent,person,switch_type = '30', \
     'A','4800','64','T','compound','0','320','1250.000','1450.000','1350.000','compound','KH','-'

# Load all-sky pointing file and select the pointings with the label for the appropriate survey:
# Replace filename with the latest if necessary.

try:
    timing_fields
except NameError:
    # labels: l=lofar; m=medium-deep; s=shallow; t=timing; g=Milky Way +/-5 in galactic latitude
    #         h=NCP that will be covered with hexagonal compound beam arrangement
    fields = Table(ascii.read('all_pointings.v4.13dec18.txt', format='fixed_width'))
    timing_fields = fields[(fields['label'] == 't')]
    weights = np.ones(len(timing_fields))  # Each timing field should be observed once.
    # weights[apertif_fields['label']=='m'] = 10     # How to modify other fields.
    timing_fields['weights'] = weights  # Add "weights" column to table.
else:
    print
    "Pointings already loaded. '> del timing_fields' if want to reset weights."

print("Number of all-sky fields are: {}".format(len(fields)))
print("Number of Timing fields are: {}".format(len(timing_fields)))

observations=atdbquery.atdbquery(obs_mode='SC4')
timing_obs=[dict(observations[i])['name'] for i in range(len(observations)) \
             if (dict(observations[i])['name'][0:2]!='3c') and (dict(observations[i])['name'][0:2]!='3C') \
             and (dict(observations[i])['name'][0:2]!='CT')]
# Adjust 'weights' field for objects that have been previously observed:
for obs in timing_obs:
    if obs in fields['name']:
        i=np.where(timing_fields['name']==obs)
        timing_fields['weights'][i]-=1

starttimeUTC = datetime.datetime.utcnow().replace(microsecond=0)  # Keep this as a datetime.time object!!
starttimeUTC = datetime.datetime(2019,1,1,8,0,0)  # Start observing at 9 am Local (8 UTC) on 01 Jan 2019.
obstimeUTC = starttimeUTC

# User supplied field (will choose nearest survey pointing from file later):
#start_pos = SkyCoord.from_name('Abell 262')
start_pos = None

# User supplied filename for the csv file of observed fields:
csv_filename = 'timing_output_file.csv'
header=['source','ra','dec','date1','time1','date2','time2','int','type','weight','beam','switch_type']

# User supplied filename for the map of observed fields:
filename = 'timing_map.png'

# Estimate the telescope starting position as on the meridian (approximately parked)
telescope_position = SkyCoord(ra=Time(obstimeUTC).sidereal_time('apparent',westerbork.lon),dec='50d00m00s')

# Do the observations: select calibrators and target fields, and write the output to a CSV file.
# Also, writes out a record of what is observed, when, and if the telescope has to wait for objects to rise.

# Create a record of the positions observed so we can plot later.
observed_pointings=[]

print("Starting observations! UTC: "+str(obstimeUTC))

if start_pos:
    sep = start_pos.separation(SkyCoord(np.array(timing_fields['hmsdms'])))
    closest_field = np.where((sep == min(sep)) & (timing_fields['weights'] != 0))[0]
else:
    closest_field = None

total_wait = 0
# Open & prepare CSV file to write parset parameters to, in format given by V.M. Moss.
with open(csv_filename,'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)

    # Always start on a calibrator
    i=1
    i,new_obstimeUTC,new_position,total_wait = do_calibration(i,starttimeUTC,telescope_position,writer,total_wait)
    obstimeUTC = new_obstimeUTC
    telescope_position = new_position
    print("UTC: "+str(obstimeUTC)+",  LST: "+ \
           str(Time(obstimeUTC).sidereal_time('apparent',westerbork.lon))+ " at end of scan.")

    # Iterate between target and calibrators for the specified amount of time & write to CSV file:
    while obstimeUTC < starttimeUTC + datetime.timedelta(days = obs_length):
#    while obstimeUTC < starttimeUTC + datetime.timedelta(hours=11):   #For 11 Dec 2018 Shakedown (19 hours)
        i+=1
        i,new_obstimeUTC,new_position,total_wait = do_target_observation(i,obstimeUTC,telescope_position,writer,total_wait,closest_field)
        closest_field = None
        obstimeUTC = new_obstimeUTC
        telescope_position = new_position
        observed_pointings.append(new_position)
        print("UTC: "+str(obstimeUTC)+",  LST: "+ \
               str(Time(obstimeUTC).sidereal_time('apparent',westerbork.lon))+ " at end of scan.")
        ask_calib = ((obstimeUTC-starttimeUTC).total_seconds()/3600.)%24 # Has it been a unit of ~24 hours since last calib?
        if (ask_calib > 22.5) | (ask_calib < 1.5):
            i+=1
            i,new_obstimeUTC,new_position,total_wait = do_calibration(i,obstimeUTC,telescope_position,writer,total_wait)
            obstimeUTC = new_obstimeUTC
            telescope_position = new_position
            print("UTC: "+str(obstimeUTC)+",  LST: "+ \
                   str(Time(obstimeUTC).sidereal_time('apparent',westerbork.lon))+ " at end of scan(s).")
        closest_field = None

print("Ending observations! UTC: "+str(obstimeUTC))
print("Total wait time: {} mins is {:3.1f}% of total.".format(total_wait,total_wait*60./(obstimeUTC-starttimeUTC).total_seconds()*100))

# Create and save a figure of all pointings selected for this survey, and which have been observed.

plt.figure(figsize=[14,10])
m = Basemap(projection='moll',lon_0=90,resolution='l',celestial=True)
m.drawparallels(np.arange(30,90,15),labels=[False,False,False,False],color='darkgray')
m.drawmeridians(np.arange(0,360,15),labels=[True,True,False,True],color='darkgray',latmax=90)
xpt_ncp,ypt_ncp=m(SkyCoord(np.array(timing_fields['hmsdms'])).ra.deg, \
                  SkyCoord(np.array(timing_fields['hmsdms'])).dec.deg)
m.plot(xpt_ncp,ypt_ncp,'o',markersize=6,label='SNS',mfc='none',color='0.3')
for i,f in enumerate(timing_fields):
    if (f['label']=='t') & (f['weights']==0):
        m.plot(xpt_ncp[i],ypt_ncp[i],'o',markersize=7,mfc='red',color='0.8')
m.plot(0,0,'o',markersize=7,label='Already in ATDB',mfc='red',color='0.8')
for o in observed_pointings:
    xpt_ncp,ypt_ncp=m(o.ra.deg,o.dec.deg)
    m.plot(xpt_ncp,ypt_ncp,'o',markersize=6,mfc='blue',color='0.3')
m.plot(xpt_ncp,ypt_ncp,'o',markersize=7,label='To be observed',mfc='blue',color='0.8')
plt.legend(loc=3)

plt.savefig(filename)