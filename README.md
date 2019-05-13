# apersched
Creates an observing schedule for Apertif imaging or timing surveys.

The programs read in a list of possible pointings.  It then calls the module atdbquery to query the ATDB for previously completed observations and removes them from the list of possible pointings.  The programs then generate an observing schedule alternating between calibrators and target pointings for a user specified length of time.  The schedule is output as a csv file, and is also accompanied with a .png image which shows which fields have been selected (as well as those which were previously observed, or previously scheduled).

## Usage:
```
python3 make_imaging_sched.py [options]
python3 make_timing_sched.py [options]
```

## Notes:
Run with -h to see the list of options.  Recommended first time operation:
```
python3 make_imaging_sched.py -a -b -o test_run
```

If you have no way to connect to ATDB, always run with the -a option.

For imaging, one may wish to change the input field list to only a list of the "first year" positions.  Default is all possible imaging positions, however, a portion of the code bounded by `##` indicates where the user can hardcode a selection.  Current version of the code includes some examples of specific field selections for week 2 of science verification, including focusing on (1) PP for MDS fields or (2) CVn for MDS fields.  These need to be uncommented out to be used.

For timing, need to modify so the choice of the next field is within a narrower range of Dec than it does now (for slewing purposes).

## Behavior:

For imaging, the program selects the lowest Dec field available, however if the Sun is closer to a target field than a user specified amount, the program selects the highest Dec field available.  At least in the beginning of the survey period, in most cases, this will give at least 20 more degrees separation.  If the Sun is within the user specified distance of a calibrator, the program will print a WARNING, but it will still schedule the calibrator. (Need to fix this for 3C286 for imaging!)

A Sun and Moon avoidance strategy has also been implemented for the timing surveys (moon is never above ~+23d).  The minimum Sun distance can be set by the user.  The Moon distance is 0.5 degrees.  (This should be increased & hardcoded to the FoV of the tied array beams).  At the moment, the program prints a warning if a pulsar calibrator is within the minimum Sun or Moon distance, but takes no other action.

## Known 'features':

For imaging, if the RA range of targets is so restricted that it takes more than 6 hours to wait for a target source to rise, then the output png plot includes the position of 3C147 and/or 3C138, however the schedule is still correct.  It's just a funky plotting thing which can be ignored.  This is rare and is only seen when the RA range in one part of the sky is less than a couple hours.
