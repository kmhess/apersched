# apersched
Creates an observing schedule for Apertif imaging or timing surveys.

The program reads in a list of possible pointings.  It then calls the module atdbquery to query the ATDB for previously completed observations and removes them from the list of possible pointings.  The program then generates an observing schedule alternating between calibrators and target pointings for a user specified length of time.  The schedule is output as a csv file, and is also accompanied with a .png image which shows which fields have been selected (as well as those which were previously observed).

## Usage:
python3 make_imaging_sched.py

python3 make_timing_sched.py

## Notes:
Need to edit the start time of the observations and the length for which you are scheduling observations by hand within the python script.

For imaging, one may wish to change the input field list to only a list of the "first year" positions.  Current default is all possible imaging positions.

For timing, need to finalize the calibration strategy within the team and modify the script as necessary.
