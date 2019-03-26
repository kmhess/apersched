# apersched
Creates an observing schedule for Apertif imaging or timing surveys.

The programs read in a list of possible pointings.  It then calls the module atdbquery to query the ATDB for previously completed observations and removes them from the list of possible pointings.  The programs then generate an observing schedule alternating between calibrators and target pointings for a user specified length of time.  The schedule is output as a csv file, and is also accompanied with a .png image which shows which fields have been selected (as well as those which were previously observed).

## Usage:
<<<<<<< HEAD
python3 make_imaging_sched.py

python3 make_timing_sched.py
=======
```
python3 make_imaging_sched.py [options]
python3 make_timing_sched.py [options]
```
>>>>>>> 2f582820c5091d2d184067bf16199fe62be19981

## Notes:
Run with -h to see the list of options.  Recommended first time operation:
```
python3 make_imaging_sched.py -a -b -o test_run
```

If you have no way to connect to ATDB, always run with the -a option.

For imaging, one may wish to change the input field list to only a list of the "first year" positions.  Current default is all possible imaging positions.

For timing, need to finalize the calibration strategy within the team and modify the script as necessary.
