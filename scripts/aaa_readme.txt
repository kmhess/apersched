These scripts originally lived in the main directory, so paths may need to be changed if any of these are run from scripts/.

footprint_check.ipynb - Jupyter notebook to visualize the input fields to make_imaging_sched.py and the previously observed fields.  Uses Basemap for plotting.

footprint_check_new.ipynb - Jupyter notebook to visualize the input fields to make_imaging_sched.py and the previously observed fields.  Uses kapteyn python package for plotting.

make_imaging_hatlas.py - Modified version of make_imaging_sched.py.  Observe only medium-deep pointings in the H-ATLAS part of the sky.  Alternate flux and pol cal each time.

make_imaging_pponly.py - Modified version of make_imaging_sched.py.  Observe only medium-deep pointings in the Perseus-Pisces part of the sky. Visit both flux and pol cal between every observation.

make_imaging_priority.py - Modified version of make_imaging_sched.py.  Prioritize very select fields at the very end of survey operations.

make_new_pointing_file.ipynb - Jupyter notebook to generate a new file of target fields with the desired labels (priority/number of visits).