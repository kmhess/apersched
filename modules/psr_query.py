#!/usr/bin/env python

import time
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from psrqpy import QueryATNF


class PsrQuery(object):
    def __init__(self):

        self.do_query()

    def do_query(self):
        """
        Query ATNF 
        """
        params = ['JNAME', 'RAJ', 'DECJ', 'S1400', 'DM', 'P0']
        condition = "S1400 > 1 && DECJ > -35"
        print("Querying ATNF pulsar catalogue")
        tstart = time.time()
        result = QueryATNF(params=params, condition=condition)
        time_taken = time.time() - tstart
        print("Query took {:.1f} seconds".format(time_taken))
        print("Found {} pulsars for {}".format(len(result.table), condition))
        self.ATNFtable = result.table
        # store coordinates as SkyCoord
        self.ATNFtable['coordinates'] = SkyCoord(self.ATNFtable['RAJ'], self.ATNFtable['DECJ'], unit=(u.hourangle, u.deg))

    def get_pulsars(self, pointing_coord):
        """
        Find which pulsars are within the Apertif FoV for each pointing
        """
        # Furthest CB center is ~1.5deg from the pointing center
        # CB radius is ~.25 deg
        # So 2 deg covers the entire field
        max_dist = 2*u.degree  # covers entire FoV
        # skycoord of this pointing
        # get separation to each pulsar
        sep = pointing_coord.separation(self.ATNFtable['coordinates'])
        in_field = sep < max_dist
        visible_pulsars = self.ATNFtable[in_field]

        self.print_results(visible_pulsars)

    @staticmethod
    def print_results(self, pulsars):
        """
        Print list of pulsars
        """
        print("Number of pulsars within field: {}".format(len(pulsars)))
        if len(pulsars) > 0:
            print("Pulsars:")
            cols = ['JNAME', 'S1400', 'DM', 'P0', 'RAJ', 'DECJ']
            print(pulsars[cols])
        print("\n")
