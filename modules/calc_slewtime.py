# calc_slewtime.py: Calculate the WSRT slewing time from one direction to another
#
# Copyright (C) 2018
# ASTRON (Netherlands Institute for Radio Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This file is part of the LOFAR software suite.
# The LOFAR software suite is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The LOFAR software suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
#
# $Id$

__author__ = "Ger van Diepen"
__date__ = "$26-sep-2018 10:47:44$"
__version__ = "0.1"

import math

def calc_slewtime (dir1, dir2, speed_factor=1.):
    """ calculate the WSRT slewing time from direction dir1 to dir2

    :param float[2] dir1: start direction (radians)
    :param float[2] dir1: target direction (radians)
    :param float speed_factor: factor to lower the maximum telescope speed
                               (in case it is very low temperatures)
    :return: slewing time (whole seconds)
    """
    dra = normalize_dra (dir2[0] - dir1[0])
    ddec = dir2[1] - dir1[1]
    # The slewing time is the maximum of the time in RA (HA) and DEC.
    return math.ceil (max (calc_settle_time (dra, speed_factor),
                           calc_settle_time (ddec, speed_factor)))


def calc_settle_time (distance, speed_factor):
    """ calculate time needed to travel the requested distance.
    The implementation is taken from the dish_lib package.
    If it is cold, the maximum speed can be lower which can be controlled
    by giving a speed_factor < 1.
    :param float distance: distance (in ra or dec) in radians
    :param float speed_factor: factor to change the maximum telescope speed
    :return: settle time (i.e., time to move the dish) in seconds
    """
    dist = abs(distance * 180 / math.pi)   # turn into degrees
    MAX_SPEED = 0.3 * speed_factor  # deg/sec
    RAMP_TIME = 5.0  # acceleration/deceleration time zero to max
    acceleration = MAX_SPEED / RAMP_TIME
    ramp_updown_distance = MAX_SPEED * RAMP_TIME
    if dist <= ramp_updown_distance:
        if dist == 0:
            settle_time = 0.
        else:
            # s = 0.5*a*t**2 formula relates distance to acc and time.
            # The distance to move consists of 2 parts (ramp_up and _down).
            # hence dist/2 = 0.5*acc*t**2, thus t = sqrt(dist/acc)
            ramp_time = math.sqrt(dist / acceleration)
            settle_time = ramp_time * 2.0  # for ramp_up and ramp_down
    else:
        settle_time = 2 * RAMP_TIME
        settle_time += (dist - ramp_updown_distance) / MAX_SPEED
    return settle_time

def normalize_dra (dra):
    """
    Normalize the RA difference between -pi and pi.
    """
    if dra > math.pi:
        dra = dra - 2*math.pi
    elif dra < -math.pi:
        dra = dra + 2*math.pi
    return dra

