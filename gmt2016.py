#!/usr/bin/python -tt

# Copyright

__author__="Ridlo W. Wibowo"

import math
import numpy as np
from sun import *
from moon import *
from util import *

year = 2016
month = 3
day = 8
UT = 23.3
delta_t = 68.2
obs_longitude = 113.92273
obs_latitude = -2.21045
obs_elevation = 1000
press = 1010 # mbar
temp = 25 # Celcius

for ut in np.linspace(UT, UT+2.5, 120):
    t = day + ut/24.0
    jd = gregorian_to_jd(year, month, t)
    x = sun_position(jd, delta_t, obs_longitude, obs_latitude, obs_elevation, press, temp)
    y = moon_position(jd, delta_t, obs_longitude, obs_latitude, obs_elevation, press, temp)
    d = math.floor(t)
    utt = degree_to_minsec((t - d)*24.0)
    time = '2016 3 '+ str(d)[:-2] + '  ' + utt[0] + ':' + utt[1] + ':' + utt[2]
    print x[0], x[1], y[0], y[1], x[2], y[2], time