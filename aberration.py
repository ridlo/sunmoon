#!/usr/bin/python -tt

# Copyright

__author__="Ridlo W. Wibowo"

import math

def aberration_correction(jme, R):
    # Aberration
    # daily variation in arcseconds of geocentric longitude of the Sun
    # in a fixed reference frame
    delta_lamb = 3548.193 + \
        118.568*math.sin(math.radians(87.5287 + 359993.7286*jme)) + \
        2.476*math.sin(math.radians(85.0561 + 719987.4571*jme)) + \
        1.376*math.sin(math.radians(27.8502 + 4452671.1152*jme)) + \
        0.119*math.sin(math.radians(73.1375 + 450368.8564*jme)) + \
        0.114*math.sin(math.radians(337.2264 + 329644.6718*jme)) + \
        0.086*math.sin(math.radians(222.54 + 659289.3436*jme)) + \
        0.078*math.sin(math.radians(162.8136 + 9224659.7915*jme)) + \
        0.054*math.sin(math.radians(82.5823 + 1079981.1857*jme)) + \
        0.052*math.sin(math.radians(171.5189 + 225184.4282*jme)) + \
        0.034*math.sin(math.radians(30.3214 + 4092677.3866*jme)) + \
        0.033*math.sin(math.radians(119.8105 + 337181.4711*jme)) + \
        0.023*math.sin(math.radians(247.5418 + 299295.6151*jme)) + \
        0.023*math.sin(math.radians(325.1526 + 315559.5560*jme)) + \
        0.021*math.sin(math.radians(155.1241 + 675553.2846*jme)) + \
        7.311*jme*math.sin(math.radians(333.4515 + 359993.7286*jme)) + \
        0.305*jme*math.sin(math.radians(330.9814 + 719987.4571*jme)) + \
        0.01*jme*math.sin(math.radians(328.5170 + 1079981.1857*jme)) + \
        0.309*jme*jme*math.sin(math.radians(241.4518 + 359993.7286*jme)) + \
        0.021*jme*jme*math.sin(math.radians(205.0482 + 719987.4571*jme)) + \
        0.004*jme*jme*math.sin(math.radians(297.861 + 4452671.1152*jme)) + \
        0.01*jme*jme*jme*math.sin(math.radians(154.7066 + 359993.7286*jme))

    aberration_corr = -0.005775518*R*delta_lamb # in arcseconds
    return aberration_corr/3600.0 # in degree