#!/usr/bin/python -tt

# Copyright

__author__="Ridlo W. Wibowo"

import math


def limit_degree(deg):
    '''limit degrees to 0 - 360 value'''
    res = deg % 360.0
    if res<0.0: res += 360.0
    return res


def limit_degree180(deg):
    '''limit degree to 0 - 180'''
    res = deg % 180.0
    if res<0.0: res += 180.0
    return res


def limit_degree180pm(deg):
    '''limit degrees to +- 180'''
    res = deg % 360.0
    if res < -180.0: res += 360.0
    elif res > 180.0: res -= 360.0
    return res


def limit_minute(minute):
    res = minute
    if (res < -20.0): res += 1440.0 # 24 hour
    elif (res > 20.0): res -= 1440.0
    return res


def limit_zero_one(val):
    res = val % 1.0
    if res < 1.0: res += 1.0
    return res


def degree_to_minsec(angle): 
    '''Convert degree or hour (number) to [deg, arcmin, arcsec] or [h, m, s] (string) '''
    absAngle = abs(angle)
    degFloor = math.floor(absAngle)
    minute = (absAngle - degFloor)*60.0
    minFloor = math.floor(minute)
    sec = (minute - minFloor)*60.0
    degFloor = str(degFloor)
    minFloor = str(minFloor) 
    sec = str(round(sec,3))
    if (angle < 0.0): degFloor = "-"+degFloor
    # return as string, remove tailing dot and zero for degFloor and minFloor
    return degFloor[:-2], minFloor[:-2], sec 


def minsec_to_degree(deg, min, sec):
    '''Convert deg, min, sec (string) to degree (float)'''
    degree = abs(float(deg)) + abs(float(minute))/60.0 + abs(float(sec))/3600.0
    if (deg[0] == '-'):
        degree = -1.0*degree
    return degree


def degree_to_hourminsec(deg):
    '''Convert deg (float) to h,m,s (string)'''
    hour = deg/15.0
    return degree_to_minsec(hour)


def poly3(a, b, c, d, x):
    '''Calculate 3rd order polynomial'''
    return a + x*(b + x*(c + x*d))


def poly4(a, b, c, d, e, x):
    '''Calculate 4th order polynomial'''
    return a + x*(b + x*(c + x*(d + x*e)))

