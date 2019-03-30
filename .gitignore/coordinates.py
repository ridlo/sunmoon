#!/usr/bin/python -tt

# Copyright

__author__="Ridlo W. Wibowo"

import math
from util import *

def true_ecliptic_obliquity(jme, del_eps):
    u = jme/10.0
    mean_obliquity = 84381.448 + u*(-4680.93 + u*(-1.55 + u*(1999.25 + u*(-51.38 + u*(-249.67 + \
        u*(  -39.05 + u*( 7.12 + u*(  27.87 + u*(  5.79 + u*2.45))))))))) # in arcseconds

    return mean_obliquity/3600.0 + del_eps # in degree


def greenwich_sidereal_time(jd, jc, del_psi, eps):
    greenwich_mean_sidereal_time = limit_degree(280.46061837 + 360.98564736629 * (jd - 2451545.0) + \
        jc*jc*(0.000387933 - jc/38710000.0))
    
    return greenwich_mean_sidereal_time + del_psi*math.cos(math.radians(eps))

def ecliptic_to_equator(lamb, beta, eps):
    '''Convert ecliptic to equatorial coordinate'''
    lamb = math.radians(lamb)
    beta = math.radians(beta)
    eps = math.radians(eps)
    alpha = math.atan2(math.sin(lamb) * math.cos(eps) - math.tan(beta) * math.sin(eps), math.cos(lamb))
    delta = math.asin(math.sin(beta) * math.cos(eps) + math.cos(beta) * math.sin(eps) * math.sin(lamb))
    return limit_degree(math.degrees(alpha)), math.degrees(delta)


# def geocentric_to_topocentric_eq():
#     pass


# def atmospheric_refraction_correction(altitude0, pressure, temperature):
#     '''Atmospheric refraction correction
#     - pressure in millibars
#     - temperature in Celcius
#     - only approximation, depend on the wavelength, here is yellow (eyes).
#     '''
    
#     threshold_refrac = -0.83337 # sun_radius + horizon atmos refract = 0.26667 + 0.5667

#     delta_alt = 0.0
#     if (altitude0 >= threshold_refrac): 
#         delta_alt = (pressure / 1010.0) * (283.0 / (273.0 + temperature)) * 1.02 / \
#         (60.0 * math.tan(math.radians(altitude0 + 10.3/(altitude0 + 5.11))))

#     if altitude0 == 90.0: # in case
#         delta_alt = 0.0

#     return delta_alt


# def topocentric_equator_to_altaz(obs_latitude, local_hour_angle_prime, delta_prime, pressure=1010, temperature=10):
#     '''Convert topocentric equatorial coordinate to topocentric horizon coordinate.'''

#     obs_latitude_rad = math.radians(obs_latitude)
#     delta_prime_rad = math.radians(delta_prime)
#     local_hour_angle_prime_rad = math.radians(local_hour_angle_prime)

#     # topocentric elevation angle without atmospheric refraction correction
#     altitude0 = math.degrees(math.asin(math.sin(obs_latitude_rad)*math.sin(delta_prime_rad) + \
#         math.cos(obs_latitude_rad)*math.cos(delta_prime_rad)*math.cos(local_hour_angle_prime)))

#     # atmospheric refraction correction
#     delta_alt = atmospheric_refraction_correction(altitude0, pressure, temperature)

#     # Topocentric elevation angle
#     altitude = altitude0 + delta_alt

#     # Topocentric azimuth angle
#     az0 = limit_degree(math.degrees(math.atan2(math.sin(local_hour_angle_prime), \
#         math.cos(local_hour_angle_prime)*math.sin(obs_latitude_rad) - math.tan(delta_prime_rad)*math.cos(obs_latitude_rad))))

#     azimuth = limit_degree(az0 + 180.0) # E of N

#     return azimuth, altitude