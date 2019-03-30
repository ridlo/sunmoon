#!/usr/bin/python -tt

# Copyright

__author__="Ridlo W. Wibowo"

import math
from subprocess import call, check_output

from util import *
from coordinates import *
from timeconv import *
from nutation import *
from aberration import *


def sun_geocentric_ecliptic(jme, jce):
    '''Get geocentric ecliptic coordinate of the Sun, from Earth's heliocentric coordinate.'''

    res = check_output(["./vsop87d_earth", str(jme)]) # Complete VSOP87 Theory
    res = res.strip().split()
    
    # Heliocentric coordinate of the Earth
    L = float(res[0])
    B = float(res[1])
    R = float(res[2])

    # Geocentric coordinate of the Sun
    phi = (L + 180.0) % 360.0 # longitude of the Sun
    beta = -B # latitude of the Sun

    # Conversion to FK5
    _lamb = math.radians(phi - 1.397*jce - 0.00031*jce*jce)
    delta_phi = -2.5091666666666667E-05 # degree
    delta_beta = 1.0877777777777778E-05*(math.cos(_lamb) - math.sin(_lamb)) # degree
    phi = phi + delta_phi
    beta = beta + delta_beta

    return phi, beta, R


def sun_eq_hor_parallax(R):
    '''Calculate Sun's equatorial horizontal parallax, in degrees'''
    return 8.794/(3600.0*R)


def sun_disk_radius(R):
    '''Calculate topocentric Sun's disk radius, in degrees'''
    return 959.63/(3600.0*R)


def sun_position(jd, delta_t, obs_longitude, obs_latitude, obs_elevation, pressure=1010, temperature=10):
    '''Calculate Sun position
    Modularity problems: del_psi, epsilon, R, local_hour_angle, etc still interconnected

    Parameter:
    - jd = Julian date, in days
    - delta_t = Delta T, in seconds
    - pressure in millibars
    - temperature in Celcius
    '''
    
    jc =julian_ephem_century(jd)
    jde = julian_ephem_day(jd, delta_t)
    jce = julian_ephem_century(jde)
    jme = julian_ephem_millennium(jce)

    # ----------------------------------------------------
    ## Geocentric ecliptic position of the Sun (phi, beta, R)
    phi, beta, R = sun_geocentric_ecliptic(jme, jce)


    # ----------------------------------------------------
    ## Convert geocentric ecliptic to geocentric equatorial (alpha, delta)
    # Apparent Sun (lambda_sun, beta, R): effect of Nutation and Aberration of light
    del_psi, del_eps = nutation_longitude_and_obliquity(jce)
    aberration_corr = aberration_correction(jme, R)
    lambda_sun = phi + del_psi + aberration_corr

    # Convert to geocentric equatorial coordinate
    epsilon = true_ecliptic_obliquity(jme, del_eps) # get true epsilon
    alpha, delta = ecliptic_to_equator(lambda_sun, beta, epsilon)

    # Calculate local hour angle, H, in degrees
    local_hour_angle = limit_degree(greenwich_sidereal_time(jd, jc, del_psi, epsilon) + obs_longitude - alpha)
    


    # ----------------------------------------------------
    ## Convert geocentric equatorial to topocentric equatorial (alpha_prime, delta_prime, local_hour_angle_prime)
    obs_latitude_rad = math.radians(obs_latitude)
    delta_rad = math.radians(delta)
    local_hour_angle_rad = math.radians(local_hour_angle)

    # Calculate equatorial horizontal parallax: is it sufficient? sin(eq_hor_parallax) = sin(8.794")/R
    eq_hor_parallax = sun_eq_hor_parallax(R) # in degrees
    eq_hor_parallax_rad = math.radians(eq_hor_parallax)

    # Earth ellipsoid correctioon
    # Calculate the term u
    u_term = math.atan(0.99664719 * math.tan(obs_latitude_rad)) # 0.99664719 = b/a = 1 - f -> (Earth's flattening)
    # Calculate the term x = rho*cos(phi_prime)
    x_term = math.cos(u_term) + obs_elevation*math.cos(obs_latitude_rad)/6378140.0
    # Calculate the term y = rho*sin(phi_prime)
    y_term = 0.99664719*math.sin(u_term) + obs_elevation*math.sin(obs_latitude_rad)/6378140.0

    # Calculate delta_alpha, parallax in the sun right ascension, Astronomical Algorithms Chap. 39
    delta_alpha_rad = math.atan2(-x_term*math.sin(eq_hor_parallax_rad)*math.sin(local_hour_angle_rad), \
        math.cos(delta_rad) - x_term*math.sin(eq_hor_parallax_rad)*math.sin(local_hour_angle_rad))
    delta_alpha = math.degrees(delta_alpha_rad)

    # Calculate topocentric sun right ascension, declination, and local hour angle
    alpha_prime = alpha + delta_alpha
    delta_prime = math.degrees(math.atan2((math.sin(delta_rad) - y_term*math.sin(eq_hor_parallax_rad))*math.cos(delta_alpha_rad), \
        math.cos(delta_rad) - x_term*math.sin(eq_hor_parallax_rad)*math.cos(local_hour_angle_rad)))
    local_hour_angle_prime = local_hour_angle - delta_alpha


    # ----------------------------------------------------
    ## Convert topocentric equatorial to topocentric altaz (azimuth, altitude)
    delta_prime_rad = math.radians(delta_prime)

    # topocentric elevation angle without atmospheric refraction correction
    altitude0 = math.degrees(math.asin(math.sin(obs_latitude_rad)*math.sin(delta_prime_rad) + \
        math.cos(obs_latitude_rad)*math.cos(delta_prime_rad)*math.cos(math.radians(local_hour_angle_prime))))

    # Atmospheric refraction correction
    # Only approximation, depend on the wavelength, here is yellow (eyes)
    delta_alt = 0.0
    if (altitude0 >= -0.83337): # sun_radius + atmos_refract = 0.26667 + 0.5667
        delta_alt = (pressure / 1010.0) * (283.0 / (273.0 + temperature)) * 1.02 / \
        (60.0 * math.tan(math.radians(altitude0 + 10.3/(altitude0 + 5.11))))

    if altitude0 == 90.0:
        delta_alt = 0.0

    # Topocentric elevation angle
    altitude = altitude0 + delta_alt

    # Topocentric azimuth angle
    az0 = limit_degree(math.degrees(math.atan2(math.sin(math.radians(local_hour_angle_prime)), \
        math.cos(math.radians(local_hour_angle_prime))*math.sin(obs_latitude_rad) - math.tan(delta_prime_rad)*math.cos(obs_latitude_rad))))

    azimuth = limit_degree(az0 + 180.0) # E of N

    sun_radius = sun_disk_radius(R)

    return azimuth, altitude, sun_radius



################################################################################

def sun_mean_longitude(jme):
    '''Calculate Sun's mean longitude.'''
    return limit_degree(280.4664567 + jme*(360007.6982779 + jme*(0.03032028 + \
        jme*(1/49931.0 + jme*(-1/15300.0 + jme*(-1/2000000.0))))))


def equation_of_time(L0, alpha, delta_psi, epsilon):
    '''Calculate the EoT in minutes'''
    return limit_minute(4 * (L0 - 0.0057183 - alpha + delta_psi * math.cos(math.radians(epsilon))))


# approx sun transit time, alpha_zero = apparent right ascension at transit
# m are times, on D, expressed as fractions of a day
def approx_sun_transit_time(alpha_zero, obs_longitude, nu):
    return (alpha_zero - obs_longitude - nu)/360.0


# find hour angle at rise or set (with h0)
# delta_zero = apparent declination, h0_prime = altitude of the body
def sun_hourangle_at_rise_set(obs_latitude, delta_zero, h0_prime):
    h0 = None # return None if circumpolar
    obs_latitude_rad = math.radians(obs_latitude)
    delta_zero_rad = math.radians(delta_zero)
    arg_hourangle = (math.sin(h0_prime) - math.sin(obs_latitude_rad)*math.sin(delta_zero_rad))/(math.cos(obs_latitude_rad)*math.cos(delta_zero_rad))

    if abs(arg_hourangle)<=1.0:
        h0 = limit_degree180(math.degrees(math.acos(arg_hourangle)))

    return h0


# rise = -1, set = +1
# m0 = approx_sun_transit_time
def approx_sun_at_rise(h0, m0):
    return m0 - h0/360.0

def approx_sun_at_set(h0, m0):
    return m0 + h0/360.0


def greenwich_sidereal_time_rts(theta0, m_rts):
    return theta0 + 360.985647*m_rts


def n_term(m_rts, delta_t):
    return m_rts + delta_t/86400.0

# alpha_or_delta = array[D-1, D, D+1]
def interpolate_3point(alpha_or_delta, n_term):
    a = alpha_or_delta[1] - alpha_or_delta[0]
    b = alpha_or_delta[2] - alpha_or_delta[1]
    if (abs(a) >= 2.0): a = limit_zero_one(a)
    if (abs(b) >= 2.0): b = limit_zero_one(b)
    return alpha_or_delta[1] + n_term * (a + b + n_term*(b - a)) / 2.0


# get new local hour angle
def get_local_hour_angle_rts(theta, obs_longitude, alpha):
    return limit_degree180pm(theta + obs_longitude - alpha) # obs_long different +- from Meus

# h_rts
def get_altitude_rts(obs_latitude, delta, local_hour_angle):
    obs_latitude_rad = math.radians(obs_latitude)
    delta_rad = math.radians(delta)
    return math.degrees(math.asin(math.sin(obs_latitude_rad)*math.sin(delta_rad) - \
        math.cos(obs_latitude_rad)*math.cos(delta_rad)*math.cos(math.radians(local_hour_angle))))


def new_transit(m_rts, local_hour_angle):
    return m_rts - local_hour_angle/360.0

def new_rise_set(m_rts, local_hour_angle, h_rts, h0_prime, delta, obs_latitude):
    obs_latitude_rad = math.radians(obs_latitude)
    delta_rad = math.radians(delta)
    return m_rts + (h_rts - h0_prime)/(360.0 * math.cos(delta_rad) \
        * math.cos(obs_latitude_rad) * math.sin(math.radians(local_hour_angle))) 


def sun_rts():
    pass
    # # calculate apparent sidereal time at greenwich
    # greenwich_sidereal_time(jd, jc, del_psi, eps)

    # # calculate geocentric coordinate of the sun: alpha delta
    # alpha, dec = ecliptic_to_equator(lambda_sun, beta, epsilon)

def sun_transit():
    pass


def sun_set():
    pass
    

def pray_times():
    pass


if __name__ == '__main__':
    #print "\nExample in Astronomical Algorithms"
    #jme = -0.00721834360027
    #jce = jme*10
    #sun_geocentric_equatorial(jme, jce)

    print "\nExample in NREL SPA"
    year = 2003
    month = 10
    day = 17
    zonetime = -7.0
    lst = 12+(30+30/60.0)/60.0
    delta_t = 67.0
    obs_longitude = -105.1786
    obs_latitude = 39.742476
    obs_elevation = 1830.14
    press = 820 # mbar
    temp = 11 # Celcius
    UT = lst_to_UT(lst, zonetime)
    jd = gregorian_to_jd(year, month, day+UT/24.0)
    print jd
    x = sun_position(jd, delta_t, obs_longitude, obs_latitude, obs_elevation, press, temp)
    print x

    # print 'Geometric lng, mean equinox of the date: ', phi, degree_to_minsec(phi)
    # print 'Apparent longitude: ', lambda_sun, degree_to_minsec(lambda_sun)
    # print 'Apparent latitude: ', beta, degree_to_minsec(beta)
    # print 'Radius Vector: ', R
    # print 'Apparent RA: ', alpha, degree_to_hourminsec(alpha)
    # print 'Apparent Dec: ', dec, degree_to_minsec(dec) 
