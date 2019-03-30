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
from sun import *

# Moon's Periodic Terms for Longitude and Distance
MLR_TERMS = \
    [
        [0,0,1,0,6288774,-20905355],
        [2,0,-1,0,1274027,-3699111],
        [2,0,0,0,658314,-2955968],
        [0,0,2,0,213618,-569925],
        [0,1,0,0,-185116,48888],
        [0,0,0,2,-114332,-3149],
        [2,0,-2,0,58793,246158],
        [2,-1,-1,0,57066,-152138],
        [2,0,1,0,53322,-170733],
        [2,-1,0,0,45758,-204586],
        [0,1,-1,0,-40923,-129620],
        [1,0,0,0,-34720,108743],
        [0,1,1,0,-30383,104755],
        [2,0,0,-2,15327,10321],
        [0,0,1,2,-12528,0],
        [0,0,1,-2,10980,79661],
        [4,0,-1,0,10675,-34782],
        [0,0,3,0,10034,-23210],
        [4,0,-2,0,8548,-21636],
        [2,1,-1,0,-7888,24208],
        [2,1,0,0,-6766,30824],
        [1,0,-1,0,-5163,-8379],
        [1,1,0,0,4987,-16675],
        [2,-1,1,0,4036,-12831],
        [2,0,2,0,3994,-10445],
        [4,0,0,0,3861,-11650],
        [2,0,-3,0,3665,14403],
        [0,1,-2,0,-2689,-7003],
        [2,0,-1,2,-2602,0],
        [2,-1,-2,0,2390,10056],
        [1,0,1,0,-2348,6322],
        [2,-2,0,0,2236,-9884],
        [0,1,2,0,-2120,5751],
        [0,2,0,0,-2069,0],
        [2,-2,-1,0,2048,-4950],
        [2,0,1,-2,-1773,4130],
        [2,0,0,2,-1595,0],
        [4,-1,-1,0,1215,-3958],
        [0,0,2,2,-1110,0],
        [3,0,-1,0,-892,3258],
        [2,1,1,0,-810,2616],
        [4,-1,-2,0,759,-1897],
        [0,2,-1,0,-713,-2117],
        [2,2,-1,0,-700,2354],
        [2,1,-2,0,691,0],
        [2,-1,0,-2,596,0],
        [4,0,1,0,549,-1423],
        [0,0,4,0,537,-1117],
        [4,-1,0,0,520,-1571],
        [1,0,-2,0,-487,-1739],
        [2,1,0,-2,-399,0],
        [0,0,2,-2,-381,-4421],
        [1,1,1,0,351,0],
        [3,0,-2,0,-340,0],
        [4,0,-3,0,330,0],
        [2,-1,2,0,327,0],
        [0,2,1,0,-323,1165],
        [1,1,-1,0,299,0],
        [2,0,3,0,294,0],
        [2,0,-1,-2,0,8752]
    ]

# Moon's Periodic Terms for Latitude
MB_TERMS = \
    [
        [0,0,0,1,5128122],
        [0,0,1,1,280602],
        [0,0,1,-1,277693],
        [2,0,0,-1,173237],
        [2,0,-1,1,55413],
        [2,0,-1,-1,46271],
        [2,0,0,1,32573],
        [0,0,2,1,17198],
        [2,0,1,-1,9266],
        [0,0,2,-1,8822],
        [2,-1,0,-1,8216],
        [2,0,-2,-1,4324],
        [2,0,1,1,4200],
        [2,1,0,-1,-3359],
        [2,-1,-1,1,2463],
        [2,-1,0,1,2211],
        [2,-1,-1,-1,2065],
        [0,1,-1,-1,-1870],
        [4,0,-1,-1,1828],
        [0,1,0,1,-1794],
        [0,0,0,3,-1749],
        [0,1,-1,1,-1565],
        [1,0,0,1,-1491],
        [0,1,1,1,-1475],
        [0,1,1,-1,-1410],
        [0,1,0,-1,-1344],
        [1,0,0,-1,-1335],
        [0,0,3,1,1107],
        [4,0,0,-1,1021],
        [4,0,-1,1,833],
        [0,0,1,-3,777],
        [4,0,-2,1,671],
        [2,0,0,-3,607],
        [2,0,2,-1,596],
        [2,-1,1,-1,491],
        [2,0,-2,1,-451],
        [0,0,3,-1,439],
        [2,0,2,1,422],
        [2,0,-3,-1,421],
        [2,1,-1,1,-366],
        [2,1,0,1,-351],
        [4,0,0,1,331],
        [2,-1,1,1,315],
        [2,-2,0,-1,302],
        [0,0,1,3,-283],
        [2,1,1,-1,-229],
        [1,1,0,-1,223],
        [1,1,0,1,223],
        [0,1,-2,-1,-220],
        [2,1,-1,-1,-220],
        [1,0,1,1,-185],
        [2,-1,-2,-1,181],
        [0,1,2,1,-177],
        [4,0,-2,-1,176],
        [4,-1,-1,-1,166],
        [1,0,1,-1,-164],
        [4,0,1,-1,132],
        [1,0,-1,-1,-119],
        [4,-1,0,-1,115],
        [2,-2,0,1,107]
    ]



# Caution! Astronomical Algorithms 1st and 2nd Edition use difference constants
def moon_mean_longitude(jce):
    '''Calculate the Moon's mean longitude, L_prime '''
    return limit_degree(poly4(218.3164477, 481267.88123421, -0.0015786, 1.0/538841, -1.0/65194000, jce))


def moon_mean_elongation(jce):
    '''Calculate the Moon's mean elongation, D '''
    return limit_degree(poly4(297.8501921, 445267.1114034, -0.0018819, 1.0/545868, -1.0/113065000, jce))


def sun_mean_anomaly(jce):
    '''Calculate the Sun's mean anomaly, M '''
    return limit_degree(poly3(357.5291092, 35999.0502909, -0.0001536, 1.0/24490000, jce))


def moon_mean_anomaly(jce):
    '''Calculate the Moon's mean anomaly, M_prime '''
    return limit_degree(poly4(134.9633964, 477198.8675055, 0.0087414, 1.0/69699, -1.0/14712000, jce))


def moon_latitude_argument(jce):
    '''Calculate the Moon's argument of latitude (mean distance of the Moon from its ascending node), F '''
    return limit_degree(poly4(93.2720950, 483202.0175233, -0.0036539, -1.0/3526000, 1.0/863310000, jce))


def moon_periodic_term_summation(d, m, m_prime, f, jce):
    '''Calculate Moon's periodic terms, 
    Jean Meeuss Algorithm, shortened version of ELP2000-82 (Chapront, 1982) Theory 
    with mean arguments from ELP2000-82/B Theory (Chapront, 1998)
    '''

    e = 1.0 - jce*(0.002516 + jce*0.0000074)

    l = 0.0
    r = 0.0
    b = 0.0

    for i in range(60):
        e_pow_LR = e**(abs(MLR_TERMS[i][1]))
        sum_LR = MLR_TERMS[i][0]*d + MLR_TERMS[i][1]*m + MLR_TERMS[i][2]*m_prime + MLR_TERMS[i][3]*f
        l += e_pow_LR*MLR_TERMS[i][4]*math.sin(math.radians(sum_LR))
        r += e_pow_LR*MLR_TERMS[i][5]*math.cos(math.radians(sum_LR))

        e_pow_B = e**(abs(MB_TERMS[i][1]))
        sum_B = MB_TERMS[i][0]*d + MB_TERMS[i][1]*m + MB_TERMS[i][2]*m_prime + MB_TERMS[i][3]*f
        b += e_pow_B*MB_TERMS[i][4]*math.sin(math.radians(sum_B))

    return l, r, b


def moon_geocentric_ecliptic(jce):
    '''Calculate the Moon's geocentric Long, Lat, and distance from the Earth's center, lambda, beta, Delta'''
    A1 = 119.75 + 131.849*jce
    A2 = 53.09 + 479264.290*jce
    A3 = 313.45 + 481266.484*jce

    L_prime = moon_mean_longitude(jce)
    D = moon_mean_elongation(jce)
    M = sun_mean_anomaly(jce)
    M_prime = moon_mean_anomaly(jce)
    F = moon_latitude_argument(jce)
    #print L_prime, D, M, M_prime, F
    
    l, r, b = moon_periodic_term_summation(D, M, M_prime, F, jce)
    #print l, r, b

    delta_l = 3958*math.sin(math.radians(A1)) + 1962*math.sin(math.radians(L_prime - F)) + 318*math.sin(math.radians(A2))
    delta_b = -2235*math.sin(math.radians(L_prime)) + 382*math.sin(math.radians(A3)) \
            + 175*math.sin(math.radians(A1 - F)) + 175*math.sin(math.radians(A1 + F)) \
            + 127*math.sin(math.radians(L_prime - M_prime)) - 115*math.sin(math.radians(L_prime + M_prime))

    lambda_prime = limit_degree(L_prime + 0.000001*(l + delta_l)) # longitude in deg
    beta = limit_degree(0.000001*(b + delta_b)) # latitude in deg
    Delta = 385000.56 + 0.001*r # distance in km

    return lambda_prime, beta, Delta


def moon_eq_hor_parallax(R):
    '''Calculate Moon's equatorial horizontal parallax, in degrees'''
    return math.degrees(math.asin(6378.14/R))


def moon_disk_radius(altitude, Delta):
    '''Calculate topocentric Moon's disk radius, in degrees'''
    return 358473400*(1 + math.sin(math.radians(altitude))*6378.14/Delta)/(3600.0*Delta)


def moon_position(jd, delta_t, obs_longitude, obs_latitude, obs_elevation, pressure=1010, temperature=10):
    '''Calculate Moon's position'''

    jc =julian_ephem_century(jd)
    jde = julian_ephem_day(jd, delta_t)
    jce = julian_ephem_century(jde)
    jme = julian_ephem_millennium(jce)

    # -----------------------------------------------------------------
    ## Geocentric ecliptic position of the the Moon
    lamb, beta, R = moon_geocentric_ecliptic(jce)


    # -----------------------------------------------------------------
    ## Convert geocentric ecliptic to geocentric equatorial (alpha, delta)
    # Apparent Moon's longitude, w/o aberration correction
    del_psi, del_eps = nutation_longitude_and_obliquity(jce)
    lambda_moon = lamb + del_psi 

    # Convert to geocentric equatorial coordinate
    epsilon = true_ecliptic_obliquity(jme, del_eps)
    alpha, delta = ecliptic_to_equator(lambda_moon, beta, epsilon)

    # Calculate local hour angle, H, in degrees
    local_hour_angle = limit_degree(greenwich_sidereal_time(jd, jc, del_psi, epsilon) + obs_longitude - alpha)


    # -----------------------------------------------------------------
    ## Convert geocentric equatorial to topocentric equatorial (alpha_prime, delta_prime, local_hour_angle_prime)
    obs_latitude_rad = math.radians(obs_latitude)
    delta_rad = math.radians(delta)
    local_hour_angle_rad = math.radians(local_hour_angle)

    # Calculate equatorial horizontal parallax
    eq_hor_parallax = moon_eq_hor_parallax(R) # in degrees # the only different
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


    # ------------------------------------------------------------------
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

    moon_radius = moon_disk_radius(altitude, R)

    return azimuth, altitude, moon_radius


if __name__ == '__main__':
    # print "\nExample in Astronomical Algorithms"
    # jce = julian_ephem_century(2448724.5)
    # moon_geocentric_equatorial(jce)

    print "\nExample in NREL SAMPA"
    year = 2009
    month = 7
    day = 22
    UT = 1 + 33/60.0
    delta_t = 66.4
    obs_longitude = 143.36167
    obs_latitude = 24.61167
    obs_elevation = 0.0
    press = 1000 # mbar
    temp = 11 # Celcius
    jd = gregorian_to_jd(year, month, day+UT/24.0)
    print(jd)
    x = moon_position(jd, delta_t, obs_longitude, obs_latitude, obs_elevation, press, temp)
    print(x)

    # print 'Geometric lng, mean equinox of the date: ', lambda_prime, degree_to_minsec(lambda_prime)
    # print 'Apparent longitude: ', lambda_moon, degree_to_minsec(lambda_moon)
    # print 'Apparent latitude: ', beta_moon, degree_to_minsec(beta_moon)
    # print 'Radius Vector: ', Delta
    # print 'Apparent RA: ', alpha, degree_to_hourminsec(alpha)
    # print 'Apparent Dec: ', dec, degree_to_minsec(dec) 