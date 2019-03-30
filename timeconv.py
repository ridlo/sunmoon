#!/usr/bin/python -tt

# Copyright

__author__="Ridlo W. Wibowo"

import math
from util import *


def gregorian_to_jd(year, month, day):
    '''Convert Gregorian Calendar (or Julian Calendar before 15 Oct 1582) to Julian day
    Reference: Astronomical Algorithms by Jean Meeus
    JD start from 1 Jan -4712 (4713 B.C.) at 12:00:00 UT
    this method valid for positive and negative years, but not for negative JD.
    '''
    if year < -4712:
        print "Value Error: this convertion algorithm only works for positive JD."
        return None

    if month <= 2:
        year -= 1
        month += 12

    A = math.floor(year / 100)

    # flip to julian calendar, B = 0
    # Gregorian start from 15 Oct 1582
    if (year < 1582):
        B = 0
    elif (year == 1582):
        if month < 10:
            B = 0
        elif month == 10:
            if day < 15:
                B = 0
            else:
                B = 2 - A + math.floor(A / 4)
        else:
            B = 2 - A + math.floor(A / 4)
    else:
        B = 2 - A + math.floor(A / 4)

    return math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (month + 1)) + day + B - 1524.5


def gregorian_to_jd_zero(year):
    '''Return JD0: January 0.0 of a given year
    Reference: Astronomical Algorithms by Jean Meeus
    ''' 
    Y = year-1
    A = math.floor(Y/100)
    return math.floor(365.25*Y) - A + math.floor(A/4) + 1721424.5


def get_mjd(year, month, day):
    '''Return Modified Julian Day: start from 1858 Nov 17, 00:00 UT
    Reference: Astronomical Algorithms by Jean Meeus
    '''
    return gregorian_to_jd(year, month, day) - 2400000.5


def julian_ephem_day(jd, delta_t):
    '''Return Julian Ephemeris Day, JDE'''
    return jd + delta_t/86400.0


def julian_century(jd):
    '''Return Julian Century, from J2000, JC'''
    return (jd - 2451545.0)/36525.0


def julian_ephem_century(jde):
    '''Return Julian Ephemeris Century, JCE'''
    return (jde - 2451545.0)/36525.0


def julian_ephem_millennium(jce):
    '''Return Julian Ephemeris Millennium, JME'''
    return jce/10.0


def jd_to_gregorian(jd):
    '''Convert JD to Gregorian calendar
    Reference: Astronomical Algorithms by Jean Meeus
    this method valid for positive and negative years, but not for negative JD.
    '''
    if jd < 0:
        print "Value Error: this convertion algorithm only works for positive JD."
        return None

    jd += 0.5
    Z = math.floor(jd)
    F = jd - Z

    if Z < 2299161:
        A = Z
    else:
        alpha = math.floor((Z - 1867216.25)/36524.25)
        A = Z + 1 + alpha - math.floor(alpha/4)

    B = A + 1524
    C = math.floor((B - 122.1)/365.25)
    D = math.floor(365.25 * C)
    E = math.floor((B - D)/30.6001)

    day = B - D - math.floor(30.6001 * E) + F
    
    if E < 14:
        month = E - 1
    else:
        month = E - 13

    if month > 2:
        year = C - 4716
    else:
        year = C - 4715

    return year, month, day


def interval_day(day1, day2):
    '''Calculate time interval between two gregorian date
    Reference: Astronomical Algorithms by Jean Meeus
    '''
    return gregorian_to_jd(day2[0], day2[1], day2[2]) - gregorian_to_jd(day1[0], day1[1], day1[2])


def get_day_from_gregorian(year, month, day):
    '''Return day of the week
    0 = Sunday, 6 = Saturday
    Reference: Astronomical Algorithms by Jean Meeus
    '''
    day = math.floor(day) # at 00:00 UT
    jd = gregorian_to_jd(year, month, day)
    jd += 1.5
    return math.fmod(jd, 7) # C modulo


def is_leapyear(year):
    '''Return Boolean, leap year (True) or not (False)
    Reference: Astronomical Algorithms by Jean Meeus
    '''
    if year < 1583: # julian calendar
        return True if (year % 4) == 0 else False
    else: # gregorian calendar
        return ((year % 4) == 0 and (year % 100) != 0) or ((year % 400) == 0)


def day_of_the_year(year, month, day):
    '''Return number of day from 1 January in that year
    e.g. 1 Jan  = 1 and 31 Dec = 365 (common year) or 366 (leap year)
    Reference: Astronomical Algorithms by Jean Meeus
    '''
    K = 1 if is_leapyear(year) else 2
    return math.floor(275*month/9) - K * math.floor((month+9)/12) + math.floor(day) - 30


def hijri_to_gregorian(H, M, D): # or Julian
    '''Convert Hijri Moslem Calendar (Lunar) to Gregorian Calendar (or Julian if year < 1583)
    meaningless for dates earlier than 622 July 16 of the Julian calendar, (1 Muharram A.H. 1)
    Reference: Astronomical Algorithms by Jean Meeus
    '''
    N = D + math.floor(29.5001*(M - 1) + 0.99)
    Q = math.floor(H/30.0)
    R = H % 30
    A = math.floor((11*R + 3)/30.0)
    W = 404*Q + 354*R + 208 + A
    Q1 = math.floor(W/1461.0)
    Q2 = math.floor(W % 1461)
    G = 621 + 4 * math.floor(7*Q + Q1)
    K = math.floor(Q2/365.2422)
    E = math.floor(365.2422*K)
    J = Q2 - E + N - 1
    X = G + K
    
    # J is the number of the day in the Julian year X
    if (J > 366) and ((X % 4) == 0):
        J -= 366
        X += 1
    
    if (J > 365) and ((X % 4) > 0):
        J -= 365
        X += 1

    # if the date is later than 1582 Oct 4
    jd = math.floor(365.25*(X - 1)) + 1721423 + J
    alpha = math.floor((jd - 1867216.25)/36524.25)
    beta = jd + 1 + alpha - math.floor(alpha/4.0)

    if jd < 2299161:
        beta = jd

    b = beta + 1524
    c = math.floor((b - 122.1)/365.25)
    d = math.floor(365.25*c)
    e = math.floor((b - d)/30.6001)

    day = b - d - math.floor(30.6001*e)
    month = (e - 1) if (e < 14) else (e - 13)
    year = (c - 4716) if (month > 2) else (c - 4715)

    leap = False
    if ((11*R + 3) % 30) > 18: leap = True

    return year, month, day, leap


def gregorian_to_hijri(X, M, D):
    '''Convert Gregorian to Hijri Moslem Calendar (Lunar)
    meaningless for dates earlier than 622 July 16 of the Julian calendar, (1 Muharram A.H. 1)
    Reference: Astronomical Algorithms by Jean Meeus 2nd edition
    '''
    if M < 3:
        X -= 1
        M += 12

    alpha = math.floor(X/100)
    beta = 2 - alpha + math.floor(alpha/4.0)
    b = math.floor(365.25*X) + math.floor(30.6001*(M + 1)) + D + 1722519 + beta
    c = math.floor((b - 122.1)/365.25)
    d = math.floor(365.25*c)
    e = math.floor((b - d)/30.6001)

    # in julian calendar
    D = b - d - math.floor(30.6001*e)
    M = (e - 1) if (e < 14) else (e - 13)
    X = (c - 4716) if (M > 2) else (c - 4715)

    W = 1 if ((X % 4) == 0) else 2
    N = math.floor(275*M/9.0) - W*math.floor((M + 9)/12.0) + D - 30
    A = X - 623
    B = math.floor(A/4.0)
    C = A % 4
    C1 = 365.25001*C
    C2 = math.floor(C1)
    if (C1 - C2) > 0.5: C2 += 1
    DD = 1461*B + 170 + C2
    Q = math.floor(DD/10631.0)
    R = DD % 10631
    J = math.floor(R/354.0)
    K = R % 354
    O = math.floor((11*J + 14)/30.0)
    H = 30*Q + J + 1
    JJ = K - O + N - 1 # the number of the day in moslem year H
    CL = H % 30
    DL = (11*CL + 3) % 30
    if DL < 19:
        JJ -= 354
        H += 1
    else:
        JJ -= 355
        H += 1

    if (JJ == 0): 
        JJ = 355
        H -= 1

    S = math.floor((JJ - 1)/29.5)
    m = 1 + S
    d = math.floor(JJ - 29.5*S)
    if JJ == 355:
        m = 12
        d = 30

    return H, m, d


def get_gregorian_time(year, month, day):
    '''Convert year, month, day (in fraction) to human language'''
    D = get_day_from_gregorian(year, month, day)

    if D == 0:
        d = 'Sunday'
    elif D == 1:
        d = 'Monday'
    elif D == 2:
        d = 'Tuesday'
    elif D == 3:
        d = 'Wednesday'
    elif D == 4:
        d = 'Thursday'
    elif D == 5:
        d = 'Friday'
    elif D == 6:
        d = 'Saturday'
    else:
        d = None

    if month == 1:
        m = 'Januari'
    elif month == 2:
        m = 'February'
    elif month == 3:
        m = 'March'
    elif month == 4:
        m = 'April'
    elif month == 5:
        m = 'May'
    elif month == 6:
        m = 'June'
    elif month == 7:
        m = 'July'
    elif month == 8:
        m = 'August'
    elif month == 9:
        m = 'September'
    elif month == 10:
        m = 'October'
    elif month == 11:
        m = 'November'
    elif month == 12:
        m = 'December'
    else:
        m = None

    date = math.floor(day)
    t = day - date
    h = t * 24.0
    hour = math.floor(h) 
    mi = (h - hour)*60.0
    minute = math.floor(mi)
    second = (mi - minute)*60.0

    # e.g. [Sunday, 13, December, 2015, 16, 49, 53]
    return d, date, m, year, hour, minute, second


def get_hijri_time(year, month, day):
    pass


def get_delta_t(year, month, day):
    '''Get Delta T in seconds
    http://maia.usno.navy.mil/
    '''
    pass


def get_UT(year, month, day, zonetime):
    dayUT = day - zonetime/24.0
    return year, month, dayUT


def lst_to_UT(lst, zonetime):
    '''local standard time in hour, zonetime in hour'''
    return lst - zonetime


def dayfrac_to_zonetime(dayfraction, timezone):
    return 24.0*limit_zero_one(dayfraction + timezone/24.0)


if __name__ == '__main__':
    # Testing
    # test convert Gregorian to JD
    print gregorian_to_jd(2000, 1, 1.5)
    print gregorian_to_jd(1999, 1, 1.0)
    print gregorian_to_jd(1987, 1, 27.0)
    print gregorian_to_jd(1987, 6, 19.5)
    print gregorian_to_jd(1988, 1, 27.0)
    print gregorian_to_jd(1988, 6, 19.5)
    print gregorian_to_jd(1900, 1, 1.0)
    print gregorian_to_jd(1600, 1, 1.0)
    print gregorian_to_jd(1600, 12, 31.0)
    print gregorian_to_jd(837, 4, 10.3) 
    print gregorian_to_jd(-123, 12, 31.0)
    print gregorian_to_jd(-122, 1, 1.0)
    print gregorian_to_jd(-1000, 7, 12.5)
    print gregorian_to_jd(-1000, 2, 29.0)
    print gregorian_to_jd(-1001, 8, 17.9)
    print gregorian_to_jd(-4712, 1, 1.5)

    # test JD to Gregorian
    print jd_to_gregorian(2436116.31)
    print jd_to_gregorian(1842713.0)
    print jd_to_gregorian(1507900.13)

    # test interval day
    print interval_day([1910, 4, 20.0], [1986, 2, 9.0])

    # test day of the week
    print get_day_from_gregorian(1954, 6, 30.0)

    # test day of the year
    print day_of_the_year(1978, 11, 14)
    print day_of_the_year(1988, 4, 22)

    # test hijri to gregorian
    print hijri_to_gregorian(1421, 1, 1)

    # test gregorian to hijri
    print gregorian_to_hijri(1991, 8, 13)

    #test get gregorian time
    print get_gregorian_time(1990, 7, 25.5)