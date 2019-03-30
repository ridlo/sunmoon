#!/usr/bin/python -tt

# Copyright

__author__="Ridlo W. Wibowo"

import math
from util import *

nutation_arg_terms = \
    [
        [0,0,0,0,1],
        [-2,0,0,2,2],
        [0,0,0,2,2],
        [0,0,0,0,2],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [-2,1,0,2,2],
        [0,0,0,2,1],
        [0,0,1,2,2],
        [-2,-1,0,2,2],
        [-2,0,1,0,0],
        [-2,0,0,2,1],
        [0,0,-1,2,2],
        [2,0,0,0,0],
        [0,0,1,0,1],
        [2,0,-1,2,2],
        [0,0,-1,0,1],
        [0,0,1,2,1],
        [-2,0,2,0,0],
        [0,0,-2,2,1],
        [2,0,0,2,2],
        [0,0,2,2,2],
        [0,0,2,0,0],
        [-2,0,1,2,2],
        [0,0,0,2,0],
        [-2,0,0,2,0],
        [0,0,-1,2,1],
        [0,2,0,0,0],
        [2,0,-1,0,1],
        [-2,2,0,2,2],
        [0,1,0,0,1],
        [-2,0,1,0,1],
        [0,-1,0,0,1],
        [0,0,2,-2,0],
        [2,0,-1,2,1],
        [2,0,1,2,2],
        [0,1,0,2,2],
        [-2,1,1,0,0],
        [0,-1,0,2,2],
        [2,0,0,2,1],
        [2,0,1,0,0],
        [-2,0,2,2,2],
        [-2,0,1,2,1],
        [2,0,-2,0,1],
        [2,0,0,0,1],
        [0,-1,1,0,0],
        [-2,-1,0,2,1],
        [-2,0,0,0,1],
        [0,0,2,2,1],
        [-2,0,2,0,1],
        [-2,1,0,2,1],
        [0,0,1,-2,0],
        [-1,0,1,0,0],
        [-2,1,0,0,0],
        [1,0,0,0,0],
        [0,0,1,2,0],
        [0,0,-2,2,2],
        [-1,-1,1,0,0],
        [0,1,1,0,0],
        [0,-1,1,2,2],
        [2,-1,-1,2,2],
        [0,0,3,2,2],
        [2,-1,0,2,2],
    ]


nutation_coeff_terms = \
    [
        [-171996,-174.2,92025,8.9],
        [-13187,-1.6,5736,-3.1],
        [-2274,-0.2,977,-0.5],
        [2062,0.2,-895,0.5],
        [1426,-3.4,54,-0.1],
        [712,0.1,-7,0],
        [-517,1.2,224,-0.6],
        [-386,-0.4,200,0],
        [-301,0,129,-0.1],
        [217,-0.5,-95,0.3],
        [-158,0,0,0],
        [129,0.1,-70,0],
        [123,0,-53,0],
        [63,0,0,0],
        [63,0.1,-33,0],
        [-59,0,26,0],
        [-58,-0.1,32,0],
        [-51,0,27,0],
        [48,0,0,0],
        [46,0,-24,0],
        [-38,0,16,0],
        [-31,0,13,0],
        [29,0,0,0],
        [29,0,-12,0],
        [26,0,0,0],
        [-22,0,0,0],
        [21,0,-10,0],
        [17,-0.1,0,0],
        [16,0,-8,0],
        [-16,0.1,7,0],
        [-15,0,9,0],
        [-13,0,7,0],
        [-12,0,6,0],
        [11,0,0,0],
        [-10,0,5,0],
        [-8,0,3,0],
        [7,0,-3,0],
        [-7,0,0,0],
        [-7,0,3,0],
        [-7,0,3,0],
        [6,0,0,0],
        [6,0,-3,0],
        [6,0,-3,0],
        [-6,0,3,0],
        [-6,0,3,0],
        [5,0,0,0],
        [-5,0,3,0],
        [-5,0,3,0],
        [-5,0,3,0],
        [4,0,0,0],
        [4,0,0,0],
        [4,0,0,0],
        [-4,0,0,0],
        [-4,0,0,0],
        [-4,0,0,0],
        [3,0,0,0],
        [-3,0,0,0],
        [-3,0,0,0],
        [-3,0,0,0],
        [-3,0,0,0],
        [-3,0,0,0],
        [-3,0,0,0],
        [-3,0,0,0],
    ]

def nutation_longitude_and_obliquity(jce):
    # Nutation
    nutation_main_terms = [0, 0, 0, 0, 0]
    nutation_main_terms[0] = poly3(297.85036, 445267.11148, -0.0019142, 1.0/189474.0, jce) # mean_elongation_moon_sun
    nutation_main_terms[1] = poly3(357.52772, 35999.05034, -0.0001603, -1.0/300000.0, jce) # mean_anomaly_sun
    nutation_main_terms[2] = poly3(134.96298, 477198.867398, 0.0086972, 1.0/56250.0, jce) # mean_anomaly_moon
    nutation_main_terms[3] = poly3(93.27191, 483202.017538, -0.0036825, 1.0/327270.0, jce) # argument_lat_moon
    nutation_main_terms[4] = poly3(125.04452, -1934.136261, 0.0020708, 1.0/450000.0, jce) # ascending_lng_moon

    sum_psi = 0.0
    sum_eps = 0.0
    for i in range(63):
        sum_xy = 0.0
        for j in range(5):
            sum_xy += nutation_main_terms[j] * nutation_arg_terms[i][j]
        sum_psi += (nutation_coeff_terms[i][0] + nutation_coeff_terms[i][1]*jce) * math.sin(math.radians(sum_xy))
        sum_eps += (nutation_coeff_terms[i][2] + nutation_coeff_terms[i][3]*jce) * math.cos(math.radians(sum_xy))

    del_psi = sum_psi/36000000.0
    del_eps = sum_eps/36000000.0

    return del_psi, del_eps # in degrees