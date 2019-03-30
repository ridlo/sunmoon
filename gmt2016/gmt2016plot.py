#!/usr/bin/python -tt

# Copyright

__author__="Ridlo W. Wibowo"

import math
import matplotlib.pyplot as plt

i = 0
with open('gmt2016.txt', 'r') as  f:
    for line in f:
        i += 1
        isi = line.strip().split()
        az_sun = float(isi[0])
        alt_sun = float(isi[1])
        az_moon = float(isi[2])
        alt_moon = float(isi[3])
        disk_sun = float(isi[4])
        disk_moon = float(isi[5])
        time = isi[6] + ' ' + isi[7] + ' ' + isi[8] + ' ' + isi[9] + ' UT'

        fig=plt.figure(1)
        ax=fig.add_subplot(1,1,1, aspect='equal')
        ax.set_xlim([az_sun-1.5, az_sun+1.5])
        ax.set_ylim([alt_sun-1.5, alt_sun+1.5])
        plt.xlabel('Azimuth ($^\circ$)')
        plt.ylabel('Altitude ($^\circ$)')
        ax.text(0.95, 0.95, time, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        sun = plt.Circle((az_sun, alt_sun), disk_sun, color='y')
        moon = plt.Circle((az_moon, alt_moon), disk_moon, color='k')
        fig = plt.gcf()
        fig.gca().add_artist(sun)
        fig.gca().add_artist(moon)
        plt.savefig('plotgmt_'+str(i)+'.png')
        plt.close()
        # fig.savefig('plotgmt_'+str(i)+'.png')



