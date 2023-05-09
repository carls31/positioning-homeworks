# -*- coding: utf-8 -*-
"""
@author: Marianna Alghisi

4th EXERCISE POS&LBS
"""

"""
Compute the coordinates of a GPS receiver using one epoch of pseudorange
observations and the given position of the satellites in view

Meaning of the variables:
    - time_rx       epoch of the observations
    - pr_C1         array of pseudo range observations
    - sat_ids       array with the list of the satellites IDs
    - xyz_sat       estimated positions of the satellites
    - dtS           estimated clock error of the satellites
    - alpha, beta   iono parameters

Useful functions:
    - topocent(X, sat) 
      Returns the elevation, azimuth and the distance of the satellite with respect to the receiver
    - iono_correction(lat, lon, A, E, GPStime, alpha, beta)
      Returns the ionospheric effect in meters by applying Klobuchar Model
    - tropo_correction(h, el)
      Returns the tropospheric effect in meters by applying Saastamoinen Model
"""
import numpy as np

c = 299792458       # speed of light

gps_week = 1906
gps_week_start = 17
year = 2016
month = 7
day = 22
hour = 2
minute = 0
second = 0

time_rx = (day - gps_week_start)*24*3600 + hour*3600 + second*60

pr_C1 = np.array([[22578312.093],
                [20984179.054],
                [22340643.025],
                [21815745.186],
                [23719962.313],
                [24558115.868],
                [25751171.706],
                [24359848.780],
                [26560055.854],
                [20547846.341],
                [25187331.212]])

sat_ids = np.transpose(np.array([28, 5, 13, 7, 20, 9, 8, 2, 21, 30, 15]))

dtS = np.array([[ 0.000538885950029137],
                [-0.000103714172891042],
                [-3.26664571204891*10**(-5)],
                [0.000440397108129438],
                [0.000425625330509237],
                [0.000171981683578018],
                [-4.36651382082638*10**(-5)],
                [0.000573964626877986],
                [-0.000528855944540131],
                [0.000141099019219313],
                [-0.000320324134333714]])

xyz_sat =  np.array([
   [2.266904303720417,   1.376019580610336,   0.242665876089085],
   [1.793483097238855,  -0.684611592059353,   1.836848309959032],
   [1.237349240389060,  -1.488069273671674,   1.810628493493055],
   [0.682633532041234,   1.366381869196001,   2.184884029980142],
   [0.141020153293916,  -1.610493243878792,   2.092127478822572],
   [0.710758426920100,   2.494566375976196,   0.565262210580487],
   [-0.670387964847087,  2.192133222131345,   1.342746581584370],
   [2.183170669948725,  -1.415437238094089,  -0.371984191760939],
   [-1.019768755765267, -1.243833666228832,   2.189467478141541],
   [1.528675973015969,   0.745824912302640,   2.042145267368744],
   [0.467596411393501,  -2.316970109165663,   1.162832980243857]])*10**7

alpha = [0.7451*10**(-8),  0.1490*10**(-7),  -0.5960*10**(-7), -0.1192*10**(-6)]
beta = [0.9216*10**5,  0.1311*10**6, -0.6554*10**5, -0.5243*10**6]