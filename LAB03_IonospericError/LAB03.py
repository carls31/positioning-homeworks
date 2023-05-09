# -*- coding: utf-8 -*-
"""
Positioning and Location Based Services
A.A. 2022/2023
3rd Exercise: Ionospheric delay computation

@author: Marianna Alghisi
"""

'''
Goals: 
   1) 4 zenithal maps of ionospheric error corrections:
   - Elevation: 90°
   - Latitude: [-80°, 80°, step = 0.5°]
   - Longitude: [-180°, 180°, step = 0.5°]
   - time: [00:00, 06:00, 12:00, 18:00]
   
   2) 2 polar maps of ionospheric error corrections:
    - Observer located in Milan
    - Elevation: [0, 90°, step = 0.5°]
    - Azimuth: [-180°, 180°, step = 0.5°]
    time: [00:00, 12:00]
'''

# Import required libraries
import numpy as np
import matplotlib.pyplot as plt
from PLBS import iono_correction as ic
import pandas as pd
import geopandas as gpd 
from LAB_variables import alpha, beta # Ionospheric correction parameters

'''
1) ZENITHAL MAPS
'''
Elev = 90
Lat  = np.arange(-80, 80, 0.5)
Lon  = np.arange(-180, 180, 0.5)
time = range(4)
#points = []
# l*m + i
GPStime = 6*3600
iono_delay = np.zeros(shape=(len(time),len(Lon)*len(Lat)))
list_Lat  = np.zeros(len(Lon)*len(Lat))
list_Lon = np.zeros(len(Lon)*len(Lat))
for j in range(0,len(Lon)):
    for k in range(0,len(Lat)):
        # points.append([k, j])
        list_Lon[j + k*len(Lon)] = Lon[j]
        list_Lat[j + k*len(Lon)] = Lat[k]

        for i in time:
            GPStime_i = GPStime*i
            iono_delay[i, j + k*len(Lon) ] = ic(Lat[k], Lon[j], 0, Elev, GPStime_i, alpha, beta)
     

# Loop on time, latitude and longitude --> compute for each point the Ionospheric delay
# TIP: store latitude, longitude and iono_delay in list objects!

# SUGGESTION FOR PLOTTING

timePrint = ['00:00', '06:00', '12:00', '18:00']

for i in range(4):
    IONO = iono_delay[i, :] #list of iono delays for considered epoch
    results = pd.DataFrame()
    
    results['longitude'] = list_Lon
    results['latitude'] = list_Lat 
    results['iono_delay'] = IONO
    print(results.head())
    gdf = gpd.GeoDataFrame(results, geometry=gpd.points_from_xy(results.longitude, results.latitude), crs = 3857)
    world = gpd.read_file('c:/Users/Utente/Documents/Python/Positioning_Lab/LAB02/world/world.shp')
    fig, ax = plt.subplots (figsize = (15,15))
    world.boundary.plot(ax=ax, color='black')
    ax.set(xlabel='Longitude', ylabel='Latitude', title='Zenithal map of ionospheric delay at '+ str(timePrint[i]))
    gdf.plot(column='iono_delay', ax = ax, marker='o', legend=True)

'''
2) POLAR MAPS
'''
# Definition of Milano's coordinates
lat_mi = 45 + 28/60 + 38.28/60**2
lon_mi = 9 + 10/60 + 53.4/60**2

# Inizialization of the parameters for the loop: time, elevation, azimuth
GPStime = 12*3600
elev = np.arange(0, 90, 0.5) #[0, 90°, step = 0.5°]
azim = np.arange(-180, 180, 0.5) 
# Loop on the parameters ---> compute for each point the Ionospheric delay
all_iono_delays=np.zeros(shape=(2,len(elev)*len(azim)))
elevation  = np.zeros(len(elev)*len(azim))
azimuth = np.zeros(len(elev)*len(azim))
for j in range(0,len(elev)):
    for k in range(0,len(azim)):
        elevation[j + k*len(elev)] = elev[j]
        azimuth[j + k*len(elev)] = azim[k]

        for i in [0,1]:
            all_iono_delays[i, j + k*len(elev) ] = ic(lat_mi,lon_mi, azim[k], elev[j], GPStime*i, alpha, beta)

# SUGGESTION FOR PLOTTING

timePrint = ['00:00', '12:00']

for i in [0,1]:
    t = timePrint[i]
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(projection='polar')
    plt.scatter(azimuth, elevation, c=all_iono_delays[i], cmap='brg',  alpha=0.75, label=all_iono_delays[i])
    ax.set_title('Ionospheric Error Polar Map for Milan Observer time = '+str(t))
    plt.colorbar(label='Ionospheric Delay')
plt.show()