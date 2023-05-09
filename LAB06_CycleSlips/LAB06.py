# -*- coding: utf-8 -*-
"""
POSITIONING & LOCATION BASED SERVICES
AA 2022/23

EX 6: Cycle Slips
@author: Marianna Alghisi

-------------------------------------------------------------------------------
Guidelines: DD analysis and cycle slip identification and repairing

Input data: contained in CycleSlipData.txt:
    - col1: Epoch [s]
    - col2: Observed DD [m]
    - col3: Approx DD [m]

Workflow:
    1) Import the data and graph the observed DDs --> plot(epoch, obsDDs)
    2) For each epoch, compute differences between observed DDs and
       approximated DDs (residual DDs) and graph them --> plot(epoch, obsDDs - approxDDs)
    3) Compute the difference between residual DDs between consecutive epochs and plot them
    4) Identify cycle slips and repair them --> see algorithm explained at lesson,
       use 19 cm for wavelenght and 0.20 cycle (3.8 cm) as threshold for cycle slips 
       identification and repairing.
    5) Graph the repaired DDs
-------------------------------------------------------------------------------
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

threshold = 0.038
wavelenght = 0.19

# 1) Import data and graph observed DDs
df = pd.read_csv('LAB06\CycleSlipsData.txt', sep="\t")
#df = pd.read_csv('LAB06\CycleSlipsDataSun.txt', sep="\t")
df.head()

epoch = df.epoch.to_numpy()
DD_obs = df.DD_obs.to_numpy()
DD_aprx = df.DD_apr.to_numpy()

# Plot example:
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(epoch, DD_obs, '-', color='blue')
ax.set(xlabel='Epoch', ylabel='DD_obs', title = 'Observed DDs [m]')

# 2) For each epoch compute differences between observed DDs and approximated DDs (residual DDs) and graph them
res = DD_obs - DD_aprx
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(epoch, res, '-', color='blue')
ax.set(xlabel='Epoch', ylabel='residual DDs', title = 'differences between observed DDs and approximated DDs [m]')

# 3) Compute differences between consecutive epochs of residual DDs and graph them
diff = []
diff = [res[i]-res[i-1] for i in range(1,len(res))]
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(epoch[1:], diff, '-', color='blue')
ax.set(xlabel='Epoch', ylabel='difference', title = 'differences between consecutive epochs of residual DDs [m]')

# 4) Identify cycle slips and repair them (hint: just one for cycle with both the actions)
for i in range(len(diff)):
    if (abs(diff[i])>threshold ):
        print("at i =",i,"diff =", diff[i])
        x = diff[i]/wavelenght
        if(wavelenght*abs(round(x)-x)<=threshold):
            for j in range(i+1,len(res)):
                DD_obs[j] = DD_obs[j] - wavelenght*round(x) 
            # then repeat from line 46

# 5) Graph the corrected DDs, the corrected residuals DDs and their differences in time
# CORRECTED Plot example:
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(epoch, DD_obs, '-', color='blue')
ax.set(xlabel='Epoch', ylabel='DD_obs', title = 'Observed DDs [m]')

# CORRECTED 2) For each epoch compute differences between observed DDs and approximated DDs (residual DDs) and graph them
res = DD_obs - DD_aprx
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(epoch, res, '-', color='blue')
ax.set(xlabel='Epoch', ylabel='residual DDs', title = 'differences between observed DDs and approximated DDs [m]')

# CORRECTED 3) Compute differences between consecutive epochs of residual DDs and graph them
diff = []
diff = [res[i]-res[i-1] for i in range(1,len(res))]
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(epoch[1:], diff, '-', color='blue')
ax.set(xlabel='Epoch', ylabel='difference', title = 'differences between consecutive epochs of residual DDs [m]')

plt.show()