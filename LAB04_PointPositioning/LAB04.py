# -*- coding: utf-8 -*-
"""
POSITIONING & LOCATION BASED SERVICES
AA 2021/2022

EX 4: Point Positioning

@author: Marianna Alghisi
"""


import numpy as np
from LAB_variables import *
import matplotlib.pyplot as plt
from PLBS import topocent, GC2Geodetic, iono_correction, tropo_correction, R_0, rad2deg

import sys 
orig_stdout = sys.stdout
f = open('C:/Users/Utente/Documents/Python/Positioning_Lab/LAB04_results.txt','w')
sys.stdout = f

# Initialization of the parameters
max_iter = 10
convergence = 0.1
X_start = np.array([[0], [0], [0], [0]])

x_new = X_start
# Least Squares iteration --> tip: use a loop!
x_hat = []
arr_b = []
arr_x = []
arr_r = []

for i in range(max_iter):
    # Initialize the structures for the Least Squares (A, b)
    b = []
    A = []

    # For each satellite compute:
        # - topocentric positions of satellite (approx distance, elevation and azimuth)
        # - tropospheric and ionospheric corrections
        # - b = LS known term vector
        # - A = LS design matrix
    x_prev = x_new
    phi_u, lambda_u, h = GC2Geodetic([x_prev[0,0],x_prev[1,0],x_prev[2,0]])
    
    arr_ro = []
    arr_el = []
    arr_cd = []
    arr_io = []

    for j in range(len(sat_ids)):

        xyz_sat_j =np.array([xyz_sat[j,:]]).T
        ro_approx, az, el = topocent(x_prev[0:3], xyz_sat_j)
        
        b.append([ ro_approx - c*dtS[j][0] + iono_correction(phi_u, lambda_u, az, el, time_rx, alpha, beta) + tropo_correction(h, el)])
        e_j = (x_prev[0:3] - xyz_sat_j)/ro_approx
        A.append([e_j[0,0], e_j[1,0], e_j[2,0], 1])
        arr_ro.append([ro_approx])
        arr_el.append([el])
        arr_io.append([iono_correction(phi_u, lambda_u, az, el, time_rx, alpha, beta) + tropo_correction(h, el)])

     
        #e_j = np.array([ (x_prev[0:3] - xyz_sat_j)/ro_approx ])
        #A.append([a for a in e_j.squeeze()])
        #A = np.array([x + [1] for x in A])
   



    # Implement LS solution and estimate the corrections to receiver pos
    arr_ro = np.array(arr_ro)
    arr_el = np.array(arr_el)
    arr_io = np.array(arr_io)
    A = np.array(A)
    b = np.array(b)
    obs_corr = pr_C1-b
    N = np.dot(np.transpose(A), A)
    est_corr = np.dot(np.linalg.inv(N), np.dot(np.transpose(A), obs_corr))

    # Estimated coord = approximate + estimated correction
    x_est = x_prev[0:3] + est_corr[0:3]
    x_new = np.array([[x_est[0,0], x_est[1,0], x_est[2,0], est_corr[3,0]/c]]).T
    
    arr_r.append(max(abs(arr_ro))[0])
    arr_x.append(max(abs(x_prev[0:3]))[0])
    arr_b.append(max(abs(b))[0])
    x_hat.append(max(abs(est_corr[0:3]))[0])
    
    
    # Check convergence of the resut, in case break the loop
    if max(abs(est_corr[0:3])) <= convergence:
        print('\n\nBreak at cycle: ', i+1)
        break
    # Check at the end that convergence didn't fail
    if i == (max_iter-1):
        print('Convergence failed')
'''    
fig, ax = plt.subplots(figsize=(10,6))
ax.plot(arr_x)
ax.set(xlabel='iteration', ylabel='x', title = 'max value of receiver coordinates x')

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(arr_r)
ax.set(xlabel='iteration', ylabel='ro', title = 'ro_approx = sqrt(x_prev[0:3] - xyz_sat_j**2)')

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(arr_b)
ax.set(xlabel='iteration', ylabel='b', title = 'max value of (ro_approx - c*dtS + iono_corr + tropo_corr)')

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(x_hat)
ax.set(xlabel='iteration', ylabel='max value of the estimated correction[0:3]', title = 'x_hat')

plt.show()
'''

# Final estimated unknowns
x_new
# LS residuals and sigma2
y = np.dot(A,x_new) + b
dy = np.abs(pr_C1 - y)    # VERIFY THAT IS AN ARRAY
sigma2 = np.dot(dy.T,dy)/(len(b) - len(x_new)) # m - n -> n. satellites - 4

# Covariance matrix of the estimated coordinates
C_xx = sigma2*np.linalg.inv(N)
# C_xx = np.linalg.inv(N)

# PDOP computeation
Q_ee = np.linalg.inv(N)[0:3,0:3]

#angles = [ phi_u, lambda_u]
angles = [rad2deg(phi_u) , rad2deg(lambda_u)]
Q_LC = np.dot(R_0(angles),np.dot(Q_ee,R_0(angles).T))  # INPUT VALUE OF R_0
#PDOP = np.sqrt(Q_LC[i,i] for i in range(len(Q_LC)))
PDOP = np.sqrt(sum(np.diag(Q_LC))  ) 
print('----------------------------------------')
print('Coordinates of the receiver:')
print('X: ', x_new[0,0], '\nY: ', x_new[1,0], '\nZ: ',x_new[2,0])
print('Clock offset of the receiver: \n', x_new[3,0])
print('PDOP: \n', PDOP)

'''Repeat with cut-off angle of 5°'''
# Initialization of the parameters
X_start = x_new # result previously obtained
cutoff = 5

# Loop over all the available satellites:
    # - Compute the elevations
    # - Discard all the satellites with elevation < 5°
new_el = []
k = -1
new_sat = []
for x in arr_el.squeeze():
    k += 1
    if x > cutoff:
        new_el.append([x])
        new_sat.append(k)
new_el = np.array(new_el)

'''
rang_sat = []
for x in range(len(arr_el)):
    if arr_el.squeeze()[x] > cutoff:
        rang_sat.append(x)

'''

# new_el = [if x < cutoff: x for x in arr_el]






















# Initialization of the parameters already done!

# Least Squares iteration --> tip: use a loop!


for i in range(max_iter):
    # Initialize the structures for the Least Squares (A, b)
    b = []
    A = []

    # For each satellite compute:
        # - topocentric positions of satellite (approx distance, elevation and azimuth)
        # - tropospheric and ionospheric corrections
        # - b = LS known term vector
        # - A = LS design matrix
    x_prev = x_new
    phi_u, lambda_u, h = GC2Geodetic([x_prev[0,0],x_prev[1,0],x_prev[2,0]])
    
    for j in new_sat:

        xyz_sat_j =np.array([xyz_sat[j,:]]).T
        ro_approx, az, el = topocent(x_prev[0:3], xyz_sat_j)
        
        b.append([ ro_approx - c*dtS[j][0] + iono_correction(phi_u, lambda_u, az, el, time_rx, alpha, beta) + tropo_correction(h, el)])
        e_j = (x_prev[0:3] - xyz_sat_j)/ro_approx
        A.append([e_j[0,0], e_j[1,0], e_j[2,0], 1])
     
   
    # Implement LS solution and estimate the corrections to receiver pos

    A = np.array(A)
    b = np.array(b)
    obs_corr = pr_C1[new_sat]-b
    N = np.dot(np.transpose(A), A)
    est_corr = np.dot(np.linalg.inv(N), np.dot(np.transpose(A), obs_corr))

    # Estimated coord = approximate + estimated correction
    x_est = x_prev[0:3] + est_corr[0:3]
    x_new = np.array([[x_est[0,0], x_est[1,0], x_est[2,0], est_corr[3,0]/c]]).T
  
    x_hat.append(max(abs(est_corr[0:3]))[0])
    
    # Check convergence of the resut, in case break the loop
    if max(abs(est_corr[0:3])) <= convergence:
        print('\n\nBreak at cycle: ', i+1)
        break
    # Check at the end that convergence didn't fail
    if i == (max_iter-1):
        print('Convergence failed')


# Final estimated unknowns
x_new
# LS residuals and sigma2
y = np.dot(A,x_new) + b
dy = np.abs(pr_C1[new_sat] - y)    # VERIFY THAT IS AN ARRAY
sigma2 = np.dot(dy.T,dy)/(len(b) - len(x_new)) # m - n -> n. satellites - 4

# Covariance matrix of the estimated coordinates
C_xx = sigma2*np.linalg.inv(N)
# C_xx = np.linalg.inv(N)

# PDOP computeation
Q_ee = np.linalg.inv(N)[0:3,0:3]

#angles = [ phi_u, lambda_u]
angles = [rad2deg(phi_u) , rad2deg(lambda_u)]
Q_LC = np.dot(R_0(angles),np.dot(Q_ee,R_0(angles).T))  # INPUT VALUE OF R_0
#PDOP = np.sqrt(Q_LC[i,i] for i in range(len(Q_LC)))
PDOP = np.sqrt(sum(np.diag(Q_LC))  ) 

# Implement LS solution for the new satellite configuration
print('----------------------------------------')
print('In view satellites: \n', len(new_sat))
print('Coordinates of the receiver:')
print('X: ', x_new[0,0], '\nY: ', x_new[1,0], '\nZ: ',x_new[2,0])
print('Clock offset of the receiver: \n', x_new[3,0])
print('PDOP: \n', PDOP)


sys.stdout = orig_stdout
f.close()