import math
import numpy as np
#np.set_printoptions(precision=3)

def deg2rad(deg): # convert from degree to radius
    return deg*np.pi/180

def sex2deg(sex): # convert from sexadegimal to degree
    return sex[0] + sex[1]/60 + sex[2]/3600 

def rad2deg(rad):
    return rad*180/np.pi

def deg2sex(deg):
	d = np.fix(deg)
	p =(deg-d)*60
	s =(p - np.fix(p))*60
	return [d, np.fix(p), s]

def R(csi,eta,alpha): # rotation matrix 
      
    R_x = np.array([[1, 0, 0],
                   [0, np.cos(csi), -np.sin(csi)],
                   [0, np.sin(csi), np.cos(csi)]])
    R_y = np.array([[np.cos(eta), 0, -np.sin(eta)],
                   [0, 1, 0],
                   [np.sin(eta), 0, np.cos(eta)]])
    R_z = np.array([[np.cos(alpha), np.sin(alpha), 0],
                   [-np.sin(alpha), np.cos(alpha), 0], 
                   [0, 0, 1]])
    return (R_z @ R_y @ R_x)

def R_LL2LC(x_LL,csi,eta,alpha): # transform from Local Level to Local Cartesian
    return R(csi,eta,alpha).T.dot(x_LL)

def R_0(Geodetic_0): # rotation matrix 
      # Geodetic_0 = [0,0,0]
    phi = deg2rad(Geodetic_0[0])
    lambd = deg2rad(Geodetic_0[1])

    return np.array([[-np.sin(lambd),  np.cos(lambd), 0],
                     [-np.sin(phi)*np.cos(lambd), -np.sin(phi)*np.sin(lambd), np.cos(phi)],
                     [ np.cos(phi)*np.cos(lambd),  np.cos(phi)*np.cos(lambd), np.sin(phi)]])
             
def R_LC2GC(x_LC,x_Geodetic_0): # transforn from Local Cartesian to Global Cartesian
    return Geodetic2GC(x_Geodetic_0) + R_0(x_Geodetic_0).T.dot(x_LC)

a = 6378137
e2 = 0.0818191908426215**2

def Geodetic2GC(Geodetic_0): # transform from Geodetic (ITRF) to Global Cartesian

    phi = deg2rad(Geodetic_0[0])
    lambd = deg2rad(Geodetic_0[1])
    h = Geodetic_0[2]

    R_N = a/np.sqrt(1 - e2*(np.sin(phi))**2)

    x = (R_N + h)*np.cos(phi)*np.cos(lambd)
    y = (R_N + h)*np.cos(phi)*np.sin(lambd)
    z = (R_N*(1-e2) + h)*np.sin(phi)

    return [x, y, z]


def GC2Geodetic(x_GC):
    
    b2 = (a**2)*(1 - e2)

    eb2 = (a**2 - b2)/(b2)
    r = np.sqrt(x_GC[0]**2 + x_GC[1]**2)
    psi = np.arctan2(x_GC[2], (r*np.sqrt(1-e2)))

    l = np.arctan2(x_GC[1], x_GC[0])
    p = np.arctan2(x_GC[2] + eb2*np.sqrt(b2)*(np.sin(psi)**3), r - e2*a*(np.cos(psi)**3))
    R_n = a/np.sqrt(1 - e2*((np.sin(p))**2))
    h = r/np.cos(p) - R_n
    return [rad2deg(p), rad2deg(l), h]

def covariance(sigma,sigma_0,csi,eta,alpha,Geodetic_0): # covariance propagation
 
    C_LL = np.array([[sigma**2,0,0],
                     [0,sigma**2,0],
                     [0,0,sigma**2]])

    C_LC = R(csi,eta,alpha).T @ C_LL @ R(csi,eta,alpha)
    
    C_GC_0 = np.array([[sigma_0**2,0,0],
                       [0,sigma_0**2,0],
                       [0,0,sigma_0**2]])
    
    C_GC = C_GC_0 + R_0(Geodetic_0).T @ C_LC @ R_0(Geodetic_0)
    
    C = R_0(Geodetic_0) @ C_GC @ R_0(Geodetic_0).T
    
    # introduce an approximation using R0 with origin coordinates instead of A, B, C ones 
    # the error made is negligible since the approximation introduced is 'small' wrt the size order of the terrestrial coordinates 
    return C

# -*- coding: utf-8 -*-
"""
Positioning & Location Based Services
AA: 2022/2023

@author: Marianna Alghisi
"""




def eccAnomaly(M, e):
    '''
    Input:
    ----------
    M : Mean anomaly.
    e : eccenticity.

    Output:
    -------
    E : Eccentricity anomaly.

    '''
    # Inizialization of the Eccentricity anomaly value:
    E = M
    
    max_iter = 12 
    # it was 10 when using only GPS (convergence was achieved at 4-6 iterations)
    # now it set to 12 because QZSS PRN 193 can take 11 iterations to converge
    
    i = 0
    dE = 1
    
    while ((dE > 10**(-12)) and (i < max_iter)):
        Etemp = E
        E = M + e*np.sin(E)
        i = i + 1
        dE = abs(E - Etemp)
    
    if i == max_iter:
        raise ValueError('WARNING: Eccentric anomaly does not converge!')
    else:
        return E
        
def topocent(X, sat):
    '''
    Parameters
    ----------
    X : np.array of the receiver position in global cartesian.
    sat : np.array of the satellite position in global cartesian.

    Returns
    -------
    ro_approx : distance between the satellite and the receiver [m].
    el : elevation of the satellite with respect the receiver [rad].
    az : azimuth of the satellite with respect the receiver [rad].

    '''
    
    delta = sat - X

    ro_approx = np.sqrt(delta[0,0]**2 + delta[1,0]**2 + delta[2,0]**2) 
    
    X_g = GC2Geodetic([X[0,0], X[1,0], X[2,0]])
    lat0 = deg2rad(X_g[0])
    lon0 = deg2rad(X_g[1])
    
    R = np.array([[-np.sin(lon0), np.cos(lon0), 0],
                  [-np.sin(lat0)*np.cos(lon0), -np.sin(lat0)*np.sin(lon0), np.cos(lat0)],
                  [ np.cos(lat0)*np.cos(lon0), np.cos(lat0)*np.sin(lon0), np.sin(lat0)]
                    ])
    
    ENU = np.dot(R, delta)
    
    E, N, U = ENU[0,0], ENU[1,0], ENU[2,0]
    
    d_h = np.sqrt(E**2 + N**2)

    if d_h < 0.1:
        az = 0
        el = 90
    else:
        az = np.arctan2(E, N)
        el = np.arctan2(U, d_h)

        az = rad2deg(az)
        el = rad2deg(el)
    
    return [ro_approx, az, el]    

def iono_correction(phi_u, lambda_u, A, E, GPStime, alpha, beta):
    '''
    Computation of the pseudorange correction due to ionospheric effect.
    Klobuchar model.
    
    Parameters
    ----------
    phi_u : latitude in degrees.
    lambda_u : longitude in degrees.
    A : azimuth in degrees.
    E : elevation in degrees.
    GPStime : time in seconds.
    alpha : list containing alpha parameters.
    beta: list containing beta parameters.

    Returns
    -------
    T_iono : ionosperic effect.

    '''

    # elevation from 0 to 90 degrees
    E = abs(E)    

    # conversion to semicircles
    phi_u = phi_u / 180
    lambda_u = lambda_u / 180
    A = A / 180
    E = E / 180
    # Speed of light
    c = 2.99792458 * 10**8

    # Psi is the Earth's central angle between the user position and the earth projection of ionospheric intersection point
    psi = 0.0137 / (E + 0.11) - 0.022       
    phi_i = phi_u + psi * math.cos(A * math.pi)

    if abs(phi_i) <= 0.416:
        phi_i = phi_i
    elif phi_i > 0.416:
        phi_i = 0.416
    elif phi_i < -0.416:
        phi_i = -0.416

    # geodetic longitude of the earth projection of the ionospheric intersection point
    lambda_i = lambda_u + psi * math.sin(A * math.pi) / math.cos(phi_i * math.pi)
    # geomagnetic latitude of the earth projection of the ionospheric intersection point
    phi_m = phi_i + 0.064 * math.cos((lambda_i - 1.617) * math.pi)
    # The local time in seconds
    t = 4.32 * 10**4 * lambda_i + GPStime  
    if t >= 86400:
        t = t - 86400
    elif t < 0:
        t = t + 86400

    # Obliquity factor
    F = 1 + 16 * math.pow((0.53 - E), 3)        
    PER = 0
    for n in range(3):
        PER = PER + beta[n] * math.pow(phi_m, n)        

    if PER < 72000:
        PER = 72000    

    x = 2 * math.pi * (t - 50400) / PER
    
    AMP = 0
    for n in range(3):
        AMP  = AMP + alpha[n] * math.pow(phi_m, n)      # the coefficients of a cubic equation representing the amplitude of the vertical delay (4 coefficients - 8 bits each)

    if AMP < 0:
        AMP = 0
    
    if abs(x) >= 1.57:
        T_iono = c *  F * 5 * 10**-9
    else:
        T_iono = c * F * ((5 * 10**-9) + AMP * (1 - (x**2/2) + (x**4 / 24)))

    return T_iono

def tropo_correction (h, el):
    '''
    Computation of the pseudorange correction due to tropospheric refraction.
    Saastamoinen model.
    
    Parameters
    ----------
    h : height f the receiver [m].
    el : elevation of the satellite with respect the receiver [degrees].

    Returns
    -------
    tropoDelay : tropospheric effect [m].

    '''
    if (h > -500) and (h < 5000):
        eta = deg2rad(el)
        #eta satellite elevation in radians
        Po = 1013.25 #mBar
        To = 291.15 #degree Kelvin
        Ho = 50/100
        ho = 0
        height = h - ho # h is the ellipsoidal height of the receiver
        Pr = Po * (1-0.0000226 * height)**5.225 #pressure
        Tr = To - 0.0065 * height #temperature
        Hr = Ho * math.exp(-0.0006396 * height)
        er = Hr * math.exp(-37.2465 + (0.213166 * Tr) - 0.000256908 * (Tr)**2) #humidity
        tropoDelay = 0.002277 / np.sin(eta) * (Pr + er * (1255/Tr + 0.05) - (np.tan(eta)) ** -2)
        return tropoDelay
    else:
        tropoDelay = 0
        return tropoDelay


def calcTrajectory(P0 , a_x, a_y, omegaz, epoch):
    # Initialize and compute velocities and delta positions in body frame

    pos = [[[P0[0,0]], [P0[1,0]]]]
    vx = [0]
    vy = [0]
    alpha = [0]
    x = pos[0][0]
    y = pos[0][1]
    dx = [0]
    dy = [0]
    ac = [0]
    for t in range(len(epoch)-1):
        dt = epoch[t+1] - epoch[t]
        vx.append(vx[t] + a_x[t+1] * dt)  
        x.append(x[t] + vx[t]*dt + 1/2*a_x[t+1]*dt**2)
        dx.append(x[t+1]-x[t])
        #dx = vx[t+1] * dt + (1/2) * a_x[t+1] * dt**2
        vy.append(vy[t] + a_y[t+1] * dt)  
        y.append(y[t] + vy[t]*dt + 1/2*a_y[t+1]*dt**2)
        dy.append(y[t+1]-y[t])

        
        ac.append(omegaz[t+1]*np.sqrt(vx[t]**2+vy[t]**2)) 
        #ac = omegaz[t+1]*vx[t]
        a_y[t+1] = a_y[t+1] - ac[t+1]
        vy[t+1] = vy[t] + a_y[t+1]*dt
        y[t+1] = y[t] + vy[t]*dt + 1/2*a_y[t+1]*dt**2
        dy[t+1] = y[t+1]-y[t]

        
        #resay = a_y[t+1] - ac
        #vy.append(vy[t] + resay*dt)
        #dy = vy[t+1] * dt + (1/2) * resay * dt**2

        alpha.append(alpha[t] + omegaz[t+1] * dt)
        R = np.matrix([ [ np.cos(alpha[t+1]), np.sin(alpha[t+1])],
                        [-np.sin(alpha[t+1]), np.cos(alpha[t+1])]])
        delta = np.array([[dx],
                          [dy]])
        delta =  R * delta
        new_pos = np.array(pos[t] + delta)
        pos.append([[new_pos[0,0]], [new_pos[1,0]]])
    return pos