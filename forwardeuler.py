# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 10:33:34 2020

@author: ethan
"""

import numpy as np
import matplotlib.pyplot as plt


# initial condition for heat equation
def utrue(x, t):
    return np.exp((-np.pi**2 * t)/4.0) * np.sin(np.pi*x*.5) + .5*np.exp(-np.pi**2 *(4.0)*t) * np.sin(2.0*np.pi*x)

def IC(x):
    return utrue(x, 0.0)

def phi0(t):
    return utrue(0.0, t)

def phi1(t):
    return utrue(1.0, t)

def forwardEuler(gridX, nt):
    # The number of grid points
    nx = gridX + 1
    
    # Spacing between grid points
    dx = 1.0 / (gridX)
    #nt = (gridX ** 2) * 2
    
    tfinal = 1.0
    dt = tfinal / nt
    
    mu = dt / dx**2
    # print(mu)

    x = np.linspace(0, 1, nx)

    # Interpolates the IC onto the mesh
    u = IC(x)
    
    unew = (np.zeros(u.shape))
    
    for j in range(nt):
        unew[0] = phi0((j+1)*dt)
        unew[-1] = phi1((j+1)*dt)
        
        for i in range(1, nx-1):
            unew[i] = u[i] + mu*(u[i-1] - 2.0*u[i] + u[i+1])
            
    
        u[:] = unew[:]
    
    #calculate true solution
    usol = utrue(x, tfinal)   
    #calculate error
    error = np.max(np.abs(u - usol))
    
    errorPlot = np.abs(u-usol)
    #plot solution
    # plt.plot(x, u, x, usol) 
    # plt.show()
    
    
    return error, errorPlot, usol

# nx = 10
# y = (nx ** 2) * 2
# print(y)
# print(forwardEuler(nx, 190)[0])
# print(forwardEuler(nx, y)[2])
# run the method at different intervals then store errors to see order of method
# for x in vectorN:
#     errorVect.append(forwardEuler(x))

# print(errorVect)
# print(mu)
