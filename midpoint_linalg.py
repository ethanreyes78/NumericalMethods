# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 14:59:21 2020

@author: ethan
"""


import numpy as np
import matplotlib.pyplot as plt


# number of intervals
errorPlot = 0
u = []
mu = 0

# initial condition for heat equation
def utrue(x, t):
    return np.exp((-np.pi**2 * t)/4.0) * np.sin(np.pi*x*.5) + .5*np.exp(-np.pi**2 *(4.0)*t) * np.sin(2.0*np.pi*x)

def IC(x):
    return utrue(x, 0.0)

def phi0(t):
    return utrue(0.0, t)

def phi1(t):
    return utrue(1.0, t)

def midpointLinAlg(gridX,nt):
    # The number of grid points
    nx = gridX + 1
    
    # Spacing between grid points
    dx = 1.0 / (nx - 1)
    
    #nt = (gridX ** 2) * 2
    
    tfinal = 1.0
    dt = tfinal / nt

    global mu    
    mu = dt / dx**2
    x = np.linspace(0, 1, nx)
    
    # Interpolates the IC onto the mesh
    global u
    u = IC(x)
    
    unew = np.zeros(u.shape)
    u_half = np.zeros(u.shape)
    
    A = np.zeros((nx, nx))
    
    for i in range(1, nx-1):
        A[i, i-1] = A[i, i+1] = 1.0
        A[i, i] = -2.0
    
    for j in range(nt):
        
        u_half[:] = u[:] + .5*mu * np.dot(A,u)
        u_half[0] = phi0((j+.5)*dt)
        u_half[-1] = phi1((j+.5)*dt)
        
        unew[:] = u[:] + mu * np.dot(A, u_half)
        unew[0] = phi0((j+1)*dt)
        unew[-1] = phi1((j+1)*dt)
        
        u[:] = unew[:]
        
        
    #calculate true solution
    usol = utrue(x, tfinal)
    
    #calculate error
    error = np.max(np.abs(u - usol))
    global errorPlot
    errorPlot = np.abs(u-usol)
    
    #graph real solution and our solution
    # plt.plot(x, u, x, usol) 
    # plt.show()
    
    return error, errorPlot

# nx = 20
# nt = (nx ** 2) * 2
# print(midpointLinAlg(nx,nt))
# for x in vectorN:
#     errorVect.append(midpointLinAlg(x))
    
# print(errorVect)
# print(mu)
