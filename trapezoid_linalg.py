# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:14:58 2020

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

def trapezoidLinAlg(gridX,nt):
    # The number of grid points
    nx = gridX + 1
    
    # Spacing between grid points
    dx = 1.0 / (nx - 1)    
    #nt = (gridX** 2) * 2
    tfinal = 1.0
    dt = tfinal / nt
    
    global mu
    mu = dt / dx**2
    
    x = np.linspace(0, 1, nx)
    # Interpolates the IC onto the mesh
    global u
    u = IC(x)
    
    unew = np.zeros(u.shape)
    u_star = np.zeros(u.shape)
    
    A = np.zeros((nx, nx))
    
    for i in range(1, nx-1):
        A[i, i-1] = A[i, i+1] = 1.0
        A[i, i] = -2.0
    
    
    for j in range(nt):
        #compute ustar and input boundary conditions
        u_star[:]  = u[:] + mu * np.dot(A,u)
        u_star[0]  = phi0((j+1)*dt)
        u_star[-1] = phi1((j+1)*dt)
        
        unew[:]  = u[:] + .5* mu * ( np.dot(A, u_star) + np.dot(A, u) )
        unew[0]  = phi0((j+1)*dt)
        unew[-1] = phi1((j+1)*dt)
        
        u[:] = unew[:]
    
    #real solution
    usol = utrue(x, tfinal)
    
    #print %error as decimal
    error = np.max(np.abs(u - usol))
    
    global errorPlot
    errorPlot = np.abs(u-usol)
    
    #graph compputed solution and real solution
    # plt.plot(x, u, x, usol)
    # plt.show()
    
    return error
# nx = 20
# nt = (nx ** 2) * 2
# print(trapezoidLinAlg(nx,nt))

# for x in vectorN:
#     errorVect.append(trapezoidLinAlg(x))

# print(errorVect)
# print(mu)