# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:56:35 2020

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

def trapezoidMethod(gridX,nt):
    # The number of grid points
    nx = gridX + 1
    
    # Spacing between grid points
    dx = 1.0 / (nx - 1)
    #nt = (gridX ** 2) * 2    
    tfinal = 1.0
    dt = tfinal / nt

    mu = dt / dx**2
    x = np.linspace(0, 1, nx)
    
    # Interpolates the IC onto the mesh
    u = IC(x)
    
    u_star = (np.zeros(u.shape))
    unew = (np.zeros(u.shape))
    
    for j in range(nt):    
        #input BC into vector
        u_star[0] = phi0(j*dt+dt)
        u_star[-1] = phi1(j*dt+dt)
        
        for i in range(1,nx-1):
            u_star[i] = u[i] + mu*( u[i-1] - 2.0*u[i] + u[i+1] )
        
        unew[0] = phi0(j*dt+dt)
        unew[-1] = phi1(j*dt+dt)
        
        for i in range(1,nx-1):
            unew[i] = u[i] + .5*mu*( u_star[i-1] - 2.0*u_star[i] + u_star[i+1] + u[i-1] -2.0*u[i] + u[i+1])
        u[:] = unew[:]
        
    #real solution
    usol = utrue(x, tfinal)
    
    #print %error as decimal
    error = np.max(np.abs(u - usol))
    errorPlot = np.abs(u-usol)
    
    #graph compputed solution and real solution
    # plt.plot(x, u, x, usol)
    # plt.show()
    
    return error, errorPlot

# nx = 20
# nt = (nx ** 2) * 2
# # print(nt)
# print(trapezoidMethod(nx,nt)[0])
