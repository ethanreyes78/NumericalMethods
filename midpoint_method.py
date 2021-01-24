# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 10:33:34 2020
@author: ethan
"""

#from scipy.sparse import diags
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

def midpointMethod(gridX,nt):
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
    u_half = (np.zeros(u.shape))
    unew   = (np.zeros(u.shape))
    
    for j in range(nt):
        u_half[0] = phi0(j*dt+dt*.5)
        u_half[-1] = phi1(j*dt+dt*.5)
        
        for i in range(1, nx-1):
            u_half[i] = u[i] + 0.5 *mu*(u[i-1] -2.0*u[i] + u[i+1])
    
        unew[0] = phi0(j*dt+dt)
        unew[-1] = phi1(j*dt+dt)
    
        for i in range(1, nx-1):
            unew[i] = u[i] + mu*(u_half[i-1] - 2.0*u_half[i] + u_half[i+1])
    
        u[:] = unew[:]
    
    #calculate true solution
    usol = utrue(x, tfinal)
    
    #calculate error
    error = np.max(np.abs(u - usol))
    errorPlot = np.abs(u-usol)
    
    #plot the solution
    plt.plot(x, u, x, usol) 
    plt.show()
        
    return error, errorPlot


nx = 80
nt = (nx ** 2) * 2
# print(midpointMethod(nx,nt)[0])
diff = np.abs( midpointMethod(nx, nt)[0] - midpointMethod(nx,1323000)[0])
print(diff)