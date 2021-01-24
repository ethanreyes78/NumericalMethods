# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 08:39:02 2020

@author: ethan
"""

import numpy as np
import matplotlib.pyplot as plt

errorPlot = 0

# initial condition for heat equation
def utrue(x, t):
    return np.exp((-np.pi**2 * t)/4.0) * np.sin(np.pi*x*.5) + .5*np.exp(-np.pi**2 *(4.0)*t) * np.sin(2.0*np.pi*x)

def IC(x):
    return utrue(x, 0.0)

def phi0(t):
    return utrue(0.0, t)

def phi1(t):
    return utrue(1.0, t)

def crankNicolSolver(gridX, nt):
    # The number of grid points
    nx = gridX + 1
    
    # Spacing between grid points
    dx = 1.0 / (nx - 1)
    
    #nt = nx 
    
    tfinal = 1.0
    dt = tfinal / nt
    
    mu = dt / dx**2


    x = np.linspace(0, 1, nx)
    
    # Interpolates the IC onto the mesh
    u    = IC(x[1:-1])  
    unew = np.zeros(u.shape)
    f    = np.zeros(u.shape)
    
    A = np.zeros((nx-2, nx-2))
    
    #builds the matrix
    A[0, 0] = 1.0+mu
    A[0, 1] = -mu*.5
    for i in range(1, nx-3):
        A[i, i-1] = A[i, i+1] = -mu*.5
        A[i, i] = 1.0+mu
    A[-1, -2] = -mu*.5
    A[-1, -1] = 1.0+mu

    for j in range(nt):
        # Difference on old solution
        f[0]  = u[0]  * (1-mu) + 0.5 * mu * (u[1]  + (phi0(j*dt) + phi0((j+1)*dt)) )
        for i in range(1, nx-3):
            f[i] = u[i] * (1 - mu) + 0.5 * mu * (u[i-1] + u[i+1])
        f[-1] = u[-1] * (1-mu) + 0.5 * mu * (u[-2] + (phi1(j*dt) + phi1((j+1)*dt)) )
        # print(u)
        unew[:] = np.linalg.solve(A, f)
        
        u[:] = unew[:]
    
    #real solution
    usol = utrue(x, tfinal)
    
    ucomputed = np.array([phi0(tfinal)] + list(u) + [phi1(tfinal)])

    #print %error as decimal
    error = np.max(np.abs(ucomputed - usol))
    
    errorPlot = np.abs(ucomputed-usol) 
    
    #graph compputed solution and real solution
    plt.plot(x, ucomputed, x, usol)
    plt.show()
    # plt.plot(x, np.abs(ucomputed-usol))
    
    
    return error, errorPlot, x

nx = 80
nt = (nx ** 2) * 2
# print(crankNicolSolver(nx, nt)[0])
diff = np.abs(crankNicolSolver(nx, 120)[0] - crankNicolSolver(nx, nt)[0])
print(diff)







