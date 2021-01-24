# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 13:00:58 2020

@author: ethan
"""
import numpy as np
import matplotlib.pyplot as plt
import makeF as F

# number of intervals
vectorN = [10,20,40,80]
errorVect = []
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

def rungeKutta(gridX):
    # The number of grid points
    nx = gridX - 1 
    
    # Spacing between grid points
    dx = 1.0 / (nx - 1)
    
    nt = (gridX ** 2) * 2
    
    tfinal = 1.0
    dt = tfinal / nt
    
    global mu
    mu = dt / dx**2
    
    x = np.linspace(0, 1, nx)
    

    # Interpolates the IC onto the mesh
    u    = IC(x)
    unew = np.zeros(u.shape)
    
    #def __init__(self, nx, phi0, phi1):
    for j in range(nt):
        rungeKutta = F.HeatEquation(gridX, phi0, phi1)
        rungeKutta(j,u,unew)
        
       #u[:] = unew[:]
        #build the matrix
    # for i in range(1, nx-1):
    #     A[i, i-1] = A[i, i+1] = 1.0
    #     A[i, i] = -2.0    
    # k1 = k2 = k3 = k4 = 0
    
        
    # for j in range(nt):
    #     #input boundary conditions
    #     unew[0]  = phi0((j+1)*dt)
    #     unew[-1] = phi1((j+1)*dt)
        
    #     for i in range (1,nx-1):
    #         k1 = mu*u[i-1] * .5
    #         k2 = mu*u[i] *.5 + k1*.5
    #         k3 = mu*u[i+1] *.5 + k2*.5
    #         k4 = mu*u[i] *.5 + k3*.5
            
    #         unew[i] = u[i] + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
            
            
    #     u[:] = unew[:]
    #real solution
    usol = utrue(x, tfinal)

    #print %error as decimal
    error = np.max(np.abs(u - usol))
    
    #graph compputed solution and real solution
    plt.plot(x, u, x, usol)
    plt.show()
    
    return error

print(rungeKutta(20))