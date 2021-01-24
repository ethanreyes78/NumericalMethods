# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 14:15:57 2020

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

def backwardEulerNoBdNodes(gridX,nt):
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
    u    = IC(x[1:-1])
    unew = np.zeros(u.shape)
    f    = np.zeros(u.shape)
    
    A = np.zeros((nx-2, nx-2))
    
    A[0, 0] = 1.0 + 2.0 * mu
    A[0, 1] = -mu
    for i in range(1, nx-3):
        A[i, i-1] = A[i, i+1] = -mu
        A[i, i] = 1.0+2.0*mu
    A[-1, -2] = -mu
    A[-1, -1] = 1.0 + 2.0 * mu

    for j in range(nt):
        f[:] = u[:]
        f[0] = f[0] + mu * phi0((j*dt+dt))
        f[-1] = f[-1] + mu * phi1((j*dt+dt))
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
    # plt.plot(x, ucomputed, x, usol)
    # plt.show()
    
    return error, errorPlot, mu
# nx = 10
# nt = (nx ** 2) * 2
# print(backwardEulerNoBdNodes(nx,100))
# for x in vectorN:
#     errorVect.append(backwardEulerNoBdNodes(x))

# print(errorVect)
# print(mu)
#Structure and Interpretation of Computer Programs (Abelson & Susan)
