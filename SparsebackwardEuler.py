# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 13:42:24 2020

@author: ethan
"""

from scipy.sparse import diags
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
# from scipy.sparse import csc_matrix
# from scipy import sparse

# number of intervals
n = 10
# The number of grid points
nx = n + 1

# Spacing between grid points
dx = 1.0 / (nx - 1)

nt = (n ** 2) * 2

tfinal = 1.0
dt = tfinal / nt

mu = dt / dx**2

x = np.linspace(0, 1, nx)

# initial condition for heat equation
def utrue(x, t):
    return np.exp((-np.pi**2 * t)/4.0) * np.sin(np.pi*x*.5) + .5*np.exp(-np.pi**2 *(4.0)*t) * np.sin(2.0*np.pi*x)

def IC(x):
    return utrue(x, 0.0)

def phi0(t):
    return utrue(0.0, t)

def phi1(t):
    return utrue(1.0, t)

# Interpolates the IC onto the mesh
u    = IC(x)
unew = np.zeros(u.shape)
f    = np.zeros(u.shape)

k = np.array([-mu * np.ones(nx-1),  np.ones(nx) + 2.0*mu,  -mu * np.ones(nx-1)])
A = diags(k,[-1,0,1], format='csc')

A[0,0] = 1.0
A[0,1] = 0.0
A[-1,-2] = 0.0
A[-1,-1] = 1.0

for j in range(nt):
    #input boundary conditions    
    
    f[0] = 0.0
    f[1: -1] = u[1: -1]
    f[-1] = np.exp((-np.pi**2 *(j+1)*dt)/4.0)
    
    unew[:] = spsolve(A,f)
    # ufake[:] = np.linalg.solve(A,f)
    
    u[:] = unew[:]
    

#real solution
usol = utrue(x, tfinal)

#print %error as decimal
print(np.max(np.abs(u - usol)))

#graph compputed solution and real solution
plt.plot(x, u, x, usol) 
plt.show()


