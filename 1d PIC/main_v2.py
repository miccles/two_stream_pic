import numpy as np
import math
import random as rnd
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve

from functions import *


# First attempt on a 1D PIC simulation for two-stream instability


dt = 0.05                           # time step
Nt = 12                             # number of time steps (set to -1 for no loops)
tfinal = Nt * dt                    # total time
tex = 10                            # time interval for data export

ng = 64                             # number of grid boxes
L = 2 * np.pi                       # domain length in dx units
dx = L / ng                         # dimensionless size of grid boxes (=1 in our units)

Np = 10                             # number of particles/cell
N = ng * Np                         # number of computational particles
n = 1                               # number density of particles in units of N/L.

#Quantities in dimaneionless units (see notes)#
qp = -1 / (Np * ng)                 #charge of computational particle normalized to the total charge of the domain
mp = 1 / (Np * ng)                  #mass of computational particle normalized to the total mass of the domain
qm = qp/mp
wp = np.sqrt(qp ** 2 * n / mp)      #plasma frequency in dimensionless units

vd = 0.2                            # drift velocity
v0 = 0.1 * vd                       # random component of velocity
                                    # final velocities of particles will be centered around +- (v0+vd)

# perturbed quantities
vp1 = 0.0                           # perturbation of velocity
xp1 = 1.0 * L / N                   # perturbation of position
mode = 1                            # perturbation mode

# initial particle loading#
xp = np.linspace(0, L - L / N, N)  # - L/N because of BCs (at x=L we have the x=0 particle)
xi = xgrid(ng, L)
vp = v0 * np.random.normal(1, 1, N)
vp = [(vp[i] + vd) * (-1) ** i for i in range(len(vp))]

# Boundary Conditions of field#
phiL = 0
phiR = 0

# Initial Rho and d matrices (Before the perturbation)#
rho = rhoavg(xp, xi, qp, dx)
d = [-r for r in rho]
# Phi vector#
phi = phifunc(d, ng)
# Ei vector #
eli = efavg(phi, ng, dx)  # Electric field at the center of cells
# Ep vector #
elp = efparticle(xp, xi, dx, eli)

# perturbation#
vp = [vp[i] + vp1 * math.sin(2 * np.pi * xp[i] * mode / L) for i in range(len(xp))]
xp = [xp[i] + xp1 * math.sin(2 * np.pi * xp[i] * mode / L) for i in range(len(xp))]


# Calculation of rhoavg, phiavg, Ei and Ep after perturbation#

# Main Computational Cycle #
if Nt >= 0:
    # Evolving v to dt/2 once#
    vp = vpup(vp, dt / 2, qm, elp)
    for i in range(Nt + 1):
        # Position Update & BCs#
        xp = bcp(xpup(xp, vp, dt), L)
        # Rho and d vectors
        rho = rhoavg(xp, xi, qp, dx)
        d = [-4 * np.pi * r for r in rho]
        # Phi vector#
        phi = phifunc(d, ng)
        # Ei vector #
        eli = efavg(phi, ng, dx)  # Electric field at the center of cells
        # Ep vector #
        elp = efparticle(xp, xi, dx, eli)
        #Velocity update#
        vp = vpup(vp, dt, qm, elp)
        plt.scatter([xp[2 * i] for i in range(int(N / 2))], [vp[2 * i] for i in range(int(N / 2))])
        plt.scatter([xp[2 * i + 1] for i in range(int(N / 2))], [vp[2 * i + 1] for i in range(int(N / 2))])
        plt.show()
    # Final evolution of vp by -dt/2 once#
    vp = vpup(vp, -dt / 2, qm, elp)

#Plotting#
plt.scatter([xp[2*i] for i in range(int(N/2))], [vp[2*i] for i in range(int(N/2))])
plt.scatter([xp[2*i + 1] for i in range(int(N/2))], [vp[2*i + 1] for i in range(int(N/2))])
plt.show()

plt.scatter(xi, phi)
plt.show()