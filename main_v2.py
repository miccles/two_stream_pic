import numpy as np
import math
import random as rnd
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve


# First attempt on a 1D PIC simulation for two-stream instability


def b0(x):
    if np.abs(x) <= 1 / 2:
        return 1
    else:
        return 0


def b1(x):
    if -1 <= x <= 1:
        return -np.abs(x) + 1
    else:
        return 0


def xgrid(ng, L):
    dx = L / ng
    x = [i * dx for i in range(ng + 1)]
    return x


def rhoavg(xp, xg, q, dx):
    global L, ng
    rho = [0] * (ng + 1)
    for i in range(ng):
        for p in range(len(xp)):
            rho[i] += (q / dx) * b1((xg[i] - xp[p]) / dx)
            if xg[ng - 1] <= xp[p] <= L and i == 0:
                rho[i] += (q / dx) * b1((xg[i] - (xp[p] - L)) / dx)
            elif 0 <= xp[p] <= xg[0] and i == ng - 1:
                rho[i] += (q / dx) * b1((xg[i] - (xp[p] + L)) / dx)
    rho [ng] = rho[0]
    return rho


def phifunc(d, ng):
    tempd = [de for de in d]
    a = [1] * (ng + 1)
    b = [-2] * (ng + 1)
    c = [1] * (ng + 1)
    b[0] = 1
    b[ng] = 1
    a[0] = 0
    a[ng] = 0
    c[0] = 0
    c[ng] = 0
    c[0] = c[0] / b[0]
    tempd[0] = tempd[0] / b[0]
    for i in range(1, ng):
        c[i] = c[i] / (b[i] - a[i] * c[i - 1])
        tempd[i] = (tempd[i] - a[i] * tempd[i - 1]) / (b[i] - a[i] * c[i - 1])
    tempd[ng] = (tempd[ng] - a[ng] * tempd[ng - 1]) / (b[ng] - a[ng] * c[ng - 1])
    x = [0] * (ng + 1)
    x[ng] = tempd[ng]
    for i in range(ng - 1, -1, -1):
        x[i] = tempd[i] - c[i] * x[i + 1]
    return x


def efavg(phi, ng, dx):
    el = [0] * (ng + 1)
    el[0] = -(phi[1] - phi[0]) / (1 * dx)
    el[ng] = -(phi[ng] - phi[ng - 1]) / (1 * dx)
    for i in range(1, ng):
        el[i] = -(phi[i + 1] - phi[i - 1]) / (2 * dx)
    return el


def efparticle(xp, xg, dx, elgrid):
    global L, N, ng
    efp = [0] * N
    for p in range(N):
        # if 0 <= xp[p] <= dx / 2:
        #     efp[p] = elgrid[0] * b1((xg[0] - xp[p]) / dx) #+ elgrid[-1] * b1((xg[-1] - xp[p]) / dx)
        # elif L - dx / 2 <= xp[p] <= L:
        #     efp[p] = elgrid[ng - 1] * b1((xg[ng - 1] - xp[p]) / dx)
        # else:
            try:
                for i in range(int(np.ceil(xp[p] / dx - 1)), int(np.floor(xp[p] / dx + 1) + 1)):
                    efp[p] += elgrid[i] * b1((xg[i] - xp[p]) / dx)
            except IndexError:
                print("Index Error for i=", i)
    return efp


def xpup(xp, vp, dt):
    newpos = [xp[i] + vp[i] * dt for i in range(len(xp))]
    return newpos


def vpup(vp, dt, qm, Ep):
    newvel = [vp[i] + qm * dt * Ep[i] for i in range(len(vp))]
    return newvel


def bcp(xp, L):
    for i in range(len(xp)):
        if xp[i] < 0:
            xp[i] += L
        elif xp[i] > L:
            xp[i] -= L
    return xp


dt = 0.05  # time step
Nt = 12 #number of time steps (set to -1 for no loops)
tfinal = Nt * dt #total time
tex = 10  # time interval for data export

ng = 64  # number of grid boxes
L = 2 * np.pi  # domain length in dx units
dx = L / ng  # dimensionless size of grid boxes (=1 in our units)

Np = 10  # number of particles/cell
N = ng * Np  # number of computational particles
n = 1  # number density of particles in units of N/L.

#Quantities in dimaneionless units (see notes)#
qp = -1 / (Np * ng) #charge of computational particle normalized to the total charge of the domain
mp = 1 / (Np * ng)  #mass of computational particle normalized to the total mass of the domain
qm = qp/mp
wp = np.sqrt(qp ** 2 * n / mp)  #plasma frequency in dimensionless units

vd = 0.2  # drift velocity
v0 = 0.1 * vd  # random component of velocity
# final velocities of particles will be centered around +- (v0+vd)

# perturbed quantities
vp1 = 0.0  # perturbation of velocity
xp1 = 1.0 * L / N  # perturbation of position
mode = 1  # perturbation mode

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