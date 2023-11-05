import numpy as np
from functions import *

dt = 0.005  # time step
Nt = 40  # number of time steps (set to -1 for no loops)
tfinal = Nt * dt  # total time
tex = 10  # time step interval for data export

ng = 64  # number of grid boxes
L = 2 * np.pi  # domain length in dx units
dx = L / ng  # dimensionless size of grid boxes
xi = xgrid(ng, L)

Np = 20  # number of particles/cell
N = ng * Np  # number of computational particles
n0 = N / L  # initial (and average) number density of particles in units of N/L.

# Quantities in dimaneionless units (see notes)#
qp = -1  # charge of computational particle
mp = 1  # mass of computational particle
qm = qp / mp
wp = np.sqrt(qp ** 2 * n0 / mp)  # plasma frequency in dimensionless units