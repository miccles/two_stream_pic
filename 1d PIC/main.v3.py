from functions import *
import numpy as np




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