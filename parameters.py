import numpy as np


##### Simulation parameters #####

# Domain parameters #
Nx = 100      # Number of cells
Lx = Nx       # Length of the domain
dx = Lx / Nx  # Size of each cell

# Particle parameters #
ppc0 = 5      # Number of particles per cell
Np = Nx * ppc0  # Number of computational particles

# Time resolution #
c = 1
eta = 0.1     # number < 1 that is used for dt
dt = eta * dx / c  # Simulation timestep

# Spatial resolution #
kappa = 2     # number > eta/2 that is used for c_omp
c_omp = kappa * dx    # Plasma skin depth in dx units

# Plasma parameters #
q_m = 1       # Charge to mass ratio
n0 = Np / Lx  # Plasma number density
omega_p = c / c_omp    # Plasma frequency
q = omega_p ** 2 / (4 * np.pi * n0 * q_m)  # Particle charge
m = q / q_m                                # Particle mass



