import numpy as np

Nx = 100      # Number of cells
ppc0 = 5      # Number of particles per cell
q_m = 1       # Charge to mass ratio
c = 1         # Speed of light

Lx = Nx       # Length of the domain
dx = Lx / Nx  # Size of each cell
Np = Nx * ppc0  # Number of computational particles
n0 = Np / Lx  # Plasma number density

eta = 0.1     # number < 1 that is used for dt
dt = eta * dx / c  # Simulation timestep

kappa = 2     # number > eta/2 that is used for c_omp
c_omp = kappa * dx    # Plasma skin depth in dx units

omega_p = c / c_omp    # Plasma frequency
q = omega_p ** 2 / (4 * np.pi * n0 * q_m)  # Particle charge
m = q / q_m                                # Particle mass



