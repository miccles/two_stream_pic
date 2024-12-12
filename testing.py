import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from functions import *
import timeit


Lap = laplacian_matrix(Nx) # Laplacian matrix
G = gradient_matrix(Nx) # Gradient matrix



# Generate initial conditions
pos, vel = generate_init_cond(Lx, Np, beam_v0, beam_dv0, beam_perp)

# Calculate average charge density
dens_avg = rho_avg(dx, Nx, pos, q)

# Calculate RHS matrix
d_matrix = -4 * np.pi * dx ** 2 * (dens_avg - q * n0)

# Solve for potential
phi = phi_sparse_solver(Lap, d_matrix)

# Calculate electric field of cells
E_cells = el_solver(G, phi, dx)

# Calculate electric field of particles
Ep = el_particles(dx, Nx, pos, E_cells)

# Calculate acceleration
acc = q * Ep / m

# Create a 4-pane plot
fig, axs = plt.subplots(2, 3, figsize=(18, 10))

# Top left: pos, vel scatter
axs[0, 0].scatter(pos, vel, s=2, color='blue', alpha=0.5)
axs[0, 0].set_title('Phase Diagram')
axs[0, 0].set_xlabel('Position')
axs[0, 0].set_ylabel('Velocity')

# Top right: dens_avg plot
axs[0, 1].plot(np.linspace(0, Nx*dx, Nx), dens_avg, color='green')
axs[0, 1].axhline(y=np.mean(dens_avg), color='red', linestyle='--')
axs[0, 1].set_title('Average Charge Density')
axs[0, 1].set_xlabel('Position')
axs[0, 1].set_ylabel('Density')


axs[0, 2].plot(np.linspace(0, Nx*dx, Nx), dens_avg - n0 * q, color='green')
axs[0, 2].axhline(y=np.mean(dens_avg - n0 * q), color='red', linestyle='--')
axs[0, 2].set_title('Average Charge Density')
axs[0, 2].set_xlabel('Position')
axs[0, 2].set_ylabel('Density')

# Bottom left: phi plot
axs[1, 0].plot(np.linspace(0, Nx*dx, Nx), phi, color='red')
axs[1, 0].set_title('Potential (phi)')
axs[1, 0].set_xlabel('Position')
axs[1, 0].set_ylabel('Potential')

# Bottom right: E_cells plot
axs[1, 1].plot(np.linspace(0, Nx*dx, Nx), E_cells, color='purple')
axs[1, 1].set_title('Electric Field of Cells')
axs[1, 1].set_xlabel('Position')
axs[1, 1].set_ylabel('Electric Field')

# Adjust layout
plt.tight_layout()
plt.show()


