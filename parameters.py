import numpy as np


##### Simulation parameters #####

# Domain parameters #
Nx = 100      # Number of cells
Lx = Nx       # Length of the domain      
dx = Lx / Nx  # Size of each cell

# Particle parameters #
ppc0 = 10      # Number of particles per cell
Np = Nx * ppc0  # Number of computational particles

# Time resolution #
c = 1
eta = 0.5     # number < 1 that is used for dt
dt = eta * dx / c  # Simulation timestep
timesteps = 100  # Number of timesteps

# Spatial resolution #
kappa = 4     # number > eta/2 that is used for c_omp
c_omp = kappa * dx    # Plasma skin depth in dx units

# Plasma parameters #
q_m = 1       # Charge to mass ratio
n0 = Np / Lx  # Plasma number density
omega_p = c / c_omp    # Plasma frequency
q = omega_p ** 2 / (4 * np.pi * n0 * q_m)  # Particle charge
m = q / q_m                                # Particle mass

# Beam parameters #
beam_v0 = 0.5 * c  # Beam velocity
beam_dv0 = 0.2 * beam_v0  # Beam velocity spread (thermal component)
beam_perp = 0.001  # Beam perturbation amplitude


# Plotting parameters #
real_time = False  # Real-time plotting


if __name__ == "__main__":
    print('--- Simulation Parameters ---')
    print(f'Number of cells: {Nx}')
    print(f'Length of the domain: {Lx}')
    print(f'Size of each cell: {dx}')
    print(f'Number of particles per cell: {ppc0}')
    print(f'Number of computational particles: {Np}')
    print(f'Simulation timestep: {dt}')
    print(f'Number of timesteps: {timesteps}')
    print(f'Plasma skin depth in cell units: {c_omp}')
    print(f'Plasma number density: {n0}')
    print(f'Plasma frequency: {omega_p}')
    print(f'Particle charge: {q}')
    print(f'Particle mass: {m}')
    print(f'Beam velocity: {beam_v0}')
    print(f'Beam velocity spread: {beam_dv0}')
    print(f'Beam perturbation amplitude: {beam_perp}')
    print(f'Charge to mass ratio: {q_m}')
    print('----------------------------')









