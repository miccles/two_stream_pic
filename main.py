import numpy as np
from functions import *
from parameters import *


Nx = 100
ppc = 10
Np = Nx * ppc
dx = 0.01
dx0 = 10 ** (-2)
q = 1

part_pos = np.random.uniform(0, Nx * dx, size=(Nx,))
part_pos = np.linspace(0, (Np - 1) * Nx * dx / Np, Np)
part_pos = np.random.normal(Nx * dx / 2, Nx * dx / 5, size=(Np,))
#part_pos = part_pos + dx0 * np.cos(2 * np.pi * part_pos / (Nx * dx))
dens_avg = rho_avg(dx, Nx, part_pos, q)
d = [-4 * np.pi * dx ** 2 * dens for dens in dens_avg]
phi = phi_tridiag_solver(d)

phi_theor = 2 * np.pi * ppc * q * np.linspace(0, Nx * dx, Nx) * (Nx * dx - np.linspace(0, Nx * dx, Nx))
plt.hist(part_pos, bins=Nx)
plt.show()


plt.plot([i for i in range(len(phi))], phi, label='Numerical solution')
plt.plot([i for i in range(len(phi))], phi_theor, label='Theoretical solution')
plt.legend(loc='best')
plt.show()

theta = np.linspace(0, 2*np.pi, len(phi))

# plt.polar(theta, dens_avg)
# plt.show()


# plt.polar(theta, phi, label='Numerical solution')
# plt.polar(theta, phi_theor, label='Theoretical solution')
# plt.legend(loc='best')
# plt.show()