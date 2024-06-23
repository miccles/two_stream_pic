import numpy as np
import matplotlib.pyplot as plt
from parameters import *


# b0 spline function #
def b0(x):
    if np.abs(x) <= 1 / 2:
        return 1
    else:
        return 0

# b1 spline function #
def b1(x):
    if np.abs(x) <= 1:
        return -np.abs(x) + 1
    else:
        return 0
    
# cell center positions #
def grid_pos(i, dx):
    return (i + 1 / 2) * dx


# cell average charge density #
def rho_avg(dx, Nx, part_pos, qp):
    Lx = Nx * dx
    dens_avg = [0] * Nx
    cell_pos = [grid_pos(i, dx) for i in range(Nx)]
    print(cell_pos)
    for i in range(len(cell_pos)):
        for p in range(len(part_pos)):
                if i == 0 and part_pos[p] >= cell_pos[Nx - 1]:
                    print(f'Particle position: {part_pos[p]}. Added density to cell {i}')
                    dens_avg[i] += qp * b1((cell_pos[i] - (part_pos[p] - Lx)) / dx)
                if i == Nx - 1 and part_pos[p] <= cell_pos[0]:
                    print(f'Particle position: {part_pos[p]}. Added density to cell {i}')
                    dens_avg[i] += qp * b1((cell_pos[i] - (part_pos[p] + Lx)) / dx)
                dens_avg[i] += qp * b1((cell_pos[i] - part_pos[p]) / dx)
    return dens_avg

#check particle pos that are exactly on the cell boundary


qp = 1
Nx = 20
ppc = 6
Np = ppc * Nx
dx = 0.1
part_pos = np.linspace(0.001, 0.999*Nx * dx, Np)
print(part_pos)

dens_avg = rho_avg(dx, Nx, part_pos, qp)
plt.plot(np.linspace(0, Nx * dx, Nx), dens_avg)
plt.show()
        