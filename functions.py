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
                    dens_avg[i] += qp * b1((cell_pos[i] - (part_pos[p] - Lx)) / dx)
                if i == Nx - 1 and part_pos[p] <= cell_pos[0]:
                    dens_avg[i] += qp * b1((cell_pos[i] - (part_pos[p] + Lx)) / dx)
                dens_avg[i] += qp * b1((cell_pos[i] - part_pos[p]) / dx)
    return dens_avg