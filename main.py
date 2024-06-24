import numpy as np





qp = 1
Nx = 20
ppc = 6
Np = ppc * Nx
dx = 0.1
part_pos = np.linspace(0, (Np - 1) * Nx * dx / Np, Np)