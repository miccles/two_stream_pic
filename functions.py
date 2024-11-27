import numpy as np
import scipy.sparse as sp

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
    for i in range(len(cell_pos)):
        for p in range(len(part_pos)):
                p_pos = part_pos[p]
                if i == 0 and part_pos[p] >= cell_pos[Nx - 1]:
                    dens_avg[i] += qp * b1((cell_pos[i] - (p_pos - Lx)) / dx)
                if i == Nx - 1 and part_pos[p] <= cell_pos[0]:
                    dens_avg[i] += qp * b1((cell_pos[i] - (p_pos + Lx)) / dx)
                dens_avg[i] += qp * b1((cell_pos[i] - p_pos) / dx)
    return dens_avg


# potential solver #
def phi_tridiag_solver(d): # Tridiagonal algorithm for non-periodic boundary conditions
    Nx = len(d)

    a = [1] * Nx
    a[0] = 0
    a[Nx - 1] = 0

    b = [-2] * Nx
    b[0] = 1        # Dirichlet boundary conditions
    b[Nx - 1] = 1   # Dirichlet boundary conditions

    c = [1] * Nx
    c[0] = 0
    c[Nx - 1] = 0

    for i in range(0, Nx):
        c[i] = c[i] / (b[i] - a[i] * c[i - 1])
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1])
    x = [0] * Nx
    x[-1] = d[-1]
    for i in range(Nx - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x



def laplacian_matrix(Nx):
    e = np.ones(Nx)
    diags = np.array([-1,0,1])
    vals  = np.vstack((e,-2*e,e))
    Lmtx = sp.spdiags(vals, diags, Nx, Nx)
    Lmtx = sp.lil_matrix(Lmtx) # tranform mtx type to modify entries
    Lmtx[0,Nx-1] = 1 # periodic boundary conditions
    Lmtx[Nx-1,0] = 1 # periodic boundary conditions
    Lmtx = sp.csr_matrix(Lmtx) # transform mtx type
    return Lmtx


def gradient_matrix(Nx):
    e = np.ones(Nx)
    diags = np.array([-1,1])
    vals  = np.vstack((-e,e))
    Gmtx = sp.spdiags(vals, diags, Nx, Nx)
    Gmtx = sp.lil_matrix(Gmtx) # tranform mtx type to modify entries
    Gmtx[0,Nx-1] = -1 # periodic boundary conditions
    Gmtx[Nx-1,0] = 1 # periodic boundary conditions
    Gmtx = sp.csr_matrix(Gmtx) # transform mtx type
    return Gmtx


# electric field solver #
def el_cell(phi, dx):
    return np.gradient(phi, dx)








