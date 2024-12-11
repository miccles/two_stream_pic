import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

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
                p_pos = part_pos[p]
                if i == 0 and part_pos[p] >= cell_pos[Nx - 1]:
                    dens_avg[i] += qp * b1((cell_pos[i] - (p_pos - Lx)) / dx)
                if i == Nx - 1 and part_pos[p] <= cell_pos[0]:
                    dens_avg[i] += qp * b1((cell_pos[i] - (p_pos + Lx)) / dx)
                dens_avg[i] += qp * b1((cell_pos[i] - p_pos) / dx)
    return dens_avg


def laplacian_matrix(Nx):
    e = np.ones(Nx)
    diags = np.array([-1,0,1])
    vals  = np.vstack((e,-2*e,e))
    Lap = sp.spdiags(vals, diags, Nx, Nx)
    Lap = sp.lil_matrix(Lap) # tranform mtx type to modify entries
    Lap[0,Nx-1] = 1 # periodic boundary conditions
    Lap[Nx-1,0] = 1 # periodic boundary conditions
    Lap = sp.csr_matrix(Lap) # transform mtx type
    return Lap


def gradient_matrix(Nx):
    e = np.ones(Nx)
    diags = np.array([-1,1])
    vals  = np.vstack((-e,e))
    Grad = sp.spdiags(vals, diags, Nx, Nx)
    Grad = sp.lil_matrix(Grad) # tranform mtx type to modify entries
    Grad[0,Nx-1] = -1 # periodic boundary conditions
    Grad[Nx-1,0] = 1 # periodic boundary conditions
    Grad = sp.csr_matrix(Grad) # transform mtx type
    return Grad


# potential solver - Tridiagonal algorithm for non-periodic boundary conditions (phiL, phiR are known) #
def phi_tridiag_solver(d):
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


# potential solver - sparse matrix solver for periodic boundary conditions #
def phi_sparse_solver(L, d):
    phi = spla.spsolve(L, d)
    return phi


# electric field solver - using the gradient matrix G #
def el_solver(G, phi, dx):
    electric_field = -G.dot(phi) / (2 * dx)
    return electric_field


# apply periodic boundary conditions #
def periodic_bc(pos, Lx):
    pos = np.mod(pos, Lx)
    return pos


# generate initial particle positions and velocities #
def generate_init_cond(Nx, Np, v0, dv0, A):
    pos_list = np.random.rand(Np, 1) * Nx # random positions (uniform distribution)
    vel_list = dv0 * np.random.randn(Np, 1) + v0 # random velocities (normal distribution)
    vel_list[int(Np / 2):] *= -1 # half of the particles have negative velocity
    vel_list *= (1 + A * np.sin(2 * np.pi * pos_list / Nx)) # add perturbation
    return pos_list, vel_list










