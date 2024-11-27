import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from functions import *

Nx = 10
dx = 0.1

# Construct matrix Lmtx to computer Laplacian (2nd derivative) 
# on mesh size Nx, separation dx
e = np.ones(Nx)
diags = np.array([-1,0,1])
vals  = np.vstack((e,-2*e,e))
Lmtx = sp.spdiags(vals, diags, Nx, Nx)
Lmtx = sp.lil_matrix(Lmtx) # tranform mtx type to modify entries
Lmtx[0,Nx-1] = 1
Lmtx[Nx-1,0] = 1
# Lmtx /= dx**2
Lmtx = sp.csr_matrix(Lmtx) # transform mtx type



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


Nx = 100
dx = 0.1
L = laplacian_matrix(Nx)
G = gradient_matrix(Nx)
rho = [1] * Nx
d = -np.random.rand(Nx)  # Example known vector d

# Solve for phi using sparse matrix solver
phi = spla.spsolve(L, d)

# Calculate the electric field using the gradient matrix G
electric_field = -G.dot(phi) / (2 * dx)


# Create subplots
fig, axs = plt.subplots(1, 2, figsize=(10, 8))

# Plot phi
axs[0].plot([i for i in range(len(phi))], phi)
axs[0].set_title('Potential (phi)')
axs[0].set_xlabel('Index')
axs[0].set_ylabel('Potential')
axs[0].set_yscale('log')

# Plot electric field
axs[1].plot([i for i in range(len(electric_field))], electric_field)
axs[1].set_title('Electric Field')
axs[1].set_xlabel('Index')
axs[1].set_ylabel('Electric Field')

# Adjust layout
plt.tight_layout()
plt.show()