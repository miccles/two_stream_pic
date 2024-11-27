import numpy as np
import matplotlib.pyplot as plt


from functions import *



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



d = [-1] * 100


phi = phi_tridiag_solver(d)
electric_field = -np.gradient(phi)


fig, axs = plt.subplots(1, 2, figsize=(10, 5))
axs[0].plot(phi)
axs[1].plot(electric_field)
plt.show()