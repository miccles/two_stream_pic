import numpy as np
import math
import random as rnd
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve

def phifunc(d, ng):
    tempd = [de for de in d]
    #tempd[0], tempd[ng - 1] = 0, 0 
    a = [1] * ng
    b = [-2] * ng
    c = [1] * ng
    a[0] = 0
    b[0] = 1
    b[ng - 1] = 1
    c[ng - 1] = 0
    c[0] = c[0] / b[0]
    tempd[0] = tempd[0] / b[0]
    for i in range(1, ng):
        c[i] = c[i] / (b[i] - a[i] * c[i - 1])
        tempd[i] = (tempd[i] - a[i] * tempd[i - 1]) / (b[i] - a[i] * c[i - 1])
    x = [0] * ng
    x[ng - 1] = tempd[ng - 1]
    for i in range(ng - 2, -1, -1):
        x[i] = tempd[i] - c[i] * x[i + 1]
    return x


def efavg(phi, ng, dx):
    el = [0] * ng
    el[0] = -(phi[1] - phi[-1]) / (1 * dx)
    el[ng - 1] = -(phi[0] - phi[ng - 2]) / (1 * dx)
    for i in range(1, ng - 1):
        el[i] = -(phi[i + 1] - phi[i - 1]) / (2 * dx)
    return el

ng = 4
dx = 0.1
rho = [1] * ng
d = [-4 * np.pi * r * dx ** 2 for r in rho]
phi = phifunc(d, ng)
el = efavg(phi, ng, dx)

print(phi)
# plt.plot([i for i in range(len(el))], phi)
# plt.show()
# plt.close()
# plt.plot([i for i in range(len(el))], el)
# plt.show()
