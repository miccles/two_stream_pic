import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from functions import *
import timeit


N = 1000
boxsize = 50
vb = 1
vth = 0.1
A = 0.01


pos, vel = generate_init_cond(boxsize, Np, vb, vth, A)


plt.scatter(pos, vel)
plt.show()

