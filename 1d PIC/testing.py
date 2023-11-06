import numpy as np
import math
import random as rnd
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve


def xgrid(ng, L):
    dx = L / ng
    x = [(i + 1 / 2) * dx for i in range(ng)]
    return x

L = 100
ng = 10

print(xgrid(ng, L))

