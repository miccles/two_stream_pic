import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from functions import *
import timeit





box_size = 10

pos = [-1, 2, 4, 8, 9, 11, 14, 18]

print(np.mod(pos, box_size))