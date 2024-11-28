import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from functions import *



xvalues = np.linspace(0, 1, 3)

phi_values = xvalues * (1 - xvalues)

random_matrix = np.random.randint(0, 11, size=(3, 3))

print(random_matrix)

print(phi_values)

print(random_matrix.dot(phi_values))
print(np.dot(random_matrix, phi_values))
print(random_matrix @ phi_values)