import numpy as np
import matplotlib.pyplot as plt


from functions import *

plt.plot(np.linspace(-2, 2, 100), [b1(x) for x in np.linspace(-2, 2, 100)])
plt.show()




