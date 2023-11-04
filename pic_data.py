import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy as sp
from scipy.signal import find_peaks


path = 'C:/Users/micha/Desktop/PIC_Plots/'
yslice = 903
i = 101

maxindex = 0
plasmoid_y = []


hf = h5py.File(path + 'flds.tot' + '.%03d' % i, 'r')
densy = np.transpose(hf['dens'], (0, 1, 2))[0]
dens = densy[yslice]
densy_avg = sum(dens) / len(dens)
maxdens = max(dens)
if maxdens > 2.5 * densy_avg:
    maxindex = np.where(dens == maxdens)[0][0]
xnum = [i for i in range(densy.shape[1])]
plt.subplot(2, 1, 1)
plt.plot(xnum, dens, label=f'slice {yslice}')
plt.plot(xnum, [densy_avg] * len(dens), '--', label='average')
plt.xlabel('x slice')
plt.ylabel('dens')
plt.title('y-Slice ' + str(yslice) + ' of density as a function of x')
densx = np.transpose(hf['dens'], (0, 2, 1))[0]
dens = densx[maxindex]
densx_avg = sum(dens) / len(dens)
for j in range(len(dens)):
    if dens[j] > 2 * densx_avg:
        plasmoid_y.append(j)
ynum = [i for i in range(densx.shape[1])]
plt.subplot(2, 1, 2)
plt.plot(ynum, dens, label=f'slice {maxindex}')
plt.plot(ynum, [densx_avg] * len(dens), '--', label='average')
plt.scatter(plasmoid_y, [dens[i] for i in plasmoid_y], color='r', marker='*')
plt.axvline(x=yslice, linestyle='--', color='r')
plt.xlabel('y slice')
plt.ylabel('density')
plt.title('max density x-Slice ' + str(maxindex) + ' as a function of y')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('C:/Users/micha/Desktop/PIC_Plots/slice.png')
plt.close()
hf.close()

peaks, _ = find_peaks(dens, height=1.5*densx_avg, threshold=0.1*densx_avg)

plt.plot(ynum, dens, label=f'slice {maxindex}')
plt.scatter(peaks, dens[peaks], color = 'r', marker = '*')
#plt.scatter(plasmoid_y, [dens[i] for i in plasmoid_y], color='r', marker='*')
plt.show()






