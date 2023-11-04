import numpy as np
import math
import pandas as pd
import os
import matplotlib.pyplot as plt
import imageio
import random as rnd
import scipy.sparse as sp
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve



pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 500)


my_path = os.path.abspath(__file__)  # Figures out the absolute path for you in case your working directory moves around.
save_results_image = 'C:/Users/micha/Desktop/PIC_Plots/images'
save_results_gif = 'C:/Users/micha/Desktop/PIC_Plots/gifs'


# First attempt on a 1D PIC simulation for two-stream instability


def b0(x):
    if np.abs(x) <= 1 / 2:
        return 1
    else:
        return 0


def b1(x):
    if -1 <= x <= 1:
        return -np.abs(x) + 1
    else:
        return 0


def xgrid(ng, L):
    dx = L / ng
    x = [i * dx for i in range(ng + 1)]
    return x


def rhoavg(xp, xg, q, dx):
    global L, ng
    rho = [0] * (ng + 1)
    for i in range(ng):
        for p in range(len(xp)):
            rho[i] += (q / dx) * b1((xg[i] - xp[p]) / dx)
            if xg[ng - 1] <= xp[p] <= L and i == 0:
                rho[i] += (q / dx) * b1((xg[i] - (xp[p] - L)) / dx)
            elif 0 <= xp[p] <= xg[0] and i == ng - 1:
                rho[i] += (q / dx) * b1((xg[i] - (xp[p] + L)) / dx)
    rho[ng] = rho[0]
    return rho


def phifunc(d, ng):
    tempd = [de for de in d]
    a = [1] * (ng + 1)
    b = [-2] * (ng + 1)
    c = [1] * (ng + 1)
    b[0] = 1
    b[ng] = 1
    a[0] = 0
    a[ng] = 0
    c[0] = 0
    c[ng] = 0
    c[0] = c[0] / b[0]
    tempd[0] = tempd[0] / b[0]
    for i in range(1, ng):
        c[i] = c[i] / (b[i] - a[i] * c[i - 1])
        tempd[i] = (tempd[i] - a[i] * tempd[i - 1]) / (b[i] - a[i] * c[i - 1])
    tempd[ng] = (tempd[ng] - a[ng] * tempd[ng - 1]) / (b[ng] - a[ng] * c[ng - 1])
    x = [0] * (ng + 1)
    x[ng] = tempd[ng]
    for i in range(ng - 1, -1, -1):
        x[i] = tempd[i] - c[i] * x[i + 1]
    return x


def efavg(phi, ng, dx):
    el = [0] * (ng + 1)
    el[0] = -(phi[1] - phi[-1]) / (2 * dx)
    el[ng] = -(phi[0] - phi[ng - 1]) / (2 * dx)
    for i in range(1, ng):
        el[i] = -(phi[i + 1] - phi[i - 1]) / (2 * dx)
    return el


def efparticle(xp, xg, dx, elgrid):
    global L, N, ng
    efp = [0] * N
    for p in range(N):
        try:
            for i in range(int(np.ceil(xp[p] / dx - 1)), int(np.floor(xp[p] / dx + 1) + 1)):
                efp[p] += elgrid[i] * b1((xg[i] - xp[p]) / dx)
        except IndexError:
            print("Index Error for i=", i)
    return efp


def xpup(xp, vp, dt):
    newpos = [xp[i] + vp[i] * dt for i in range(len(xp))]
    return newpos


def vpup(vp, dt, qm, Ep):
    newvel = [vp[i] + qm * dt * Ep[i] for i in range(len(vp))]
    return newvel


def bcp(xp, L):
    for i in range(len(xp)):
        if xp[i] < 0 or xp[i] > L:
            xp[i] -= (xp[i] // L) * L
    return xp


def byenoise(list):
    for i in range(len(list)):
        if np.abs(list[i]) < 10 ** (-10):
            list[i] = 0
    return list


def plotfunc(x, y, type):
    if type == "scatter":
        plt.scatter(x, y)
    elif type == "plot":
        plt.plot(x, y)
    plt.show()


def get_frame(xp, vp, phi, elp, t):
    global dt, N, xi
    xphi = pd.DataFrame({'x': xi, 'y': phi})
    xep = pd.DataFrame({'x' : xp, 'y': elp})
    xphi.sort_values('x', inplace=True)
    xep.sort_values('x', inplace=True)
    fig = plt.figure(figsize=(6, 6))
    plt.scatter([xp[2 * i] for i in range(int(N / 2))], [vp[2 * i] for i in range(int(N / 2))])
    plt.scatter([xp[2 * i + 1] for i in range(int(N / 2))], [vp[2 * i + 1] for i in range(int(N / 2))])
    plt.xlim([0, 2 * np.pi])
    plt.xlabel('x', fontsize=14)
    plt.ylim([-1.2 * max(vp), 1.2 * max(vp)])
    plt.ylabel('v', fontsize=14)
    plt.title(f'Phase Space at t = {t * dt}', fontsize=14)
    plt.savefig(save_results_image + f'/phase_{t}.png', dpi=300, transparent=False, facecolor='white')
    plt.close()
    plt.plot(xphi['x'], xphi['y'])
    plt.xlabel('x')
    plt.ylabel(r'$\phi$')
    plt.title(f'Electric potential at t = {t * dt}', fontsize=14)
    plt.savefig(save_results_image + f'/phi_{t}.png', dpi=300, transparent=False, facecolor='white')
    plt.close()
    plt.plot(xep['x'], xep['y'])
    plt.xlabel('x')
    plt.ylabel(r'$E_{p}$')
    plt.title(f'Electric field at t = {t * dt}', fontsize=14)
    plt.savefig(save_results_image + f'/ep_{t}.png', dpi=300, transparent=False, facecolor='white')
    plt.close()


def velhist(vp, t):
    plt.hist(vp, bins=15)
    plt.xlabel("v")
    plt.title(f'Velocity, t = {t * dt}', fontsize=14)
    plt.savefig(save_results_image + f'/vhist_{t}.png', dpi=300, transparent=False, facecolor='white')
    plt.close()


def kinetic(vp):
    global mp
    kin = 0
    for i in range(len(vp)):
        kin += 0.5 * mp * vp[i] * vp[i]
    return kin


def potential(xp):
    pot = 0
    global qp, N
    for i in range(N):
        for j in range(i + 1, N):
            pot += qp * qp / np.abs(xp[i] - xp[j])
    return pot

