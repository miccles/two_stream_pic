import numpy as np


def b0(x):
    if np.abs(x) <= 1 / 2:
        return 1
    else:
        return 0


def b1(x):
    x = np.abs(x)
    if x <= 1:
        return -x + 1
    else:
        return 0

## fix, maybe i + 1 / 2 instead of i? For the cell centers
def xgrid(ng, L):
    dx = L / ng
    x = [i * dx for i in range(ng + 1)]
    return x

##fix
def rhoavg(xp, xg, q, dx):
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
    a = [1] * ng
    b = [-2] * ng
    c = [1] * ng
    a[0] = 0
    a[ng - 1] = 0
    b[0] = 1
    b[ng - 1] = 1
    c[0] = 0
    c[ng - 1] = 0
    c[0] = c[0] / b[0]
    tempd[0] = d[0] / b[0]
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

## fix
def efparticle(xp, xg, dx, elgrid):
    global L, N, ng
    efp = [0] * N
    for p in range(N):
        # if 0 <= xp[p] <= dx / 2:
        #     efp[p] = elgrid[0] * b1((xg[0] - xp[p]) / dx) #+ elgrid[-1] * b1((xg[-1] - xp[p]) / dx)
        # elif L - dx / 2 <= xp[p] <= L:
        #     efp[p] = elgrid[ng - 1] * b1((xg[ng - 1] - xp[p]) / dx)
        # else:
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
        if xp[i] < 0:
            xp[i] += L
        elif xp[i] > L:
            xp[i] -= L
    return xp