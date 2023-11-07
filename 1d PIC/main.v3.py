from functions import *
import numpy as np
import matplotlib.pyplot as plt
import math



# First attempt on a 1D PIC simulation for two-stream instability


dt = 0.05                           # time step
Nt = 0                              # number of time steps (set to -1 for no loops)
tfinal = Nt * dt                    # total time
tex = -1                            # step interval for data export

ng = 64                             # number of grid boxes
L = 2 * np.pi                       # domain length in dx units
dx = L / ng                         # dimensionless size of grid boxes

ppc = 10                            # number of particles/cell
N = ng * ppc                        # number of computational particles
n = N / L                           # number density of particles

#Quantities in dimaneionless units (see notes)#
qp = -1 / (ppc * ng)                # charge of computational particle normalized to the total charge of the domain
mp = 1 / (ppc * ng)                 # mass of computational particle normalized to the total mass of the domain
qm = qp / mp                        # charge, mass ratio
wp = np.sqrt(qp ** 2 * n / mp)      # plasma frequency in dimensionless units

vd = 0.2                            # drift velocity
v0 = 0.1 * vd                       # random component of velocity
                                    # final velocities of particles will be centered around +- (v0+vd)

# perturbed quantities
vp1 = 0.0                           # perturbation of velocity
xp1 = 1.0 * L / N                   # perturbation of position
mode = 1                            # perturbation mode

# initial particle loading #
xp = np.linspace(0, L - L / N, N)  # - L/N because of BCs (at x=L we have the x=0 particle)
xi = xgrid(ng, L)
vp = v0 * np.random.normal(1, 1, N)
vp = [(vp[i] + vd) * (-1) ** i for i in range(len(vp))]


# Boundary Conditions of field #
phiL = 0
phiR = 0

# Initial Rho and d matrices (Before the perturbation) #
rho = rhoavg(xp, xi, ng, L, qp)
d = [-r for r in rho]
phi = phifunc(d, ng)                    # Cell averaged scalar potential
eli = efavg(phi, ng, dx)                # Electric field at the center of cells
elp = efparticle(xp, xi, dx, eli, N)       # Electric field of particles

# perturbation#
vp = [vp[i] + vp1 * np.sin(2 * np.pi * xp[i] * mode / L) for i in range(len(xp))]
xp = [xp[i] + xp1 * np.sin(2 * np.pi * xp[i] * mode / L) for i in range(len(xp))]

fig = plt.figure(figsize=(5, 4), dpi=80)
plt.cla()
plt.scatter([xp[2*i] for i in range(int(N/2))], [vp[2*i] for i in range(int(N/2))])
plt.scatter([xp[2*i + 1] for i in range(int(N/2))], [vp[2*i + 1] for i in range(int(N/2))])
plt.axis([0, L, -3, 3])
plt.pause(0.001)


# Calculation of rhoavg, phiavg, Ei and Ep after perturbation#

t = 0
for i in range(Nt):
    vp = vpup(vp, elp, dt / 2., qm)  # Half push of vp
    xp = bcp(xpup(xp, vp, dt), L)    # Push of xp + BCs
    rho = rhoavg(xp, xi, ng, L, qp)
    d = [- r for r in rho]
    phi = phifunc(d, ng)
    eli = efavg(phi, ng, dx)
    elp = efparticle(xp, xi, dx, eli, N)
    vp = vpup(vp, elp, dt / 2., qm)  # Half push of vp
    t += dt
    # Plotting #
    plt.cla()
    plt.scatter([xp[2*i] for i in range(int(N/2))], [vp[2*i] for i in range(int(N/2))])
    plt.scatter([xp[2*i + 1] for i in range(int(N/2))], [vp[2*i + 1] for i in range(int(N/2))])
    plt.axis([0, L, -3, 3])
    plt.pause(0.001)
plt.xlabel('x')
plt.ylabel('v')
plt.show()

# # Main Computational Cycle #
# if Nt >= 0:
#     # Evolving v to dt/2 once#
#     vp = vpup(vp, dt / 2, qm, elp)
#     for i in range(Nt + 1):
#         # Position Update & BCs#
#         xp = bcp(xpup(xp, vp, dt), L)
#         # Rho and d vectors
#         rho = rhoavg(xp, xi, qp, dx)
#         d = [-4 * np.pi * r for r in rho]
#         # Phi vector#
#         phi = phifunc(d, ng)
#         # Ei vector #
#         eli = efavg(phi, ng, dx)  # Electric field at the center of cells
#         # Ep vector #
#         elp = efparticle(xp, xi, dx, eli)
#         #Velocity update#
#         vp = vpup(vp, dt, qm, elp)
#         #plt.scatter([xp[2 * i] for i in range(int(N / 2))], [vp[2 * i] for i in range(int(N / 2))])
#         #plt.scatter([xp[2 * i + 1] for i in range(int(N / 2))], [vp[2 * i + 1] for i in range(int(N / 2))])
#         #plt.show()
#     # Final evolution of vp by -dt/2 once#
#     vp = vpup(vp, -dt / 2, qm, elp)

# #Plotting#
# plt.scatter([xp[2*i] for i in range(int(N/2))], [vp[2*i] for i in range(int(N/2))])
# plt.scatter([xp[2*i + 1] for i in range(int(N/2))], [vp[2*i + 1] for i in range(int(N/2))])
# plt.show()

# plt.scatter(xi, phi)
# plt.show()


# def rhoavg(xp, xg, q, dx):
#     global L, ng
#     rho = [0] * (ng + 1)
#     for i in range(ng):
#         for p in range(len(xp)):
#             rho[i] += (q / dx) * b1((xg[i] - xp[p]) / dx)
#             if xg[ng - 1] <= xp[p] <= L and i == 0:
#                 rho[i] += (q / dx) * b1((xg[i] - (xp[p] - L)) / dx)
#             elif 0 <= xp[p] <= xg[0] and i == ng - 1:
#                 rho[i] += (q / dx) * b1((xg[i] - (xp[p] + L)) / dx)
#     rho[ng] = rho[0]
#     return rho


# def phifunc(d, ng):
#     tempd = [de for de in d]
#     a = [1] * (ng + 1)
#     b = [-2] * (ng + 1)
#     c = [1] * (ng + 1)
#     b[0] = 1
#     b[ng] = 1
#     a[0] = 0
#     a[ng] = 0
#     c[0] = 0
#     c[ng] = 0
#     c[0] = c[0] / b[0]
#     tempd[0] = tempd[0] / b[0]
#     for i in range(1, ng):
#         c[i] = c[i] / (b[i] - a[i] * c[i - 1])
#         tempd[i] = (tempd[i] - a[i] * tempd[i - 1]) / (b[i] - a[i] * c[i - 1])
#     tempd[ng] = (tempd[ng] - a[ng] * tempd[ng - 1]) / (b[ng] - a[ng] * c[ng - 1])
#     x = [0] * (ng + 1)
#     x[ng] = tempd[ng]
#     for i in range(ng - 1, -1, -1):
#         x[i] = tempd[i] - c[i] * x[i + 1]
#     return x


# def efavg(phi, ng, dx):
#     el = [0] * (ng + 1)
#     el[0] = -(phi[1] - phi[-1]) / (2 * dx)
#     el[ng] = -(phi[0] - phi[ng - 1]) / (2 * dx)
#     for i in range(1, ng):
#         el[i] = -(phi[i + 1] - phi[i - 1]) / (2 * dx)
#     return el


# def efparticle(xp, xg, dx, elgrid):
#     global L, N, ng
#     efp = [0] * N
#     for p in range(N):
#         try:
#             for i in range(int(np.ceil(xp[p] / dx - 1)), int(np.floor(xp[p] / dx + 1) + 1)):
#                 efp[p] += elgrid[i] * b1((xg[i] - xp[p]) / dx)
#         except IndexError:
#             print("Index Error for i=", i)
#     return efp


# def xpup(xp, vp, dt):
#     newpos = [xp[i] + vp[i] * dt for i in range(len(xp))]
#     return newpos


# def vpup(vp, dt, qm, Ep):
#     newvel = [vp[i] + qm * dt * Ep[i] for i in range(len(vp))]
#     return newvel


# def bcp(xp, L):
#     for i in range(len(xp)):
#         if xp[i] < 0 or xp[i] > L:
#             xp[i] -= (xp[i] // L) * L
#     return xp


# def byenoise(list):
#     for i in range(len(list)):
#         if np.abs(list[i]) < 10 ** (-10):
#             list[i] = 0
#     return list


# def plotfunc(x, y, type):
#     if type == "scatter":
#         plt.scatter(x, y)
#     elif type == "plot":
#         plt.plot(x, y)
#     plt.show()


# def get_frame(xp, vp, phi, elp, t):
#     global dt, N, xi
#     xphi = pd.DataFrame({'x': xi, 'y': phi})
#     xep = pd.DataFrame({'x' : xp, 'y': elp})
#     xphi.sort_values('x', inplace=True)
#     xep.sort_values('x', inplace=True)
#     fig = plt.figure(figsize=(6, 6))
#     plt.scatter([xp[2 * i] for i in range(int(N / 2))], [vp[2 * i] for i in range(int(N / 2))])
#     plt.scatter([xp[2 * i + 1] for i in range(int(N / 2))], [vp[2 * i + 1] for i in range(int(N / 2))])
#     plt.xlim([0, 2 * np.pi])
#     plt.xlabel('x', fontsize=14)
#     plt.ylim([-1.2 * max(vp), 1.2 * max(vp)])
#     plt.ylabel('v', fontsize=14)
#     plt.title(f'Phase Space at t = {t * dt}', fontsize=14)
#     plt.savefig(save_results_image + f'/phase_{t}.png', dpi=300, transparent=False, facecolor='white')
#     plt.close()
#     plt.plot(xphi['x'], xphi['y'])
#     plt.xlabel('x')
#     plt.ylabel(r'$\phi$')
#     plt.title(f'Electric potential at t = {t * dt}', fontsize=14)
#     plt.savefig(save_results_image + f'/phi_{t}.png', dpi=300, transparent=False, facecolor='white')
#     plt.close()
#     plt.plot(xep['x'], xep['y'])
#     plt.xlabel('x')
#     plt.ylabel(r'$E_{p}$')
#     plt.title(f'Electric field at t = {t * dt}', fontsize=14)
#     plt.savefig(save_results_image + f'/ep_{t}.png', dpi=300, transparent=False, facecolor='white')
#     plt.close()


# def velhist(vp, t):
#     plt.hist(vp, bins=15)
#     plt.xlabel("v")
#     plt.title(f'Velocity, t = {t * dt}', fontsize=14)
#     plt.savefig(save_results_image + f'/vhist_{t}.png', dpi=300, transparent=False, facecolor='white')
#     plt.close()


# def kinetic(vp):
#     global mp
#     kin = 0
#     for i in range(len(vp)):
#         kin += 0.5 * mp * vp[i] * vp[i]
#     return kin


# def potential(xp):
#     pot = 0
#     global qp, N
#     for i in range(N):
#         for j in range(i + 1, N):
#             pot += qp * qp / np.abs(xp[i] - xp[j])
#     return pot


# dt = 0.005  # time step
# Nt = 40  # number of time steps (set to -1 for no loops)
# tfinal = Nt * dt  # total time
# tex = 10  # time step interval for data export

# ng = 64  # number of grid boxes
# L = 2 * np.pi  # domain length in dx units
# dx = L / ng  # dimensionless size of grid boxes
# xi = xgrid(ng, L)

# Np = 20  # number of particles/cell
# N = ng * Np  # number of computational particles
# n0 = N / L  # initial (and average) number density of particles in units of N/L.

# # Quantities in dimaneionless units (see notes)#
# qp = -1  # charge of computational particle
# mp = 1  # mass of computational particle
# qm = qp / mp
# wp = np.sqrt(qp ** 2 * n0 / mp)  # plasma frequency in dimensionless units


# # Stability Check#
# cond1 = wp * dt < 2
# cond2 = dt < dx
# if not cond1 or not cond2:
#     print("Warning: Stability issues")

# vd = L / (ng * dt)  # drift velocity
# v0 = 0.1 * vd  # random component of velocity
# # final velocities of particles will be centered around +- (v0+vd)

# # perturbed quantities
# vp1 = 0.0 * vd  # perturbation of velocity
# xp1 = 1.0 * L / N  # perturbation of position
# mode = 1  # perturbation mode

# # initial particle loading#
# xp = np.linspace(0, L - L / N, N)  # - L/N because of BCs (at x=L we have the x=0 particle)
# vp = v0 * np.random.normal(1, 1, N)
# vp = [(vp[i] + vd) * (-1) ** i for i in range(len(vp))]


# # Boundary Conditions of field#
# phiL = 0
# phiR = 0


# # perturbation#
# vp = [vp[i] + vp1 * math.sin(2 * np.pi * xp[i] * mode / L) for i in range(len(xp))]
# xp = bcp([xp[i] + xp1 * math.sin(2 * np.pi * xp[i] * mode / L) for i in range(len(xp))], L)

# kinp = [0] * (Nt + 1)
# potp = [0] * (Nt + 1)

# # Rho and d matrices (After the perturbation)#
# rho = rhoavg(xp, xi, qp, dx)
# d = [-(r - n0 * qp) for r in rho]
# # Phi vector#
# phi = phifunc(d, ng)
# # Ei vector #
# eli = efavg(phi, ng, dx)  # Electric field at the center of cells
# # Ep vector #
# elp = efparticle(xp, xi, dx, eli)


# # Calculation of rhoavg, phiavg, Ei and Ep after perturbation#

# # Main Computational Cycle #
# tcount = 0
# kinp[0] = kinetic(vp)
# potp[0] = potential(xp)
# get_frame(xp, vp, phi, elp, tcount)
# velhist(vp, tcount)
# if Nt >= 0:
#     for i in range(Nt):
#         tcount += 1
#         # Evolving v to dt/2 #
#         vp = vpup(vp, dt / 2, qm, elp)
#         # Position Update & BCs#
#         xp = bcp(xpup(xp, vp, dt), L)
#         # Rho and d vectors
#         rho = rhoavg(xp, xi, qp, dx)
#         d = [-(r - n0 * qp) for r in rho]
#         # Phi vector#
#         phi = phifunc(d, ng)
#         # Ei vector #
#         eli = efavg(phi, ng, dx)  # Electric field at the center of cells
#         # Ep vector #
#         elp = efparticle(xp, xi, dx, eli)
#         # Velocity update#
#         vp = vpup(vp, dt / 2, qm, elp)
#         kinp[i + 1] = kinetic(vp)
#         potp[i + 1] = potential(xp)
#         get_frame(xp, vp, phi, elp, tcount)
#         velhist(vp, tcount)
#         print("Timestep", tcount,  "of", Nt)


# frames_phase = []
# frames_phi = []
# frames_ep = []
# frames_hist = []
# for t in range(Nt + 1):
#     image_phase = imageio.v2.imread(save_results_image + f'/phase_{t}.png')
#     image_phi = imageio.v2.imread(save_results_image + f'/phi_{t}.png')
#     image_ep = imageio.v2.imread(save_results_image + f'/ep_{t}.png')
#     image_hist = imageio.v2.imread(save_results_image + f'/vhist_{t}.png')
#     frames_phase.append(image_phase)
#     frames_phi.append(image_phi)
#     frames_ep.append(image_ep)
#     frames_hist.append(image_hist)
# imageio.mimsave(save_results_gif + '/phase.gif', frames_phase, fps=2)
# imageio.mimsave(save_results_gif + '/phi.gif', frames_phi, fps=2)
# imageio.mimsave(save_results_gif + '/ep.gif', frames_ep, fps=2)
# imageio.mimsave(save_results_gif + '/vhist.gif', frames_hist, fps=2)

# plt.plot([i * dt for i in range(Nt + 1)], [kinp[i]/kinp[0] for i in range(len(kinp))], label='K/K0')
# plt.plot([i * dt for i in range(Nt + 1)], [potp[i]/potp[0] for i in range(len(potp))], label='U/U0')
# plt.plot([i * dt for i in range(Nt + 1)], [(kinp[i] + potp[i])/(kinp[0]) for i in range(len(kinp))], label='K+U/K0')
# plt.xlabel('t')
# plt.ylabel('Normalized Energies')
# plt.legend(loc='best')
# plt.savefig(save_results_image + '/energy.png', dpi=300, transparent=False, facecolor='white')
# # for t in range(Nt + 1):
# #     os.remove(save_results_image + f'/img_{t}.png')