import numpy as np
from functions import *
from parameters import *
import matplotlib.pyplot as plt



def main():
    Lap = laplacian_matrix(Nx) # Laplacian matrix
    G = gradient_matrix(Nx) # Gradient matrix

    # Generate initial conditions
    pos, vel = generate_init_cond(Lx, Np, beam_v0, beam_dv0, beam_perp)

    # Calculate initial acceleration
    acc = find_acc(pos, Nx, dx, q, m)

    fig = plt.figure(figsize=(5, 4), dpi=80)

    for i in range(timesteps):
        pos, vel, acc = leapfrog(pos, vel, acc, dt)
        pop1_pos, pop1_vel = pos[:int(Np/2)], vel[:int(Np/2)]
        pop2_pos, pop2_vel = pos[int(Np/2):], vel[int(Np/2):]
        plt.cla()
        plt.scatter(pop1_pos, pop1_vel, s=.4, color='blue', alpha=0.5)
        plt.scatter(pop2_pos, pop2_vel, s=.4, color='red', alpha=0.5)
        plt.pause(0.001)
    
    plt.show()






if __name__ == "__main__":
    main()

