import numpy as np
import matplotlib.pyplot as plt
from functions import *
from parameters import *

def main():
    # Generate initial conditions
    pos, vel = generate_init_cond(Lx, Np, beam_v0, beam_dv0, beam_perp)

    # Calculate initial acceleration
    acc = find_acc(pos, Nx, dx, q, m)

    # Set up the plot
    fig, ax = plt.subplots(figsize=(8, 6))
    scatter1 = ax.scatter([], [], s=2, color='blue', alpha=0.5)
    scatter2 = ax.scatter([], [], s=2, color='red', alpha=0.5)
    ax.set_xlim(0, Lx)
    ax.set_ylim(-3 * beam_v0, 3 * beam_v0)
    ax.set_xlabel('x')
    ax.set_ylabel('v')
    plt.ion()
    plt.show()

    # Simulation loop
    for i in range(timesteps):
        pos, vel, acc = leapfrog(pos, vel, acc, dt)
        pop1_pos, pop1_vel = pos[:int(Np/2)], vel[:int(Np/2)]
        pop2_pos, pop2_vel = pos[int(Np/2):], vel[int(Np/2):]

        # Update scatter plot data
        scatter1.set_offsets(np.c_[pop1_pos, pop1_vel])
        scatter2.set_offsets(np.c_[pop2_pos, pop2_vel])
        plt.pause(0.001)

    # Save the final plot
    plt.ioff()
    plt.savefig('pic.png', dpi=240)
    plt.show()

if __name__ == "__main__":
    main()
