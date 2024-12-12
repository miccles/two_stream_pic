import os
import numpy as np
import matplotlib.pyplot as plt
from functions import *
from parameters import *

def main():
    # Generate initial conditions
    pos, vel = generate_init_cond(Lx, Np, beam_v0, beam_dv0, beam_perp)
    # Calculate initial acceleration
    dens, phi, el, acc = find_acc(pos, Nx, dx, q, m, n0)
    
    # Create a 4-pane plot
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    
    # Top left: pos, vel scatter
    scatter1 = axs[0, 0].scatter(pos[:int(Np/2)], vel[:int(Np/2)], s=2, color='blue', alpha=0.5)
    scatter2 = axs[0, 0].scatter(pos[int(Np/2):], vel[int(Np/2):], s=2, color='red', alpha=0.5)
    axs[0, 0].set_title('Phase Diagram')
    axs[0, 0].set_xlabel('Position')
    axs[0, 0].set_ylabel('Velocity')
    
    # Top right: dens_avg plot
    dens_line, = axs[0, 1].plot(np.linspace(0, Nx*dx, Nx), dens - n0 * q, color='green')
    dens_avg_line = axs[0, 1].axhline(y=np.mean(dens - q * n0), color='red', linestyle='--')
    axs[0, 1].set_title('Average Charge Density')
    axs[0, 1].set_xlabel('Position')
    axs[0, 1].set_ylabel('Density')
    
    # Bottom left: phi plot
    phi_line, = axs[1, 0].plot(np.linspace(0, Nx*dx, Nx), phi, color='red')
    axs[1, 0].set_title('Potential (phi)')
    axs[1, 0].set_xlabel('Position')
    axs[1, 0].set_ylabel('Potential')
    
    # Bottom right: E_cells plot
    E_cells_line, = axs[1, 1].plot(np.linspace(0, Nx*dx, Nx), el, color='purple')
    axs[1, 1].set_title('Electric Field of Cells')
    axs[1, 1].set_xlabel('Position')
    axs[1, 1].set_ylabel('Electric Field')
    
    plt.tight_layout()
    
    if real_time:
        plt.ion()
        plt.show()

    # Directory to save plots if real_time is False
    if not real_time:
        if not os.path.exists('plots'):
            os.makedirs('plots')

    # Simulation loop
    for i in range(timesteps):
        pos, vel, dens, phi, el, acc = leapfrog(pos, vel, acc, dt, n0)
        pop1_pos, pop1_vel = pos[:int(Np/2)], vel[:int(Np/2)]
        pop2_pos, pop2_vel = pos[int(Np/2):], vel[int(Np/2):]

        # Update scatter plot data
        scatter1.set_offsets(np.c_[pop1_pos, pop1_vel])
        scatter2.set_offsets(np.c_[pop2_pos, pop2_vel])
        
        # Update dens_avg plot data
        dens_line.set_ydata(dens - n0 * q)
        dens_avg_line.set_ydata([np.mean(dens - q * n0)] * 2)
        
        # Update phi plot data
        phi_line.set_ydata(phi)
        
        # Update E_cells plot data
        E_cells_line.set_ydata(el)
        
        if real_time:
            plt.pause(0.001)
        else:
            plt.savefig(f'plots/frame_{i:04d}.png')

    if real_time:
        plt.ioff()
        plt.savefig('pic.png', dpi=240)
        plt.show()
    else:
        # Create a movie from the saved frames
        os.system('ffmpeg -r 10 -i plots/frame_%04d.png -vcodec libx264 -y movie.mp4')

if __name__ == "__main__":
    main()