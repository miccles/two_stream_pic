import numpy as np
from functions import *
from parameters import *



def main():
    Lap = laplacian_matrix(Nx) # Laplacian matrix
    G = gradient_matrix(Nx) # Gradient matrix

    # Generate initial conditions
    pos, vel = generate_init_cond(Lx, Np, beam_v0, beam_dv0, beam_perp)

    






if __name__ == "__main__":
    main()

