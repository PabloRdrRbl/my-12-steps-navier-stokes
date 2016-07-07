import numpy as np
import matplotlib.pyplot as plt


def linearconv(nx, sigma):
    dx = 2 / (nx - 1)  # Spatial domain has 2 units length
    nt = 25  # Number of time steps
    c = 1  # Wave speed

    # Courant number, CFL condition
    # CFL < 1
    dt = sigma * dx

    # Initial conditions
    u = np.ones(nx)
    u[0.5 / dx: 1 / dx + 1] = 2  # u = 2 between 0.5 and 1

    # Printing u
    plt.plot(np.linspace(0, 2, nx), u)

    un = np.empty(nx)  # Initialize a temporary array

    # Solution
    for n in range(nt):
        un = u.copy()
        for i in range(1, nx):
            u[i] = un[i] - c * dt / dx * (un[i] - un[i - 1])

    # Printing the new u
    plt.plot(np.linspace(0, 2, nx), u)

    plt.title('Using $nx = %d$ and $\sigma = %1.1f$' % (nx, sigma))

    plt.savefig('image_output/u_profile_sigma_%d_nx_%d.png' % (sigma * 10, nx))
    # plt.show()
    # Cleaning the figure
    # [http://stackoverflow.com/questions/8213522/matplotlib-
    # clearing-a-plot-when-to-use-cla-clf-or-close]
    plt.clf()

if __name__ == '__main__':
    for sigma in [0.5, 1, 1.5]:
        for nx in [41, 51, 61, 71, 81, 91, 341]:
            linearconv(nx, sigma)
