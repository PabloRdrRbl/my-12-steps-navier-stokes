import numpy as np
import matplotlib.pyplot as plt


def linearconv(nx):
    dx = 2 / (nx - 1)  # Spatial domain has 2 units length
    nt = 25  # Number of time steps
    sigma = 0.2
    nu = 0.3  # Viscosity

    # Courant number, CFL condition
    # CFL < 1
    dt = sigma * dx**2 / nu

    # Initial conditions
    u = np.ones(nx)
    u[0.5 / dx: 1 / dx + 1] = 2  # u = 2 between 0.5 and 1

    # Printing u
    plt.plot(np.linspace(0, 2, nx), u)

    un = np.ones(nx)  # Initialize a temporary array

    # Solution
    for n in range(nt):
        un = u.copy()
        for i in range(1, nx - 1):
            u[i] = (un[i] + nu * dt / dx**2 *
                    (un[i + 1] - 2 * un[i] + un[i - 1]))

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
    linearconv(41)
