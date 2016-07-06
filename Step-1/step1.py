import numpy as np
import matplotlib.pyplot as plt


def linearconv(nx):
    dx = 2 / (nx - 1)  # Spatial domain has 2 units length
    nt = 25  # Number of time steps
    dt = 0.025  # Amount of time each time steps takes
    c = 1  # Wave speed

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

    plt.savefig('u_profile.png')
    plt.show()


if __name__ == '__main__':
    for nx in [41, 51, 61, 71, 81, 91]:
        linearconv(nx)
