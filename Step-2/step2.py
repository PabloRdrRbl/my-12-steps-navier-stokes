import numpy as np
import matplotlib.pyplot as plt


def linearconv(nx):
    # Grid parameters
    nx = 41
    dx = 2. / (nx - 1)
    nt = 20
    dt = 0.025

    # Initial conditions
    u = np.ones(nx)
    u[0.5 / dx: 1 / dx + 1] = 2

    # Printing u
    plt.plot(np.linspace(0, 2, nx), u)

    un = np.empty(nx)  # Temprary array to hold the time iteration

    # Solution
    for i in range(nt):
        un = u.copy()
        for i in range(1, nx):  # We are not taking indexes 0 and -1

            # given the boundary condition
            u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1])

    # Printing the new u
    plt.plot(np.linspace(0, 2, nx), u)

    plt.savefig('u_profile.png')
    plt.show()


if __name__ == '__main__':
    for nx in [41, 51, 61, 71, 81, 91, 141]:
        linearconv(nx)
