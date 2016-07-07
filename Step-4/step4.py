import numpy as np
import sympy
from sympy.utilities.lambdify import lambdify
import matplotlib.pyplot as plt
from datetime import datetime


def unidim_burgles(nx=101, nt=100, nu0=0.07, t0=0):
    # Initial conditions for this problems
    x, nu, t = sympy.symbols('x nu t')

    # u at the initial condition depend on the derivative of phi
    phi = sympy.exp(-(x - 4 * t)**2 / (4 * nu * (t + 1))) + \
        sympy.exp(-(x - 4 * t - 2 * np.pi)**2 / (4 * nu * (t + 1)))

    # Derive phi to get the initial condition at fot u
    phiprime = phi.diff(x)

    # Compute u as initial condition
    u = -2 * nu * (phiprime / phi) + 4

    # Lambdify
    ufunc = lambdify((t, x, nu), u)

    # Grid parameters
    # nx and nt are function args
    dx = 2 * np.pi / (nx - 1)
    nu = nu0  # Viscosity. Because nu is used by sympy
    dt = dx * nu

    x = np.linspace(0, 2 * np.pi, nx)
    un = np.empty(nx)
    t = t0  # Because t is used by sympy

    u = np.asarray([ufunc(t, x0, nu) for x0 in x])

    # Set up new figure
    plt.figure(figsize=(11, 7), dpi=100)

    # Ploting the initial problem
    plt.plot(x, u, marker='o', lw=2, label='Initial')

    # Solving equation
    for n in range(nt):
        un = u.copy()
        for i in range(1, nx - 1):
            u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + nu * dt / dx**2 *\
                (un[i + 1] - 2 * un[i] + un[i - 1])
            u[0] = un[0] - un[0] * dt / dx * (un[0] - un[-2]) + nu * dt / dx**2 *\
                (un[1] - 2 * un[0] + un[-2])
            u[-1] = un[-1] - un[-1] * dt / dx * (un[-1] - un[-2]) + nu * dt / dx**2 *\
                (un[0] - 2 * un[-1] + un[-2])

    # Adding the analytical solution
    u_analytical = np.asanyarray([ufunc(nt * dt, xi, nu) for xi in x])

    # Ploting both solutions and initial conditions
    plt.plot(x, u, marker='o', lw=2, label='Computational')
    plt.plot(x, u_analytical, lw=1.5, label='Analytical')

    # Figure parameters
    plt.xlim([0, 2 * np.pi])
    plt.ylim([0, 10])
    plt.legend()

    plt.title(r'$nx = %d$, $nt = %d$, $\nu = %1.2f$, $t = %d$'
              % (nx, nt, nu, t))

    # Figure output settings
    plt.savefig('image_output/solution' +
                str(datetime.now().microsecond) + '.png')
    # plt.show()

    # Closing the figure
    # [http://stackoverflow.com/questions/8213522/matplotlib-
    # clearing-a-plot-when-to-use-cla-clf-or-close]
    plt.close()


if __name__ == '__main__':
    for nx in [101, 201, 301]:
        for nt in [100, 200, 500]:
            for nu in [0.07, 0.1]:
                for t in [0, 2, 10]:
                    unidim_burgles(nx, nt, nu, t)
