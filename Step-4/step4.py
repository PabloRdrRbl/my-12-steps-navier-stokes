import numpy as np
import sympy
from sympy.utilities.lambdify import lambdify
import matplotlib.pyplot as plt


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
nx = 101
nt = 100
dx = 2 * np.pi / (nx - 1)
nu = 0.07  # Viscosity
dt = dx * nu

x = np.linspace(0, 2 * np.pi, nx)
un = np.empty(nx)
t = 0

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

plt.title('Burgles Equation solution')

# Figure output settings
plt.savefig('image_output/solution.png')
# plt.show()

# Cleaning the figure
# [http://stackoverflow.com/questions/8213522/matplotlib-
# clearing-a-plot-when-to-use-cla-clf-or-close]
plt.clf()
