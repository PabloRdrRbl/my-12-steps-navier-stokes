import numpy as np
import matplotlib.pyplot as plt


# plt.ion()  # Turn interactive mode on for plotting

# Grid parameters
nx = 41
nt = 6
sigma = 0.8  # (adim)
c = 1
dx = 0.05

# CFL condition
dt = sigma * dx / c

x = np.arange(0, 2 + dx, dx)
u = np.zeros(nx)
un = np.zeros(nx)
uzero = np.zeros(nx)

# Inital condition (witch hat)
for i in range(nx):
    if 0.9 <= x[i] and x[i] <= 1:
        u[i] = 10 * (x[i] - 0.9)
    if 1.0 <= x[i] and x[i] <= 1.1:
        u[i] = 10 * (1.1 - x[i])

uzero = u.copy()

# Plotting setup
plt.figure()
line, = plt.plot(x, uzero, 'k--')
plt.axis([0, 2, -2, 2])
plt.xlabel('x')
plt.ylabel('u')

for it in range(nt):
    un = u.copy()

    for i in range(1, nx - 1):
        # Maintainig always the same boundary condition (see range())
        # Uncomment as necessary:

        # BD in x
        u[i] = un[i] - c * dt / dx * (un[i] - un[i - 1])

        # CD in x:
#        u[i] = un[i] - c * dt / dx * (un[i + 1] - un[i - 1])

    line.set_ydata(u)

    plt.draw()
    plt.show()
