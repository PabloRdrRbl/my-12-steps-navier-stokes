import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

# Grid parameters
nx = 81
ny = 81
nt = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (nx - 1)
sigma = 0.2  # CFL condition (adim)
dt = sigma * dx

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

u = np.ones((nx, ny))

# Initial conditions
u[0.5 / dy: 1 / dy + 1, 0.5 / dx: 1 / dx + 1] = 2  # Hat function

# Plot initial conditions
fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(x, y)
ax.plot_wireframe(X, Y, u, color='blue')

# Iterating in two dimensions
for n in range(nt + 1):
    un = u.copy()
    row, col = u.shape
    u[1:, 1:] = (un[1:, 1:] - (c * dt / dx * (un[1:, 1:] - un[1:, :-1])) -
                 (c * dt / dy * (un[1:, 1:] - un[:-1, 1:])))
    # Boundary conditions
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1


ax.plot_wireframe(X, Y, u, color='blue')

plt.show()
