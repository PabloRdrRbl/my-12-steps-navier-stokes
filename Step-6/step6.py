import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

nx = 101
ny = 101
nt = 80
dx = 2 / (nx - 1)  # We are braking 2 in 100 parts
dy = 2 / (nx - 1)
sigma = 0.2
dt = sigma * dx

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

u = np.ones((ny, nx))
v = np.ones((ny, nx))

un = u.copy()
vn = v.copy()

# Initial conditions
u[.5 / dy:1 / dy + 1, .5 / dx:1 / dx + 1] = 2
v[.5 / dy:1 / dy + 1, .5 / dx:1 / dx + 1] = 2

# Solving
for n in range(nt + 1):
    un = u.copy()
    vn = v.copy()

    ku1 = un[1:, 1:]
    ku2 = un[1:, 1:] * dt / dx
    ku3 = un[1:, 1:] - un[1:, :-1]
    ku4 = dt / dy * (un[1:, 1:] - un[:-1, 1:])

    kv1 = vn[1:, 1:]
    kv2 = vn[1:, 1:] * dt / dx
    kv3 = vn[1:, 1:] - vn[1:, :-1]
    kv4 = dt / dy * (vn[1:, 1:] - vn[:-1, 1:])

    u[1:, 1:] = ku1 - (ku2 * ku3) - vn[1:, 1:] * ku4
    v[1:, 1:] = kv1 - (kv2 * kv3) - vn[1:, 1:] * kv4

    # Boundary conditions
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1

# Plotings
from matplotlib import cm

fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)

ax.plot_surface(X, Y, u, cmap=cm.coolwarm)
plt.show()

fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)

ax.plot_surface(X, Y, v, cmap=cm.coolwarm)
plt.show()
