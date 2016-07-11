from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt


# Grid parameters
nx = 41
ny = 41
nt = 120
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = 0.0009
nu = 0.01
dt = sigma * dx * dy / nu

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

u = np.ones((ny, nx))
v = np.ones((ny, nx))
un = np.ones((ny, nx))
vn = np.ones((ny, nx))
comb = np.ones((ny, nx))

# Initial conditions
u[0.5 / dy:1 / dy + 1, 0.5 / dx:1 / dx + 1] = 2
v[0.5 / dy:1 / dy + 1, 0.5 / dx:1 / dx + 1] = 2

# Initial plot
fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
wire1 = ax.plot_wireframe(X, Y, u, cmap=cm.coolwarm)
wire2 = ax.plot_wireframe(X, Y, v, cmap=cm.coolwarm)
ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.set_zlim(1, 2.2)

plt.show()

for n in range(nt + 1):
    un = u.copy()
    vn = v.copy()

    ku1 = un[1:-1, 1:-1]
    ku2 = un[1:-1, 1:-1] * (un[1:-1, 1:-1] - un[1:-1, 0:-2])
    ku3 = vn[1:-1, 1:-1] * (un[1:-1, 1:-1] - un[0:-2, 1:-1])
    ku4 = (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2])
    ku5 = (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])

    u[1:-1, 1:-1] = (ku1 - dt / dx * ku2 - dt / dy * ku3 +
                     nu * dt / dx**2 * ku4 + nu * dt / dy**2 * ku5)

    kv1 = vn[1:-1, 1:-1]
    kv2 = un[1:-1, 1:-1] * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2])
    kv3 = vn[1:-1, 1:-1] * (vn[1:-1, 1:-1] - vn[0:-2, 1:-1])
    kv4 = vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]
    kv5 = vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1]

    v[1:-1, 1:-1] = (kv1 - dt / dx * kv2 - dt / dy * kv3 +
                     nu * dt / dx**2 * kv4 + nu * dt / dy**2 * kv5)

    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1


fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
wire1 = ax.plot_wireframe(X, Y, u)
wire2 = ax.plot_wireframe(X, Y, v)

ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.set_zlim(1, 2.2)

plt.show()
