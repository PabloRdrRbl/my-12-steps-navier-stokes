import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


# Grid parameters
nx = 31
ny = 31
nt = 17
nu = 0.05
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = 0.25
dt = sigma * dx * dy / nu

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

u = np.ones((ny, nx))
un = u.copy()

# Initial conditions
u[0.5 / dy:1 / dy + 1, 0.5 / dx:1 / dx + 1] = 2

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.set_zlim(1, 2.5)

plt.savefig('gif/photo-001.png')


def diffuse(nt):
    """
    Solve
    """
    u[0.5 / dy:1 / dy + 1, 0.5 / dx:1 / dx + 1] = 2

    for n in range(nt + 1):
        un = u.copy()

        k1 = un[1:-1, 1:-1]
        k2 = (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2])
        k3 = (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])

        u[1:-1, 1:-1] = (k1 + nu * dt / dx**2 * k2 + nu * dt / dy**2 * k3)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, u, rstride=1, cstride=1,
                           cmap=cm.coolwarm, linewidth=0, antialiased=True)

    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2)
    ax.set_zlim(1, 2.5)


if __name__ == '__main__':
    for i, nt in enumerate(range(0, 90, 2)):
        diffuse(nt)
        plt.savefig('gif/photo-%03d.png' % (i + 2))
