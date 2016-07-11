import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D


plt.ion()  # Allows interactiva plotings

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
un = u.copy()

# Initial conditions
u[0.5 / dy: 1 / dy + 1, 0.5 / dx: 1 / dx + 1] = 2  # Hat function

# Plot initial conditions
fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u[:])

plt.draw()
plt.pause(5)
