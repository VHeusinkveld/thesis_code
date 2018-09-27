import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Phi = np.arange(0, 2*np.pi, 2*np.pi/100)
Theta = np.arange(0, np.pi, np.pi/100) # Polar angle

[theta, phi] = np.meshgrid(Theta, Phi)

nx = np.sin(theta)*np.cos(phi)
ny = np.sin(theta)*np.sin(phi)
nz = np.cos(theta)

W = 1

vx = 0
vy = 0
vz = 0

Wx = np.sign(nx)*W*nx**2
Wy = np.sign(ny)*W*ny**2
Wz = np.sign(nz)*W*nz**2

vxn = np.sign(vx**2 + Wx)*np.sqrt(np.abs(vx**2 + Wx))
vyn = np.sign(vy**2 + Wy)*np.sqrt(np.abs(vy**2 + Wy))
vzn = np.sign(vz**2 + Wz)*np.sqrt(np.abs(vz**2 + Wz))


Ekin_proj = vxn**2 + vyn**2 + vzn**2

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(vxn, vyn, vzn)

ax.set_xlim(np.min(vxn), np.max(vxn))
ax.set_ylim(np.min(vyn), np.max(vyn))
ax.set_zlim(np.min(vzn), np.max(vzn))


ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()