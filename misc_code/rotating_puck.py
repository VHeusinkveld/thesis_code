# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 11:46:13 2018

@author: vince
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

theta = np.pi/4
phi = 0

n = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)], dtype=float)

R = 5
L = 1
x0 = 5
y0 = 5
z0 = 5

x = np.linspace(0, 10, 20) - x0
y = np.linspace(0, 10, 20) - y0
z = np.linspace(0, 10, 20) - z0

[X, Y, Z] = np.meshgrid(x, y, z)


Plane_up = -(n[0]*X + n[1]*Y + n[2]*Z) - L/2 < 0
Plane_down = -(n[0]*X + n[1]*Y + n[2]*Z) + L/2 > 0

Sphere = X**2 + Y**2 + Z**2 < R**2

Puck = Plane_up*Plane_down*Sphere

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.voxels(Puck, edgecolor='k')

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

plt.show()
