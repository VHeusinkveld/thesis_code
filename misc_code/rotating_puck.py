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
phi = np.pi/4

n =  np.array([np.sin(theta)*np.cos(-phi + np.pi/2), np.sin(theta)*np.sin(-phi + np.pi/2), np.cos(theta)], dtype=float)
npu = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)], dtype=float)


R = 25/2
L = 2
x0 = 25/2
y0 = 25/2
z0 = 25/2

x = np.linspace(0, 25, 25) - x0
y = np.linspace(0, 25, 25) - y0
z = np.linspace(0, 25, 25) - z0

[X, Y, Z] = np.meshgrid(x, y, z)


Plane_up = (n[0]*X + n[1]*Y + n[2]*Z) - L/2 < 0
Plane_down = (n[0]*X + n[1]*Y + n[2]*Z) + L/2 > 0

Sphere = X**2 + Y**2 + Z**2 < R**2

Puck = Plane_up*Plane_down*Sphere

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.voxels(Puck, edgecolor='k')
ax.quiver(x0,  y0, z0, npu[0], npu[1], npu[2], length=6.25, color='y')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.set_aspect('equal', 'box')

plt.show()
