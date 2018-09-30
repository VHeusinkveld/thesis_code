import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Phi = np.arange(0, 2.1*np.pi, 2*np.pi/30)
Theta = np.arange(0, 1.1*np.pi, np.pi/30) # Polar angle

[theta, phi] = np.meshgrid(Theta, Phi)

n = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)], dtype = float)
v = np.array([0, 0, 0], dtype  = float)

W = 3

vn = 0*n
for i in range(3):
    wsgn = np.sign(n[i]*v[i]) + (n[i]*v[i] == 0)*np.sign(n[i])
    w = np.array([wsgn*W*n[i]**2],  dtype = float)
    vsgn =  1*(v[i] >= 0) * ((v[i]**2 + w) > 0) + \
           -1*(v[i] >= 0) * ((v[i]**2 + w) < 0) + \
            1*(v[i] <  0) * ((v[i]**2 + w) < 0) + \
           -1*(v[i] <  0) * ((v[i]**2 + w) > 0) 
    
    vn[i] = vsgn*np.sqrt(np.abs(v[i]**2 + w))

#Ekin_proj = vxn**2 + vyn**2 + vzn**2
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(vn[0], vn[1], vn[2])
ax.quiver(0, 0, 0, v[0], v[1], v[2], color='g')

ax.set_xlim(np.min(vn[0]), np.max(vn[0]))
ax.set_ylim(np.min(vn[1]), np.max(vn[1]))
ax.set_zlim(np.min(vn[2]), np.max(vn[2]))


ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()