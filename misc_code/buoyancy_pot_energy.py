# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 10:56:08 2018

@author: vince
"""

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt


p0 = 1e5
Rsp = 287

z = np.arange(0, 40, 0.1, dtype=float)
T = 273 + 0.5*z

rho = p0/(Rsp*T)
f = plt.figure(figsize=(10, 5))
ax1 = plt.subplot(121)
ax2 = plt.subplot(122)

ax1.plot(rho, z)
ax1.set_xlabel('rho [kg/m3]')
ax1.set_ylabel('z [m]')

ax2.plot(T, z)
ax2.set_xlabel('T [K]')
ax2.set_ylabel('z [m]')