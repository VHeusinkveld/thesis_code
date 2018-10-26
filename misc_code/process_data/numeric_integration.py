# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 23:06:31 2018

@author: vince
"""

import numpy as np

P = 10
W = 0.3
m = 5

dt = 0.0001

x = 0
v = 0
E = 0

while x < W:
    x += v*dt
    E += P*dt
    v = np.sqrt(2*E/m)
    
v_calc = (3*P*W/m)**(1/3)    
print(x, v, v_calc)


