# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 12:43:15 2018

@author: vince
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data_files = []

data_dir = os.getcwd() + "\\data_new"
for filename in os.listdir(data_dir):
    data_files.append(filename)
#%%
fig = plt.figure()
ax = plt.axes(xlim=(0,0.1), ylim=(0,5))
ln, = ax.plot([],[])

def init():
    ln.set_data([], [])
    return ln,


def update(i):
    f = open(".\\data_new\\" +  data_files[i], "r")
    y = []
    b = []
    
    for line in f:
        y.append(line.split()[0])
        b.append(line.split()[1])
        
    y = np.asarray(y[2:]).astype(np.float)
    b = np.asarray(b[2:]).astype(np.float)
    ln.set_data(b, y)
    
    f.close()
    return ln,
    
ani = FuncAnimation(fig, update, init_func=init, 
                    frames = len(data_files), 
                    interval=100,
                    blit=True)
plt.show()