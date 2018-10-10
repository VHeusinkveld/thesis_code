# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 12:43:15 2018

@author: vince
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=1, metadata=dict(artist='Me'), bitrate=1800)

data_files = []

data_dir = os.getcwd() + "\\data"
for filename in os.listdir(data_dir):
    data_files.append(filename)
#%%
fig = plt.figure()
ax = plt.axes(xlim=(0,0.6), ylim=(0,6))
ln, = ax.plot([],[])

def init():
    ln.set_data([], [])
    return ln,


def update(i):
    f = open(".\\data\\" +  data_files[i], "r")
    y = []
    b = []
    
    for line in f:
        y.append(line.split()[0])
        b.append(line.split()[1])
        
    y = np.asarray(y[1:]).astype(np.float)
    b = np.asarray(b[1:]).astype(np.float)
    ln.set_data(b, y)

    f.close()
    return ln,
    
ani = animation.FuncAnimation(fig, update, 
                              init_func=init, 
                              frames = len(data_files))
ani.save('im.mp4', writer=writer)