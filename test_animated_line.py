#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 23:51:49 2019

@author: meyerhof
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.close("all")

fig, ax = plt.subplots()

x = np.arange(0, 2*np.pi, 0.01)
line, = ax.plot(x, x)

def pp(): 
    print("perro")

def iterat(n): 
    num = 5
    n = num +n
    while num < n: 
        pp()
        yield num 
        num += 1
        


def init():  # only required for blitting to give a clean slate.
    line.set_ydata([np.nan] * len(x))
    return line,


def animate(i):
    line.set_ydata((x + i / 200))  # update the data.
    return line,


ani = animation.FuncAnimation(fig, animate, init_func=init, interval=2, blit=True, save_count=50)

# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# from matplotlib.animation import FFMpegWriter
# writer = FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)

plt.show()
