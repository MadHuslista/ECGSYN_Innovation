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

gen = 0

def pp(i): 
    print("perro")
    gen = i
    return gen , 2
    #global gen
    #gen = i

def iterat(n): 
    num = 5
    n = num +n
    while num < n: 
        d = pp(num)
        print(d)
        yield d
        num += 1    

def animate(h):
    line.set_ydata((x + h[0] / 200))  # update the data.
    return line,


ani = animation.FuncAnimation(fig, animate, frames=iterat(i), interval=2, blit=True, save_count=50)


plt.show()
