#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 23:51:49 2019

@author: meyerhof
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button, RadioButtons
plt.close("all")

#Anim
fig, ax = plt.subplots()
x = np.arange(0, 2*np.pi, 0.01)
line, = ax.plot(x, x)

#Slider 
fig_s, ax_s = plt.subplots()
fig_s.subplots_adjust(left=0.25, bottom=0.25)
axcolor = 'lightgoldenrodyellow'
axfreq = fig_s.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
sfreq = Slider(axfreq, 'Freq', 1, 1000, valinit=0)

gen = 0

def update(val):
    global gen
    gen = sfreq.val
    print(gen)
    #l.set_ydata(amp*np.sin(2*np.pi*freq*t))
    #fig.canvas.draw_idle()
    #print(next(iterat()))

sfreq.on_changed(update)


def pp(i): 
    print("perro")
    global gen
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
