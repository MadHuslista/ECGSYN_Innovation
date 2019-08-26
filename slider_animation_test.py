#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 16:17:00 2019

@author: madhuslista
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider


fig, ax = plt.subplots()

x1 = np.arange(0, np.pi, 0.1)
x2 = np.arange(np.pi, 2*np.pi, 0.1)
y = np.sin(x1)

line, = ax.plot(x1, y)



#ax.margins(x=0)

axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor='lightgoldenrodyellow')

samp = Slider(axamp, 'Amp', 0.1, 5.0, valinit=1)


#def update(val):
#    amp = samp.val
#    y = amp*np.sin(x)
#    line.set_ydata(y)
#    fig.canvas.draw_idle()

#samp.on_changed(update)



def animate(i, data):
#    samp.on_changed(update)
    line.set_data(data, np.sin(data + i / 100))  # update the data.
    return line,


animation.FuncAnimation(fig, animate, frames =len(x1), fargs = (x1,), interval=2, blit=True, repeat=0)

animation.FuncAnimation(fig1, animate, frames =len(x2), fargs = (x2,), interval=2, blit=True, repeat=0)



