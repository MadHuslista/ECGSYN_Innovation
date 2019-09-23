#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 18:58:29 2019

@author: madhuslista
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
plt.close("all")

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
t = np.arange(0.0, 1.0, 0.001)
a0 = 5
f0 = 3
delta_f = 5.0

s = a0 * np.sin(2 * np.pi * f0 * t)

l, = plt.plot(t, s, lw=2)

ax.margins(x=0)

axcolor = 'lightgoldenrodyellow'

axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

sfreq = Slider(axfreq, 'Freq', 0.1, 30.0, valinit=f0, valstep=delta_f)
samp = Slider(axamp, 'Amp', 0.1, 10.0, valinit=a0)


def update(val):
    amp = samp.val
    freq = sfreq.val
    l.set_ydata(amp*np.sin(2*np.pi*freq*t))
    fig.canvas.draw_idle()
    print(next(iterat()))

def pp(arg): 
    print("perro",arg)

n = 1012012
def iterat(): 
    while True: 
        pp(n)
        yield n 
"""
Con esto demostré que: 
    dentro de un iterador puedo llamar una función, pasarle un argumento nuevo y obtener eso p. 

Necesito: 
    Que el update altere un valor global 
    Y que luego el iterador oucpe ese valor global en el cálculo de los nuevos valores de la función.
    
Problema: 
    Qué tanto demorará? 
    Probar.
"""
sfreq.on_changed(update)
samp.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sfreq.reset()
    samp.reset()
button.on_clicked(reset)

rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)


def colorfunc(label):
    l.set_color(label)
    fig.canvas.draw_idle()
    
radio.on_clicked(colorfunc)

plt.show()
