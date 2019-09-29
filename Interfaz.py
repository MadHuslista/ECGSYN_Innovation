#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:43:59 2019

@author: meyerhof
"""

"""
Llevar todo el Model  a una función. 
Aquí armar luego, varias sliders. 
Cada slider para cada parámetro 
    - Aquí armar el circulito de las posiciones. 
    - Quizás armar una gaussiana genérica para setear amplitud y duración 
Luego, crear un botón de "GENERAR", 
    - Aquí en el update llamar a todo el Model con todo lo seteado
""" 



from scipy.integrate import ode
import math as m
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3

import matplotlib.pyplot as plt 
from matplotlib import animation

from rr_gen import RR_gen
from din_fun import dinamic_function
from model_func import model
import variables_func as varfun

from datetime import datetime
plt.close("all")



"""
####################### 0.- Obtención de parámetros ####################################
"""

param_gener     = varfun.param_gener
param_Artf      = varfun.param_Artf
param_HVR       = varfun.param_HVR
theta_vals      = varfun.theta_vals
a_vals          = varfun.a_vals
b_vals          = varfun.b_vals
y0              = varfun.y0

"""
####################### 1.- Cálculo de datos ####################################
"""

x_val, y_val, z_val, t = model(param_gener, param_Artf, param_HVR, theta_vals, a_vals, b_vals, y0)


"""
####################### 2.- Animacióin 2D ####################################
"""

fig_2d, ax_2d = plt.subplots()
#Agregar grid reglamentaria del papel al gráfico 

mtr = 2 #Monitor Time Range
FPS = param_gener[5]
dt  = param_gener[4]

data_2d = [t, z_val]

sign, = ax_2d.plot([],[],'g')
signr, = ax_2d.plot([],[],'g')

ax_2d.set_xlim([0,mtr])
ax_2d.set_xlabel('t [s]')

ax_2d.set_aspect(0.4)

ax_2d.xaxis.grid(True, which='major', lw= 1.5)
ax_2d.xaxis.grid(True, which='minor', lw= 0.5)

ax_2d.set_ylim(Amp_ECG*-0.15,Amp_ECG*1.09)
ax_2d.set_ylabel('V [mV]')

ax_2d.set_yticks(np.arange(Amp_ECG*-0.15,Amp_ECG*1.09, step=0.5), minor = False)                
ax_2d.set_yticks(np.arange(Amp_ECG*-0.15,Amp_ECG*1.09, step=0.1), minor = True)
ax_2d.yaxis.grid(True, which='major', lw= 1.5)
ax_2d.yaxis.grid(True, which='minor', lw= 0.5)


xdata1, ydata1 = [], []
xdata2, ydata2 = [], []

FI = 1 / FPS    #Frame Interval 
DpF = FI/ dt    #Datos por frame


def init():                     #Sin esta función también funciona. Documentación sugiere que es más eficiente. No lo sé
    ax_2d.set_ylim(Amp_ECG*-0.15,Amp_ECG*1.09)
    ax_2d.set_xlim(0, mtr)
    del xdata1[:]
    del ydata1[:]
    del xdata2[:]
    del ydata2[:]
    sign, = ax_2d.plot([],[])
    signr,= ax_2d.plot([],[])
    return sign, signr


def ecg_beat(num, data, sign, signr, hrmean, dt, mtr, DpF):
    
    
    
    t = data[0]
    z = data[1]
    #Posible mejora: Usar el argumento 'Frames' para pasar la data. Ahora, para cada frame, le paso la lista completa de datos. Mucho 
    
    time_gap = 0.01   #Separación entre la nueva señal y la anterior. En [s]
    
    data_gap = time_gap/dt    #
    growth_cursor = int(round(num*DpF - int(data_gap/2)))
    decrease_cursor = int(round(num*DpF + int(data_gap/2)))
    
    #print()
    
    xmin, xmax = ax_2d.get_xlim()
    pos_inf = int(xmin/dt)
    pos_sup = int(xmax/dt)

    if pos_sup > len(t)-1:
        pos_sup = len(t)-1
    
    
    
    xdata1 = t[pos_inf:growth_cursor]
    ydata1 = z[pos_inf:growth_cursor]
    
    
    if growth_cursor*dt < mtr:
        xdata2 = []
        ydata2 = []
    else:
        xdata2 = t[decrease_cursor:pos_sup]
        ydata2 = z[decrease_cursor-int(mtr/dt):pos_sup-int(mtr/dt)]

        
    
    if num <= 1: 
        ax_2d.set_xlim(0,mtr)
        
        ax_2d.set_xticks(np.arange(0,mtr, step=0.2), minor = False)                
        ax_2d.set_xticks(np.arange(0,mtr, step=0.04), minor = True)

        
        ax_2d.figure.canvas.draw()
        
    elif growth_cursor*dt > xmax: 
        ax_2d.set_xlim(xmin+mtr,xmax+mtr)                

        ax_2d.set_xticks(np.arange(xmin+mtr,xmax+mtr, step=0.2), minor = False)                
        ax_2d.set_xticks(np.arange(xmin+mtr,xmax+mtr, step=0.04), minor = True)

        ax_2d.figure.canvas.draw()
 
    sign.set_data(xdata1,ydata1)
    signr.set_data(xdata2,ydata2)  
 
    return sign, signr
    


ani_2d = animation.FuncAnimation(fig_2d,ecg_beat, frames = round(len(psoln)/DpF), init_func=init, fargs = (data_2d,sign,signr,hrmean,dt, mtr, DpF), interval=FI*1000, blit=1)
plt.show()