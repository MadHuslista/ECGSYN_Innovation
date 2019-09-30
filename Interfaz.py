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
from matplotlib.widgets import Slider, Button, RadioButtons

from rr_gen import RR_gen
from din_fun import dinamic_function
from model_func import model
import variables_func as varfun
import Slider_Interfaz as slid

from datetime import datetime
plt.close("all")



"""
####################### 0.- Inicialización de Parámetros ####################################
#######################       (Variables Globales)       ####################################
"""

param_gener     = varfun.param_gener
param_Artf      = varfun.param_Artf
last_Artf       = param_Artf
param_HVR       = varfun.param_HVR
theta_vals      = varfun.theta_vals
a_vals          = varfun.a_vals
b_vals          = varfun.b_vals
y0              = varfun.y0


"""
####################### 2.2.- Elementos Slider ####################################
"""

def update_Artf(event):
    global param_Artf
    global Flag
    param_Artf[0] = slid.s_Anoise.val
    param_Artf[1] = slid.s_Hznoise.val
    param_Artf[2] = slid.s_AHznoise.val
    Flag = True
slid.fig_Artf.show()
slid.sim_Artf.on_clicked(update_Artf)



def update_gen(event):
    global param_gener
    global Flag
    param_gener[0] = slid.s_hrmean.val      #necesita generación completa de la señal
    param_gener[1] = slid.s_resp.val
    param_gener[2] = slid.s_Amp_ECG.val     #Revisar máximos graph
    param_gener[3] = slid.s_n.val
    param_gener[4] = 1/(10**slid.s_dt.val)  #Revisar FPs abajo. Genera problema con el FPS si no se actualiza hacia abajo también
    param_gener[5] = slid.s_FPS.val
    Flag = True
slid.fig_gen.show()    
slid.sim_gen.on_clicked(update_gen)


"""
####################### 2.1- Generador para Animación ####################################
"""

Flag = False

def generator(dpf): 
    i = 0
    data = []
    
    global param_gener
    global Flag
    
    
    x_val, y_val, z_val, t = model(param_gener, param_Artf, param_HVR, theta_vals, a_vals, b_vals, y0)
    z_m = z_val 
    while True: 
        actual_point = int(dpf*i)
        if Flag:
            print("yas1")
            x_val, y_val, z_m, t = model(param_gener, param_Artf, param_HVR, theta_vals, a_vals, b_vals, y0)
            print("yas2")
            print(len(z_val), len(z_m))
            z_val = z_m
            #z_val[actual_point:] = z_m[actual_point:]
            print("yas3")
            Flag = False
        data = [x_val, y_val, z_val, t, i]
        yield data
        
        i += 1
        n_frames = round(len(t)/DpF)
        if i+1 >= n_frames: 
            i = 0
            z_val = z_m


"""
####################### 3.- Animacióin 2D ####################################
"""

fig_2d, ax_2d = plt.subplots()
#Agregar grid reglamentaria del papel al gráfico 

mtr = 7 #Monitor Time Range
hrmean = param_gener[0]
Amp_ECG = param_gener[2]
FPS = param_gener[5]
dt  = param_gener[4]


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


def ecg_beat(data, sign, signr, hrmean, dt, mtr, DpF):
    
    
    num = data[4]
    t = data[3]
    z = data[2]
    #t = data[0]
    #z = data[1]
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

ani_2d = animation.FuncAnimation(fig_2d,ecg_beat, frames = generator(DpF), init_func=init, fargs = (sign,signr,hrmean,dt, mtr, DpF), interval=FI*1000, blit=1)
plt.show()