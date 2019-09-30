#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 29 15:00:05 2019

@author: meyerhof
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

from datetime import datetime
plt.close("all")


"""
####################### 0.- PARÁMETROS GENERALES ####################################
"""
axcolor = 'lightgoldenrodyellow'

#Ventana para control de Parámetros Generales
fig_gen = plt.figure(figsize=[10,2])

ax_hrmean       = fig_gen.add_axes([0.25, 0.15, 0.65, 0.05], facecolor=axcolor)
ax_Resp_by_min  = fig_gen.add_axes([0.25, 0.25, 0.65, 0.05], facecolor=axcolor)
ax_Amp_ECG      = fig_gen.add_axes([0.25, 0.35, 0.65, 0.05], facecolor=axcolor)
ax_n            = fig_gen.add_axes([0.25, 0.45, 0.65, 0.05], facecolor=axcolor)
ax_dt           = fig_gen.add_axes([0.25, 0.55, 0.65, 0.05], facecolor=axcolor)
ax_FPS          = fig_gen.add_axes([0.25, 0.65, 0.65, 0.05], facecolor=axcolor)
ax_sim_gen      = fig_gen.add_axes([0.8, 0.020, 0.1, 0.1], facecolor=axcolor)


s_hrmean    = Slider(ax_hrmean, 'Frecuencia Cardíaca Promedio', 20, 200, valinit=60, valstep=1)
s_resp      = Slider(ax_Resp_by_min, 'Frecuencia Respiratoria Promedio', 0, 70, valinit=15, valstep=1)
s_Amp_ECG   = Slider(ax_Amp_ECG, 'Amplitud Máxima ECG', 0, 10, valinit=1.7, valstep=0.1)
s_n         = Slider(ax_n, 'Pulsaciones Simuladas', 1, 100, valinit=2, valstep=1)
s_dt        = Slider(ax_dt, 'Frecuencia de Muestreo 10^x',0,5, valinit=3, valstep=1)
s_FPS       = Slider(ax_FPS, 'Cuadros por Segundo', 1, 100, valinit=40, valstep=1)
sim_gen     = Button(ax_sim_gen, 'Simular', color=axcolor, hovercolor='0.975')

#Ventana para control de Parámetros Artefactos
fig_Artf = plt.figure(figsize=[8,1])

ax_Anoise       = fig_Artf.add_axes([0.25, 0.25, 0.65, 0.1], facecolor=axcolor)
ax_Hz_noise     = fig_Artf.add_axes([0.25, 0.5, 0.65, 0.1], facecolor=axcolor)
ax_AHznoise     = fig_Artf.add_axes([0.25, 0.75, 0.65, 0.1], facecolor=axcolor)
ax_sim_Artf     = fig_Artf.add_axes([0.8, 0.020, 0.1, 0.20], facecolor=axcolor)

s_Anoise    = Slider(ax_Anoise, 'Amplitud Ruido Aleatorio', 0, 1, valinit=0.15, valstep=0.01)
s_Hznoise   = Slider(ax_Hz_noise, 'Frecuencia Interferencia', 0, 100, valinit=50, valstep=1)
s_AHznoise  = Slider(ax_AHznoise, 'Amplitud Interferencia', 0, 1, valinit=0.15, valstep=0.01)
sim_Artf    = Button(ax_sim_Artf, 'Simular', color=axcolor, hovercolor='0.975')

#Ventana para control de Parámetros HVR
fig_HVR = plt.figure(figsize=[8,1])

ax_hrstd    = fig_HVR.add_axes([0.25, 0.1, 0.65, 0.1], facecolor=axcolor)
ax_c1       = fig_HVR.add_axes([0.25, 0.3, 0.65, 0.1], facecolor=axcolor)
ax_c2       = fig_HVR.add_axes([0.25, 0.5, 0.65, 0.1], facecolor=axcolor)
ax_f1       = fig_HVR.add_axes([0.25, 0.7, 0.65, 0.1], facecolor=axcolor)
ax_f2       = fig_HVR.add_axes([0.25, 0.8, 0.65, 0.1], facecolor=axcolor)
ax_sim_HVR  = fig_HVR.add_axes([0.8, 0.020, 0.1, 0.20], facecolor=axcolor)

s_hrstd     = Slider(ax_hrstd, 'Desv.Est. Frec. Cardíaca', 0, 1, valinit=0., valstep=0.01)
s_c1        = Slider(ax_c1, 'Desv.Est. Onda Mayer', 0, 0.5, valinit=0.01, valstep=0.01)
s_c2        = Slider(ax_c2, 'Desv.Est. Onda RSA', 0, 0.5, valinit=0.15, valstep=0.01)
s_f1        = Slider(ax_f1, 'Frecuencia Central Onda Mayer', 0, 0.5, valinit=0.1, valstep=0.01)
s_f2        = Slider(ax_f2, 'Frecuencia Central Onda RSA', 0, 0.5, valinit=0.25, valstep=0.01)
sim_HVR     = Button(ax_sim_HVR, 'Simular', color=axcolor, hovercolor='0.975')

#Ventana para control de Parámetros Theta
fig_theta = plt.figure(figsize=[8,4])
ax_circle = plt.subplot(frame_on=False)
plt.axis('off')
circle = plt.Circle((0,0), radius = 1, fill = False)
ax_circle.add_patch(circle)
ax_circle.axis('scaled')
ax_circle.set_position([0.3,0.02,0.9,0.9])

ax_tP       = fig_theta.add_axes([0.05, 0.7, 0.4, 0.05], facecolor=axcolor)
ax_tQ       = fig_theta.add_axes([0.05, 0.6, 0.4, 0.05], facecolor=axcolor)
ax_tR       = fig_theta.add_axes([0.05, 0.5, 0.4, 0.05], facecolor=axcolor)
ax_tS       = fig_theta.add_axes([0.05, 0.4, 0.4, 0.05], facecolor=axcolor)
ax_tTd      = fig_theta.add_axes([0.05, 0.3, 0.4, 0.05], facecolor=axcolor)  
ax_tTu      = fig_theta.add_axes([0.05, 0.2, 0.4, 0.05], facecolor=axcolor)
ax_sim_th   = fig_theta.add_axes([0.05, 0.1, 0.1, 0.05], facecolor=axcolor)

s_tP   = Slider(ax_tP, 'P',     -1, 1, valinit=-1/3, valstep=0.01)
s_tQ   = Slider(ax_tQ, 'Q',     -1, 1, valinit=-1/12, valstep=0.01)
s_tR   = Slider(ax_tR, 'R',     -1, 1, valinit=0.00, valstep=0.01)
s_tS   = Slider(ax_tS, 'S',     -1, 1, valinit=1/12, valstep=0.01)
s_tTd  = Slider(ax_tTd, 'Td',   -1, 1, valinit=(5/9 - 1/60), valstep=0.01)
s_tTu  = Slider(ax_tTu, 'Tu',   -1, 1, valinit=5/9, valstep=0.01)
sim_th = Button(ax_sim_th, 'Simular', color=axcolor, hovercolor='0.975')

p_P  = ax_circle.scatter(np.cos(-1/3*m.pi), np.sin(-1/3*m.pi))
p_Q  = ax_circle.scatter(np.cos(-1/12*m.pi), np.sin(-1/12*m.pi))
p_R  = ax_circle.scatter(np.cos(0), np.sin(0))
p_S  = ax_circle.scatter(np.cos(1/12*m.pi), np.sin(1/12*m.pi))
p_Td = ax_circle.scatter(np.cos((5/9 - 1/60)*m.pi), np.sin((5/9 - 1/60)*m.pi))
p_Tu = ax_circle.scatter(np.cos(5/9*m.pi), np.sin(5/9*m.pi))

an_P  = ax_circle.annotate('P', (1.1*np.cos(-1/3*m.pi), 1.1*np.sin(-1/3*m.pi)))
an_Q  = ax_circle.annotate('Q', (1.1*np.cos(-1/12*m.pi), 1.1*np.sin(-1/12*m.pi)))
an_R  = ax_circle.annotate('R', (1.1*np.cos(0), 1.1*np.sin(0)))
an_S  = ax_circle.annotate('S', (1.1*np.cos(1/12*m.pi), 1.1*np.sin(1/12*m.pi)))
an_Td = ax_circle.annotate('Td', (1.1*np.cos((5/9 - 1/60)*m.pi), 1.1*np.sin((5/9 - 1/60)*m.pi)))
an_Tu = ax_circle.annotate('Tu', (1.1*np.cos(5/9*m.pi), 1.1*np.sin(5/9*m.pi)))


def update_circle(val):
    hr_factor = np.sqrt(s_hrmean.val/60)
    tP      = s_tP.val*m.pi      *np.sqrt(hr_factor)
    tQ      = s_tQ.val*m.pi      *hr_factor
    tR      = s_tR.val*m.pi
    tS      = s_tS.val*m.pi      *hr_factor
    tTd     = s_tTd.val*m.pi     *np.sqrt(hr_factor)
    tTu     = s_tTu.val*m.pi     *np.sqrt(hr_factor)
    
    p_P.set_offsets( [np.cos(tP),np.sin(tP)] )
    p_Q.set_offsets( [np.cos(tQ),np.sin(tQ)] )
    p_R.set_offsets( [np.cos(tR),np.sin(tR)] )
    p_S.set_offsets( [np.cos(tS),np.sin(tS)] )
    p_Td.set_offsets( [np.cos(tTd),np.sin(tTd)] )
    p_Tu.set_offsets( [np.cos(tTu),np.sin(tTu)] )
    
    an_P.set_position((1.1*np.cos(tP), 1.1*np.sin(tP)))
    an_Q.set_position((1.1*np.cos(tQ), 1.1*np.sin(tQ)))
    an_R.set_position((1.1*np.cos(tR), 1.1*np.sin(tR)))
    an_S.set_position((1.1*np.cos(tS), 1.1*np.sin(tS)))
    an_Td.set_position((1.12*np.cos(tTd), 1.12*np.sin(tTd)))
    an_Tu.set_position((1.12*np.cos(tTu), 1.12*np.sin(tTu)))
    fig_theta.canvas.draw_idle()

s_hrmean.on_changed(update_circle)
s_tP.on_changed(update_circle)
s_tQ.on_changed(update_circle)
s_tR.on_changed(update_circle)
s_tS.on_changed(update_circle)
s_tTd.on_changed(update_circle)
s_tTu.on_changed(update_circle)