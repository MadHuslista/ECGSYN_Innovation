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
#param_HVR   = [hrstd, c1, c2, f1, f2]