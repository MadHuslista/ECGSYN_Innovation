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
fig_gen = plt.figure(figsize=[11,4])


#Slider HR_Mean
ax_hrmean       = fig_gen.add_axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
ax_Resp_by_min  = fig_gen.add_axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)
ax_Amp_ECG      = fig_gen.add_axes([0.25, 0.3, 0.65, 0.03], facecolor=axcolor)
ax_n            = fig_gen.add_axes([0.25, 0.4, 0.65, 0.03], facecolor=axcolor)
ax_dt           = fig_gen.add_axes([0.25, 0.5, 0.65, 0.03], facecolor=axcolor)
ax_FPS          = fig_gen.add_axes([0.25, 0.6, 0.65, 0.03], facecolor=axcolor)

Amp_ECG = 1.7                       #Amplitud Máxima ECG
n = 10                              #Cantidad de Pulsaciones simuladas
dt = 0.001                           # En segundos
FPS = 40


s_hrmean    = Slider(ax_hrmean, 'Frecuencia Cardíaca Promedio', 20, 200, valinit=60, valstep=1)
s_resp      =  Slider(ax_Resp_by_min, 'Frecuencia Respiratoria Promedio', 0, 70, valinit=15, valstep=1)
s_Amp_ECG   = Slider(ax_Amp_ECG, 'Amplitud Máxima ECG', 0, 10, valinit=1.7, valstep=0.1)
s_n         = Slider(ax_n, 'Pulsaciones Simuladas', 20, 200, valinit=60, valstep=1)
s_dt        = Slider(ax_dt, 'Frecuencia de Muestreo 10^x',0,5, valinit=3, valstep=1)
s_FPS       = Slider(ax_FPS, 'Cuadros por Segundo', 25, 100, valinit=40, valstep=1)




