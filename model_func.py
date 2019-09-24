#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 19:53:03 2019

@author: meyerhof
"""


from scipy.integrate import ode
import math as m
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3

import matplotlib.pyplot as plt 
from matplotlib import animation

from rr_gen import RR_gen
from din_fun import dinamic_function

from datetime import datetime
plt.close("all")

"""
####################### 0.- PARÁMETROS DE CONFIGURACIÓN ####################################
"""

#Parámetros Generales

hrmean = 60                         #Frecuencia Cardíaca
Resp_by_min = 15                    #Frecuencia Respiratoria
Amp_ECG = 1.7                       #Amplitud Máxima ECG
n = 8                              #Cantidad de Pulsaciones simuladas
dt = 0.001                           # En segundos
FPS = 40

#Control de Artefactos
Anoise = 0.15                       #Amplitud del Ruido Aleatorio
Hz_Noise = 50                       #Frecuencia de la Interferencia
Hz_Anoise = 0.05                    #Amplitud de la Interferencia


#Variabilidad del Pulso Cardíaco
hrstd = 0                           #Desviación Estándar de la Frecuencia Cardíaca
c1 = 2*m.pi*0.01                    #Desviación Estándar Onda Mayer
c2 = 2*m.pi*0.01                    #Desviación Estándar Onda RSA
f1 = 0.1*2*m.pi                     #Frecuencia Central Onda Mayer
f2 = 0.25*2*m.pi                    #Frecuencia Central Onda RSA


#Posición angular de cada Peak
theta_P = -(1/3)*m.pi 
theta_Q = -(1/12)*m.pi 
theta_R = 0
theta_S = (1/12)*m.pi 
theta_Td = ((5/9)-(1/60))*m.pi 
theta_Tu = (5/9)*m.pi 

#Determina el alto de cada peak
a_P = 0.8
a_Q = -5
a_R = 30
a_S = -7.5
a_Td = 0.5
a_Tu = 0.75


#Determina la duración de cada peak 
b_P = 0.2 
b_Q = 0.1
b_R = 0.1
b_S = 0.1 
b_Td = 0.4 
b_Tu = 0.2 


#Valores Iniciales y empaquetamiento                                     
X0 = 1
Y0 = 0
Z0 = 0.04

#Empaquetamiento 
param_gener = [hrmean, Resp_by_min, Amp_ECG, n, dt, FPS]
param_Artf  = [Anoise, Hz_Noise, Hz_Anoise]
param_HVR   = [hrstd, c1, c2, f1, f2]

theta_vals  = [theta_P, theta_Q, theta_R, theta_S, theta_Td, theta_Tu]
a_vals      = [a_P, a_Q, a_R, a_S, a_Td, a_Tu]
b_vasl      = [b_P, b_Q, b_R, b_S, b_Td, b_Tu]
y0 = [X0, Y0, Z0] 



def model(param_g, param_A, param_H, theta, a, b, y0): 
    
    """#### Desempaquetado ### """
    
    hrmean, Resp_by_min, Amp_ECG, n, dt, FPS    = param_g
    Anoise, Hz_Noise, Hz_Anoise                 = param_A
    hrstd, c1, c2, f1, f2                       = param_H
    
    theta_P, theta_Q, theta_R, theta_S, theta_Td, theta_Tu  = theta
    a_P, a_Q, a_R, a_S, a_Td, a_Tu                          = a
    b_P, b_Q, b_R, b_S, b_Td, b_Tu                          = b
    
    
    """#### CÁLCULO DE VALORES ### """
    
    hr_factor = np.sqrt(hrmean/60)            #Factor que permite la adaptabilidad de las posiciones al ritmo cardíaco 
    fresp = Resp_by_min/60      

    #Posición angular de cada Peak
    theta_P = theta_P   * np.sqrt(hr_factor)
    theta_Q = theta_Q   * hr_factor
    theta_R = theta_R
    theta_S = theta_S   * hr_factor
    theta_Td = theta_Td * np.sqrt(hr_factor)
    theta_Tu = theta_Tu * np.sqrt(hr_factor)
    
    #Determina el alto de cada peak
    a_P = a_P
    a_Q = a_Q
    a_R = a_R
    a_S = a_S
    a_Td = a_Td *(hr_factor**(2.5))
    a_Tu = a_Tu  *(hr_factor**(2.5))
    
    
    #Determina la duración de cada peak 
    b_P = b_P * hr_factor
    b_Q = b_Q * hr_factor
    b_R = b_R * hr_factor
    b_S = b_S * hr_factor
    b_Td = b_Td * (1/hr_factor)
    b_Tu = b_Tu * hr_factor
    
    """#### Creación Tacograma ### """
    
    rr_times = RR_gen(f1, f2, c1, c2, hrmean, hrstd, n)

    rr_axis = []
    cumulative_time = 0

    for i in rr_times:
        cumulative_time += i 
        rr_axis.append(cumulative_time)
    
    t = np.arange(0, rr_axis[-1], dt)       #rr_axis[-1] representa al último elemento de rr_axis

    RR = rr_times[0]                        #Esta definición está aquí para poder iniciar el empaquetamiento para el ODE. 
                                            #Si bien 'RR' actúa como constante, como a lo largo del tiempo debe ser actualizada, aquí sería como especificar otro valor inicial 

    params = [theta_P, theta_Q, theta_R, theta_S, theta_Td, theta_Tu, a_P, a_Q, a_R, a_S, a_Td, a_Tu, b_P, b_Q, b_R, b_S, b_Td, b_Tu, RR, fresp]
    
    
    """#### Utilización del solver ### """

    
    solver = ode(dinamic_function)          #Creación de una instancia del ODE, con la 'dinamic_function' como función llamable
    solver.set_integrator("lsoda")          #Se setea el método de integración
    solver.set_initial_value(y0)            #Se setean los valores iniciales para X, Y, Z
    solver.set_f_params(params)             #Se setean los parametros para la función llamable
    
    
    pos = 0
    psoln = []                                                  
    for i in range(len(t)):         #Este for permite que cada vez que el parámetro RR, se actualice cada vez que el t alcance al siguiente intervalo RR. Para esto, la serie de intervalos de rr_time se pasaron a una escala temporal en rr_axis. Y cada vez que t alcanza al siguiente rr_axis, rr_time[pos] se actualiza al siguiente valor. 
        if t[i] > rr_axis[pos]:
            params[-2] = rr_times[pos+1]
            solver.set_f_params(params)        
            pos = pos+1    
            
        solver.integrate(solver.t+dt)
        psoln.append(solver.y)

    """#### Escalamiento y Ruido ### """
        
    psoln_transp = np.array(psoln).T
    z = psoln_transp[2]
    
    
    zmin = min(z)
    zmax = max(z)
    zrange = zmax - zmin              #Aquí se obtiene el rango máximo de z
    z = (z - zmin)*(Amp_ECG)/zrange #-0.4    #Aquí cada dato, es escalado en proporción zrange:1.6 con una regla de 3 simple =>  Zrange/(Z- zmin) = 1.6 / x ; donde x es el nuevo valor de z
    
    white_noise = 2*np.random.rand(len(z), 1) -1    #Aquí el np.random.rand() genera ruido aleatorio de distribución uniforme, entre [0,1]. Luego al multiplicar por 2, el rango queda en [0,2], y finalmente al restar en uno, queda [-1,1] => Conclusión: Ruido aleatorio entre -1 y 1
    for i in range(len(z)):
        z[i] = z[i] + Anoise*white_noise[i]         #Aquí el ruido aleatorio entre [-1,1] se escala a la magnitud deseada del ruido (Anoise) y se suma a cada valor de z[i]
        
    noise = np.sin(2*np.pi*t*Hz_Noise)
    z = z + Hz_Anoise*noise
    
    
    x_values = np.array(psoln).T[0]
    y_values = np.array(psoln).T[1]
    z_values = z
    
    return x_values, y_values, z_values, t
    
    
#model(param_gener, param_Artf, param_HVR, theta_vals, a_vals, b_vasl, y0)