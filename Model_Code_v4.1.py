# -*- coding: utf-8 -*-


from scipy.integrate import ode
import math as m
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3

import matplotlib.pyplot as plt 
from matplotlib import animation

from rr_gen import RR_gen
from din_fun import dinamic_function

from datetime import datetime


"""
####################### 0.- PARÁMETROS DE CONFIGURACIÓN ####################################
"""

#Parámetros Tacograma

Resp_by_min = 15
hrmean = 60
hrstd = 1
sfrr = 1
c1 = 2*m.pi*0.01
c2 = 2*m.pi*0.01
f1 = 0.1*2*m.pi
f2 = 0.25*2*m.pi




Anoise = 0.15                             #Additive uniformly distributed measurement noise [0 mV]
Hz_Noise = 50
Hz_Anoise = 0.05

dt = 0.01
n = 20

"""
########################### 1.- CREACIÓN DEL TACOGRAMA ########################### 
"""

rr_times = RR_gen(f1, f2, c1, c2, hrmean, hrstd, n)

rr_axis = []
cumulative_time = 0

for i in rr_times:
    cumulative_time += i 
    rr_axis.append(cumulative_time)


"""
########################### 2.- DEFINICIÓN DE PARÁMETROS Y EMPAQUETAMIENTO ########################### 
"""

hr_factor = np.sqrt(hrmean/60)            #Factor que permite la adaptabilidad de las posiciones al ritmo cardíaco 
fresp = Resp_by_min/60      

#Posición angular de cada Peak
theta_P = -(1/3)*m.pi * np.sqrt(hr_factor)
theta_Q = -(1/12)*m.pi * hr_factor
theta_R = 0
theta_S = (1/12)*m.pi *hr_factor
theta_Td = ((5/9)-(1/60))*m.pi * np.sqrt(hr_factor)
theta_Tu = (5/9)*m.pi * np.sqrt(hr_factor)

#Determina el alto de cada peak
a_P = 0.8
a_Q = -5
a_R = 30
a_S = -7.5
a_Td = 0.5*(hr_factor**(2.5))
a_Tu = 0.75*(hr_factor**(2.5))


#Determina la duración de cada peak 
b_P = 0.2 * hr_factor
b_Q = 0.1 * hr_factor
b_R = 0.1 * hr_factor
b_S = 0.1 * hr_factor
b_Td = 0.4 * (1/hr_factor)
b_Tu = 0.2 * hr_factor

RR = rr_times[0]                        #Esta definición está aquí para poder iniciar el empaquetamiento para el ODE. 
                                        #Si bien 'RR' actúa como constante, como a lo largo del tiempo debe ser actualizada, aquí sería como especificar otro valor inicial 

params = [theta_P, theta_Q, theta_R, theta_S, theta_Td, theta_Tu, a_P, a_Q, a_R, a_S, a_Td, a_Tu, b_P, b_Q, b_R, b_S, b_Td, b_Tu, RR, fresp]

"""Valores Iniciales y empaquetamiento"""

                                        #Son estos los vaores iniciales para el ODE, porque los ciclos parten en el R-peak; el cual está ubicado a 0° en el círculo unitario (i.e X0 = 1, Y0 = 0), y a 0.04 en el eje Z (serán mV?) Con otros valores, la misma secuencia tiende a estos valores iniciales
X0 = 1
Y0 = 0
Z0 = 0.04
y0 = [X0, Y0, Z0]                       #Empaquetamiento de los valores iniciales en una sóla variable. 

"""Construcción del step size para la integración numérica"""
#dt = 0.01                              #Comentado porque se presenta arriba 
t = np.arange(0, rr_axis[-1], dt)       #rr_axis[-1] representa al último elemento de rr_axis



"""
########################### 3.- UTILIZACIÓN DEL ODE SOLVER ########################### 
"""

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


"""
######################## 4.- ESCALAMIENTO y RUIDO ############
"""

psoln_transp = np.array(psoln).T
z = psoln_transp[2]


zmin = min(z)
zmax = max(z)
zrange = zmax - zmin              #Aquí se obtiene el rango máximo de z
z = (z - zmin)*(1.6)/zrange #-0.4    #Aquí cada dato, es escalado en proporción zrange:1.6 con una regla de 3 simple =>  Zrange/(Z- zmin) = 1.6 / x ; donde x es el nuevo valor de z

white_noise = 2*np.random.rand(len(z), 1) -1    #Aquí el np.random.rand() genera ruido aleatorio de distribución uniforme, entre [0,1]. Luego al multiplicar por 2, el rango queda en [0,2], y finalmente al restar en uno, queda [-1,1] => Conclusión: Ruido aleatorio entre -1 y 1
for i in range(len(z)):
    z[i] = z[i] + Anoise*white_noise[i]         #Aquí el ruido aleatorio entre [-1,1] se escala a la magnitud deseada del ruido (Anoise) y se suma a cada valor de z[i]
    
noise = np.sin(2*np.pi*t*Hz_Noise)
print(max(noise), min(noise))
z = z + Hz_Anoise*noise


"""
####################### 5.- GRAFICACÍON CON MATPLOTLIB ###################
"""


x_values = np.array(psoln).T[0]
y_values = np.array(psoln).T[1]
z_values = z

"""Gráfico 2D (t, Z)"""
plt.figure()
plt.plot(t, z_values)
plt.xlabel('time')
plt.ylabel('z')
plt.show()



#Gráfico 3D (X, Y, Z)
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot(x_values, y_values, z_values)
#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')
#plt.show()

"""
###################### 6.- ANIMACIÓN MATPLOTLIB 2D ##############################
"""
fig_2d, ax_2d = plt.subplots()
#Agregar grid reglamentaria del papel al gráfico 

mtr = 5 #Monitor Time Range

data_2d = [t, z_values]

sign, = ax_2d.plot([],[],'g')
signr, = ax_2d.plot([],[],'g')

ax_2d.set_xlim([0,mtr])
ax_2d.set_xlabel('X')

ax_2d.set_ylim(-0.25,1.75)
ax_2d.set_ylabel('Y')

xdata1, ydata1 = [], []
xdata2, ydata2 = [], []


def init():
    ax.set_ylim(-0.25,1.75)
    ax.set_xlim(0, 5)
    del xdata1[:]
    del ydata1[:]
    del xdata2[:]
    del ydata2[:]
    sign, = ax_2d.plot([],[])
    signr,= ax_2d.plot([],[])
    return sign,


def ecg_beat(num, data, sign, signr, hrmean, dt, mtr):
    
    t = data[0]
    z = data[1]

    xmin, xmax = ax_2d.get_xlim()
    pos_inf = int(xmin/dt)
    pos_sup = int(xmax/dt)
    gap = 2

    xdata1 = t[pos_inf:num-gap]
    ydata1 = z[pos_inf:num-gap]
    if num*dt < mtr:
        xdata2 = []
        ydata2 = []
    else:
        xdata2 = t[num+gap:pos_sup]
        ydata2 = z[num+gap-int(mtr/dt):pos_sup-int(mtr/dt)]
    
    if num <= 1: 
        ax_2d.set_xlim(0,mtr)
        ax_2d.figure.canvas.draw()
        
    elif dt*num > xmax: 
        ax_2d.set_xlim(xmin+mtr,xmax+mtr)
        ax_2d.figure.canvas.draw()
 
    sign.set_data(xdata1,ydata1)
    signr.set_data(xdata2,ydata2)  
 
    return sign, signr
    

ani_2d = animation.FuncAnimation(fig_2d,ecg_beat, frames = len(psoln), fargs = (data_2d,sign, signr,hrmean,dt, mtr), interval=1000*dt, blit=1)
plt.show()


"""
###################### 7.- ANIMACIÓN MATPLOTLIB 3D ##############################
"""

fig = plt.figure()
ax = p3.Axes3D(fig)

data = np.array(psoln).T
data[2] = z_values

line, = ax.plot([],[],[])
liner, = ax.plot([],[],[], 'r')

# Setting the axes properties
ax.set_xlim3d([-1.0, 1.0])
ax.set_xlabel('X')

ax.set_ylim3d([-1.0, 1.0])
ax.set_ylabel('Y')

ax.set_zlim3d([-0.5, 1.5])
ax.set_zlabel('Z')


def update(num, data, line, liner, hrmean, dt):
    
    rr_interval = 60/(hrmean*dt)
    segunda_vuelta = rr_interval*2
    
    
    rr_interval = int(rr_interval)    
    segunda_vuelta = int(segunda_vuelta)
    
    if num > segunda_vuelta :
        line.set_data(data[:2, (num-segunda_vuelta):(num-rr_interval+1)])
        line.set_3d_properties(data[2, (num-segunda_vuelta):(num-rr_interval+1)])
        
        liner.set_data(data[:2, (num-rr_interval):num])
        liner.set_3d_properties(data[2, (num-rr_interval):num])
        
                       
        return line,liner


    else:
        liner.set_data(data[:2, 0:num])
        liner.set_3d_properties(data[2, 0:num])
        
        return line,liner
    



ani = animation.FuncAnimation(fig, update, frames=len(psoln), fargs=(data, line, liner, hrmean, dt), interval=1000*dt, blit=1)

plt.show()




"""
###################### ANIMACIÓN MATPLOTLIB 2D ##############################
"""


"""REDISEÑO

OBJETIVOS: 
    + Visualización completa para un número especificado de ciclos y parámetros 
            5.-
	- Visualización animada 2D actualizabley que responde a interacción
            Ya está animada
    - Visualización animada 3D actualizable y que responde a interacción

PROBLEMAS A SOLUCIONAR: 
    - Primero graficar animación 2D tal que se vea igual que la 3D anterior
    - Luego, para cada vez que se actualice el input:
            - Generar un nuevo set de RR_gen
            - Generar un nuevo set de puntos de la función.


"""


"""
PROBLEMAS A SOLUCIONAR 

+ Preparar Input serie temporal (Tuve que pasar del ODEINT a la clase ODE y hacerlo como 'manual')
+ Preparar Input de parámetros.  
+ Velocidad de Animación. Debería ser 60 por min. hasta ahora no lo da. 
+ Abstraer la posición de la sobreposición (El problema era que las abstracciones me estaban quedando en float, pero las listas sólo reconocen int. Había que sólo transformar el tipo)
+ Revisar el parámetro 'n' 
- el tema de la velocidad de sampleo. 
+ Ordenar el tema de la vairable Z que recibe el ruido. Quizás hacer dos psol. psol sólo que sería la señal sin nada y psol_noise que sería la señal con ruido blanco. 
+ Cachar como agregar un ruido de 50 Hz!. 
+ Animar ECG 2D! 
- FUTURO cachar como colocar patologías! 
"""
