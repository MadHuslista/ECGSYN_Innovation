# -*- coding: utf-8 -*-


from scipy.integrate import ode
import scipy as sp
import math as m
import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3

import matplotlib.pyplot as plt 
from matplotlib import animation

"""
###################################################################################
####################### PARÁMETROS DE CONFIGURACIÓN ####################################
###################################################################################

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




Anoise = 0.025 #Additive uniformly distributed measurement noise [0 mV]
Hz_Noise = 50
Hz_Anoise = 0.05

dt = 0.01
n = 100



"""
###################################################################################
####################### CREACIÓN DE TACOGRAMA ####################################
###################################################################################
"""
"""
f = open('rr_list.dat')
lines =  f.readlines()
rr_times = []
for i in lines:
    rr_times.append(float(i))
f.close()
"""
#rr_times = [1, 2, 3, 4, 5]


rrstd = 60*(hrstd)/(hrmean*hrmean)
rrmean = 60/hrmean
sfrr = 1
theta1 = 0.5
theta2 = 1
df = sfrr/n
w = np.arange(0,n,1)*2*m.pi*df
def s1(x):
    return(((theta1)*m.exp(-0.5*((x-f1)/c1)**2))/m.sqrt(2*m.pi*c1**2))
def s2(y):
    return(((theta2)*m.exp(-0.5*((y-f2)/c2)**2))/m.sqrt(2*m.pi*c2**2))
sf = []
for i in w:
    suma = s1(i)+s2(i)
    sf.append(suma)
plt.plot(sf)
sf0 = []
for i in range(0,int(n/2)):
    piv = sf[i]
    sf0.append(piv)
for i in range(int((n/2)-1),-1,-1):
    piv = sf[i]
    sf0.append(piv)
sf1 = []
for i in sf0:
    piv = (sfrr/2)*m.sqrt(i)
    sf1.append(piv)

var = np.random.rand(int(n/2+1))
ph0 = []
for i in range(int((n/2)-1)):
    piv= 1*var[i]*2*m.pi
    ph0.append(piv)
ph = [0]
for i in range(len(ph0)):
    piv = ph0[i]
    ph.append(piv)
ph.append(0)
for i in range(len(ph0)-1,-1,-1):
    piv = -ph0[i]
    ph.append(piv)
sfw = []
for i in range(len(ph)):
    piv = sf1[i]*np.exp(1j*ph[i])
    sfw.append(piv)
x = (1/n)*sp.real(sp.ifft(sfw))
xstd = np.std(x)
ratio = rrstd / xstd
rr = []
for i in range(len(x)):
    piv = rrmean + x[i]*ratio
    rr.append(piv)

rr_times = rr       #SERIE DE INTERVALOS RR

#Time Tools

rr_axis = []
cumulative_time = 0

for i in rr_times:
    cumulative_time += i 
    rr_axis.append(cumulative_time)

"""
###################################################################################
##############################  RESOLUCIÓN DE ECUACIONES DIFERENCIALES#############
###################################################################################
"""

"""Definición de Parámetros y empaquetamiento"""

hr_factor = np.sqrt(hrmean/60)  #Factor que permite la adaptabilidad de las posiciones al ritmo cardíaco 
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

RR = rr_times[0]    #Esta definición está aquí para poder iniciar el empaquetamiento para el ODE. 
                    #Si bien 'RR' actúa como constante, como a lo largo del tiempo debe ser actualizada, aquí sería como especificar otro valor inicial 

params = [theta_P, theta_Q, theta_R, theta_S, theta_Td, theta_Tu, a_P, a_Q, a_R, a_S, a_Td, a_Tu, b_P, b_Q, b_R, b_S, b_Td, b_Tu, RR, fresp]

"""Valores Iniciales y empaquetamiento"""

#Son estos los vaores iniciales para el ODE, porque los ciclos parten en el R-peak; el cual está ubicado a 0° en el círculo unitario (i.e X0 = 1, Y0 = 0), y a 0.04 en el eje Z (serán mV?) Con otros valores, la misma secuencia tiende a estos valores iniciales
X0 = 1
Y0 = 0
Z0 = 0.04
y0 = [X0, Y0, Z0]       #Empaquetamiento de los valores iniciales en una sóla variable. 

"""Construcción del step size para la integración numérica"""
#dt = 0.01      #Comentado porque se presenta arriba 
t = np.arange(0, rr_axis[-1], dt)      #rr_axis[-1] representa al último elemento de rr_axis


"""
Argumentos func para ODEINT. 
Función Python que retorna una lista de valores correspondientes a las n 
funciones correspondientes al modelo dinámico, a un tiempo t. 
"""

def dinamic_function(t, y, params):
    
    X, Y, Z = y     #Desempaquetamiento de los valores de las funciones involucradas en un tiempo 't'
    theta_P, theta_Q, theta_R, theta_S, theta_Td, theta_Tu, a_P, a_Q, a_R, a_S, \
    a_Td, a_Tu, b_P, b_Q, b_R, b_S, b_Td, b_Tu, RR, fresp = params #Desempaquetamiento de los parámetros involucrados en cada función derivada
    
    
    #Variables utilitarias para mantener el código legible y ordenado
    alfa = 1 - m.sqrt(X**2 + Y**2)
    w = (2*m.pi)/RR      # w = 2pi/RR . Se asume RR = 1 seg
    #print(w, RR)
    theta = m.atan2(Y, X)
    #print(theta)
    
    delta_theta_P = np.fmod((theta - theta_P),(2*m.pi))
    delta_theta_Q = np.fmod((theta - theta_Q),(2*m.pi))
    delta_theta_R = np.fmod((theta - theta_R),(2*m.pi))
    delta_theta_S = np.fmod((theta - theta_S),(2*m.pi))
    delta_theta_Td = np.fmod((theta - theta_Td),(2*m.pi))
    delta_theta_Tu = np.fmod((theta - theta_Tu),(2*m.pi))


    Zo = 0.005*m.sin(2*m.pi*(fresp)*t)
    #Zo = 0
    
    sumatoria_P = a_P * delta_theta_P * m.exp(  ( - (delta_theta_P**2)/(2 * b_P**2)   )  ) 
    sumatoria_Q = a_Q * delta_theta_Q * m.exp(  ( - (delta_theta_Q**2)/(2 * b_Q**2)   )  ) 
    sumatoria_R = a_R * delta_theta_R * m.exp(  ( - (delta_theta_R**2)/(2 * b_R**2)   )  ) 
    sumatoria_S = a_S * delta_theta_S * m.exp(  ( - (delta_theta_S**2)/(2 * b_S**2)   )  ) 
    sumatoria_Td = a_Td * delta_theta_Td * m.exp(  ( - (delta_theta_Td**2)/(2 * b_Td**2)   )  ) 
    sumatoria_Tu = a_Tu * delta_theta_Tu * m.exp(  ( - (delta_theta_Tu**2)/(2 * b_Tu**2)   )  )
    
    derivs = [
            
          alfa*X - w*Y,                                                                         #Eq para X'
          alfa*Y + w*X,                                                                         #Eq para Y'
          -(sumatoria_P + sumatoria_Q + sumatoria_R + sumatoria_S + sumatoria_Td + sumatoria_Tu) - (Z - Zo)     #Eq para Z'
            
            ]
    
    #print('z ', derivs[2])
    
    return derivs


"""Utilización del ODE solver"""

#psoln = ode(dinamic_function, y0, t, args=(params,))       ODEINT


solver = ode(dinamic_function)          #Creación de una instancia del ODE, con la 'dinamic_function' como función llamable
solver.set_integrator("lsoda")          #Se setea el método de integración
solver.set_initial_value(y0)            #Se setean los valores iniciales para X, Y, Z
solver.set_f_params(params)             #Se setean los parametros para la función llamable
"""
#Método más lento
for i in range(len(t)):
    for h in range(len(rr_axis)):
        if t[i] < rr_axis[h]:
            params[-1] = rr_times[h]
            #print(params[-1])
            break
        elif t[i] > rr_axis[h]:
            #break
            continue
    solver.set_f_params(params)        
    solver.integrate(solver.t+dt)
    psoln.append(solver.y)
"""


pos = 0
psoln = []                                                  
for i in range(len(t)):         #Este for permite que cada vez que el parámetro RR, se actualice cada vez que el t alcance al siguiente intervalo RR. Para esto, la serie de intervalos de rr_time se pasaron a una escala temporal en rr_axis. Y cada vez que t alcanza al siguiente rr_axis, rr_time[pos] se actualiza al siguiente valor. 
    if t[i] > rr_axis[pos]:
        params[-2] = rr_times[pos+1]
        solver.set_f_params(params)        
        pos = pos+1    
        
    solver.integrate(solver.t+dt)
    psoln.append(solver.y)

#print(psoln)

"""
######################## ESCALAMIENTO y RUIDO ############
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
#######################GRAFICACÍON CON MATPLOTLIB###################
"""


x_values = np.array(psoln).T[0]
y_values = np.array(psoln).T[1]
z_values = z

#Gráfico 3D (X, Y, Z)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x_values, y_values, z_values)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()


#Gráfico 2D (t, Z)

plt.figure()
plt.plot(t, z_values)
plt.xlabel('time')
plt.ylabel('z')
plt.show()



"""
###################### ANIMACIÓN MATPLOTLIB 3D ##############################
"""

#def update(num, data, line):
#    line.set_data(data[:2, :num])
#    line.set_3d_properties(data[2, :num])





def update(num, data, line, liner, pointer, hrmean, dt):
    
    rr_interval = 60/(hrmean*dt)
    segunda_vuelta = rr_interval*2
    
    
    rr_interval = int(rr_interval)    
    segunda_vuelta = int(segunda_vuelta)
    
    if num > segunda_vuelta :
        line.set_data(data[:2, (num-segunda_vuelta):(num-rr_interval+1)])
        line.set_3d_properties(data[2, (num-segunda_vuelta):(num-rr_interval+1)])
        
        liner.set_data(data[:2, (num-rr_interval):num])
        liner.set_3d_properties(data[2, (num-rr_interval):num])
        
#        pointer.set_data(data[:2, num-1])
#        pointer.set_3d_properties(data[2, num-1])
    else:
        liner.set_data(data[:2, 0:num])
        liner.set_3d_properties(data[2, 0:num])
        
#        pointer.set_data(data[:2, num-1])
#        pointer.set_3d_properties(data[2, num-1])
    

fig = plt.figure()
ax = p3.Axes3D(fig)

data = np.array(psoln).T
data[2] = z_values
line, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1])
liner, = ax.plot(data[0, 0:1], data[1, 0:1], data[2, 0:1], 'r')
pointer, = ax.plot(data[0, 1:1], data[1, 1:1], data[2, 1:1], 'rD')

# Setting the axes properties
ax.set_xlim3d([-1.0, 1.0])
ax.set_xlabel('X')

ax.set_ylim3d([-1.0, 1.0])
ax.set_ylabel('Y')

ax.set_zlim3d([-0.5, 1.5])
ax.set_zlabel('Z')

ani = animation.FuncAnimation(fig, update, frames=len(psoln), fargs=(data, line, liner, pointer, hrmean, dt), interval=1000*dt, blit=False)

plt.show()




"""
###################### ANIMACIÓN MATPLOTLIB 2D ##############################
"""

"""REDISEÑO

OBJETIVOS: 
    - Visualización completa para un número especificado de ciclos y parámetros 
	- Visualización animada 2D actualizabley que responde a interacción
    - Visualización animada 3D actualizable y que responde a interacción


"""


"""
PROBLEMAS A SOLUCIONAR 

+ Preparar Input serie temporal (Tuve que pasar del ODEINT a la clase ODE y hacerlo como 'manual')
+ Preparar Input de parámetros.  
- Velocidad de Animación. Debería ser 60 por min. hasta ahora no lo da. 
+ Abstraer la posición de la sobreposición (El problema era que las abstracciones me estaban quedando en float, pero las listas sólo reconocen int. Había que sólo transformar el tipo)
+ Revisar el parámetro 'n' 
- el tema de la velocidad de sampleo. 
+ Ordenar el tema de la vairable Z que recibe el ruido. Quizás hacer dos psol. psol sólo que sería la señal sin nada y psol_noise que sería la señal con ruido blanco. 
+ Cachar como agregar un ruido de 50 Hz!. 
- Animar ECG 2D! 
- FUTURO cachar como colocar patologías! 
"""
