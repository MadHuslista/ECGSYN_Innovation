#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 21:53:37 2019

@author: madhuslista
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# radius of the circle
R = 1

def circle(a, b, r):
    # (a,b): the center of the circle
    # r: the radius of the circle
    # T: The number of the segments
    T = 100
    x, y = [0]*T, [0]*T
    for i,theta in enumerate(np.linspace(0,2*np.pi,T)):
        x[i] = a + r*np.cos(theta)
        y[i] = b + r*np.sin(theta)
    return x, y

# Calculate the cycloid line
thetas = np.linspace(0,4*np.pi,100)
cycloid_x = R*(thetas-np.sin(thetas))
cycloid_y = R*(1-np.cos(thetas))
cycloid_c = R*thetas

fig = plt.figure()

lns = []
#trans = plt.axes().transAxes
for i in range(len(thetas)):
    x,y = circle(cycloid_c[i], R, R)
#    ln1, = plt.plot(x, y, 'g-', lw=2)
    ln2, = plt.plot(cycloid_x[:i+1] ,cycloid_y[:i+1], 'r-', lw=2)
#    ln3, = plt.plot(cycloid_x[i], cycloid_y[i], 'bo', markersize=4)
#    ln4, = plt.plot([cycloid_c[i], cycloid_x[i]], [R,cycloid_y[i]], 'y-', lw=2)
#    tx1  = plt.text(0.05, 0.8, r'$\theta$ = %.2f $\pi$' % (thetas[i]/np.pi), transform=trans)
#    lns.append([ln1,ln2,ln3,ln4,tx1])
    lns.append([ln2])
plt.xlim(0,15)
plt.ylim(0,3)
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.axes().set_aspect('equal')

ani = animation.ArtistAnimation(fig, lns, interval=2)