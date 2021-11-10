# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 22:35:36 2021

@author: rr3
"""

import numpy as np
import matplotlib.pyplot as plt

plot_data = np.loadtxt('delta_lamba_vs_B_trial_values.csv', delimiter=',', skiprows=1)
#actual data
delta_left = plot_data[:,0]
delta_right = plot_data[:,1]
leftB = plot_data[:,2]
rightB = plot_data[:,3]

#simulated linear data
left_x = np.linspace(-300,0,1000)
right_x = np.linspace(0,300,1000)
left_y = ((-1)*(left_x))/20.2
right_y = (right_x)/20.2

#plot actual data
delta_left = (-1)*delta_left
plt.plot(delta_left,leftB,'x')
plt.plot(delta_right,rightB,'o')
#plot simualted data
plt.plot(left_x,left_y,'y--')
plt.plot(right_x,right_y,'c--')
plt.xlabel('delta_lambda (Angstroms)')
plt.ylabel('B Field Strength (MG)')
plt.grid()
plt.show()

