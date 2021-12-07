# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 16:51:17 2021

@author: rr3
"""

import numpy as np
import matplotlib.pyplot as plt

#%%
#simulated linear data
#left_x = np.linspace(-300,0,1000)
#right_x = np.linspace(0,300,1000)
#left_y = ((-1)*(left_x))/20.2
#right_y = (right_x)/20.2
#plt.figure()
#plt.plot(left_x,left_y,'y--')
#plt.plot(right_x,right_y,'c--')
#plt.xlabel('delta_lambda (Angstroms)')
#plt.ylabel('B Field Strength (MG)')
#plt.grid()
#plt.show()
#%%
#simulated Preston linear split and quadratic splits for right component
B = np.linspace(0,1E6,50000)
a = 2.011385057E-25/((1E-10)**2)
left_linear = ((-a)*(B))
right_linear = a*B
#plt.figure()
#plt.plot(left_linear,B,'y--')
#plt.plot(right_linear,B,'k--')
#plt.grid()
q = 1.737366235E-33/((1E-10)**2)
ml1 = 1
ml2 = 2
B_sq = np.square(B)
right_quad_1 = (-2)*q*B_sq
right_quad_2 = (-5)*q*B_sq
left_quad_1 = (2)*q*B_sq
left_quad_2 = (5)*q*B_sq
centre_quad = q*B_sq
#plt.plot(right_quad_1,B,'r--')
#plt.plot(right_quad_2,B,'m--')
#plt.plot(left_quad_1,B,'c--')
#plt.plot(left_quad_2,B,'y--')
#plt.xlabel('delta_lambda (m)')
#plt.ylabel('B Field Strength (G)')
#plt.grid()
#plt.show()
plt.figure()
plt.plot(right_linear,B,'k')
plt.plot(left_linear,B,'k')
plt.axvline(0, color='k')
plt.plot(right_linear+right_quad_1,B,'r--')
plt.plot(right_linear+right_quad_2,B,'m--')
plt.plot(right_linear+left_quad_1,B,'y--')
plt.plot(right_linear+left_quad_2,B,'c--')
plt.plot(right_linear+centre_quad,B,'b--')
plt.plot(right_quad_1,B,'r--')
plt.plot(right_quad_2,B,'m--')
plt.plot(left_quad_1,B,'y--')
plt.plot(left_quad_2,B,'c--')
plt.plot(centre_quad,B,'b--')
plt.plot(left_linear+right_quad_1,B,'r--')
plt.plot(left_linear+right_quad_2,B,'m--')
plt.plot(left_linear+left_quad_1,B,'y--')
plt.plot(left_linear+left_quad_2,B,'c--')
plt.plot(left_linear+centre_quad,B,'b--')
plt.xlabel(r'$\Delta\lambda (\AA)$')
plt.ylabel('B Field Strength (G)')
plt.grid()
plt.show()
