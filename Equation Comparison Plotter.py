# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:28:33 2021

@author: rr3
"""
import numpy as np
import matplotlib.pyplot as plt

#simulated Preston linear split and quadratic splits for central component
B = np.linspace(0,20,50)
B_sq = np.square(B)
zeroline = np.zeros(len(B))
a = 20.11385057
#%%original Preston quadratic coefficient
preston_original = -0.173737153
preston_original_right_quad_1 = (2)*preston_original*B_sq
preston_original_right_quad_2 = (5)*preston_original*B_sq
preston_original_left_quad_1 = (-2)*preston_original*B_sq
preston_original_left_quad_2 = (-5)*preston_original*B_sq
preston_original_centre_quad = (preston_original)*B_sq
plt.figure()
plt.xlabel(r'$\Delta\lambda (\AA)$')
plt.ylabel('B Field Strength (MG)')
plt.plot(zeroline, B, color='k')
plt.plot(preston_original_right_quad_1,B,'r--')
plt.plot(preston_original_right_quad_2,B,'m--')
plt.plot(preston_original_left_quad_1,B,'y--')
plt.plot(preston_original_left_quad_2,B,'c--')
plt.plot(preston_original_centre_quad,B,'b--')
#%%corrected Preston quadratic coefficient
preston_corrected = -1.927477339
preston_corrected_right_quad_1 = (2)*preston_corrected*B_sq
preston_corrected_right_quad_2 = (5)*preston_corrected*B_sq
preston_corrected_left_quad_1 = (-2)*preston_corrected*B_sq
preston_corrected_left_quad_2 = (-5)*preston_corrected*B_sq
preston_corrected_centre_quad = (preston_corrected)*B_sq
plt.figure()
plt.xlabel(r'$\Delta\lambda (\AA)$')
plt.ylabel('B Field Strength (MG)')
plt.plot(zeroline, B, color='k')
plt.plot(preston_corrected_right_quad_1,B,'r--')
plt.plot(preston_corrected_right_quad_2,B,'m--')
plt.plot(preston_corrected_left_quad_1,B,'y--')
plt.plot(preston_corrected_left_quad_2,B,'c--')
plt.plot(preston_corrected_centre_quad,B,'b--')
#%%original Jenkins and Segre quadratic coefficient
JS_original = -0.057912384
JS_original_right_quad_1 = (2)*JS_original*B_sq
JS_original_right_quad_2 = (5)*JS_original*B_sq
JS_original_left_quad_1 = (-2)*JS_original*B_sq
JS_original_left_quad_2 = (-5)*JS_original*B_sq
JS_original_centre_quad = (JS_original)*B_sq
plt.figure()
plt.xlabel(r'$\Delta\lambda (\AA)$')
plt.ylabel('B Field Strength (MG)')
plt.plot(zeroline, B, color='k')
plt.plot(JS_original_right_quad_1,B,'r--')
plt.plot(JS_original_right_quad_2,B,'m--')
plt.plot(JS_original_left_quad_1,B,'y--')
plt.plot(JS_original_left_quad_2,B,'c--')
plt.plot(JS_original_centre_quad,B,'b--')
#%%corrected Jenkins and Segre quadratic coefficient
JS_corrected = -121.4066842
JS_corrected_right_quad_1 = (2)*JS_corrected*B_sq
JS_corrected_right_quad_2 = (5)*JS_corrected*B_sq
JS_corrected_left_quad_1 = (-2)*JS_corrected*B_sq
JS_corrected_left_quad_2 = (-5)*JS_corrected*B_sq
JS_corrected_centre_quad = (JS_corrected)*B_sq
plt.figure()
plt.xlabel(r'$\Delta\lambda (\AA)$')
plt.ylabel('B Field Strength (MG)')
plt.plot(zeroline, B, color='k')
plt.plot(JS_corrected_right_quad_1,B,'r--')
plt.plot(JS_corrected_right_quad_2,B,'m--')
plt.plot(JS_corrected_left_quad_1,B,'y--')
plt.plot(JS_corrected_left_quad_2,B,'c--')
plt.plot(JS_corrected_centre_quad,B,'b--')
#%%modified Jenkins and Segre quadratic coefficient
JS_mod = -6.419204046E-9
JS_mod_right_quad_1 = (2)*JS_mod*B_sq
JS_mod_right_quad_2 = (5)*JS_mod*B_sq
JS_mod_left_quad_1 = (-2)*JS_mod*B_sq
JS_mod_left_quad_2 = (-5)*JS_mod*B_sq
JS_mod_centre_quad = (JS_mod)*B_sq
plt.figure()
plt.xlabel(r'$\Delta\lambda (\AA)$')
plt.ylabel('B Field Strength (MG)')
plt.plot(zeroline, B, color='k')
plt.plot(JS_mod_right_quad_1,B,'r--')
plt.plot(JS_mod_right_quad_2,B,'m--')
plt.plot(JS_mod_left_quad_1,B,'y--')
plt.plot(JS_mod_left_quad_2,B,'c--')
plt.plot(JS_mod_centre_quad,B,'b--')