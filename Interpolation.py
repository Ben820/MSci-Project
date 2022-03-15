# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 19:12:32 2022

@author: 44743
"""
"""
Snippets of code to perform individual cubic splines on different transitions
Interpolate over B-field values to find wavelengths and then match with target
wavelength values found from a Gaussian fit 
RUN 'hydrogen_transitions_v2.py' BEFORE RUNNING THIS OTHERWISE VARIABLE NAMES WILL NOT BE DEFINED
"""
import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as spi
from scipy.misc import derivative
import lmfit as lm
import pylab
import glob
import os
import pandas as pd
import re
import astropy.io.fits as fits

#%%
start_B = 3E8 #starting B-field for spline 
end_B = 9E8 #end B-field for spline 
#target_lambda_1 = 7455 # transition 12a 
#target_lambda_2 = 6962 # transition 11a
#target_lambda_3 = 5035 # 16b  
#target_lambda_4 = 5457 # 7a
#target_lambda_5 = 4627 # 5a, 6a, 5b 
# tar_lam = [12a, 11a, 16b, 7a, (5a,6a,5b)]
tar_lam = [12, 11, 7, 5] # list of transition values
tar_lam_beta = 16 # leave for now - gonna be tricky to implement
target_lambda = [7455, 6962, 5457, 4627] # list of wavelengths (corresponding to transitions!)

""" Two different (or multiple different) transitions for a given lambda is okay 
provided both are properly included in tar_lam and target_lambda 
See example below """
# tar_lam = [4, 6, 8, 3, 9, 12]
# target_lambda = [X, Y, A, A, A]
#%%
B_output_list = []
for i in range(len(target_lambda)):
#    if tar_lam[i] >= 0:                # USE IF AKSHAY INCORPORATES H BETA (b_beta_list)
#        index = tar_lam[i]-1
#    else:
#        tar_lam
    index = tar_lam[i]-1
    start_trans = int(np.where(b_alpha_list[index] == min(b_alpha_list[index], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
    end_trans = int(np.where(b_alpha_list[index] == min(b_alpha_list[index], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
#    print(start_trans, end_trans)
    cs1 = spi.CubicSpline(b_alpha_list[index][start_trans-1:end_trans+1], \
                          wavelength_alpha_list[index][start_trans-1:end_trans+1])
    target_B_range = np.arange(0.5E8, 8.1E8, 1E4) #find target_lambda from the Gaussian fit to the peak
    B_output_list.append(cs1(target_B_range))
##%%
"""
Chi-squared calculations
"""
B_sq_list = []
B_chi_list = []
for i in range(len(B_output_list)):
    B_sq = np.square((B_output_list[i]-target_lambda[i]))
    B_sq_list.append(B_sq)
    B_chi_list.append((1/target_lambda[i])*(B_sq))

chi_tot = sum(B_chi_list)
# If there are multiple transitions per wavelength value then need to tailor transitions 
# To tailor chi_tot you need to do this manually - since may be many different transition combinations
#e.g. 
#chi_tot_2 = B_8574_chi + B_6516_8_chi + B_6318_12_chi + B_4517_16_chi + B_4199_11_chi

chi_min = chi_tot.min()
print(chi_min)
chi_min_idx = int((np.where(chi_tot == chi_tot.min()))[0])

#chi_min_2 = chi_tot_2.min()
#print(chi_min_2)
#chi_min_idx_2 = int((np.where(chi_tot_2 == chi_tot_2.min()))[0])

#print(chi_min < chi_min_2)

chi_B_guess = target_B_range[chi_min_idx]
#chi_B_guess_2 = target_B_range[chi_min_idx_2]

print(chi_B_guess/1E6)
#print(chi_B_guess_2/1E6)

#%%
"""
Chi-squared confidence intervals
Minimum chi-squared + 2.30 gives 68.3% of data
Minimum chi-squared + 6.17 gives 95.45% of data
Minimum chi-squared + 11.8 gives 99.73% of data
"""
confidence_range_idx = np.where((chi_tot <= (chi_min + 2.30)))[0]
B_chi_range = np.array([target_B_range[i] for i in confidence_range_idx])
B_chi_range = (1E-6)*B_chi_range
B_chi_lower = min(B_chi_range)
B_chi_upper = max(B_chi_range)
print(B_chi_lower, B_chi_upper)
##%%
"""
Chi-squared plotting
"""
x_arr = np.arange(0, len(chi_tot))
chi_con = [chi_tot[i] for i in confidence_range_idx]
plt.figure()
plt.grid()
plt.plot(x_arr,chi_tot)
plt.plot(confidence_range_idx, chi_con,linewidth=5)
plt.plot(chi_min_idx,chi_min,'x',markersize=10,linewidth=20)
plt.show()



#%%
# =============================================================================
# """
# Error propagation
# """
# 
# deriv_8574 = derivative(cs1,8531)
# err_8574 = ((29.92949067)**2)*((deriv_8574)**2)
# print(deriv_8574, err_8574)
# 
# deriv_8574 = derivative(cs1,8531)
# err_8574 = ((21.95414907)**2)*((deriv_8574)**2)
# print(deriv_8574, err_8574)
# 
# # =============================================================================
# # deriv_8574 = derivative(cs1,8531)
# # err_8574 = ((29.92949067)**2)*((deriv_8574)**2)
# # print(deriv_8574, err_8574)
# # 
# # deriv_8574 = derivative(cs1,8531)
# # err_8574 = ((29.92949067)**2)*((deriv_8574)**2)
# # print(deriv_8574, err_8574)
# # 
# # deriv_8574 = derivative(cs1,8531)
# # err_8574 = ((29.92949067)**2)*((deriv_8574)**2)
# # print(deriv_8574, err_8574)
# # =============================================================================
# 
# deriv_8574 = derivative(cs1,8531)
# err_8574 = ((29.92949067)**2)*((deriv_8574)**2)
# print(deriv_8574, err_8574)
# 
# prop_error = np.sqrt(err_8574)
# =============================================================================