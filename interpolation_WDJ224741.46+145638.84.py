# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 18:50:28 2022

@author: rr3

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
start_B = 2E8 #starting B-field for spline 
end_B = 4E8 #end B-field for spline 
target_lambda_1 = 8455
target_lambda_2 = 6976
target_lambda_3 = 7128
target_lambda_4 = 5871
target_lambda_5 = 4490
target_lambda_6 = 4141
#%%
"""
CS for the 13th Halpha transition (lambda=8455)
"""
start_8455 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_8455 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_8455, end_8455)
cs1 = spi.CubicSpline(b_alpha_list[12][start_8455-1:end_8455+1], wavelength_alpha_list[12][start_8455-1:end_8455+1])
target_B_range = np.arange(start_B, end_B, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_8455 = cs1(target_B_range)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_8455 = int(np.where(B_output_8455 == min(B_output_8455, key=lambda x:abs(x-target_lambda_1)))[0])
B_val_8455 = target_B_range[B_idx_8455]
print(B_idx_8455)
print(B_val_8455/1E6)
idx_range_8455 = np.where((B_output_8455 > target_lambda_1-15) & (B_output_8455 < target_lambda_1+15))
B_filtered_8455 = np.take(target_B_range,idx_range_8455)
B_filtered_8455 = (1E-6)*(B_filtered_8455)
B_sort_8455 = np.sort(B_filtered_8455)
B_sort_8455 = B_sort_8455[0]
#%%
"""
CS for the 12th Halpha transition (lambda=6976)
"""
start_6976_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_6976_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_6976_12, end_6976_12)
cs2_12 = spi.CubicSpline(b_alpha_list[11][start_6976_12-1:end_6976_12+1], wavelength_alpha_list[11][start_6976_12-1:end_6976_12+1])
target_B_range = np.arange(start_B, end_B, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_6976_12 = cs2_12(target_B_range)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_6976_12 = int(np.where(B_output_6976_12 == min(B_output_6976_12, key=lambda x:abs(x-target_lambda_2)))[0])
B_val_6976_12 = target_B_range[B_idx_6976_12]
print(B_idx_6976_12)
print(B_val_6976_12/1E6)
idx_range_6976_12 = np.where((B_output_6976_12 > target_lambda_2-15) & (B_output_6976_12 < target_lambda_2+15))
B_filtered_6976_12 = np.take(target_B_range,idx_range_6976_12)
B_filtered_6976_12 = (1E-6)*(B_filtered_6976_12)
B_sort_6976_12 = np.sort(B_filtered_6976_12)
B_sort_6976_12 = B_sort_6976_12[0]
#%%
"""
CS for the 11th Halpha transition (lambda=7128)
"""
start_7128_11 = int(np.where(b_alpha_list[10] == min(b_alpha_list[10], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_7128_11 = int(np.where(b_alpha_list[10] == min(b_alpha_list[10], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_7128_11, end_7128_11)
cs1 = spi.CubicSpline(b_alpha_list[10][start_7128_11-1:end_7128_11+1], wavelength_alpha_list[10][start_7128_11-1:end_7128_11+1])
target_B_range = np.arange(start_B, end_B, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_7128_11 = cs1(target_B_range)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_7128_11 = int(np.where(B_output_7128_11 == min(B_output_7128_11, key=lambda x:abs(x-target_lambda_3)))[0])
B_val_7128_11 = target_B_range[B_idx_7128_11]
print(B_idx_7128_11)
print(B_val_7128_11/1E6)
idx_range_7128_11 = np.where((B_output_7128_11 > target_lambda_3-15) & (B_output_7128_11 < target_lambda_3+15))
B_filtered_7128_11 = np.take(target_B_range,idx_range_7128_11)
B_filtered_7128_11 = (1E-6)*(B_filtered_7128_11)
B_sort_7128_11 = np.sort(B_filtered_7128_11)
B_sort_7128_11 = B_sort_7128_11[0]
#%%
"""
CS for the 8th Halpha transition (lambda=5871)
"""
start_5871_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_5871_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_5871_8, end_5871_8)
cs2_8 = spi.CubicSpline(b_alpha_list[7][start_5871_8-1:end_5871_8+1], wavelength_alpha_list[7][start_5871_8-1:end_5871_8+1])
target_B_range = np.arange(start_B, end_B, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_5871_8 = cs2_8(target_B_range)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_5871_8 = int(np.where(B_output_5871_8 == min(B_output_5871_8, key=lambda x:abs(x-target_lambda_4)))[0])
B_val_5871_8 = target_B_range[B_idx_5871_8]
print(B_idx_5871_8)
print(B_val_5871_8/1E6)
idx_range_5871_8 = np.where((B_output_5871_8 > target_lambda_4-15) & (B_output_5871_8 < target_lambda_4+15))
B_filtered_5871_8 = np.take(target_B_range,idx_range_5871_8)
B_filtered_5871_8 = (1E-6)*(B_filtered_5871_8)
B_sort_5871_8 = np.sort(B_filtered_5871_8)
B_sort_5871_8 = B_sort_5871_8[0]
#%%
"""
CS for the 16th Hbeta transition (lambda=4490)
"""
start_4490_16 = int(np.where(b_beta_list[15] == min(b_beta_list[15], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4490_16 = int(np.where(b_beta_list[15] == min(b_beta_list[15], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4490_16, end_4490_16)
cs3_16 = spi.CubicSpline(b_beta_list[15][start_4490_16-1:end_4490_16+1], wavelength_beta_list[15][start_4490_16-1:end_4490_16+1])
target_B_range = np.arange(start_B, end_B, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_4490_16 = cs3_16(target_B_range)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_4490_16 = int(np.where(B_output_4490_16 == min(B_output_4490_16, key=lambda x:abs(x-target_lambda_5)))[0])
B_val_4490_16 = target_B_range[B_idx_4490_16]
print(B_idx_4490_16)
print(B_val_4490_16/1E6)
idx_range_4490_16 = np.where((B_output_4490_16 > target_lambda_5-15) & (B_output_4490_16 < target_lambda_5+15))
B_filtered_4490_16 = np.take(target_B_range,idx_range_4490_16)
B_filtered_4490_16 = (1E-6)*(B_filtered_4490_16)
B_sort_4490_16 = np.sort(B_filtered_4490_16)
B_sort_4490_16 = B_sort_4490_16[0]
#%%
"""
CS for the 11th Hbeta transition (lambda=4141)
"""
start_4141_11 = int(np.where(b_beta_list[10] == min(b_beta_list[10], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4141_11 = int(np.where(b_beta_list[10] == min(b_beta_list[10], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4141_11, end_4141_11)
cs3_8 = spi.CubicSpline(b_beta_list[10][start_4141_11-1:end_4141_11+1], wavelength_beta_list[10][start_4141_11-1:end_4141_11+1])
target_B_range = np.arange(start_B, end_B, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_4141_11 = cs3_8(target_B_range)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_4141_11 = int(np.where(B_output_4141_11 == min(B_output_4141_11, key=lambda x:abs(x-target_lambda_6)))[0])
B_val_4141_11 = target_B_range[B_idx_4141_11]
print(B_idx_4141_11)
print(B_val_4141_11/1E6)
idx_range_4141_11 = np.where((B_output_4141_11 > target_lambda_6-15) & (B_output_4141_11 < target_lambda_6+15))
B_filtered_4141_11 = np.take(target_B_range,idx_range_4141_11)
B_filtered_4141_11 = (1E-6)*(B_filtered_4141_11)
B_sort_4141_11 = np.sort(B_filtered_4141_11)
B_sort_4141_11 = B_sort_4141_11[0]

#%%
"""
Chi-squared calculations
"""
B_8455_sq = np.square((B_output_8455-target_lambda_1))
B_8455_chi = (1/target_lambda_1)*(B_8455_sq)

B_6976_12_sq = np.square((B_output_6976_12-target_lambda_2))
B_6976_12_chi = (1/target_lambda_2)*(B_6976_12_sq)

B_7128_11_sq = np.square((B_output_7128_11-target_lambda_3))
B_7128_11_chi = (1/target_lambda_3)*(B_7128_11_sq)

B_5871_8_sq = np.square((B_output_5871_8-target_lambda_4))
B_5871_8_chi = (1/target_lambda_4)*(B_5871_8_sq)

B_4490_16_sq = np.square((B_output_4490_16-target_lambda_5))
B_4490_16_chi = (1/target_lambda_5)*(B_4490_16_sq)

B_4141_11_sq = np.square((B_output_4141_11-target_lambda_6))
B_4141_11_chi = (1/target_lambda_6)*(B_4141_11_sq)

chi_tot = B_8455_chi + B_6976_12_chi + B_7128_11_chi + B_4490_16_chi + B_4141_11_chi
chi_tot_2 =B_8455_chi + B_6976_12_chi + B_7128_11_chi + B_4490_16_chi + B_4141_11_chi + B_5871_8_chi

chi_min = chi_tot.min()
print(chi_min)
chi_min_idx = int((np.where(chi_tot == chi_tot.min()))[0])

chi_min_2 = chi_tot_2.min()
print(chi_min_2)
chi_min_idx_2 = int((np.where(chi_tot_2 == chi_tot_2.min()))[0])

print(chi_min < chi_min_2)

chi_B_guess = target_B_range[chi_min_idx]
chi_B_guess_2 = target_B_range[chi_min_idx_2]

print(chi_B_guess/1E6)
print(chi_B_guess_2/1E6)
#%%
"""
Chi-squared confidence intervals

Minimum chi-squared + 2.30 gives 68.3% of data
Minimum chi-squared + 6.17 gives 95.45% of data
Minimum chi-squared + 11.8 gives 99.73% of data
"""
confidence_range_idx = np.where((chi_tot <= (chi_min + 1)))[0]
B_chi_range = np.array([target_B_range[i] for i in confidence_range_idx])
B_chi_range = (1E-6)*B_chi_range
B_chi_lower = min(B_chi_range)
B_chi_upper = max(B_chi_range)
print(B_chi_lower, B_chi_upper)
#%%
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
# deriv_8455 = derivative(cs1,8531)
# err_8455 = ((29.92949067)**2)*((deriv_8455)**2)
# print(deriv_8455, err_8455)
# 
# deriv_8455 = derivative(cs1,8531)
# err_8455 = ((21.95414907)**2)*((deriv_8455)**2)
# print(deriv_8455, err_8455)
# 
# # =============================================================================
# # deriv_8455 = derivative(cs1,8531)
# # err_8455 = ((29.92949067)**2)*((deriv_8455)**2)
# # print(deriv_8455, err_8455)
# # 
# # deriv_8455 = derivative(cs1,8531)
# # err_8455 = ((29.92949067)**2)*((deriv_8455)**2)
# # print(deriv_8455, err_8455)
# # 
# # deriv_8455 = derivative(cs1,8531)
# # err_8455 = ((29.92949067)**2)*((deriv_8455)**2)
# # print(deriv_8455, err_8455)
# # =============================================================================
# 
# deriv_8455 = derivative(cs1,8531)
# err_8455 = ((29.92949067)**2)*((deriv_8455)**2)
# print(deriv_8455, err_8455)
# 
# prop_error = np.sqrt(err_8455)
# =============================================================================
