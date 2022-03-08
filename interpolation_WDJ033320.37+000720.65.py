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
start_B = 3E8 #starting B-field for spline 
end_B = 8.1E8 #end B-field for spline 
target_lambda_1 = 8531
target_lambda_2 = 6640
target_lambda_3 = 6120
#%%
"""
Original interpolation that has wavelength as the horizontal axis and B as the vertical axis
"""
# =============================================================================
# cs = spi.CubicSpline(wavelength_alpha_list[7][start_idx-1:end_idx+1], b_alpha_list[7][start_idx-1:end_idx+1])
#  #find target_lambda from the Gaussian fit to the peak
# print(cs(target_lambda_1)/1E6)
# =============================================================================
#%%
"""
CS for the 13th Halpha transition (lambda=8530)
"""
start_8531 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_8531 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_8531, end_8531)
cs1 = spi.CubicSpline(b_alpha_list[12][start_8531-1:end_8531+1], wavelength_alpha_list[12][start_8531-1:end_8531+1])
target_B_range = np.arange(3E8, 8.1E8, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_8531 = cs1(target_B_range)
B_8531_sq = np.square((B_output_8531-target_lambda_1))
B_8531_chi = (1/target_lambda_1)*(B_8531_sq)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_8531 = int(np.where(B_output_8531 == min(B_output_8531, key=lambda x:abs(x-target_lambda_1)))[0])
B_val_8531 = target_B_range[B_idx_8531]
print(B_idx_8531)
print(B_val_8531/1E6)
idx_range_8531 = np.where((B_output_8531 > target_lambda_1-15) & (B_output_8531 < target_lambda_1+15))
B_filtered_8531 = np.take(target_B_range,idx_range_8531)
B_filtered_8531 = (1E-6)*(B_filtered_8531)
B_sort_8531 = np.sort(B_filtered_8531)
B_sort_8531 = B_sort_8531[0]
#%%
"""
CS for the 12th Halpha transition (lambda=6640)
"""
start_6640_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_6640_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_6640_12, end_6640_12)
cs2_12 = spi.CubicSpline(b_alpha_list[11][start_6640_12-1:end_6640_12+1], wavelength_alpha_list[11][start_6640_12-1:end_6640_12+1])
target_B_range = np.arange(3E8, 8.1E8, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_6640_12 = cs2_12(target_B_range)
B_6640_12_sq = np.square((B_output_6640_12-target_lambda_2))
B_6640_12_chi = (1/target_lambda_2)*(B_6640_12_sq)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_6640_12 = int(np.where(B_output_6640_12 == min(B_output_6640_12, key=lambda x:abs(x-target_lambda_2)))[0])
B_val_6640_12 = target_B_range[B_idx_6640_12]
print(B_idx_6640_12)
print(B_val_6640_12/1E6)
idx_range_6640_12 = np.where((B_output_6640_12 > target_lambda_2-15) & (B_output_6640_12 < target_lambda_2+15))
B_filtered_6640_12 = np.take(target_B_range,idx_range_6640_12)
B_filtered_6640_12 = (1E-6)*(B_filtered_6640_12)
B_sort_6640_12 = np.sort(B_filtered_6640_12)
B_sort_6640_12 = B_sort_6640_12[0]
#%%
"""
CS for the 8th Halpha transition (lambda=6640)
"""
start_6640_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_6640_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_6640_8, end_6640_8)
cs2_8 = spi.CubicSpline(b_alpha_list[7][start_6640_8-1:end_6640_8+1], wavelength_alpha_list[7][start_6640_8-1:end_6640_8+1])
target_B_range = np.arange(3E8, 8.1E8, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_6640_8 = cs2_8(target_B_range)
B_6640_8_sq = np.square((B_output_6640_8-target_lambda_2))
B_6640_8_chi = (1/target_lambda_2)*(B_6640_8_sq)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_6640_8 = int(np.where(B_output_6640_8 == min(B_output_6640_8, key=lambda x:abs(x-target_lambda_2)))[0])
B_val_6640_8 = target_B_range[B_idx_6640_8]
print(B_idx_6640_8)
print(B_val_6640_8/1E6)
idx_range_6640_8 = np.where((B_output_6640_8 > target_lambda_2-15) & (B_output_6640_8 < target_lambda_2+15))
B_filtered_6640_8 = np.take(target_B_range,idx_range_6640_8)
B_filtered_6640_8 = (1E-6)*(B_filtered_6640_8)
B_sort_6640_8 = np.sort(B_filtered_6640_8)
B_sort_6640_8 = B_sort_6640_8[0]
#%%
"""
CS for the 12th Halpha transition (lambda=6120)
"""
start_6120_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_6120_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_6120_12, end_6120_12)
cs3_12 = spi.CubicSpline(b_alpha_list[11][start_6120_12-1:end_6120_12+1], wavelength_alpha_list[11][start_6120_12-1:end_6120_12+1])
target_B_range = np.arange(3E8, 8.1E8, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_6120_12 = cs3_12(target_B_range)
B_6120_12_sq = np.square((B_output_6120_12-target_lambda_3))
B_6120_12_chi = (1/target_lambda_3)*(B_6120_12_sq)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_6120_12 = int(np.where(B_output_6120_12 == min(B_output_6120_12, key=lambda x:abs(x-target_lambda_3)))[0])
B_val_6120_12 = target_B_range[B_idx_6120_12]
print(B_idx_6120_12)
print(B_val_6120_12/1E6)
idx_range_6120_12 = np.where((B_output_6120_12 > target_lambda_3-15) & (B_output_6120_12 < target_lambda_3+15))
B_filtered_6120_12 = np.take(target_B_range,idx_range_6120_12)
B_filtered_6120_12 = (1E-6)*(B_filtered_6120_12)
B_sort_6120_12 = np.sort(B_filtered_6120_12)
B_sort_6120_12 = B_sort_6120_12[0]
#%%
"""
CS for the 8th Halpha transition (lambda=6120)
"""
start_6120_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_6120_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_6120_8, end_6120_8)
cs3_8 = spi.CubicSpline(b_alpha_list[7][start_6120_8-1:end_6120_8+1], wavelength_alpha_list[7][start_6120_8-1:end_6120_8+1])
target_B_range = np.arange(3E8, 8.1E8, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_6120_8 = cs3_8(target_B_range)
B_6120_8_sq = np.square((B_output_6120_8-target_lambda_3))
B_6120_8_chi = (1/target_lambda_3)*(B_6120_8_sq)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_6120_8 = int(np.where(B_output_6120_8 == min(B_output_6120_8, key=lambda x:abs(x-target_lambda_3)))[0])
B_val_6120_8 = target_B_range[B_idx_6120_8]
print(B_idx_6120_8)
print(B_val_6120_8/1E6)
idx_range_6120_8 = np.where((B_output_6120_8 > target_lambda_3-15) & (B_output_6120_8 < target_lambda_3+15))
B_filtered_6120_8 = np.take(target_B_range,idx_range_6120_8)
B_filtered_6120_8 = (1E-6)*(B_filtered_6120_8)
B_sort_6120_8 = np.sort(B_filtered_6120_8)
B_sort_6120_8 = B_sort_6120_8[0]
#%%
"""
Chi-squared calculations
"""
chi_tot = B_8531_chi + B_6640_12_chi + B_6120_8_chi
chi_tot_2 = B_8531_chi + B_6640_8_chi + B_6120_12_chi

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
confidence_range_idx = np.where((chi_tot <= (chi_min + 2.30)))[0]
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
# deriv_8531 = derivative(cs1,8531)
# err_8531 = ((29.92949067)**2)*((deriv_8531)**2)
# print(deriv_8531, err_8531)
# 
# deriv_8531 = derivative(cs1,8531)
# err_8531 = ((21.95414907)**2)*((deriv_8531)**2)
# print(deriv_8531, err_8531)
# 
# # =============================================================================
# # deriv_8531 = derivative(cs1,8531)
# # err_8531 = ((29.92949067)**2)*((deriv_8531)**2)
# # print(deriv_8531, err_8531)
# # 
# # deriv_8531 = derivative(cs1,8531)
# # err_8531 = ((29.92949067)**2)*((deriv_8531)**2)
# # print(deriv_8531, err_8531)
# # 
# # deriv_8531 = derivative(cs1,8531)
# # err_8531 = ((29.92949067)**2)*((deriv_8531)**2)
# # print(deriv_8531, err_8531)
# # =============================================================================
# 
# deriv_8531 = derivative(cs1,8531)
# err_8531 = ((29.92949067)**2)*((deriv_8531)**2)
# print(deriv_8531, err_8531)
# 
# prop_error = np.sqrt(err_8531)
# =============================================================================
