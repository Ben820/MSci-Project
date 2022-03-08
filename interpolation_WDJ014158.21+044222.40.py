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
start_B = 2.5E8 #starting B-field for spline 
end_B = 8.1E8 #end B-field for spline 
target_lambda_1 = 8540
target_lambda_2 = 6900
target_lambda_3 = 6037
target_lambda_4 = 4460
target_lambda_5 = 4145

#%%
"""
CS for the 13th Halpha transition (lambda=8530)
"""
start_8540 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_8540 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_8540, end_8540)
cs1 = spi.CubicSpline(b_alpha_list[12][start_8540-1:end_8540+1], wavelength_alpha_list[12][start_8540-1:end_8540+1])
target_B_range = np.arange(3E8, 8.1E8, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_8540 = cs1(target_B_range)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_8540 = int(np.where(B_output_8540 == min(B_output_8540, key=lambda x:abs(x-target_lambda_1)))[0])
B_val_8540 = target_B_range[B_idx_8540]
print(B_idx_8540)
print(B_val_8540/1E6)
idx_range_8540 = np.where((B_output_8540 > target_lambda_1-27.5) & (B_output_8540 < target_lambda_1+27.5))
B_filtered_8540 = np.take(target_B_range,idx_range_8540)
B_filtered_8540 = (1E-6)*(B_filtered_8540)
B_sort_8540 = np.sort(B_filtered_8540)
B_sort_8540 = B_sort_8540[0]
#%%
# =============================================================================
# """
# CS for the 12th Halpha transition (lambda=6640)
# """
# start_6640_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
# end_6640_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
# print(start_6640_12, end_6640_12)
# cs2_12 = spi.CubicSpline(b_alpha_list[11][start_6640_12-1:end_6640_12+1], wavelength_alpha_list[11][start_6640_12-1:end_6640_12+1])
# target_B_range = np.arange(3E8, 8.1E8, 1E4) #find target_lambda from the Gaussian fit to the peak
# B_output_6640_12 = cs2_12(target_B_range)
# #plt.scatter(target_B_range, B_output,marker='+')
# B_idx_6640_12 = int(np.where(B_output_6640_12 == min(B_output_6640_12, key=lambda x:abs(x-target_lambda_2)))[0])
# B_val_6640_12 = target_B_range[B_idx_6640_12]
# print(B_idx_6640_12)
# print(B_val_6640_12/1E6)
# idx_range_6640_12 = np.where((B_output_6640_12 > target_lambda_2-15) & (B_output_6640_12 < target_lambda_2+15))
# B_filtered_6640_12 = np.take(target_B_range,idx_range_6640_12)
# B_filtered_6640_12 = (1E-6)*(B_filtered_6640_12)
# B_sort_6640_12 = np.sort(B_filtered_6640_12)
# B_sort_6640_12 = B_sort_6640_12[0]
# =============================================================================
#%%
"""
CS for the 8th Halpha transition (lambda=6037)
"""
start_6037_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_6037_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_6037_8, end_6037_8)
cs2_8 = spi.CubicSpline(b_alpha_list[7][start_6037_8-1:end_6037_8+1], wavelength_alpha_list[7][start_6037_8-1:end_6037_8+1])
target_B_range = np.arange(3E8, 8.1E8, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_6037_8 = cs2_8(target_B_range)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_6037_8 = int(np.where(B_output_6037_8 == min(B_output_6037_8, key=lambda x:abs(x-target_lambda_3)))[0])
B_val_6037_8 = target_B_range[B_idx_6037_8]
print(B_idx_6037_8)
print(B_val_6037_8/1E6)
idx_range_6037_8 = np.where((B_output_6037_8 > target_lambda_3-58) & (B_output_6037_8 < target_lambda_3+58))
B_filtered_6037_8 = np.take(target_B_range,idx_range_6037_8)
B_filtered_6037_8 = (1E-6)*(B_filtered_6037_8)
B_sort_6037_8 = np.sort(B_filtered_6037_8)
B_sort_6037_8 = B_sort_6037_8[0]
#%%
"""
CS for the 16th Hbeta transition (lambda=4460)
"""
start_4460_16 = int(np.where(b_beta_list[15] == min(b_beta_list[15], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4460_16 = int(np.where(b_beta_list[15] == min(b_beta_list[15], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4460_16, end_4460_16)
cs3_16 = spi.CubicSpline(b_beta_list[15][start_4460_16-1:end_4460_16+1], wavelength_beta_list[15][start_4460_16-1:end_4460_16+1])
target_B_range = np.arange(3E8, 8.1E8, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_4460_16 = cs3_16(target_B_range)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_4460_16 = int(np.where(B_output_4460_16 == min(B_output_4460_16, key=lambda x:abs(x-target_lambda_4)))[0])
B_val_4460_16 = target_B_range[B_idx_4460_16]
print(B_idx_4460_16)
print(B_val_4460_16/1E6)
idx_range_4460_16 = np.where((B_output_4460_16 > target_lambda_4-20) & (B_output_4460_16 < target_lambda_4+20))
B_filtered_4460_16 = np.take(target_B_range,idx_range_4460_16)
B_filtered_4460_16 = (1E-6)*(B_filtered_4460_16)
B_sort_4460_16 = np.sort(B_filtered_4460_16)
B_sort_4460_16 = B_sort_4460_16[0]
#%%
"""
CS for the 11th Hbeta transition (lambda=4145)
"""
start_4145_11 = int(np.where(b_beta_list[10] == min(b_beta_list[10], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4145_11 = int(np.where(b_beta_list[10] == min(b_beta_list[10], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4145_11, end_4145_11)
cs3_8 = spi.CubicSpline(b_beta_list[10][start_4145_11-1:end_4145_11+1], wavelength_beta_list[10][start_4145_11-1:end_4145_11+1])
target_B_range = np.arange(3E8, 8.1E8, 1E4) #find target_lambda from the Gaussian fit to the peak
B_output_4145_11 = cs3_8(target_B_range)
#plt.scatter(target_B_range, B_output,marker='+')
B_idx_4145_11 = int(np.where(B_output_4145_11 == min(B_output_4145_11, key=lambda x:abs(x-target_lambda_5)))[0])
B_val_4145_11 = target_B_range[B_idx_4145_11]
print(B_idx_4145_11)
print(B_val_4145_11/1E6)
idx_range_4145_11 = np.where((B_output_4145_11 > target_lambda_5-15) & (B_output_4145_11 < target_lambda_5+15))
B_filtered_4145_11 = np.take(target_B_range,idx_range_4145_11)
B_filtered_4145_11 = (1E-6)*(B_filtered_4145_11)
B_sort_4145_11 = np.sort(B_filtered_4145_11)
B_sort_4145_11 = B_sort_4145_11[0]
#%%
"""
Chi-squared calculations 
"""
B_8540_sq = np.square((B_output_8540-target_lambda_1))
B_8540_chi = (1/target_lambda_1)*(B_8540_sq)

B_6037_8_sq = np.square((B_output_6037_8-target_lambda_3))
B_6037_8_chi = (1/target_lambda_3)*(B_6037_8_sq)

B_4460_16_sq = np.square((B_output_4460_16-target_lambda_4))
B_4460_16_chi = (1/target_lambda_4)*(B_4460_16_sq)

B_4145_11_sq = np.square((B_output_4145_11-target_lambda_5))
B_4145_11_chi = (1/target_lambda_5)*(B_4145_11_sq)

chi_tot = B_8540_chi + B_6037_8_chi + B_4460_16_chi + B_4145_11_chi

chi_min = chi_tot.min()
print(chi_min)
chi_min_idx = int((np.where(chi_tot == chi_tot.min()))[0])

chi_B_guess = target_B_range[chi_min_idx]

print(chi_B_guess/1E6)

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
# deriv_8540 = derivative(cs1,8531)
# err_8540 = ((29.92949067)**2)*((deriv_8540)**2)
# print(deriv_8540, err_8540)
# 
# deriv_8540 = derivative(cs1,8531)
# err_8540 = ((21.95414907)**2)*((deriv_8540)**2)
# print(deriv_8540, err_8540)
# 
# # =============================================================================
# # deriv_8540 = derivative(cs1,8531)
# # err_8540 = ((29.92949067)**2)*((deriv_8540)**2)
# # print(deriv_8540, err_8540)
# # 
# # deriv_8540 = derivative(cs1,8531)
# # err_8540 = ((29.92949067)**2)*((deriv_8540)**2)
# # print(deriv_8540, err_8540)
# # 
# # deriv_8540 = derivative(cs1,8531)
# # err_8540 = ((29.92949067)**2)*((deriv_8540)**2)
# # print(deriv_8540, err_8540)
# # =============================================================================
# 
# deriv_8540 = derivative(cs1,8531)
# err_8540 = ((29.92949067)**2)*((deriv_8540)**2)
# print(deriv_8540, err_8540)
# 
# prop_error = np.sqrt(err_8540)
# =============================================================================
