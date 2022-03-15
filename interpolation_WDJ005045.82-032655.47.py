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
start_B = 2E8 #starting B-field for spline 1
end_B = 5E8 #end B-field for spline 1
target_B_range = np.arange(start_B, end_B, 1E5) #find target_lambda from the Gaussian fit to the peak
target_lambda_1 = 8560
target_lambda_2 = 6978
target_lambda_3 = 5946
target_lambda_4 = 4145
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
CS for the 13th Halpha transition (lambda=8491)
"""
start_8560_13 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_8560_13 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_8560_13, end_8560_13)
cs1 = spi.CubicSpline(b_alpha_list[12][start_8560_13-1:end_8560_13+1], wavelength_alpha_list[12][start_8560_13-1:end_8560_13+1])
B_output_8560_13 = cs1(target_B_range)
B_8560_13_sq = np.square((B_output_8560_13-target_lambda_1))
B_8560_13_chi = (1/target_lambda_1)*(B_8560_13_sq)
#%%
"""
CS for the 12th Halpha transition (lambda=6978)
"""
start_6978_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_6978_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_6978_12, end_6978_12)
cs2 = spi.CubicSpline(b_alpha_list[11][start_6978_12-1:end_6978_12+1], wavelength_alpha_list[11][start_6978_12-1:end_6978_12+1])
B_output_6978_12 = cs2(target_B_range)
B_6978_12_sq = np.square((B_output_6978_12-target_lambda_2))
B_6978_12_chi = (1/target_lambda_2)*(B_6978_12_sq)
#%%
"""
CS for the 8th Halpha transition (lambda=5946)
"""
start_5946_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_5946_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_5946_8, end_5946_8)
cs3 = spi.CubicSpline(b_alpha_list[7][start_5946_8-1:end_5946_8+1], wavelength_alpha_list[7][start_5946_8-1:end_5946_8+1])
B_output_5946_8 = cs3(target_B_range)
B_5946_8_sq = np.square((B_output_5946_8-target_lambda_3))
B_5946_8_chi = (1/target_lambda_3)*(B_5946_8_sq)

#%%
"""
CS for the 11th Hbeta transition (lambda=4145)
"""
start_4145_11 = int(np.where(b_beta_list[10] == min(b_beta_list[10], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4145_11 = int(np.where(b_beta_list[10] == min(b_beta_list[10], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4145_11, end_4145_11)
cs4 = spi.CubicSpline(b_beta_list[10][start_4145_11-1:end_4145_11+1], wavelength_beta_list[10][start_4145_11-1:end_4145_11+1])
B_output_4145_11 = cs4(target_B_range)
B_4145_11_sq = np.square((B_output_4145_11-target_lambda_4))
B_4145_11_chi = (1/target_lambda_4)*(B_4145_11_sq)
#%%
"""
Chi-squared calculations
"""
chi_tot = B_8560_13_chi + B_6978_12_chi + B_5946_8_chi + B_4145_11_chi 
chi_tot_2 = B_8560_13_chi + B_5946_8_chi + B_4145_11_chi 
chi_tot_3 = B_6978_12_chi + B_5946_8_chi + B_4145_11_chi 

chi_min = chi_tot.min()
red_chi_min = chi_min/4
print(chi_min)
print(red_chi_min)
chi_min_idx = int((np.where(chi_tot == chi_tot.min()))[0])

chi_min_2 = chi_tot_2.min()
red_chi_min_2 = chi_min_2/3
print(chi_min_2)
print(red_chi_min_2)
chi_min_idx_2 = int((np.where(chi_tot_2 == chi_tot_2.min()))[0])

chi_min_3 = chi_tot_3.min()
red_chi_min_3 = chi_min_3/3
print(chi_min_3)
print(red_chi_min_3)
chi_min_idx_3 = int((np.where(chi_tot_3 == chi_tot_3.min()))[0])
 
print(chi_min < chi_min_2)
print(chi_min_2 < chi_min_3)
print(red_chi_min < red_chi_min_2)
print(red_chi_min_2 < red_chi_min_3)
# =============================================================================
# chi_min_4 = chi_tot_4.min()
# red_chi_min_4 = chi_min_4/3
# print(chi_min_4)
# print(red_chi_min_4)
# chi_min_idx_4 = int((np.where(chi_tot_4 == chi_tot_4.min()))[0])
# 
# =============================================================================

# =============================================================================
# print(chi_min_2 < chi_min_4)
# print(red_chi_min_2 < red_chi_min_4)
# =============================================================================

chi_B_guess = target_B_range[chi_min_idx]
chi_B_guess_2 = target_B_range[chi_min_idx_2]
chi_B_guess_3 = target_B_range[chi_min_idx_3]

# =============================================================================
# chi_B_guess_4 = target_B_range_1[chi_min_idx_4]
# =============================================================================

print(chi_B_guess/1E6)
print(chi_B_guess_2/1E6)
print(chi_B_guess_3/1E6)
# =============================================================================
# print(chi_B_guess_4/1E6)
# =============================================================================
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
plt.plot(confidence_range_idx, chi_con, linewidth=5)
plt.plot(chi_min_idx, chi_min, 'x', markersize=10, linewidth=20)
plt.show()
