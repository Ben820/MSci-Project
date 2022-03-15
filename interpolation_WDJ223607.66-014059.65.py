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
start_B_1 = 3.5E8 #starting B-field for spline 1
end_B_1 = 4.4E8 #end B-field for spline 1
start_B_2 = 2E8 #starting B-field for spline 2
end_B_2 = 3.25E8 #end B-field for spline 2
target_B_range_1 = np.arange(start_B_1, end_B_1, 1E4) #find target_lambda from the Gaussian fit to the peak
target_B_range_2 = np.arange(start_B_2, end_B_2, 1E4) #find target_lambda from the Gaussian fit to the peak
target_lambda_1 = 8491
target_lambda_2 = 7186
target_lambda_3 = 6892
target_lambda_4 = 4797
target_lambda_5 = 4468
target_lambda_6 = 4140
target_lambda_7 = 5870
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
start_8491 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_8491 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_8491, end_8491)
cs1 = spi.CubicSpline(b_alpha_list[12][start_8491-1:end_8491+1], wavelength_alpha_list[12][start_8491-1:end_8491+1])
B_output_8491 = cs1(target_B_range_1)
B_8491_sq = np.square((B_output_8491-target_lambda_1))
B_8491_chi = (1/target_lambda_1)*(B_8491_sq)
#%%
"""
CS for the 12th Halpha transition (lambda=7186)
"""
start_7186_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_7186_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_7186_12, end_7186_12)
cs2 = spi.CubicSpline(b_alpha_list[11][start_7186_12-1:end_7186_12+1], wavelength_alpha_list[11][start_7186_12-1:end_7186_12+1])
B_output_7186_12 = cs2(target_B_range_2)
B_7186_12_sq = np.square((B_output_7186_12-target_lambda_2))
B_7186_12_chi = (1/target_lambda_2)*(B_7186_12_sq)
#%%
"""
CS for the 11th Halpha transition (lambda=6892)
"""
start_6892_11 = int(np.where(b_alpha_list[10] == min(b_alpha_list[10], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_6892_11 = int(np.where(b_alpha_list[10] == min(b_alpha_list[10], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_6892_11, end_6892_11)
cs3 = spi.CubicSpline(b_alpha_list[10][start_6892_11-1:end_6892_11+1], wavelength_alpha_list[10][start_6892_11-1:end_6892_11+1])
B_output_6892_11 = cs3(target_B_range_2)
B_6892_11_sq = np.square((B_output_6892_11-target_lambda_3))
B_6892_11_chi = (1/target_lambda_3)*(B_6892_11_sq)
#%%
"""
CS for the 12th Halpha transition (lambda=6892)
"""
start_6892_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_6892_12 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_6892_12, end_6892_12)
cs4 = spi.CubicSpline(b_alpha_list[11][start_6892_12-1:end_6892_12+1], wavelength_alpha_list[11][start_6892_12-1:end_6892_12+1])
B_output_6892_12 = cs4(target_B_range_1)
B_6892_12_sq = np.square((B_output_6892_12-target_lambda_3))
B_6892_12_chi = (1/target_lambda_3)*(B_6892_12_sq)
#%%
"""
CS for the 11th Halpha transition (lambda=7186)
"""
start_7186_11 = int(np.where(b_alpha_list[10] == min(b_alpha_list[10], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_7186_11 = int(np.where(b_alpha_list[10] == min(b_alpha_list[10], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_7186_11, end_7186_11)
cs5 = spi.CubicSpline(b_alpha_list[10][start_7186_11-1:end_7186_11+1], wavelength_alpha_list[10][start_7186_11-1:end_7186_11+1])
B_output_7186_11 = cs5(target_B_range_1)
B_7186_11_sq = np.square((B_output_7186_11-target_lambda_2))
B_7186_11_chi = (1/target_lambda_2)*(B_7186_11_sq)
#%%
"""
CS for the 9th Halpha transition (lambda=4468)
"""
start_4468_9 = int(np.where(b_alpha_list[8] == min(b_alpha_list[8], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4468_9 = int(np.where(b_alpha_list[8] == min(b_alpha_list[8], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4468_9, end_4468_9)
cs6 = spi.CubicSpline(b_alpha_list[8][start_4468_9-1:end_4468_9+1], wavelength_alpha_list[8][start_4468_9-1:end_4468_9+1])
B_output_4468_9 = cs6(target_B_range_1)
B_4468_9_sq = np.square((B_output_4468_9-target_lambda_5))
B_4468_9_chi = (1/target_lambda_5)*(B_4468_9_sq)
#%%
"""
CS for the 7th Halpha transition (lambda=4140)
"""
start_4140_7 = int(np.where(b_alpha_list[6] == min(b_alpha_list[6], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4140_7 = int(np.where(b_alpha_list[6] == min(b_alpha_list[6], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4140_7, end_4140_7)
cs7 = spi.CubicSpline(b_alpha_list[6][start_4140_7-1:end_4140_7+1], wavelength_alpha_list[6][start_4140_7-1:end_4140_7+1])
B_output_4140_7 = cs7(target_B_range_1)
B_4140_7_sq = np.square((B_output_4140_7-target_lambda_6))
B_4140_7_chi = (1/target_lambda_6)*(B_4140_7_sq)
#%%
"""
CS for the 9th Halpha transition (lambda=4797)
"""
start_4797_9 = int(np.where(b_alpha_list[8] == min(b_alpha_list[8], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4797_9 = int(np.where(b_alpha_list[8] == min(b_alpha_list[8], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4797_9, end_4797_9)
cs8 = spi.CubicSpline(b_alpha_list[8][start_4797_9-1:end_4797_9+1], wavelength_alpha_list[8][start_4797_9-1:end_4797_9+1])
B_output_4797_9 = cs8(target_B_range_2)
B_4797_9_sq = np.square((B_output_4797_9-target_lambda_4))
B_4797_9_chi = (1/target_lambda_4)*(B_4797_9_sq)
#%%
"""
CS for the 7th Halpha transition (lambda=4468)
"""
start_4468_7 = int(np.where(b_alpha_list[6] == min(b_alpha_list[6], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4468_7 = int(np.where(b_alpha_list[6] == min(b_alpha_list[6], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4468_7, end_4468_7)
cs9 = spi.CubicSpline(b_alpha_list[6][start_4468_7-1:end_4468_7+1], wavelength_alpha_list[6][start_4468_7-1:end_4468_7+1])
B_output_4468_7 = cs9(target_B_range_2)
B_4468_7_sq = np.square((B_output_4468_7-target_lambda_5))
B_4468_7_chi = (1/target_lambda_5)*(B_4468_7_sq)
#%%
"""
CS for the 11th Hbeta transition (lambda=4140)
"""
start_4140_11 = int(np.where(b_beta_list[10] == min(b_beta_list[10], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4140_11 = int(np.where(b_beta_list[10] == min(b_beta_list[10], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4140_11, end_4140_11)
cs10 = spi.CubicSpline(b_beta_list[10][start_4140_11-1:end_4140_11+1], wavelength_beta_list[3][start_4140_11-1:end_4140_11+1])
B_output_4140_11 = cs10(target_B_range_2)
B_4140_11_sq = np.square((B_output_4140_11-target_lambda_6))
B_4140_11_chi = (1/target_lambda_6)*(B_4140_11_sq)
#%%
"""
Chi-squared calculations
"""
chi_tot = B_8491_chi + B_6892_12_chi + B_7186_11_chi + B_4468_9_chi + B_4140_7_chi 
chi_tot_2 = B_7186_12_chi + B_6892_11_chi + B_4797_9_chi + B_4468_7_chi + B_4140_11_chi
chi_tot_3 = B_6892_12_chi + B_7186_11_chi + B_4468_9_chi + B_4140_7_chi 
chi_tot_4 = B_8491_chi + B_4468_9_chi + B_4140_7_chi  

chi_min = chi_tot.min()
red_chi_min = chi_min/5
print(chi_min)
print(red_chi_min)
chi_min_idx = int((np.where(chi_tot == chi_tot.min()))[0])

chi_min_2 = chi_tot_2.min()
red_chi_min_2 = chi_min_2/5
print(chi_min_2)
print(red_chi_min_2)
chi_min_idx_2 = int((np.where(chi_tot_2 == chi_tot_2.min()))[0])

chi_min_3 = chi_tot_3.min()
red_chi_min_3 = chi_min_3/4
print(chi_min_3)
print(red_chi_min_3)
chi_min_idx_3 = int((np.where(chi_tot_3 == chi_tot_3.min()))[0])

chi_min_4 = chi_tot_4.min()
red_chi_min_4 = chi_min_4/3
print(chi_min_4)
print(red_chi_min_4)
chi_min_idx_4 = int((np.where(chi_tot_4 == chi_tot_4.min()))[0])


print(chi_min < chi_min_2)
print(chi_min_2 < chi_min_3)
print(red_chi_min < red_chi_min_2)
print(red_chi_min_2 < red_chi_min_3)
print(chi_min_2 < chi_min_4)
print(red_chi_min_2 < red_chi_min_4)

chi_B_guess = target_B_range_1[chi_min_idx]
chi_B_guess_2 = target_B_range_2[chi_min_idx_2]
chi_B_guess_3 = target_B_range_1[chi_min_idx_3]
chi_B_guess_4 = target_B_range_1[chi_min_idx_4]

print(chi_B_guess/1E6)
print(chi_B_guess_2/1E6)
print(chi_B_guess_3/1E6)
print(chi_B_guess_4/1E6)
#%%
"""
Chi-squared confidence intervals

Minimum chi-squared + 2.30 gives 68.3% of data
Minimum chi-squared + 6.17 gives 95.45% of data
Minimum chi-squared + 11.8 gives 99.73% of data
"""
confidence_range_idx = np.where((chi_tot_2 <= (chi_min_2 + 2.30)))[0]
B_chi_range = np.array([target_B_range_2[i] for i in confidence_range_idx])
B_chi_range = (1E-6)*B_chi_range
B_chi_lower = min(B_chi_range)
B_chi_upper = max(B_chi_range)
print(B_chi_lower, B_chi_upper)
#%%
"""
Chi-squared plotting
"""
x_arr = np.arange(0, len(chi_tot_2))
chi_con = [chi_tot_2[i] for i in confidence_range_idx]
plt.figure()
plt.grid()
plt.plot(x_arr,chi_tot_2)
plt.plot(confidence_range_idx, chi_con, linewidth=5)
plt.plot(chi_min_idx_2, chi_min_2, 'x', markersize=10, linewidth=20)
plt.show()
