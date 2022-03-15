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
start_B = 0.5E8 #starting B-field for spline 
end_B = 1.5E8 #end B-field for spline 
target_B_range = np.arange(start_B, end_B, 1E4) #find target_lambda from the Gaussian fit to the peak
target_lambda_1 = 7800
target_lambda_2 = 7443
target_lambda_3 = 7030
target_lambda_4 = 6163
target_lambda_5 = 5087
target_lambda_6 = 4717
target_lambda_7 = 4490
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
CS for the 13th Halpha transition (lambda=7800)
"""
start_7800 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_7800 = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_7800, end_7800)
cs1 = spi.CubicSpline(b_alpha_list[12][start_7800-1:end_7800+1], wavelength_alpha_list[12][start_7800-1:end_7800+1])
B_output_7800 = cs1(target_B_range)
B_7800_sq = np.square((B_output_7800-target_lambda_1))
B_7800_chi = (1/target_lambda_1)*(B_7800_sq)
#%%
"""
CS for the 12th Halpha transition (lambda=7443)
"""
start_7443 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_7443 = int(np.where(b_alpha_list[11] == min(b_alpha_list[11], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_7443, end_7443)
cs2 = spi.CubicSpline(b_alpha_list[11][start_7443-1:end_7443+1], wavelength_alpha_list[11][start_7443-1:end_7443+1])
B_output_7443 = cs2(target_B_range)
B_7443_sq = np.square((B_output_7443-target_lambda_2))
B_7443_chi = (1/target_lambda_2)*(B_7443_sq)
#%%
"""
CS for the 11th Halpha transition (lambda=7030)
"""
start_7030 = int(np.where(b_alpha_list[10] == min(b_alpha_list[10], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_7030 = int(np.where(b_alpha_list[10] == min(b_alpha_list[10], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_7030, end_7030)
cs3 = spi.CubicSpline(b_alpha_list[10][start_7030-1:end_7030+1], wavelength_alpha_list[10][start_7030-1:end_7030+1])
B_output_7030 = cs3(target_B_range)
B_7030_sq = np.square((B_output_7030-target_lambda_3))
B_7030_chi = (1/target_lambda_3)*(B_7030_sq)
#%%
"""
CS for the 9th Halpha transition (lambda=6163)
"""
start_6163 = int(np.where(b_alpha_list[8] == min(b_alpha_list[8], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_6163 = int(np.where(b_alpha_list[8] == min(b_alpha_list[8], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_6163, end_6163)
cs4 = spi.CubicSpline(b_alpha_list[8][start_6163-1:end_6163+1], wavelength_alpha_list[8][start_6163-1:end_6163+1])
B_output_6163 = cs4(target_B_range)
B_6163_sq = np.square((B_output_6163-target_lambda_4))
B_6163_chi = (1/target_lambda_4)*(B_6163_sq)
#%%
"""
CS for the 6th Halpha transition (lambda=5087)
"""
start_5087 = int(np.where(b_alpha_list[5] == min(b_alpha_list[5], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_5087 = int(np.where(b_alpha_list[5] == min(b_alpha_list[5], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_5087, end_5087)
cs5 = spi.CubicSpline(b_alpha_list[5][start_5087-1:end_5087+1], wavelength_alpha_list[5][start_5087-1:end_5087+1])
B_output_5087 = cs5(target_B_range)
B_5087_sq = np.square((B_output_5087-target_lambda_5))
B_5087_chi = (1/target_lambda_5)*(B_5087_sq)
#%%
"""
CS for the 4th Halpha transition (lambda=4717)
"""
start_4717 = int(np.where(b_alpha_list[3] == min(b_alpha_list[3], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4717 = int(np.where(b_alpha_list[3] == min(b_alpha_list[3], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4717, end_4717)
cs6 = spi.CubicSpline(b_alpha_list[3][start_4717-1:end_4717+1], wavelength_alpha_list[3][start_4717-1:end_4717+1])
B_output_4717 = cs6(target_B_range)
B_4717_sq = np.square((B_output_4717-target_lambda_6))
B_4717_chi = (1/target_lambda_6)*(B_4717_sq)
#%%
"""
CS for the 3th Halpha transition (lambda=4490)
"""
start_4490 = int(np.where(b_alpha_list[2] == min(b_alpha_list[2], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_4490 = int(np.where(b_alpha_list[2] == min(b_alpha_list[2], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_4490, end_4490)
cs7 = spi.CubicSpline(b_alpha_list[2][start_4490-1:end_4490+1], wavelength_alpha_list[12][start_4490-1:end_4490+1])
B_output_4490 = cs7(target_B_range)
B_4490_sq = np.square((B_output_4490-target_lambda_7))
B_4490_chi = (1/target_lambda_7)*(B_4490_sq)
#%%
"""
Chi-squared calculations
"""
chi_tot = B_7800_chi + B_7443_chi + B_7030_chi + B_6163_chi + B_5087_chi + B_4490_chi
chi_tot_2 = B_7800_chi + B_7443_chi + B_7030_chi + B_6163_chi + B_5087_chi + B_4717_chi + B_4490_chi
chi_tot_3 = B_7800_chi + B_7030_chi + B_6163_chi + B_5087_chi + B_4717_chi
chi_tot_4 = B_7800_chi + B_6163_chi + B_5087_chi + B_4717_chi
chi_tot_5 = B_7800_chi + + B_7443_chi + B_6163_chi + B_5087_chi + B_4717_chi + B_4490_chi

chi_min = chi_tot.min()
print(chi_min)
chi_min_idx = int((np.where(chi_tot == chi_tot.min()))[0])

chi_min_2 = chi_tot_2.min()
print(chi_min_2)
chi_min_idx_2 = int((np.where(chi_tot_2 == chi_tot_2.min()))[0])

chi_min_3 = chi_tot_3.min()
print(chi_min_3)
chi_min_idx_3 = int((np.where(chi_tot_3 == chi_tot_3.min()))[0])

chi_min_4 = chi_tot_4.min()
print(chi_min_4)
chi_min_idx_4 = int((np.where(chi_tot_4 == chi_tot_4.min()))[0])

chi_min_5 = chi_tot_5.min()
print(chi_min_5)
chi_min_idx_5 = int((np.where(chi_tot_5 == chi_tot_5.min()))[0])

print(chi_min < chi_min_2)
print(chi_min_3 < chi_min_4)

chi_B_guess = target_B_range[chi_min_idx]
chi_B_guess_2 = target_B_range[chi_min_idx_2]
chi_B_guess_3 = target_B_range[chi_min_idx_3]
chi_B_guess_4 = target_B_range[chi_min_idx_4]
chi_B_guess_5 = target_B_range[chi_min_idx_5]

print(chi_B_guess/1E6)
print(chi_B_guess_2/1E6)
print(chi_B_guess_3/1E6)
print(chi_B_guess_4/1E6)
print(chi_B_guess_5/1E6)


#%%
"""
Chi-squared confidence intervals

Minimum chi-squared + 2.30 gives 68.3% of data
Minimum chi-squared + 6.17 gives 95.45% of data
Minimum chi-squared + 11.8 gives 99.73% of data
"""
confidence_range_idx = np.where((chi_tot_4 <= (chi_min_4 + 2.30)))[0]
B_chi_range = np.array([target_B_range[i] for i in confidence_range_idx])
B_chi_range = (1E-6)*B_chi_range
B_chi_lower = min(B_chi_range)
B_chi_upper = max(B_chi_range)
print(B_chi_lower, B_chi_upper)
#%%
"""
Chi-squared plotting
"""
x_arr = np.arange(0, len(chi_tot_4))
chi_con = [chi_tot_4[i] for i in confidence_range_idx]
plt.figure()
plt.grid()
plt.plot(x_arr,chi_tot_4)
plt.plot(confidence_range_idx, chi_con, linewidth=5)
plt.plot(chi_min_idx_4, chi_min_4, 'x', markersize=10, linewidth=20)
plt.show()
