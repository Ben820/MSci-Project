# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 14:58:25 2022

@author: rr3
"""
import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as spi
import lmfit as lm
import pylab
import glob
import os
import pandas as pd
import re
import astropy.io.fits as fits

#%%
column_names = ['NumDataPoints', 'Bfield', 'Wavelength', 'Final State']
transitions_df = pd.read_csv(r'C:\Local\Akshay Laptop Backup 29Dec2019\Imperial\MSci Project\DESI\hydrogencolumns.csv', skiprows=1, names=column_names)
notnan = transitions_df['NumDataPoints'].notnull()
#%%
transitions_alpha = [['2p-1','3s0'], ['2s0','3p+1'],['2p+1','3d+2'],['2p0','3d+1'], \
               ['2p-1','3d0'],['2p0','3s0'],['2p-1','3d-1'],['2s0','3p0'], \
               ['2p0','3d0'],['2p+1','3s0'],['2s0','3p-1'],['2p-1','3d-2'], \
               ['2p0','3d-1'],['2p+1','3d0']] #H-alpha transitions
alpha_list = []

for trans in transitions_alpha:
    alpha_list.append(transitions_df[(transitions_df['Bfield']  == trans[0]) & \
                 (transitions_df['Final State'] == trans[1])].index.tolist())
    #indices of the H-alpha transitions (start of the list)
    
num_list=[]
for alpha in alpha_list:
    num_list.append( np.int(transitions_df.iloc[ alpha[0] ].iloc[0]) )
    #Number of data points in each transition
    
alpha_list = [i[0] for i in alpha_list]
alpha_array = np.array(alpha_list)
num_alpha_array = np.array(num_list) 
end_alpha_array = alpha_array + num_alpha_array #end index of data describing a Halpha transition

wavelength_alpha_list = []
b_alpha_list = []


for i in range(len(alpha_array)):
    wavelength_alpha_list.append((transitions_df.loc[alpha_array[i]+2:end_alpha_array[i],'Wavelength']).to_numpy().astype(float))
    b_alpha_list.append((transitions_df.loc[alpha_array[i]+2:end_alpha_array[i],'Bfield']).to_numpy().astype(float))
    #wavelength and B-field arrays for all the Halpha transitions
#%%
transitions_beta = [['2p-1','4s0'],['2s0','4p+1'],['2p+1','4d+2'],['2p0','4d+1'], \
               ['2s0','4f+1'],['2p-1','4d0'],['2p0','4s0'],['2p-1','4d-1'],['2s0','4p0'], \
               ['2p0','4d0'],['2s0','4f0'],['2p+1','4s0'],['2s0','4p-1'],['2p-1','4d-2'], \
               ['2p0','4d-1'],['2s0','4f-1'],['2p+1','3d0']] #H-beta transitions
beta_list = []

for trans in transitions_beta:
    beta_list.append(transitions_df[(transitions_df['Bfield']  == trans[0]) & \
                 (transitions_df['Final State'] == trans[1])].index.tolist())
    #indices of the H-beta transitions (start of the list)
    
num_beta_list=[]
for beta in beta_list:
    num_beta_list.append( np.int(transitions_df.iloc[ beta[0] ].iloc[0]) )
    #Number of data points in each transition
    
beta_list = [i[0] for i in beta_list]
beta_array = np.array(beta_list)
num_beta_array = np.array(num_beta_list) 
end_beta_array = beta_array + num_beta_array #end index of data describing a H-beta transition

wavelength_beta_list = []
b_beta_list = []


for i in range(len(beta_array)):
    wavelength_beta_list.append((transitions_df.loc[beta_array[i]+2:end_beta_array[i],'Wavelength']).to_numpy().astype(float))
    b_beta_list.append((transitions_df.loc[beta_array[i]+2:end_beta_array[i],'Bfield']).to_numpy().astype(float))
    #wavelength and B-field arrays for all the H-beta transitions
#%%   
labels_alpha = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n']
labels_beta = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17']
labels_combined = labels_alpha + labels_beta   

plt.figure()
plt.grid()
plt.xlabel('Wavelength $[\AA]$')
plt.ylabel('B [G]')
for i in range(len(alpha_array)):
    plt.plot(wavelength_alpha_list[i],b_alpha_list[i], '--')
for i in range(len(alpha_array)):
    plt.plot(wavelength_alpha_list[i],b_alpha_list[i],'x')
# =============================================================================
# for i in range(len(beta_array)):
#     plt.plot(wavelength_beta_list[i],b_beta_list[i])
# for i in range(len(beta_array)):
#     plt.plot(wavelength_beta_list[i],b_beta_list[i],'x')
# =============================================================================
    
plt.legend(labels=labels_alpha)
#plt.legend(labels=labels_beta)
#plt.legend(labels=labels_combined)
#plt.yscale('log')
plt.show()

#%%
plt.figure()
plt.grid()
plt.xlabel('B [G]')
plt.ylabel('Wavelength $[\AA]$')
for i in range(len(alpha_array)):
    plt.plot(b_alpha_list[i],wavelength_alpha_list[i],'--')
for i in range(len(alpha_array)):
    plt.plot(b_alpha_list[i],wavelength_alpha_list[i],'x')
# =============================================================================
# for i in range(len(beta_array)):
#     plt.plot(b_beta_list[i],wavelength_beta_list[i])
# for i in range(len(beta_array)):
#     plt.plot(b_beta_list[i],wavelength_beta_list[i],'x')
# =============================================================================
    
plt.legend(labels=labels_alpha)
#plt.legend(labels=labels_beta)
#plt.legend(labels=labels_combined)
#plt.yscale('log')
plt.show()

#%%
# =============================================================================
# start_B = 3E8 #starting B-field for spline 
# end_B = 8.1E8 #end B-field for spline 
# start_idx = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-start_B)))[0])
# end_idx = int(np.where(b_alpha_list[12] == min(b_alpha_list[12], key=lambda x:abs(x-end_B)))[0])
# print(start_idx, end_idx)
# target_lambda_1 = 8530
# target_lambda_2 = 6640
# target_lambda_3 = 6120
# =============================================================================
#%%
# =============================================================================
# cs = spi.CubicSpline(wavelength_alpha_list[7][start_idx-1:end_idx+1], b_alpha_list[7][start_idx-1:end_idx+1])
#  #find target_lambda from the Gaussian fit to the peak
# print(cs(target_lambda)/1E6)
# =============================================================================
#%%
# =============================================================================
# """
# CS for the 13th transition (lambda=8530)
# """
# cs1 = spi.CubicSpline(b_alpha_list[12][start_idx-1:end_idx+1], wavelength_alpha_list[12][start_idx-1:end_idx+1])
# target_B_range = np.arange(3E8, 8.1E8, 1E5) #find target_lambda from the Gaussian fit to the peak
# B_output = cs2(target_B_range)
# plt.scatter(target_B_range, B_output,marker='+')
# B_idx = int(np.where(B_output == min(B_output, key=lambda x:abs(x-target_lambda_1)))[0])
# B_val = target_B_range[B_idx]
# print(B_idx)
# print(B_val/1E6)
# #%%
# """
# CS for the 13th transition
# """
# cs1 = spi.CubicSpline(b_alpha_list[12][start_idx-1:end_idx+1], wavelength_alpha_list[12][start_idx-1:end_idx+1])
# target_B_range = np.arange(3E8, 8.1E8, 1E5) #find target_lambda from the Gaussian fit to the peak
# B_output = cs2(target_B_range)
# plt.scatter(target_B_range, B_output,marker='+')
# B_idx = int(np.where(B_output == min(B_output, key=lambda x:abs(x-target_lambda_1)))[0])
# B_val = target_B_range[B_idx]
# print(B_idx)
# print(B_val/1E6)
# =============================================================================
#%%
# =============================================================================
# wavelength_flip = []
# b_flip = []
# 
# for i in range(len(wavelength_alpha_list)):
#     wavelength_flip.append(np.flip(wavelength_alpha_list[i]))
#     b_flip.append(np.flip(b_alpha_list[i]))
# start_B = 3.5E8 #starting B-field for spline --> HIGHER VALUE 
# end_B = 7.5E8 #end B-field for spline --> LOWER VALUE
# start_idx = int(np.where(b_flip[0] == min(b_flip[0], key=lambda x:abs(x-start_B)))[0])
# end_idx = int(np.where(b_flip[0] == min(b_flip[0], key=lambda x:abs(x-end_B)))[0])
# print(start_idx, end_idx)
# plt.figure()
# plt.grid()
# for i in range(len(alpha_array)):
#     plt.plot(wavelength_flip[i],b_flip[i])
# plt.legend(labels=labels)
# plt.show()
# #%%    
# cs2 = spi.CubicSpline(wavelength_alpha_list[12][start_idx-3:end_idx+3], b_alpha_list[12][start_idx-3:end_idx+3])
# target_lambda = 8530 #find target_lambda from the Gaussian fit to the peak
# print(cs2(target_lambda))
# =============================================================================
#%%
# =============================================================================
# plt.figure()
# plt.plot(wavelength_alpha_list[12][start_idx-1:end_idx+1], b_alpha_list[12][start_idx-1:end_idx+1],'x')
# #%%
# idx_range = np.where((B_output > 8526) & (B_output < 8536))
# #min_idx = np.min(idx_range)
# #max_idx = np.max(idx_range)
# #print(min_idx,max_idx)
# B_filtered = np.take(target_B_range,idx_range)
# B_filtered_MG = (1E-6)*(B_filtered)
# B_sort = np.sort(B_filtered_MG)
# B_sort = B_sort[0]
# 
# =============================================================================
