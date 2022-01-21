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
transitions = [['2p-1','3s0'], ['2s0','3p+1'],['2p+1','3d+2'],['2p0','3d+1'], \
               ['2p-1','3d0'],['2p0','3s0'],['2p-1','3d-1'],['2s0','3p0'], \
               ['2p0','3d0'],['2p+1','3s0'],['2s0','3p-1'],['2p-1','3d-2'], \
               ['2p0','3d-1'],['2p+1','3d0']] #Halpha transitions
alpha_list = []

for trans in transitions:
    alpha_list.append(transitions_df[(transitions_df['Bfield']  == trans[0]) & \
                 (transitions_df['Final State'] == trans[1])].index.tolist())
    #indices of the H-alpha transitions (start of the list)
    
num_list=[]
for alpha in alpha_list:
    num_list.append( np.int(transitions_df.iloc[ alpha[0] ].iloc[0]) )
    #Number of data points in each transition
    
alpha_list = [i[0] for i in alpha_list]
alpha_array = np.array(alpha_list)
num_array = np.array(num_list) 
end_array = alpha_array + num_array #end index of data describing a Halpha transition

wavelength_list = []
b_list = []


for i in range(len(alpha_array)):
    wavelength_list.append((transitions_df.loc[alpha_array[i]+2:end_array[i],'Wavelength']).to_numpy().astype(float))
    b_list.append((transitions_df.loc[alpha_array[i]+2:end_array[i],'Bfield']).to_numpy().astype(float))
    #wavelength and B-field arrays for all the Halpha transitions
    
labels=['a','b','c','d','e','f','g','h','i','j','k','l','m','n']
plt.figure()
plt.grid()
plt.xlabel('Wavelength $[\AA]$')
plt.ylabel('B [G]')
for i in range(len(alpha_array)):
    plt.plot(wavelength_list[i],b_list[i])
plt.legend(labels=labels)
plt.yscale('log')
plt.show()


