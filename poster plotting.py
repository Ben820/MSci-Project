# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 19:03:39 2022

@author: rr3
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
"""
Retrieve Gaia data and retrieve WD data
"""
gaia_data = pd.read_csv('Gaia catalogue.csv')
IM_fields_data = pd.read_csv('IM fields input table.csv')
linear_fields_data = pd.read_csv('Low-field Linear Systems_poster.csv')
linear_just_fields = pd.read_csv('Low-field Linear Systems_poster_2.csv')
high_fields_data = pd.read_csv('HF and extra.csv')
IM_just_fields = pd.read_csv('IM fields input table_2.csv')
frames = [linear_just_fields, IM_just_fields, high_fields_data]
df_concat = pd.concat(frames)
#%%
"""
Cell that plots a histogram of MWD mass vs MWD population
"""
masses = gaia_data[gaia_data.columns[1]]
masses = masses.dropna()
masses_array = masses.to_numpy()
plt.figure()
plt.xlabel('MWD Mass' + ' (' '$M_{\odot}$' + ')')
plt.ylabel('Number of MWDs')
plt.hist(masses_array, bins=42, color='slategray')
#%%
"""
Cell that plots a histogram of MWD field vs number
"""
linear_x_array = linear_fields_data[linear_fields_data.columns[2]].to_numpy()
IM_x_array = IM_fields_data[IM_fields_data.columns[1]].to_numpy()
HF_x_array = high_fields_data[high_fields_data.columns[1]].to_numpy()
combined_fields = np.concatenate((linear_x_array, IM_x_array, HF_x_array)) 
plt.figure()
plt.xlabel('MWD Field Strength' + ' (MG)')
plt.ylabel('Number of MWDs')
plt.hist(combined_fields, bins=24, color='slategrey')
#%%
linear_names = linear_fields_data[linear_fields_data.columns[0]].to_list()
IM_names = IM_fields_data[IM_fields_data.columns[0]].to_list()
#%%
combined_names = linear_names + IM_names
combined_names = [ '{} '.format(x) for x in combined_names ]
#%%
trial_names = ['DESI_WDJ034513.72-111452.15_bin0p2.dat ', 'DESI_WDJ234338.90+165031.50_bin0p2.dat ']
df_combined = gaia_data[gaia_data['WD name'].isin(combined_names)]
df_trial = gaia_data[gaia_data['WD name'].isin(trial_names)]
#%%
df_notnan = df_combined[df_combined.T_eff.notnull()]
temps = df_combined[df_combined.columns[4]].to_numpy()
names_temps = df_combined[df_combined.columns[0]].to_list()
#names_temps = [ '{} '.format(x) for x in names_temps ]
#%%
temps_sort = np.sort(temps)
names_temps_sort = np.sort(names_temps)
#%%
df_combined_sort = df_combined.sort_values('WD name')
df_concat_sort = df_concat.sort_values('WD name')
combined_temps_sort = df_combined_sort[df_combined_sort.columns[4]].to_numpy()
combined_fields_sort = df_concat_sort[df_concat_sort.columns[1]].to_numpy()
#%%
plt.figure()
plt.grid()
plt.scatter(combined_fields_sort, combined_temps_sort)
plt.xlabel('MWD Field Strength' + ' (MG)')
plt.ylabel('Effective Temperature' + ' (K)')
plt.xscale('log')
plt.yscale('log')
plt.show()