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

#!!! NOTE ONLY Top two and BOTTOM TWO CELLS ARE ACTUALLY USEFUL AT THIS POINT!!! 
#%%
"""
Retrieve Gaia data and retrieve WD data
"""
gaia_data = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Gaia catalogue.csv')
IM_fields_data = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\IM fields input table.csv')
linear_fields_data = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Low-field Linear Systems_poster.csv')
linear_just_fields = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Low-field Linear Systems_poster_2.csv')
high_fields_data = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\HF and extra.csv')
IM_just_fields = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\IM fields input table_2.csv')
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
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.xlabel('Mass' + ' (' '$M_{\odot}$' + ')', fontsize = "15", labelpad = 10)
plt.xlabel(r'Mass $ [\mathrm{M_{\odot}}] $', fontsize = "15", labelpad = 10)

plt.ylabel('Number of MWDs', fontsize = "15", labelpad = 10)
plt.hist(masses_array, bins=32, color='tomato') # ADC9FF

#top=0.97,
#bottom=0.14,
#left=0.155,
#right=0.97,
#hspace=0.2,
#wspace=0.2
#plt.savefig("216massdist_poster", dpi = 1000)

plt.show()
#%%
"""
Cell that plots a histogram of MWD field vs number
"""
linear_x_array = linear_fields_data[linear_fields_data.columns[2]].to_numpy()
IM_x_array = IM_fields_data[IM_fields_data.columns[1]].to_numpy()
HF_x_array = high_fields_data[high_fields_data.columns[1]].to_numpy()
combined_fields = np.concatenate((linear_x_array, IM_x_array, HF_x_array)) 

hist, bins, _ = plt.hist(combined_fields, bins=17)

plt.figure()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Magnetic Field Strength' + ' (MG)', fontsize = "15", labelpad = 10)
plt.ylabel('Number of MWDs', fontsize = "15", labelpad = 10)
#plt.hist(combined_fields, bins=32, color='tomato')

#hist, bins, _ = plt.hist(combined_fields, bins=32)

# histogram on log scale. 
# Use non-equal bin sizes, such that they look equal on log scale.
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
#plt.subplot(212)
plt.hist(combined_fields, bins=logbins, color='tomato')

plt.xscale('log')

#top=0.97,
#bottom=0.14,
#left=0.135,
#right=0.97,
#hspace=0.2,
#wspace=0.2

#plt.savefig("Cropped Final", dpi = 1000)

plt.show()

#%%
linear_names = linear_fields_data[linear_fields_data.columns[0]].to_list()
IM_names = IM_fields_data[IM_fields_data.columns[0]].to_list()
##%%
combined_names = linear_names + IM_names
combined_names = [ '{} '.format(x) for x in combined_names ]
##%%
trial_names = ['DESI_WDJ034513.72-111452.15_bin0p2.dat ', 'DESI_WDJ234338.90+165031.50_bin0p2.dat ']
df_combined = gaia_data[gaia_data['WD name'].isin(combined_names)]
df_trial = gaia_data[gaia_data['WD name'].isin(trial_names)]
##%%
df_notnan = df_combined[df_combined.T_eff.notnull()]
temps = df_combined[df_combined.columns[4]].to_numpy()
names_temps = df_combined[df_combined.columns[0]].to_list()
#names_temps = [ '{} '.format(x) for x in names_temps ]
##%%
temps_sort = np.sort(temps)
names_temps_sort = np.sort(names_temps)
##%%
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
#%%

Bhist = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\B hist csv.csv')
B_values = abs(Bhist[Bhist.columns[1]]).to_numpy()

#hist, bins, _ = plt.hist(B_values, bins=10)

plt.figure()

#plt.hist(combined_fields, bins=32, color='tomato')

#hist, bins, _ = plt.hist(combined_fields, bins=32)

# histogram on log scale. 
# Use non-equal bin sizes, such that they look equal on log scale.
#logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
#plt.subplot(212)
plt.hist(B_values, bins=50, color='tomato')
#%%
#fields_array = abs(B_values.to_numpy())
#%%
plt.figure()
plt.xlabel('MWD Field Strength' + ' (MG)')
plt.ylabel('Number of MWDs')
#plt.hist(fields_array,bins=100)

# histogram on linear scale
plt.subplot(211)
hist, bins, _ = plt.hist(B_values,bins=30)

# histogram on log scale. 
# Use non-equal bin sizes, such that they look equal on log scale.
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
plt.subplot(212)
plt.hist(B_values, bins=logbins)
plt.xscale('log')
plt.show() 
#%%
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Magnetic Field Strength' + ' (MG)', fontsize = "13", labelpad = 10)
plt.ylabel('Number of MWDs', fontsize = "13", labelpad = 10)
plt.xscale('log')
plt.show()

#%%
""" Reads in data for next cell IMPORTANT """
Bhist = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\B hist csv.csv')
B_values = abs(Bhist[Bhist.columns[1]]).to_numpy() # prev sample (for viva) of 188 DA WDs

B_values = datalist[0] # 192 (ALL DA WDs)
#%%
""" (Over)plot for magnetic field histograms (literature and our data) """

def plot_loghist(x, bins, colour, label, alpha):
    hist, bins = np.histogram(x, bins=bins)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(x, bins=logbins, color = colour, alpha = alpha, label = label) # , density = True
    plt.legend()
    plt.xscale('log')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('Magnetic Field Strength' + ' (MG)', fontsize = "13", labelpad = 10)
    plt.ylabel('Number of Magnetic White Dwarfs', fontsize = "13", labelpad = 10)
    plt.xscale('log')
    plt.show()


plt.figure()
#plot_loghist(dataF, 18, 'navy', 'Ferrario et al. (2020)', 1) # bins = 17
plot_loghist(dataFnew, 20, 'navy', 'Ferrario et al. (2020)', 1) # bins = 17
plot_loghist(np.array(B_values), 17, 'tomato', 'Data', 0.8) # bins = 10
#plt.savefig(f'B field hist overplot.png', dpi = 1000, bbox_inches = 'tight')

#plot_loghist(dataF, 17, 'blue', 'Ferrario et al.') # bins = 17

# bins 10, 17 match heights 


