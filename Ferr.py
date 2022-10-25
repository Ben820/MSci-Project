# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 15:41:26 2022

@author: 44743
"""
#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
Ferr_data = pd.read_csv(r'C:\Users\44743\Downloads\Default Dataset (2).csv')

#%%
plt.figure()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Magnetic Field Strength' + ' (MG)', fontsize = "15", labelpad = 10)
plt.ylabel('Number of MWDs', fontsize = "15", labelpad = 10)
#plt.hist(combined_fields, bins=32, color='tomato')

#hist, bins, _ = plt.hist(combined_fields, bins=32)


x_array = Ferr_data[Ferr_data.columns[0]].to_numpy()
y_array = Ferr_data[Ferr_data.columns[1]].to_numpy()

#hist, bins, _ = plt.hist(y_array, bins=56)

# histogram on log scale. 
# Use non-equal bin sizes, such that they look equal on log scale.
#logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
#plt.subplot(212)
plt.hist(x_array, bins=10, color='tomato')

#plt.xscale('log')
#%%
""" Magnetic field distribution - literature (Ferrario 2020) """
dataFerr = np.array([0.0125, 0.0225, 0.0425, 0.1075, 0.23, 0.425, 1.075, 23, 4.25, 7.75, 13, \
                     22, 42.5, 77.5, 130, 355, 775])
dataFerr = np.array([0.012, 0.02, 0.04, 0.2, 0.3, 0.5, 1.2, 2.2, 4.5, 8, 14, 25, 45, 80, 140, 300, 800]) #20,18
#dataFerr = np.array([0.012, 0.02, 0.04, 0.1, 0.2, 0.4, 1.0, 2, 4, 7, 13, 22, 42, 75, 130, 300, 700]) #17
rangeFerr = [0.01, 0.015, 0.03, 0.055, 0.16, 0.3, 0.55, 1.6, 3, 5.5, 10, 16, 30, 55, 100, 160, 550, 1000]

multplier = np.array([2, 0, 4, 3, 7, 3, 11, 43, 28, 37, 43, 21, 22, 10, 3, 5, 7])
dataF = []
# append this  values x number of times

#OLD WAY (VIVA)
for i in range(len(multplier)):
    for j in range(multplier[i]):
        dataF.append(dataFerr[i])

dataFnew = []
#NEW WAY (BETTER FOR BINNING)
for i in range(len(multplier)):
    for j in range(multplier[i]):
        dataFnew.append(np.random.uniform(rangeFerr[i], rangeFerr[i+1]))
        
##%%
def plot_loghist(x, bins):
  hist, bins = np.histogram(x, bins=bins)
  logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
  plt.hist(x, bins=logbins)
  plt.xscale('log')

plt.figure()
#plot_loghist(np.array(B_values), 10)
plot_loghist(dataF, 18)

#%%
""" Mass distribution for MAGNETIC WDs - Ferrario 2020 literature """
MassMwdFerr = np.arange(0.375, 1.375, 0.05)
multiplier = np.array([2,3,2,7,3,11,12,9,19,16,8,11,8,4,13,9,11,9,1,1])
massMWD = []

for i in range(len(multiplier)):
    for j in range(multiplier[i]):
        massMWD.append(MassMwdFerr[i])
##%%
#def plot_loghist(x, bins):
#  hist, bins = np.histogram(x, bins=bins)
#  logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
#  plt.hist(x, bins=logbins)
#  plt.xscale('log')

plt.figure()

plt.bar(MassMwdFerr, multiplier, width = 0.05, color = 'navy', label = 'Ferrario et al. (2020)')

# THIS PLOTS THE DATA FOR MAGNETIC WD MASS DISTRIBUTION - need to run poster plotting first for masses_array
#plt.hist(masses_array, bins=32, color='tomato', label = 'Data') # ADC9FF

#plot_loghist(np.array(B_values), 10)
#plot_loghist(massMWD, 11)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(r'Mass $ (\mathrm{M_{\odot}}) $', fontsize = "15", labelpad = 10)
plt.ylabel('Number of Magnetic White Dwarfs', fontsize = "15", labelpad = 10)
plt.legend()
#plt.savefig(f'Mass field hist data.png', dpi = 1000, bbox_inches = 'tight')

plt.show()

#%%

""" Mass distribution for NON MAGNETIC WDs - Kilic 2018 literature """

MassNonwdFerr = np.array([0.22, 0.27, 0.3, 0.33, 0.36, 0.39, 0.42, 0.46, 0.50, 0.54, 0.57, 0.61, \
                         0.64, 0.67, 0.69, 0.72, 0.75, 0.79, 0.82, 0.84, 0.87, 0.9, 0.93, 0.96, 0.99, \
                         1.02, 1.05, 1.08, 1.11, 1.14, 1.17, 1.2, 1.25, 1.29, 1.32])
multiplier = np.array([2,4,8,6,5,0,6,2,10,33,76,85,55,31,14,18,15,22,15,17,19,21,12,14,9,5,1,10,4,7])
massNonWD = []

for i in range(len(multiplier)):
    for j in range(multiplier[i]):
        massNonWD.append(MassNonwdFerr[i])
##%%
def plot_loghist(x, bins, colour, label, alpha):
    hist, bins = np.histogram(x, bins=bins)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(x, bins=logbins, color = colour, alpha = alpha, label = label) # , density = True
    #  plt.xscale('log')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel(r'Mass $ (\mathrm{M_{\odot}}) $', fontsize = "15", labelpad = 10)
    plt.ylabel('Number of Non-Magnetic White Dwarfs', fontsize = "14", labelpad = 10)
    plt.legend()

  

plt.figure()
plot_loghist(massNonWD, 21, 'navy', 'Kilic et al. (2018)', 1) # bins = 17

#plt.bar(MassMwdFerr, multiplier, width = 0.05, color = 'navy', label = 'Ferrario et al. (2020)')

#plt.hist(masses_array, bins=32, color='tomato', label = 'Data') # ADC9FF

#plot_loghist(np.array(B_values), 10)
#plot_loghist(massMWD, 11)
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)
#plt.xlabel(r'Mass $ (\mathrm{M_{\odot}}) $', fontsize = "15", labelpad = 10)
#plt.ylabel('Number of Magnetic White Dwarfs', fontsize = "15", labelpad = 10)
#plt.legend()
#plt.savefig(f'non mag Mass hist lit data.png', dpi = 1000, bbox_inches = 'tight')

plt.show()


#%%
















