# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 19:16:05 2022

@author: 44743
"""
import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm
import pickle
from pathlib import Path
import pandas as pd
import astropy.io.fits as fits
from scipy.stats import gaussian_kde
from scipy import stats
from matplotlib import colors
#%%
Bhist = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\B hist csv.csv')
B_values = abs(Bhist[Bhist.columns[1]]).to_numpy() # prev sample (for viva) of 188 DA WDs

B_values = datalist[0] # 192 (ALL DA WDs)
#%%
""" Plot for magnetic field histograms (JUST our data) """

def plot_loghist(x, bins):
    hist, bins = np.histogram(x, bins=bins)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(x, bins=logbins, color = 'darkblue', alpha = 1) # , density = True
#    plt.hist(y, bins=logbins, color = 'crimson', alpha = 0.7) # , density = True
#    plt.legend()
    plt.xscale('log')
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.xlabel('Magnetic Field Strength' + ' (MG)', fontsize = "13")#, labelpad = 10)
    plt.ylabel('Number of MWDs', fontsize = "13")#, labelpad = 10)
    plt.xscale('log')
    plt.show()

plt.figure()
plot_loghist(datalist[0], 17) # bins = 17, 21

#fossillower = np.zeros((1, 100))+10
#fossilupper = np.zeros((1, 100))+1000
#plt.plot(fossillower[0], np.linspace(0, 25 , num = len(fossillower[0])), '--', 'cyan', \
#         linewidth = 2, label = 'Field strength\ncut-off')
#plt.plot(fossilupper[0], np.linspace(0, 25 , num = len(fossilupper[0])), '--', 'cyan', \
#         linewidth = 2, label = 'Field strength\ncut-off')
#plt.axvspan(np.linspace(10, 1000, num = len(fossillower[0])), fossillower[0], fossilupper[0], color = 'red')

#plt.axvline(10, ls = '--', color = 'crimson', lw = 1)
#plt.axvline(1000, ls = '--', color = 'crimson', lw = 1)
#plt.axvspan(10, 1000, color = 'tomato', alpha = 0.2, label = 'Fossil field regime')

plt.axvline(0.1, ls = '--', color = 'green', lw = 1)
plt.axvline(1, ls = '--', color = 'green', lw = 1)
plt.axvspan(0.1, 1, color = 'green', alpha = 0.2, label = 'Crystallisation dynamo')

plt.axvline(1, ls = '--', color = 'firebrick', lw = 1)
plt.axvline(10, ls = '--', color = 'firebrick', lw = 1)
plt.axvspan(1, 10, color = 'firebrick', alpha = 0.2, label = 'Crystallisation and\nrotation dynamo')

plt.axvline(10, ls = '--', color = 'purple', lw = 1)
plt.axvline(100, ls = '--', color = 'purple', lw = 1)
plt.axvspan(10, 100, color = 'purple', alpha = 0.2, label = 'Crystallisation and\nrotation dynamo')


#plt.legend(loc = 'upper right')
plt.xlim(0.1, 1200)

#plt.savefig(f'Bfieldhistcrystshade.pdf', bbox_inches = 'tight')

#%%
""" Overplot for magnetic field histograms (literature and our data) """

def plot_loghist(x, y, bins):
    hist, bins = np.histogram(np.hstack((x,y)), bins=bins)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(x, bins=logbins, color = 'darkblue', alpha = 1, label = 'Ferrario et al. (2020)') # , density = True
    plt.hist(y, bins=logbins, color = 'tomato', alpha = 0.8, label = 'Sample of 192\nDA MWDs (this work)') # , density = True
    plt.legend()
    plt.xscale('log')
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.xlabel('Magnetic Field Strength' + ' (MG)', fontsize = "13")#, labelpad = 10)
    plt.ylabel('Number of MWDs', fontsize = "13")#, labelpad = 10)
    plt.xscale('log')
    plt.show()

plt.figure()
plot_loghist(dataF, datalist[0], 18) # bins = 18

#plt.axvline(10, ls = '--', color = 'crimson', lw = 1)
#plt.axvline(1000, ls = '--', color = 'crimson', lw = 1)
#plt.axvspan(10, 1000, color = 'tomato', alpha = 0.2, label = 'Fossil field regime')
#
#plt.axvline(0.1, ls = '--', color = 'green', lw = 1)
#plt.axvline(1, ls = '--', color = 'green', lw = 1)
#plt.axvspan(0.1, 1, color = 'green', alpha = 0.2, label = 'Crystallisation dynamo')
#
#plt.axvline(1, ls = '--', color = 'blue', lw = 1)
#plt.axvline(10, ls = '--', color = 'blue', lw = 1)
#plt.axvspan(1, 10, color = 'blue', alpha = 0.2, label = 'Crystallisation and\nrotation dynamo')
#
#
#plt.legend(loc = 'upper right')
#plt.savefig(f'Bfieldoverplot 2.pdf', bbox_inches = 'tight')

#%%
""" Statistical Tests; Kolmogorov Smirnov Test """
from scipy import stats

#for j in range(len(g)):
#    KS = stats.ks_2samp(dataFnew, datalist[0])
#    #KS = stats.kstest(LogBin[j][1], Rand_deg_dist())
#    print("KS Test m =", g[j])
#    print(KS)
KS = stats.ks_2samp(massMWDnew, datalist[2])
print(KS)
#%%
""" Plot for mass distribution histograms (JUST our data) """
plt.figure()

plt.hist(datalist[2], bins=25, color='darkblue') # ADC9FF # bins = 25

plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
#plt.xlabel('Mass' + ' (' '$M_{\odot}$' + ')', fontsize = "15", labelpad = 10)
plt.xlabel(r'Mass $ (\mathrm{M_{\odot}}) $', fontsize = "13")#, labelpad = 10)
plt.ylabel('Number of MWDs', fontsize = "13")#, labelpad = 10)

plt.show()

#plt.savefig(f'Masshist.pdf', bbox_inches = 'tight')

#%%
#def plot_hist(x, y, bins):
##    hist, bins = np.histogram(np.hstack((x,y)), bins=bins)
##    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
#    plt.hist(x, bins=bins, color = 'darkblue', alpha = 1, label = 'Ferrario et al. (2020)') # , density = True
#    plt.hist(y, bins=bins, color = 'tomato', alpha = 0.7) # , density = True
#    plt.legend()
##    plt.xscale('log')
#    plt.xticks(fontsize=11)
#    plt.yticks(fontsize=11)
#    plt.xlabel('Magnetic Field Strength' + ' (MG)', fontsize = "13", labelpad = 10)
#    plt.ylabel('Number of MWDs', fontsize = "13", labelpad = 10)
#    plt.xscale('log')
#    plt.show()

MassMwdFerr = np.arange(0.375, 1.375, 0.05)
multiplier = np.array([2,3,2,7,3,11,12,9,19,16,8,11,8,4,13,9,11,9,1,1])
rangeMWDFerr = np.arange(0.35, 1.4, 0.05)
massMWD = []

for i in range(len(multiplier)):
    for j in range(multiplier[i]):
        massMWD.append(MassMwdFerr[i])

massMWDnew = []
for i in range(len(multiplier)):
    for j in range(multiplier[i]):
        massMWDnew.append(np.random.uniform(rangeMWDFerr[i], rangeMWDFerr[i+1]))
#%%
""" (Over)plot for mass distribution histograms (literature and our data) """
Mass = np.array(datalist[2])
Mass = Mass[np.logical_not(np.isnan(Mass))]

massMWDnewA = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\massMWDnew csv.csv',\
                           skiprows = 0)
massMWDnewB = []
[massMWDnewB.append(massMWDnewA['data'][i]) for i in range(0, len(massMWDnewA))]

bins = np.histogram(np.hstack((massMWDnew, Mass)), bins=25)[1] #get the bin edges

plt.figure()

plt.hist(massMWDnew, bins, color = 'tomato', label = 'Ferrario et al.\n(2020)')

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.xlabel(r'Mass $ (\mathrm{M_{\odot}}) $', fontsize = "16")#, labelpad = 10)
plt.ylabel('Number of MWDs', fontsize = "16")#, labelpad = 10)
plt.legend(loc = 'upper left', prop={'size': 15})
#plt.figure()
#plt.hist(Mass, bins)
plt.xscale('linear')

#plt.savefig(f'MasshistFerrario 2.pdf', bbox_inches = 'tight')

##%%
##plt.figure()
#plt.bar(MassMwdFerr, multiplier, width = 0.05, color = 'navy', label = 'Ferrario et al. (2020)')
#plt.hist(MassMwdFerr, 20, color = '#61c9a8')


#%%
plt.figure()
Mass = np.array(datalist[2])
Mass = Mass[np.logical_not(np.isnan(Mass))]
plot_hist(massMWDnew, Mass, 32)
##%%
#def plot_loghist(x, bins):
#  hist, bins = np.histogram(x, bins=bins)
#  logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
#  plt.hist(x, bins=logbins)
#  plt.xscale('log')

#plt.figure()

#plt.bar(MassMwdFerr, multiplier, width = 0.05, color = 'navy', label = 'Ferrario et al. (2020)')
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
def plot_hist(x, bins, colour, label, alpha):
    hist, bins = np.histogram(x, bins=bins)
#    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(x, bins=bins, color = colour, alpha = alpha, label = label) # , density = True
    #  plt.xscale('log')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel(r'Mass $ (\mathrm{M_{\odot}}) $', fontsize = "16")#, labelpad = 10)
    plt.ylabel('Number of WDs', fontsize = "16")#, labelpad = 10)
    plt.legend(prop={'size': 15})

  

plt.figure()
plot_hist(massNonWD, 32, 'darkblue', 'Kilic et al. (2018)', 1) # bins = 17 bins = 32! 
plt.show()
#plt.savefig(f'nonMWDMasshistKilic 2.pdf', bbox_inches = 'tight')
