# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 17:49:05 2022

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
Gaia_data_216 = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\AlldataDA csv.csv', \
                        skiprows = 0)#, unpack = True)#names = columns)

Gaia_data_216 = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\UROP Year 4\Gaia Data new systems csv.csv', \
                        skiprows = 0)#, unpack = True)#names = columns)

Gaia_data_216 = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\UROP Year 4\AlldataDA - UROP inc.csv', \
                        skiprows = 0)#, unpack = True)#names = columns)

# Can use Alldata but remove all types other than 1 (DA)
bprp_216 = []
absG_216 = []

[bprp_216.append(Gaia_data_216['BPRP'][i]) for i in range(0, len(Gaia_data_216))]
[absG_216.append(Gaia_data_216['G_abs'][i]) for i in range(0, len(Gaia_data_216))]

#mass = []
#bfield = []
#mass_err = []
#bfield_err = []
key = ['Bfield', 'eBfield', 'Mass', 'eMass', 'T_eff', 'eT_eff', 'BPRP', 'G_abs']
datalist = [[] for i in range(0,8)]

for j in range(len(datalist)):
    [datalist[j].append(Gaia_data_216[key[j]][i]) for i in range(0, len(Gaia_data_216))]

for i in range(len(datalist[0])):
    datalist[0][i] = abs(datalist[0][i])

datalistarray = [[] for i in range(0,8)]
for i in range(len(datalist)):
    datalistarray[i] = (np.array(datalist[i]))

#%%
# Calculating sp.stats.pearsonr (cannot do with nans inside)
Mass = np.array(datalist[2])
eMass = np.array(datalist[3])

Bfield = np.array(datalist[0])
eBfield = np.array(datalist[1])
##%%
Bfield = Bfield[np.logical_not(np.isnan(Mass))] # remove nans related to mass nans first
Mass = Mass[np.logical_not(np.isnan(Mass))] # remove mass nans (cannot remove B nans after since mass nans now removed)
Mass = Mass[np.logical_not(np.isnan(Bfield))] # remove mass nans (cannot remove B nans after since mass nans now removed)
Bfield = Bfield[np.logical_not(np.isnan(Bfield))] # remove nans related to mass nans first

eBfield = eBfield[np.logical_not(np.isnan(eMass))] # remove nans related to mass nans first
eMass = eMass[np.logical_not(np.isnan(eMass))] # remove mass nans (cannot remove B nans after since mass nans now removed)
eMass = eMass[np.logical_not(np.isnan(eBfield))] # remove mass nans (cannot remove B nans after since mass nans now removed)
eBfield = eBfield[np.logical_not(np.isnan(eBfield))] # remove nans related to mass nans first

#%%

plt.figure()
plt.errorbar(Mass, Bfield, yerr = eBfield, xerr = eMass, color='darkblue', fmt ='x', alpha = 0.2)#, color='darkblue')#, ms = 1)
plt.plot(Mass, Bfield, 'x', color='darkblue')
plt.yscale('log')
plt.show()

arr2d = np.zeros((len(Bfield),2))
arr2d[:,0] = np.array(Mass)
arr2d[:,1] = np.array(Bfield)
#
A = zip(np.corrcoef(arr2d))

sp.stats.pearsonr(Mass, Bfield)
#%%
""" B field vs Mass """
plt.figure()
plt.errorbar(datalistarray[2], datalistarray[0], yerr = datalistarray[1], xerr = datalistarray[3], \
             color='darkblue', fmt ='x', alpha = 0.2)#, color='darkblue')#, ms = 1)
plt.plot(datalistarray[2], datalistarray[0], 'x', color='darkblue')

plt.xlabel(r'Mass $ (\mathrm{M_{\odot}}) $', size = "13")
plt.ylabel("Magnetic field (MG)", size = "13")
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)

plt.yscale('log')
plt.show()
#%%
def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)

    # now determine nice limits by hand:
#    binwidth = 0.25
#    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
#    lim = (int(xymax/binwidth) + 1) * binwidth
#
#    bins = np.arange(-lim, lim + binwidth, binwidth)
    
    
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation='horizontal')
    ax.set_yscale('log')

#%%
# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005


rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a square Figure
fig = plt.figure(figsize=(8, 8))

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# use the previously defined function
scatter_hist(Mass, Bfield, ax, ax_histx, ax_histy)

plt.show()



#%%
plt.figure()
plt.plot(datalist[4], datalist[6], 'x')
plt.show()


#%%
import statistics
a = np.where(np.array(datalist[2]) < 0)[0]
b = []
bc = []
for i in range(len(a)):
    b.append(datalist[0][a[i]])
    bc.append(datalist[2][a[i]])

b = np.array(b)
bc = np.array(bc)

#print(statistics.median(b))
#print(statistics.median(bc))

#%%
aa = np.where((bc) < 10)[0]
bb = []
for i in range(len(aa)):
    bb.append(bc[aa[i]])


print(statistics.median(bb))








