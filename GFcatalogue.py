# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 10:12:16 2021

@author: rr3
"""
"""
This script works with the Gentile Fusilo (GF) catalogue, extracting quantities for our DESI analysis,
and processing the data so we can present values such as effective temperature and WD surface gravity
to accompany our determined B-field values for the DESI MWDs
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

#%%
"""
Part 1:
    Open and read the GF FITS file
    Extract quantities from the FITS file
"""
hdulist = fits.open("{:s}final_catalogue_v6.fits".format("C:\Local\Akshay Laptop Backup 29Dec2019\Imperial\MSci Project\\"))

GFdata = hdulist[1].data 
GFdata0 = hdulist[0].data

header = hdulist[1].header
header0 = hdulist[0].header

#Extract quantities from the overarching GF catalogue
T_eff = GFdata['Teff'] #WD effective temperature
log_g = GFdata['log_g'] #WD surface gravity
WD_name = GFdata['WDJ_name'] #GF identifier for the WD
G_abs = GFdata['absG'] #the absolute Gaia magnitude for the WD
BPRP = GFdata['bp_rp'] #the difference between the BP and RP filters used when observing
parallax = GFdata['parallax'] #Gaia parallax of source
parallax_err = GFdata['parallax_error'] #error on the parralax of the source
SN = GFdata['S/N'] #signal to noise
bp_flux = GFdata['phot_bp_mean_flux']
rp_flux = GFdata['phot_rp_mean_flux']
bp_flux_err = GFdata['phot_bp_mean_flux_error']
rp_flux_err = GFdata['phot_rp_mean_flux_error']

#%%
"""
Part 2:
    Compile all the quantities which will be used to characterise the DESI MWDs into one Pandas dataframe
    Process the data to remove NaN values (wherever there is a NaN for T_eff there is a NaN for log_g)
    Tweak the formatting of the data to make it more suitable for our analysis
"""
names = ["G_abs", "BP-RP", "Parallax", "T_eff", "log_g", "SN"]
info_df = pd.DataFrame(np.array([G_abs, BPRP, parallax, T_eff, log_g, SN]).transpose(), WD_name, columns=names)

nans = info_df[info_df['T_eff'].isna()]

dataset = info_df.dropna().reset_index()
#%%
dataset = dataset.rename(columns={"index": "WDJ Name"}) #relabel the relic column header from the original dataframe
#dataset.columns = dataset.columns.str.rstrip() #remove all white spaces in the strings
#%%
info2 = ['WDJ014202.20-003332.01 ', 
'WDJ075804.43+205957.86 '
]

info3 = ['WDJ014202.20-003332.01', 
'WDJ075804.43+205957.86'
]
#
df2 = dataset[dataset['WDJ Name'].isin(info2)]
df3 = dataset[dataset['WDJ Name'].isin(info3)]
print(df2)
print(df3)

#%%
"""
Part 3: Plotting script for MWD HRD 
"""
dist_pc = np.reciprocal(parallax/1000)
idx = np.where(dist_pc < 100)[0]
#%%
dist2 = [dist_pc[i] for i in idx]
#dist2 = dist_pc.copy()
#%%
BPRP_sel = [BPRP[i] for i in idx]
G_sel = [G_abs[i] for i in idx]
BPRP_sel2 = BPRP_sel.copy()
BPRP_neg = (-1)*np.array(BPRP_sel2)
##%%
#plt.scatter(BPRP_sel,G_sel,s=0.1)
#%%
plt.figure()
plt.gca().invert_yaxis()
plt.xlabel('BP-RP')
plt.ylabel('G_abs')
#plt.gca().invert_xaxis()
plt.plot(BPRP_sel,G_sel,'o', markersize=0.25)
#%%
filename = "DESI_WDJ172329.14+540755.79_bin0p2.dat"
data = np.genfromtxt(f'{filename}', delimiter=' ')
wavelength = data[:,0]
flux = data[:,1]
error = data[:,2]
plt.figure("Whole spectrum")
plt.plot(wavelength,flux, label = f"{filename}")
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Flux", size = "15")
plt.grid()
plt.legend()
plt.show()
#%%
info_trial = ['WDJ172329.14+540755.79 ']
df_trial = dataset[dataset['WDJ Name'].isin(info_trial)]
print(df_trial)
#%%
T_val = df_trial.iloc[0]['T_eff']
log_g_val = df_trial.iloc[0]['log_g']
G_abs_val = df_trial.iloc[0]['G_abs']
parallax_val = df_trial.iloc[0]['Parallax']
SN_val = df_trial.iloc[0]['SN']
BPRP_val = df_trial.iloc[0]['BP-RP']
name_val = info_trial[0]
print(T_val)
#%%
textstr = "\n".join((
    'Name: %s' % (name_val),
    'G = %2f' % (G_abs_val, ),
    'S/N = %2f' % (SN_val, ),
    'Parallax = %2f' % (parallax_val, ),
    'T_eff = %2f' % (T_val, ),
    'log_g = %2f' % (log_g_val, )))
#%%
f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
a0.set_xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
a0.set_ylabel("Flux", size = "15")
a0.plot(wavelength, flux, label = f"{filename}")
a1.invert_yaxis()
a1.set_xlabel('BP-RP')
a1.set_ylabel('G_abs')
a1.plot(BPRP_sel,G_sel,'o', markersize=0.25)
a1.plot(BPRP_val,G_abs_val,'o',color='red',markersize=5)
props = dict(boxstyle='round', alpha=0.5)
a1.text(0.05, 0.05, textstr, transform=a1.transAxes, fontsize=9,
        verticalalignment='bottom', bbox=props)
#f.tight_layout()
