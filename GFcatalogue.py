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

data = hdulist[1].data 
data0 = hdulist[0].data

header = hdulist[1].header
header0 = hdulist[0].header

#Extract quantities from the overarching GF catalogue
T_eff = data['Teff'] #WD effective temperature
log_g = data['log_g'] #WD surface gravity
WD_name = data['WDJ_name'] #GF identifier for the WD
G_abs = data['absG'] #the absolute Gaia magnitude for the WD
BPRP = data['bp_rp'] #the difference between the BP and RP filters used when observing
parallax = data['parallax'] #Gaia parallax of source
parallax_err = data['parallax_error'] #error on the parralax of the source
bp_flux = data['phot_bp_mean_flux']
rp_flux = data['phot_rp_mean_flux']
bp_flux_err = data['phot_bp_mean_flux_error']
rp_flux_err = data['phot_rp_mean_flux_error']
#%%
"""
Part 2:
    Compile all the quantities which will be used to characterise the DESI MWDs into one Pandas dataframe
    Process the data to remove NaN values (wherever there is a NaN for T_eff there is a NaN for log_g)
    Tweak the formatting of the data to make it more suitable for our analysis
"""
names = ["G_abs", "BP-RP", "Parallax", "T_eff", "log_g"]
info_df = pd.DataFrame(np.array([G_abs, BPRP, parallax, T_eff, log_g]).transpose(), WD_name, columns=names)

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
#%%
plt.
plt.scatter(BPRP_sel,G_sel,s=0.1)
#%%
plt.figure()
plt.gca().invert_yaxis()
#plt.gca().invert_xaxis()
plt.plot(BPRP_sel,G_sel,'x', markersize=0.4)
#%%
#i=WDJ_col[WDJ_col[1]=='WDJ005212.26+135302.04 ']
#print(i)
#print(type(header))
#print(header)