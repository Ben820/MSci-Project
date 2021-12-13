# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 16:49:46 2021

@author: rr3
"""

import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm
import pandas as pd
#%%
""" Part 1: Load data
Notes: data - MWD spectrum
       wavelength - x-values 
       flux - y-values
"""
#load data and sort into appropriate variables
#Ideal system (first to fit)
#DESI_WDJ110344.93+510048.61

#Unfriendly systems 
#DESI_WDJ073615.92+403335.18
#DESI_WDJ082302.41+334534.27
#DESI_WDJ112328.50+095619.35
#DESI_WDJ112926.23+493931.89
#DESI_WDJ165200.65+411029.96 # Quadratic regime 
filename = "DESI_WDJ172329.14+540755.79_bin0p2.dat"
data = np.genfromtxt(f'{filename}', delimiter=' ')
#DESI_WDJ110344.93+510048.61_bin0p2
#DESI_WDJ112328.50+095619.35_bin0p2
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
#%% # THIS IS THE PART SEPARATING THE TWO CELLS ------------------------------
""" Part 2: Performs cuts on the data to isolate the H-alpha region
Notes: start/start_Ha - beginning of cut
       last/last_Ha - end of cut 
       masked_Ha_flux - y-values included within the cut
       masked_Ha_reg - x-values included within the cut
"""
#user input for start and end regions
#start = int(input("Startpoint for H-alpha cut:"))
#last = int(input("Endpoint for H-alpha cut:"))

#define the initial cut region (manual values are inputted 
#into the x:abs(x-XXXXXXXX) section of the np.where command)

start=6200
last=7000
start_Ha = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-start)))[0])
end_Ha = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-last)))[0])

masked_Ha_flux = flux[start_Ha:end_Ha]
masked_Ha_reg = wavelength[start_Ha:end_Ha]

plt.figure("Continuum with absorption feature")
plt.grid()
plt.plot(masked_Ha_reg,masked_Ha_flux,'x')
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Flux", size = "15")
plt.show()
##%%
""" Part 3: Removes the absorption region to isolate the continuum
Notes: abs_start and abs_end are the absorption analogues of start_Ha and end_Ha
       dip_start and dip_end mark the start/end regions of the actual absorption feature
       masked_flux - the y-values with the absorption feature removed
       masked_wavelength - the x-values with the absorption feature removed
"""
#find where the absorption feature starts and ends (manual values are 
#inputted into the x:abs(x-XXXXXXXX) section of the np.where command)
#abs_start = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-6400)))[0])
#abs_end = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-6800)))[0])
#
#abs_flux = flux[abs_start:abs_end]
#abs_reg = wavelength[abs_start:abs_end]
#
#plt.figure("Absoprtion feature region +")
#plt.grid()
#plt.plot(abs_reg,abs_flux)
#plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
#plt.ylabel("Flux", size = "15")
#plt.show()

#find the wavelength values where the feature should be cut (manual values are 
#inputted into the x:abs(x-XXXXXXXX) section of the np.where command)
dip_start = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-6400)))[0])
dip_end = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-6800)))[0])

masked_flux = list(masked_Ha_flux)
masked_wavelength = list(masked_Ha_reg)
#try and remove unncessary features from the spectrum (specific to spectrum,
#COMMENT OUT THIS SECTION IF NOT NEEDED)
del masked_flux[dip_start:dip_end]
del masked_wavelength[dip_start:dip_end]
#del masked_flux[505:517]
#del masked_wavelength[505:517]

plt.figure("Spectra with absorption feature removed")
plt.grid()
plt.plot(masked_wavelength,masked_flux,'x')
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Flux", size = "15")
plt.show()