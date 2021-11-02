# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 13:31:08 2021

@author: rr3

Directories: 
    C:\Local\Akshay Laptop Backup 29Dec2019\Imperial\MSci Project\DESI\DESI_DxH\DESI_DxH
    
Method: 
    Part 1 - Loads data
    Part 2 - Performs cuts on the data to isolate the H-alpha region 
    Part 3 - Removes the absorption region to isolate the continuum
    Part 4 - Fits a polynomial to the continuum and normalises the spectrum
    
"""

""" Part 1: Load data
Notes: data - MWD spectrum
       wavelength - x-values 
       flux - y-values
"""

import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

#laod data and sort into appropriate variables
data = np.genfromtxt('DESI_WDJ110344.93+510048.61_bin0p2.dat', delimiter=' ')
wavelength = data[:,0]
flux = data[:,1]
error = data[:,2]

plt.grid()
plt.plot(wavelength,flux,'x')
#%%
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
start_Ha = np.where(wavelength == min(wavelength, key=lambda x:abs(x-start)))
end_Ha = np.where(wavelength == min(wavelength, key=lambda x:abs(x-last)))
start_Ha = start_Ha[0] 
end_Ha = end_Ha[0]
start_Ha = int(start_Ha)
end_Ha = int(end_Ha)

masked_Ha_flux = flux[start_Ha:end_Ha]
masked_Ha_reg = wavelength[start_Ha:end_Ha]


#%%
plt.figure()
plt.grid()
plt.plot(masked_Ha_reg,masked_Ha_flux,'x')
#%%
""" Part 3: Removes the absorption region to isolate the continuum
Notes: abs_start and abs_end are the absorption analogues of start_Ha and end_Ha
       dip_start and dip_end mark the start/end regions of the actual absorption feature
       masked_flux - the y-values with the absorption feature removed
       masked_wavelength - the x-values with the absorption feature removed
"""
#find where the absorption feature starts and ends (manual values are 
#inputted into the x:abs(x-XXXXXXXX) section of the np.where command)
abs_start = np.where(wavelength == min(wavelength, key=lambda x:abs(x-6400)))
abs_end = np.where(wavelength == min(wavelength, key=lambda x:abs(x-6800)))
abs_start = abs_start[0]
abs_end = abs_end[0]
abs_start = int(abs_start)
abs_end = int(abs_end)
#%%
abs_flux = flux[abs_start:abs_end]
abs_reg = wavelength[abs_start:abs_end]
plt.figure()
plt.grid()
plt.plot(abs_reg,abs_flux)
#%%
#find the wavelength values where the feature should be cut (manual values are 
#inputted into the x:abs(x-XXXXXXXX) section of the np.where command)
dip_start = np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-6500)))
dip_end = np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-6650)))
#%%
dip_start = dip_start[0]
dip_end = dip_end[0]
dip_start = int(dip_start)
dip_end = int(dip_end)
#%%
masked_flux = list(masked_Ha_flux)
masked_wavelength = list(masked_Ha_reg)
#try and remove unncessary features from the spectrum (specific to spectrum,
#COMMENT OUT THIS SECTION IF NOT NEEDED)
del masked_flux[dip_start:dip_end]
del masked_wavelength[dip_start:dip_end]
del masked_flux[505:517]
del masked_wavelength[505:517]
plt.figure()
plt.grid()
plt.plot(masked_wavelength,masked_flux,'x')
#%%
#avg_start_array = masked_flux[0:10]
#avg_end_array = masked_flux[-11:-1]
#avg_start_val = np.mean(masked_flux)
#avg_end_val = np.mean(masked_flux)
#%%
""" Part 4: Fits a polynomial to the continuum and normalises the spectrum
Notes: Poly_3o is a third-order polynomial function which fits the continuum using a least-squares method
      
"""
def Poly_3o(x, a, b, c, d):
    y = a*x**3 + b*x**2 + c*x + d
    return y

x_array = masked_Ha_reg
y_array = masked_Ha_flux

p0 = np.array([1, 1, 1, 1]) #fit parameters
p, cov = opt.curve_fit(Poly_3o, masked_wavelength, masked_flux, p0) # do not change
new_y = Poly_3o(x_array, p[0], p[1], p[2], p[3])
plt.figure()
plt.grid()
#plt.yticks([0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30])
plt.plot(x_array, Poly_3o(x_array, p[0], p[1], p[2], p[3])/new_y, zorder=4,color = 'red', label = "Poly")
plt.plot(x_array, y_array/new_y) #plot the normalised spectrum