# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 16:49:50 2021

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

#regime 3
begin = 5600
finish = 7600
start1 = 5800
end1 = 6200
start2 = 6400
end2 = 6600
start3 = 6800
end3 = 7200


start_Ha = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-begin)))[0])
end_Ha = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-finish)))[0])
masked_Ha_flux = flux[start_Ha:end_Ha]
masked_Ha_reg = wavelength[start_Ha:end_Ha]
dip_start1 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-start1)))[0])
dip_end1 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-end1)))[0])
dip_start2 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-start2)))[0])
dip_end2 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-end2)))[0])
dip_start3 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-start3)))[0])
dip_end3 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-end3)))[0])
masked_flux = list(masked_Ha_flux)
masked_wavelength = list(masked_Ha_reg)


#del masked_flux_3[dip_start1_3:dip_end1_3, dip_start2_3:dip_end2_3, dip_start3_3:dip_end3_3]
#del masked_wavelength_3[dip_start1_3:dip_end1_3, dip_start2_3:dip_end2_3, dip_start3_3:dip_end3_3]

#%%
flxframe = pd.DataFrame(masked_flux)
flxframe.loc[dip_start1:dip_end1] = np.nan
flxframe.loc[dip_start2:dip_end2] = np.nan
flxframe.loc[dip_start3:dip_end3] = np.nan
fframe = flxframe.copy()
ff = fframe.dropna()
ff_arr = ff.to_numpy()
wavframe = pd.DataFrame(masked_wavelength)
wavframe.loc[dip_start1:dip_end1] = np.nan
wavframe.loc[dip_start2:dip_end2] = np.nan
wavframe.loc[dip_start3:dip_end3] = np.nan
wframe = wavframe.copy()
wf = wframe.dropna()
wf_arr = wf.to_numpy()
""" Part 4: Fits a polynomial to the continuum and normalises the spectrum
Notes: Poly_3o is a third-order polynomial function which fits the continuum using a least-squares method
      
"""
ff_arr = np.array(ff_arr)
wf_arr = np.array(wf_arr)
flux_list = ff_arr[:,0]
wav_list = wf_arr[:,0]
#%%
plt.figure()
plt.plot(wav_list,flux_list,'x')
#%%
def Poly_3o(x, a, b, c, d):
    y = a*x**3 + b*x**2 + c*x + d
    return y

p0 = np.array([1, 1, 1, 1]) #fit parameters
p, cov = opt.curve_fit(Poly_3o, wav_list,flux_list, p0) # do not change
Continuum = Poly_3o(masked_Ha_reg, p[0], p[1], p[2], p[3]) # use masked_Ha_reg since more data points in array

#for c in zip(p, np.sqrt(np.diag(cov))):# zips root of diag of cov matrix with related value in curve fit
#    print('%.15f pm %.4g' % (c[0], c[1]))# prints value and uncertainty, f is decimal places and G is sig figs

norm_spectra = masked_Ha_flux/Continuum

plt.figure("Normalised absorption feature")
plt.grid()
#plt.yticks([0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30])
#plt.plot(masked_Ha_reg, Poly_3o(masked_Ha_reg, p[0], p[1], p[2], p[3])/Continuum, \
#         zorder=4,color = 'red', label = "Poly")
plt.plot(wav_list,flux_list,'x')
plt.plot(masked_Ha_reg, Continuum)
#%%
plt.figure()
plt.plot(masked_Ha_reg, norm_spectra, label = f"{filename}") #plot the normalised spectrum
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Normalised Flux", size = "15")
plt.legend()
plt.show()