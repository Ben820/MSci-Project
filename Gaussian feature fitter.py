# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 18:09:09 2022
@author: rr3
"""

import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm
import pylab
import glob
import pandas as pd 
#%%
""" Part 1: Load data
Notes: data - MWD spectrum
       wavelength - x-values 
       flux - y-values
"""
#load data and sort into appropriate variables
filename = "DESI_WDJ224741.46+145638.84_bin0p2.dat"
filename = "DESI_WDJ023420.63+264801.74_bin0p2.dat"
filename = "DESI_WDJ012930.35+155106.95_bin0p2.dat"
filename = "DESI_WDJ034513.72-111452.15_bin0p2.dat"


data = np.genfromtxt(f'{filename}', delimiter=' ')

wavelength = data[:,0]
flux = data[:,1]
error = data[:,2]

plt.figure()
plt.errorbar(wavelength,flux, yerr = error ,label = f"{filename}", fmt ='')#,linewidth=1)
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Flux", size = "15")

#plt.xlim(3660, 9300)
#plt.ylim(0,140)
plt.grid()
plt.legend()
plt.show()
##%%
#scaled_flux = 330E6*norm_spectra
##wavelength[:] = [(wavelength-100) for number in wavelength]
#plt.figure()
#plt.plot(masked_Ha_reg, scaled_flux, linewidth=1)
#
#plt.figure()
#plt.plot(wavelength_alpha_list[12],b_alpha_list[12],'--',linewidth=2.5)
#plt.plot(wavelength_alpha_list[11],b_alpha_list[11],'--',linewidth=2.5)
#plt.plot(wavelength_alpha_list[10],b_alpha_list[10],'--',linewidth=2.5)
#plt.plot(wavelength_alpha_list[7],b_alpha_list[7],'--',linewidth=2.5)
#plt.plot(wavelength_alpha_list[12],b_alpha_list[12],'x')
#plt.plot(wavelength_alpha_list[11],b_alpha_list[11],'x')
#plt.plot(wavelength_alpha_list[10],b_alpha_list[10],'x')
#plt.plot(wavelength_alpha_list[7],b_alpha_list[7],'x')
#plt.plot(wavelength_beta_list[15],b_beta_list[15],'--',linewidth=2.5)
#plt.plot(wavelength_beta_list[10],b_beta_list[10],'--',linewidth=2.5)
#plt.plot(wavelength_beta_list[15],b_beta_list[15],'x')
#plt.plot(wavelength_beta_list[10],b_beta_list[10],'x')
#plt.plot(masked_Ha_reg, scaled_flux, linewidth=1,color='darkblue')

#%%
""" Part 2: Performs cuts on the data to isolate the H-alpha region
Notes: start/start_Ha - beginning of cut
       last/last_Ha - end of cut 
       masked_Ha_flux - y-values included within the cut
       Gauss_Ha_reg - x-values included within the cut
begin/ finish define the whole region including the triplet feature 
startx/endx define the specific region to be cut out (the absorption feature) """

begin_cut = 6650
finish_cut = 7250

start_Gauss = 6800
end_Gauss = 7065

p0_Gauss = [6900,50,3] # p0 for Gaussian 


start_cut = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-begin_cut)))[0])
end_cut = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-finish_cut)))[0])

Gauss_Ha_flux = flux[start_cut:end_cut]
Gauss_Ha_reg = wavelength[start_cut:end_cut]
Gauss_Ha_err = error[start_cut:end_cut]

dip_start = int(np.where(Gauss_Ha_reg == min(Gauss_Ha_reg, key=lambda x:abs(x-start_Gauss)))[0])
dip_end = int(np.where(Gauss_Ha_reg == min(Gauss_Ha_reg, key=lambda x:abs(x-end_Gauss)))[0])

masked_flux = list(Gauss_Ha_flux)
masked_wavelength = list(Gauss_Ha_reg)
masked_err = list(Gauss_Ha_err)

##%%
flxframe = pd.DataFrame(masked_flux)
flxframe.loc[dip_start:dip_end] = np.nan
flxframe = (flxframe.dropna()).to_numpy()

wavframe = pd.DataFrame(masked_wavelength)
wavframe.loc[dip_start:dip_end] = np.nan
wavframe = (wavframe.dropna()).to_numpy()

flxframe = np.array(flxframe)
wavframe = np.array(wavframe)
flux_list_Gauss = flxframe[:,0]
wav_list_Gauss = wavframe[:,0]

# =============================================================================
# plt.figure()
# plt.plot(wav_list_Gauss,flux_list_Gauss,'x')
# plt.show()
# =============================================================================

""" Part 4: Fits a polynomial to the continuum and normalises the spectrum
Notes: Poly_3o is a third-order polynomial function which fits the continuum using a least-squares method
      
"""

###%%
#plt.figure()
#plt.plot(wav_list,flux_list,'x')
#plt.show()
##%%
def Poly_3o(x, a, b, c, d):
    y = a*x**3 + b*x**2 + c*x + d
    return y

p0 = np.array([1, 1, 1, 1]) #fit parameters
p, cov = opt.curve_fit(Poly_3o, wav_list_Gauss,flux_list_Gauss, p0) # do not change
Continuum = Poly_3o(Gauss_Ha_reg, p[0], p[1], p[2], p[3]) # use Gauss_Ha_reg since more data points in array

#for c in zip(p, np.sqrt(np.diag(cov))):# zips root of diag of cov matrix with related value in curve fit
#    print('%.15f pm %.4g' % (c[0], c[1]))# prints value and uncertainty, f is decimal places and G is sig figs

Gauss_spectra = Gauss_Ha_flux/Continuum

#plt.figure()
#plt.grid()
##plt.yticks([0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30])
#plt.plot(Gauss_Ha_reg, Poly_3o(Gauss_Ha_reg, p[0], p[1], p[2], p[3])/Continuum, \
#      zorder=4,color = 'red', label = "Poly")
#plt.plot(wav_list_Gauss,flux_list_Gauss,'x')
#plt.plot(Gauss_Ha_reg, Continuum)
#plt.show()
##%%
# =============================================================================
# plt.figure()
# plt.plot(Gauss_Ha_reg, Gauss_spectra, label = f"{filename}") #plot the normalised spectrum
# plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
# plt.ylabel("Normalised Flux", size = "15")
# plt.legend()
# plt.show()
# =============================================================================
##%%
# =============================================================================
# plt.figure()
# plt.grid()
# plt.plot(masked_wavelength,masked_flux)
# =============================================================================
##%%
# =============================================================================
#plt.figure()
#plt.plot(Gauss_Ha_reg, Gauss_spectra, label = f"{filename}") #plot the normalised spectrum
#plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
#plt.ylabel("Normalised Flux", size = "15")
#plt.legend()
#plt.show()
# =============================================================================
##%%
""" Gaussian Profile 
Notes: 
"""
xp_feature = []
yp_feature = []
err_feature = []
for a in range(len(Gauss_Ha_reg)):
    if Gauss_Ha_reg[a] > (start_Gauss) and Gauss_Ha_reg[a] < (end_Gauss):
        xp_feature.append(Gauss_Ha_reg[a])
        yp_feature.append(Gauss_spectra[a])
        err_feature.append(masked_err[a])
xp_feature = np.array(xp_feature)
yp_feature = np.array(yp_feature)
err_feature = np.array(err_feature)

def Gaussian(x,mu,sig,A):
    gaus = (-(A)*((1/np.sqrt((2*sp.pi)*sig))*(sp.exp(-(x-mu)**2/(2*sig**2)))))+1
    return gaus

#p0 = [7454,100,3]
#p0 = [6980,100,3] # p0 for Gaussian 


popt_Gauss, cov_Gauss = opt.curve_fit(Gaussian, xp_feature, yp_feature, p0_Gauss, sigma = err_feature)

for c in zip(popt_Gauss, np.sqrt(np.diag(cov_Gauss))):
    print("%.8f pm %.3g" % (c[0], c[1]))
##%%
plt.figure()
plt.grid()
plt.plot(xp_feature,yp_feature)
plt.plot(xp_feature,Gaussian(xp_feature,*popt_Gauss))
plt.show()
##%%
G_array = [popt_Gauss[0],popt_Gauss[1],popt_Gauss[2]]
print(G_array)




