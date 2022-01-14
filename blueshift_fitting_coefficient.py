# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 14:18:46 2022

@author: rr3
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 22:03:21 2021

@author: rr3
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 13:31:08 2021
Directories: 
    Users\44743\Documents\Imperial Year 4\MSci Project\DESI_DxH\DESI_DxH
    
Method: 
    Part 1 - Loads data
    Part 2 - Performs cuts on the data to isolate the H-alpha region 
    Part 3 - Removes the absorption region to isolate the continuum
    Part 4 - Fits a polynomial to the continuum and normalises the spectrum
    
"""
import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pandas as pd
import lmfit as lm
import pickle
from pathlib import Path

#%%
""" Part 1: Load data
Notes: data - MWD spectrum
       wavelength - x-values 
       flux - y-values
"""
#load data and sort into appropriate variables
SourceID = 'DESI_WDJ165200.65+411029.96'
ParamID = 'Params_'+SourceID
CovID = 'Cov_'+SourceID
ErrorID = 'Errors_'+SourceID
file = SourceID+'_bin0p2'
filename = file+'.dat'
imgfile = file+'.png'
resfile = 'Residuals_'+imgfile
#%%
data = np.genfromtxt(filename, delimiter=' ')
#DESI_WDJ110344.93+510048.61_bin0p2
#DESI_WDJ112328.50+095619.35_bin0p2
wavelength = data[:,0]
flux = data[:,1]
error = data[:,2]

plt.figure("Whole spectrum")
plt.plot(wavelength,flux)
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Flux", size = "15")
plt.grid()
plt.show()
#%%
"""
New normalisation code that cuts out each individual absorption feature
"""

begin = 5500
finish = 7500
start1 = 5800
end1 = 6200
start2 = 6450
end2 = 6600
start3 = 6850
end3 = 7150


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
ff_arr = np.array(ff_arr)
wf_arr = np.array(wf_arr)
flux_list = ff_arr[:,0]
wav_list = wf_arr[:,0]
plt.figure()
plt.plot(wav_list,flux_list,'x')

#%%
""" Part 4: Fits a polynomial to the continuum and normalises the spectrum
Notes: Poly_3o is a third-order polynomial function which fits the continuum using a least-squares method
      
"""
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
plt.figure()
plt.plot(masked_Ha_reg, norm_spectra, label = f"{filename}") #plot the normalised spectrum
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Normalised Flux", size = "15")
plt.legend()
plt.show()
#%%
""" Triplet Voigt Profile 
Notes: 
"""
xp_triplet = []
yp_triplet =[]
for a in range(len(masked_Ha_reg)):
    if masked_Ha_reg[a] > begin and masked_Ha_reg[a] < finish:
        xp_triplet.append(masked_Ha_reg[a])
        yp_triplet.append(norm_spectra[a])
xp_triplet = np.array(xp_triplet)
yp_triplet = np.array(yp_triplet)

#%%
""" Ha Data; Lorentzian fit _3Lorentzian; Run after clipping data xp_triplet yp_triplet """

# DEFINE ANOTHER FUNCTION HERE 
# Zeeman splitting - returns delta_lambda 
def delta_lambda(B):
    A = np.square(6562.8/6564)
    y = 20.2*A*B
    return y 

def delta_lambda2(B):
    return (4.67*10**-7)*(np.square(6562.8))*B*(1*10**6)

def _3Lorentzian(x, lambda0, B, amp1, wid1, amp2, wid2, Q):
#    lambda0 = 6562.8
   
    A = np.square(lambda0/6564)
    C = np.square(lambda0/4101)
    delt_lam = 20.2*A*B
    delt_quad = (-1)*Q*(np.square(B))
    lambda_minus = lambda0 - delt_lam
    lambda_plus = lambda0 + delt_lam   
    lambda_pi = delt_quad
    lambda_sigma = 2*delt_quad
    
    #without quad blueshift
#    return -((amp1*wid1**2/((x-lambda_minus)**2+wid1**2)) +\
#            (amp2*wid2**2/((x-lambda0)**2+wid2**2)) +\
#                (amp1*wid1**2/((x-lambda_plus)**2+wid1**2)))+1
    
    #with quad blueshift
#    return -(((amp1*wid1**2/((x-lambda_minus-lambda_sigma)**2+wid1**2))) +\
#            ((amp2*wid2**2/((x-lambda0-lambda_pi)**2+wid2**2))) +\
#            ((amp1*wid1**2/((x-lambda_plus-lambda_sigma)**2+wid1**2))))+1
    
    #with quad blueshift just for pi
    return -(((amp1*wid1**2/((x-lambda_minus)**2+wid1**2))) +\
            ((amp2*wid2**2/((x-lambda0-lambda_pi)**2+wid2**2))) +\
            ((amp1*wid1**2/((x-lambda_plus)**2+wid1**2))))+1

#def _3Lorentzian(x, amp1, cen1, wid1, amp2,cen2,wid2, amp3,cen3,wid3):
#    return -((amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
#            (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
#                (amp3*wid3**2/((x-cen3)**2+wid3**2)))+1

popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, \
                                            p0=[6515, 30, 3, 5, 5, 5, 0.1])
#popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, \
#                                            p0=[0.8, 6525, 10, 0.8, 6565, 10, 0.8, 6600, 10]) 
#parameters from old Lorentzian
# Here p0 comes from varying the parameters from a overlaid Voigt function (not a fitted function)
# p0 are guesstimate parameters from a previous overlaid function 

for c in zip(popt_3lorentz, np.sqrt(np.diag(cov_3lorentz))):
    print("%.8f pm %.3g" % (c[0], c[1]))
#%%
plt.figure(f"Lorentzian fit {filename}")
plt.plot(xp_triplet,yp_triplet,'x', label = "WD Ha data")
plt.plot(xp_triplet, _3Lorentzian(xp_triplet, *popt_3lorentz), label = "Lorentzian c_fit")
plt.xlabel("Wavelength, $[\AA]$" , size = "15")
plt.ylabel("Normalised flux", size = "15")
plt.grid()
plt.legend()
#plt.show()
#plt.savefig('Linear Halpha Lorentzian B Fit Plots/'+imgfile)
plt.show()
##%%
""" Lorentzian Residuals Plot """
Residuals= _3Lorentzian(xp_triplet, *popt_3lorentz)-yp_triplet
print("Residual sum of squares = ", sum(np.square(Residuals)))

plt.figure(f"Residuals Lorentz Fit {filename}")
plt.plot(xp_triplet, Residuals, linewidth=2, label = "Lorentzian fit Residuals")
plt.xlabel("Wavelength, $[\AA]$", size = "15")
plt.ylabel("Residuals", size = "15")
plt.legend()
plt.grid()
#plt.savefig("Gaussian Residuals", dpi = 1000)
#plt.savefig('Linear Halpha Lorentzian B Fit Residuals/'+resfile)
plt.show()

