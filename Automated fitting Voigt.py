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
import lmfit as lm
import pylab
##%%
""" Part 1: Load data

Notes: data - MWD spectrum
       wavelength - x-values 
       flux - y-values
"""
#load data and sort into appropriate variables
filename = "DESI_WDJ110344.93+510048.61_bin0p2.dat"
data = np.genfromtxt(f'{filename}', delimiter=' ')
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
abs_start = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-6400)))[0])
abs_end = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-6800)))[0])

abs_flux = flux[abs_start:abs_end]
abs_reg = wavelength[abs_start:abs_end]

plt.figure("Absoprtion feature region +")
plt.grid()
plt.plot(abs_reg,abs_flux)
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Flux", size = "15")
plt.show()

#find the wavelength values where the feature should be cut (manual values are 
#inputted into the x:abs(x-XXXXXXXX) section of the np.where command)
dip_start = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-6500)))[0])
dip_end = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-6650)))[0])

masked_flux = list(masked_Ha_flux)
masked_wavelength = list(masked_Ha_reg)
#try and remove unncessary features from the spectrum (specific to spectrum,
#COMMENT OUT THIS SECTION IF NOT NEEDED)
del masked_flux[dip_start:dip_end]
del masked_wavelength[dip_start:dip_end]
del masked_flux[505:517]
del masked_wavelength[505:517]

plt.figure("Spectra with absorption feature removed")
plt.grid()
plt.plot(masked_wavelength,masked_flux,'x')
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Flux", size = "15")
plt.show()

#avg_start_array = masked_flux[0:10]
#avg_end_array = masked_flux[-11:-1]
#avg_start_val = np.mean(masked_flux)
#avg_end_val = np.mean(masked_flux)
##%%
""" Part 4: Fits a polynomial to the continuum and normalises the spectrum

Notes: Poly_3o is a third-order polynomial function which fits the continuum using a least-squares method
      
"""
def Poly_3o(x, a, b, c, d):
    y = a*x**3 + b*x**2 + c*x + d
    return y

p0 = np.array([1, 1, 1, 1]) #fit parameters
p, cov = opt.curve_fit(Poly_3o, masked_wavelength, masked_flux, p0) # do not change
Continuum = Poly_3o(masked_Ha_reg, p[0], p[1], p[2], p[3]) # use masked_Ha_reg since more data points in array

for c in zip(p, np.sqrt(np.diag(cov))):# zips root of diag of cov matrix with related value in curve fit
    print('%.15f pm %.4g' % (c[0], c[1]))# prints value and uncertainty, f is decimal places and G is sig figs

norm_spectra = masked_Ha_flux/Continuum

plt.figure("Normalised absorption feature")
plt.grid()
#plt.yticks([0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30])
plt.plot(masked_Ha_reg, Poly_3o(masked_Ha_reg, p[0], p[1], p[2], p[3])/Continuum, \
         zorder=4,color = 'red', label = "Poly")
plt.plot(masked_Ha_reg, norm_spectra) #plot the normalised spectrum
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Normalised Flux", size = "15")
plt.show()

#%%
""" Triplet Voigt Profile 

Notes: 
"""
xp_triplet = []
yp_triplet =[]
for a in range(len(masked_Ha_reg)):
    if masked_Ha_reg[a] > 6400 and masked_Ha_reg[a] < 6750:
        xp_triplet.append(masked_Ha_reg[a])
        yp_triplet.append(norm_spectra[a])
xp_triplet = np.array(xp_triplet)
yp_triplet = np.array(yp_triplet)

##%% VOIGT TRIPLET
""" Initial Voigt Triplet Fitting attempt; VoigtNew """

def VoigtNew(x, lambda0, B, A1, sigma1, gamma1, A2, sigma2, gamma2):
    #    lambda0 = 6562.8
    A = np.square(lambda0/6564)
    delt_lam = 20.2*A*B
    lambda_minus = lambda0 - delt_lam
    lambda_plus = lambda0 + delt_lam     

    return -((lm.models.voigt(x, A1, lambda_minus, sigma1, gamma1)) \
             + (lm.models.voigt(x, A2, lambda0, sigma2, gamma2)) \
             + (lm.models.voigt(x, A1 , lambda_plus, sigma1, gamma1)))+1
    
# parameters for new Voigt (8 params)

Data = {}
Res_list = []
B_list = []
lambda0_list = []

for i in range(1,50):
    p0 = np.array([6562.8, i, 5, 10, 1.85032, .5, 10, 1.85032])
    
    try:
        popt_VoigtNew, cov_VoigtNew = opt.curve_fit(VoigtNew, xp_triplet, yp_triplet, p0)
    except:
        pass
    
#    for c in zip(popt_VoigtNew, np.sqrt(np.diag(cov_VoigtNew))):
#        print("%.8f pm %.3g" % (c[0], c[1]))
    Residuals= VoigtNew(xp_triplet, *popt_VoigtNew)-yp_triplet
    Res_list.append(sum(np.square(Residuals)))
    B_list.append(popt_VoigtNew[1])
    lambda0_list.append(popt_VoigtNew[0])
    
Data[1] = [Res_list, B_list, lambda0_list]

index = np.where(Res_list == np.amin(Res_list))[0][0]

# The determined B value is BACK INTO curvefit to re-estimate the final parameters of the fit 
popt_VoigtNew, cov_VoigtNew = opt.curve_fit(VoigtNew, xp_triplet, yp_triplet, \
                                            p0=[6562.8, B_list[index], 5, 10, 1.85032, .5, 10, 1.85032])

for c in zip(popt_VoigtNew, np.sqrt(np.diag(cov_VoigtNew))):
    print("%.8f pm %.3g" % (c[0], c[1]))

Residuals= VoigtNew(xp_triplet, *popt_VoigtNew)-yp_triplet
print("Voigt Residual sum of squares = ", sum(np.square(Residuals)))

##%% PLotting Cell
fig, axs = plt.subplots(2, gridspec_kw={'height_ratios': [2.5, 1]})
fig.suptitle(f"Voigt fit {filename}", fontsize = "13")
axs[0].plot(xp_triplet, yp_triplet,'x', label = "WD Ha data")
axs[0].plot(xp_triplet, VoigtNew(xp_triplet, *popt_VoigtNew), linewidth=2, color = "red", \
                              label = "Voigt c_fit")
for ax in axs.flat:
    ax.set_xlabel('Wavelength $[\AA]$', fontsize = "13")
    ax.set_ylabel('Normalised flux', fontsize = "13")
for ax in axs.flat:
    ax.label_outer()
axs[1].set_ylabel('Flux residuals')
axs[0].grid()
#""" Voigt Residuals Plot """
axs[1].plot(xp_triplet, Residuals, linewidth=2)#, label = "Lorentzian fit Residuals")
axs[0].legend()
#axs[1].legend()
axs[1].grid()
#plt.savefig("Gaussian Residuals", dpi = 1000)
plt.show()













































