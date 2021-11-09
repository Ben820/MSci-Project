# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 15:01:17 2021

@author: Ben Amroota
"""
"""
Reading data from DA white dwarfs
- start date: 10/10/21 
Directories: 
    Documents Imperial Year 4 MSci Project DA_models DA_models
    Documents\Imperial Year 4\MSci Project\DESI_DxH\DESI_DxH
    

Method: 
    Part 1 - Removes absorption featrure from continuum
    Part 2 - Curvefits the remaining continuum to get smooth continuum function
    Part 3 - Divides the original data spectrum by continuum function; returns normalised spectra
    Part 4 - Fit Voigt profiles to absorption features in normalised spectra 
"""
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm

##DA_models\DA_models
#data = np.genfromtxt('database1.dat', delimiter=' ')
#stdwave = np.arange(3000,11000.1, 0.1)


##DESI_DxH\DESI_DxH
data, stdwave, error = np.loadtxt('DESI_WDJ073615.92+403335.18_bin0p2.dat', unpack = True)

#%%
""" Part 1: Removes absorption featrure from continuum

Notes: reg = x region
       flux = y continuum 
       masked_Ha = region occupied by absorption feature (masked from continuum)
       data spectra == original data set 
"""
#"""
## Hardcoded wavelength region 
#stdwave_mod = np.around(stdwave,1)
#start_Ha = np.where(stdwave_mod==6200)
#end_Ha = np.where(stdwave_mod==7000)
#start_Ha = start_Ha[0]
#end_Ha = end_Ha[0]
#start_Ha = int(start_Ha)
#end_Ha = int(end_Ha)
#"""

# Selects the region occupied by absorption feature and ISOLATES from data spectra (e.g. for later removal) 
## User input for start and end regions
start = 6200#int(input("Startpoint for H-alpha cut:"))
last = 7200#int(input("Endpoint for H-alpha cut:"))
stdwave_mod = np.around(stdwave,1)
start_Ha = int(np.where(stdwave_mod==start)[0])
end_Ha = int(np.where(stdwave_mod==last)[0])
masked_Ha_flux = data[start_Ha:end_Ha]
masked_Ha_reg = stdwave_mod[start_Ha:end_Ha]

##%%
plt.figure()
plt.plot(masked_Ha_reg,masked_Ha_flux,'x')
plt.grid()
plt.show()

# Extract linear region where there is guaranteed continuum
mid_Ha = start_Ha+100
cont_flux = data[start_Ha:mid_Ha]
cont_reg = stdwave_mod[start_Ha:mid_Ha]

#plt.figure()
#plt.plot(cont_reg,cont_flux)
#plt.grid()
#plt.show()

#%%
""" Part 1b """
# Gradient - Used in defining continuum without absorption
startpoint_reg = cont_reg[1]
startpoint_flux = cont_flux[1]
endpoint_reg = cont_reg[-2]
endpoint_flux = cont_flux[-2]
gradient=(endpoint_flux-startpoint_flux)/(endpoint_reg-startpoint_reg)
print(gradient)
#%%
""" Part 1c """
# Further clips region to isolate absorption feature 
abs_start = int(np.where(masked_Ha_reg==6400)[0])
abs_end = int(np.where(masked_Ha_reg==6800)[0])
abs_flux = masked_Ha_flux[abs_start:abs_end]
abs_reg = masked_Ha_reg[abs_start:abs_end]

#plt.figure()
#plt.plot(abs_reg,abs_flux,'x')
#plt.grid()
#plt.show()
#%%
""" Part 1d """
# good_val defines the continuum without the absorption feature; masked_Ha is removed from data spectra
good_vals=[]
bad_vals=[]
for i in range(len(masked_Ha_flux)-2):
    if ((masked_Ha_flux[i+1]-masked_Ha_flux[i])/(masked_Ha_reg[i+1]-masked_Ha_reg[i]))<0 and abs((masked_Ha_flux[i+1]-masked_Ha_flux[i])/(masked_Ha_reg[i+1]-masked_Ha_reg[i])) < (-2)*gradient:
        good_vals.append(i)
    else:
        bad_vals.append(i)
x_goodvals = [masked_Ha_reg[j] for j in good_vals]
y_goodvals = [masked_Ha_flux[j] for j in good_vals]
# bad_vals are optional - not used atm 
x_badvals = [masked_Ha_reg[j] for j in bad_vals]
y_badvals = [masked_Ha_flux[j] for j in bad_vals]

#plt.figure()
#plt.plot(x_goodvals,y_goodvals,'x',color='blue')
##plt.plot(x_badvals,y_badvals,'o',color='red')
#plt.grid()
#plt.show()
#%% 
""" Part 2: Curvefits the remaining continuum to get smooth continuum function

Notes: norm_spectra = normalised spectra (with absorption features)
"""
# Curvefit Polynomial to Unnormalised spectra 
def Poly_3o(x, a, b, c, d):
    y = a*x**3 + b*x**2 + c*x + d
    return y
p0 = np.array([1, 1, 1, 1])
p, cov = opt.curve_fit(Poly_3o, x_goodvals, y_goodvals, p0) # do not change
Continuum = Poly_3o(masked_Ha_reg, p[0], p[1], p[2], p[3]) # use masked_Ha_reg since more data points in array

for c in zip(p, np.sqrt(np.diag(cov))):# zips root of diag of cov matrix with related value in curve fit
    print('%.15f pm %.4g' % (c[0], c[1]))# prints value and uncertainty, f is decimal places and G is sig figs

#plt.figure()
#plt.plot(x_goodvals, y_goodvals, 'x')
#plt.plot(masked_Ha_reg, Continuum, linewidth=2, color = "red", label = "Continuum fit")
#plt.xlabel('Wavelength [$nm$]')
#plt.ylabel("Amplitude")
#plt.legend()
#plt.grid()
#plt.show()
#%%
""" Part 3: Divides the original data spectrum by continuum function; returns normalised spectra 

Notes:
"""
norm_spectra = masked_Ha_flux/Continuum

#plt.figure()
#plt.plot(masked_Ha_reg, norm_spectra, 'x')
#plt.show()
##%%
""" Part 4: Fit Voigt profiles to absorption features in normalised spectra 

Notes:
"""
xp = []
yp =[]
for a in range(len(masked_Ha_reg)):
    if masked_Ha_reg[a] > 6000 and masked_Ha_reg[a] < 7200:
        xp.append(masked_Ha_reg[a])
        yp.append(norm_spectra[a])

xp = np.array(xp)
yp = np.array(yp)

# Curvefit Voigt to Normalised clipped spectra (only the absorption feature)
p0 = np.array([25.355, 6570, -2.6, 0.85032])
def voigt_fit(x, a, b, c, d):
    y = -lm.models.voigt(x, a, b, c, d)+1
    return y
p, cov = opt.curve_fit(voigt_fit, xp, yp, p0, maxfev = 100000)

for c in zip(p, np.sqrt(np.diag(cov))):# zips root of diag of cov matrix with related value in curve fit
    print('%.15f pm %.4g' % (c[0], c[1]))# prints value and uncertainty, f is decimal places and G is sig figs

plt.figure()
plt.plot(masked_Ha_reg, norm_spectra, 'x')
plt.plot(masked_Ha_reg, voigt_fit(masked_Ha_reg, p[0], p[1], p[2], p[3]), linewidth=2, color = "red", label = "Voigt fit")
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Frequency")
plt.legend()
plt.grid()
plt.show()
#%%
# Magnetic dipole 
# infinites in covariance matrix 
# field strength configurations 
#%%
""" Triplet Voigt Profile 

Notes: 
"""
xp_triplet = []
yp_triplet =[]
for a in range(len(masked_Ha_reg)):
    if masked_Ha_reg[a] > 6000 and masked_Ha_reg[a] < 7200:
        xp_triplet.append(masked_Ha_reg[a])
        yp_triplet.append(norm_spectra[a])
xp_triplet = np.array(xp_triplet)
yp_triplet = np.array(yp_triplet)

##%% VOIGT TRIPLET
#plt.plot(x, y, color = "red", label = "Corrected Hg \n doublet")
p0Vt1 = np.array([25.355, 6470, -2.6, 0.85032])
p0Vt2= np.array([25.355, 6570, -2.6, 0.85032])
p0Vt3= np.array([25.355, 6670, -2.6, 0.85032])
#pylab.plot(xg, G(xg, alpha), ls=':', label='Gaussian')
#pylab.plot(xg, L(xg, gamma), ls='--', label='Lorentzian')
#lm.models.voigt()
def VoigtNew(x, A1, centre1, sigma1, gamma1, A2 , centre2, sigma2, gamma2, A3 , centre3, sigma3, gamma3):
    return (lm.models.voigt(x, A1, centre1, sigma1, gamma1)) + (lm.models.voigt(x, A2, centre2, sigma2, gamma2)) + (lm.models.voigt(x, A3 , centre3, sigma3, gamma3))
pVt1, covVt1 = opt.curve_fit(lm.models.voigt, xp_triplet, yp_triplet, p0Vt1)
pVt2, covVt2 = opt.curve_fit(lm.models.voigt, xp_triplet, yp_triplet, p0Vt2)
pVt3, covVt3 = opt.curve_fit(lm.models.voigt, xp_triplet, yp_triplet, p0Vt3)
x2 = np.arange(574,583,0.000005)
#plt.plot(xp, yp, color = "blue")

plt.plot(x2, VoigtNew(x2, pVt1[0], pVt1[1], pVt1[2], pVt1[3], pVt2[0], pVt2[1], pVt2[2], pVt2[3], pVt3[0], pVt3[1], pVt3[2], pVt3[3]), linewidth=2, color = "royalblue", label = "Voigt fit")

plt.legend()
#plt.xlim(575,582)
#plt.ylim(-0.01,0.4)
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Frequency")
plt.grid()
plt.show()
#plt.savefig("Voigt fit doublet", dpi = 1000)
for c in zip(pVt1, np.sqrt(np.diag(covVt1))):
    print("%.8f pm %.3g" % (c[0], c[1]))
for c in zip(pVt2, np.sqrt(np.diag(covVt2))):
    print("%.8f pm %.3g" % (c[0], c[1]))
for c in zip(pVt3, np.sqrt(np.diag(covVt3))):
    print("%.8f pm %.3g" % (c[0], c[1]))

#plt.savefig("Combined fit doublet reshuff", dpi = 1000)



#%%
plt.figure()
x_test = np.arange(6000,7200,0.01)
y_test = VoigtNew(x_test, p0Vt1[0], p0Vt1[1], p0Vt1[2], p0Vt1[3], p0Vt2[0], p0Vt2[1], p0Vt2[2], p0Vt2[3], p0Vt3[0], p0Vt3[1], p0Vt3[2], p0Vt3[3])
plt.plot(x_test, y_test)
plt.show()




































































