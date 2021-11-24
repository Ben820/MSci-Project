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
import glob
#%%
""" Part 1: Load data

Notes: data - MWD spectrum
       wavelength - x-values 
       flux - y-values
"""
#datasets = np.loadtxt('Prelim set of Linear WDs.csv',skiprows = 1, delimiter = ',', unpack = True)
#files = [file for file in glob.glob('C:\\Users\44743\Documents\Imperial Year 4\MSci Project\DESI_DxH\DESI_DxH')]
#for file_name in files:
#    with io.open(file_name, 'rb') as image_file:
#        content = image_file.read()

#load data and sort into appropriate variables
#Ideal system (first to fit)
#DESI_WDJ110344.93+510048.61

#Unfriendly systems 
#DESI_WDJ073615.92+403335.18
#DESI_WDJ082302.41+334534.27
#DESI_WDJ112328.50+095619.35
#DESI_WDJ112926.23+493931.89
#DESI_WDJ165200.65+411029.96 # Quadratic regime 
filename = "DESI_WDJ153349.03+005916.21_bin0p2.dat"
data = np.genfromtxt(f'{filename}', delimiter=' ')
#DESI_WDJ110344.93+510048.61_bin0p2
#DESI_WDJ112328.50+095619.35_bin0p2
wavelength = data[:,0]
flux = data[:,1]
error = data[:,2]

plt.figure("Whole spectrum")
plt.errorbar(wavelength,flux, yerr = error ,label = f"{filename}", fmt ='')
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Flux", size = "15")
plt.grid()
plt.legend()
plt.show()
##%% # THIS IS THE PART SEPARATING THE TWO CELLS ------------------------------
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
masked_err = error[start_Ha:end_Ha]

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
masked_error = list(masked_err)
#try and remove unncessary features from the spectrum (specific to spectrum,
#COMMENT OUT THIS SECTION IF NOT NEEDED)
del masked_flux[dip_start:dip_end]
del masked_wavelength[dip_start:dip_end]
del masked_error[dip_start:dip_end]
del masked_flux[505:517]
del masked_wavelength[505:517]
del masked_error[505:517]

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

#for c in zip(p, np.sqrt(np.diag(cov))):# zips root of diag of cov matrix with related value in curve fit
#    print('%.15f pm %.4g' % (c[0], c[1]))# prints value and uncertainty, f is decimal places and G is sig figs

norm_spectra = masked_Ha_flux/Continuum

plt.figure("Normalised absorption feature")
plt.grid()
#plt.yticks([0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30])
#plt.plot(masked_Ha_reg, Poly_3o(masked_Ha_reg, p[0], p[1], p[2], p[3])/Continuum, \
#         zorder=4,color = 'red', label = "Poly")
plt.plot(masked_Ha_reg, norm_spectra, label = f"{filename}") #plot the normalised spectrum
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Normalised Flux", size = "15")
plt.legend()
plt.show()

#%%
""" Triplet Lorentzian Profile 

Notes: 
"""
xp_triplet = []
yp_triplet = []
err_triplet = []
for a in range(len(masked_Ha_reg)):
    if masked_Ha_reg[a] > 6400 and masked_Ha_reg[a] < 6750:
        xp_triplet.append(masked_Ha_reg[a])
        yp_triplet.append(norm_spectra[a])
        err_triplet.append(masked_err[a])
xp_triplet = np.array(xp_triplet)
yp_triplet = np.array(yp_triplet)
err_triplet = np.array(err_triplet)

""" Ha Data; Lorentzian fit _3Lorentzian; Run after clipping data xp_triplet yp_triplet """

# DEFINE ANOTHER FUNCTION HERE 
# Zeeman splitting - returns delta_lambda 
def delta_lambda(B):
    A = np.square(6562.8/6564)
    y = 20.2*A*B
    return y 

def delta_lambda2(B):
    return (4.67*10**-7)*(np.square(6562.8))*B*(1*10**6)

def _3Lorentzian(x, lambda0, B, amp1, wid1, amp2, wid2):
#    lambda0 = 6562.8
    A = np.square(lambda0/6564)
    delt_lam = 20.2*A*B
    lambda_minus = lambda0 - delt_lam
    lambda_plus = lambda0 + delt_lam     
    return -((amp1*wid1**2/((x-lambda_minus)**2+wid1**2)) +\
            (amp2*wid2**2/((x-lambda0)**2+wid2**2)) +\
                (amp1*wid1**2/((x-lambda_plus)**2+wid1**2)))+1

Data = {}
Res_list = []
B_list = []
lambda0_list = []

rangeBval = np.arange(0,50.5,0.5)
for i in rangeBval:
    p0 = np.array([6562.8, i, 0.8, 10, 0.8, 10])
    
    try:
        popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, p0, sigma = err_triplet)
    except:
        pass
    
#    popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, p0)
#    for c in zip(popt_VoigtNew, np.sqrt(np.diag(cov_VoigtNew))):
#        print("%.8f pm %.3g" % (c[0], c[1]))
    Residuals= _3Lorentzian(xp_triplet, *popt_3lorentz)-yp_triplet
    Res_list.append(sum(np.square(Residuals)))
    B_list.append(popt_3lorentz[1])
    lambda0_list.append(popt_3lorentz[0])
    
Data[1] = [Res_list, B_list, lambda0_list]

index = np.where(Res_list == np.amin(Res_list))[0][0]

# Method 
# The determined B value is BACK INTO curvefit to re-estimate the final parameters of the fit 
popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, \
                                            p0=[6562.8, B_list[index], 0.8, 10, 0.8, 10], \
                                            sigma = err_triplet)#, \
                                            #bounds = ((-np.inf, 0, -np.inf, -np.inf, -np.inf, -np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))

## ALTERNATIVE METHOD: Just take first curvefit estimate - prevents over fitting and finding a 
# wrong local minima in least squares
#p0 = np.array([6562.8, rangeBval[2], 0.8, 10, 0.8, 10])
#popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, p0, sigma = err_triplet)

## ALTERNATIVE METHOD 2: (Supercedes Method):
# Puts the list of returned B values from first curvefit into another curvefit, scans through and 
# looks to see which has the lowest residuals (on the resulting curvefit); 
# selects the corresponding new index for that fit, and then uses that index to identify the correct
# B value from B_list (the first! returned set of B values) (NOTE THIS IS THE SECOND TIME IT IS CURVEFITTED)
Res_list2 = []
B_list2 = []
lambda0_list2 = []

for k in range(len(B_list)):
    p0 = np.array([6562.8, B_list[k], 0.8, 10, 0.8, 10])
    
    try:
        popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, p0, sigma = err_triplet)
    except:
        pass
    Residuals2= _3Lorentzian(xp_triplet, *popt_3lorentz)-yp_triplet
    Res_list2.append(sum(np.square(Residuals2)))
    B_list2.append(popt_3lorentz[1])
    lambda0_list2.append(popt_3lorentz[0])

index2 = np.where(Res_list2 == np.amin(Res_list2))[0][0]

# Take the minimum residual of both residual lists combined - thus doesnt double iterate 
# if not requied/ beneficial 
A = Res_list+Res_list2
index3 = np.where(A == np.amin(A))[0][0]

# The determined B value is BACK INTO curvefit to re-estimate the final parameters of the fit 
# NOTE THIS IS STILL ONLY THE SECOND CURVEFIT - NOTE USE OF B_list NOT! B_list2
popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, \
                                            p0=[6562.8, B_list[index3], 0.8, 10, 0.8, 10], \
                                            sigma = err_triplet)#, \
                                            #bounds = ((-np.inf, 0, -np.inf, -np.inf, -np.inf, -np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))

for c in zip(popt_3lorentz, np.sqrt(np.diag(cov_3lorentz))):
    print("%.8f pm %.3g" % (c[0], c[1]))

Residuals= _3Lorentzian(xp_triplet, *popt_3lorentz)-yp_triplet
print("Lorentzian Residual sum of squares = ", sum(np.square(Residuals)))

##%% PLotting Cell
fig, axs = plt.subplots(2, gridspec_kw={'height_ratios': [2.5, 1]})
fig.suptitle(f"Lorentzian fit {filename} \n B = {popt_3lorentz[1]} +/- {np.sqrt(np.diag(cov_3lorentz))[1]}", \
             fontsize = "13")
axs[0].plot(xp_triplet, yp_triplet,'x', label = "WD Ha data")
axs[0].plot(xp_triplet, _3Lorentzian(xp_triplet, *popt_3lorentz), linewidth=2, color = "orange", \
                              label = "Lorentzian c_fit")
#axs[0].text(r"B = {popt_3lorentz}")
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











































