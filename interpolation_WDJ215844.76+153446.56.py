# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 22:46:33 2022

@author: rr3
"""

import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as spi
import lmfit as lm
import pylab
import glob
import os
import pandas as pd 
import re


#%%
""" Part 1: Load data

Notes: data - MWD spectrum
       wavelength - x-values 
       flux - y-values
"""
#load data and sort into appropriate variables
filename = "DESI_WDJ215844.76+153446.56_bin0p2.dat"
data = np.genfromtxt(f'{filename}', delimiter=' ')

wavelength = data[:,0]
flux = data[:,1]
error = data[:,2]

#wavelength = wavelength[0:28251]
#flux = flux[0:28251]
#error = error[0:28251]

plt.figure(f'{filename}')
plt.errorbar(wavelength,flux, yerr = error ,label = f"{filename}", fmt ='')
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Flux", size = "15")
#plt.xlim(5000, 9000)
#plt.ylim(0,20)
plt.grid()
plt.legend()
plt.show()
#%%
""" Part 2: Performs cuts on the data to isolate the H-alpha region

Notes: start/start_Ha_Gauss - beginning of cut
       last/last_Ha - end of cut 
       masked_Gauss_flux - y-values included within the cut
       masked_Gauss_reg - x-values included within the cut

begin/ finish define the whole region including the triplet feature 
startx/endx define the specific region to be cut out (the absorption feature) """

begin = 5000
finish = 8000
start_Gauss_1 = 5109
end_Gauss_1 = 5142
start_Gauss_2 = 6158
end_Gauss_2 = 6548
start_Gauss_3 = 7062
end_Gauss_3 = 7106
			
start_Ha_Gauss = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-begin)))[0])
end_Ha_Gauss = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-finish)))[0])
masked_Gauss_flux = flux[start_Ha_Gauss:end_Ha_Gauss]
masked_Gauss_reg = wavelength[start_Ha_Gauss:end_Ha_Gauss]
masked_Gauss_err = error[start_Ha_Gauss:end_Ha_Gauss]

dip_start1 = int(np.where(masked_Gauss_reg == min(masked_Gauss_reg, key=lambda x:abs(x-start_Gauss_1)))[0])
dip_end1 = int(np.where(masked_Gauss_reg == min(masked_Gauss_reg, key=lambda x:abs(x-end_Gauss_1)))[0])
dip_start2 = int(np.where(masked_Gauss_reg == min(masked_Gauss_reg, key=lambda x:abs(x-start_Gauss_2)))[0])
dip_end2 = int(np.where(masked_Gauss_reg == min(masked_Gauss_reg, key=lambda x:abs(x-end_Gauss_2)))[0])
dip_start3 = int(np.where(masked_Gauss_reg == min(masked_Gauss_reg, key=lambda x:abs(x-start_Gauss_3)))[0])
dip_end3 = int(np.where(masked_Gauss_reg == min(masked_Gauss_reg, key=lambda x:abs(x-end_Gauss_3)))[0])

masked_flux = list(masked_Gauss_flux)
masked_wavelength = list(masked_Gauss_reg)
masked_err = list(masked_Gauss_err)

##%%
flxframe = pd.DataFrame(masked_flux)
flxframe.loc[dip_start1:dip_end1] = np.nan
flxframe.loc[dip_start2:dip_end2] = np.nan
flxframe.loc[dip_start3:dip_end3] = np.nan
flxframe = (flxframe.dropna()).to_numpy()

wavframe = pd.DataFrame(masked_wavelength)
wavframe.loc[dip_start1:dip_end1] = np.nan
wavframe.loc[dip_start2:dip_end2] = np.nan
wavframe.loc[dip_start3:dip_end3] = np.nan
wavframe = (wavframe.dropna()).to_numpy()

flxframe = np.array(flxframe)
wavframe = np.array(wavframe)
Gauss_flux_list = flxframe[:,0]
Gauss_wav_list = wavframe[:,0]

# =============================================================================
# plt.figure()
# plt.plot(Gauss_wav_list,Gauss_flux_list,'x')
# plt.show()
# =============================================================================

""" Part 4: Fits a polynomial to the continuum and normalises the spectrum
Notes: Poly_3o is a third-order polynomial function which fits the continuum using a least-squares method
      
"""

###%%
#plt.figure()
#plt.plot(Gauss_wav_list,Gauss_flux_list,'x')
#plt.show()
##%%
def Poly_3o(x, a, b, c, d):
    y = a*x**3 + b*x**2 + c*x + d
    return y

p0 = np.array([1, 1, 1, 1]) #fit parameters
p, cov = opt.curve_fit(Poly_3o, Gauss_wav_list,Gauss_flux_list, p0) # do not change
Continuum = Poly_3o(masked_Gauss_reg, p[0], p[1], p[2], p[3]) # use masked_Gauss_reg since more data points in array

#for c in zip(p, np.sqrt(np.diag(cov))):# zips root of diag of cov matrix with related value in curve fit
#    print('%.15f pm %.4g' % (c[0], c[1]))# prints value and uncertainty, f is decimal places and G is sig figs

norm_Gauss_spectra = masked_Gauss_flux/Continuum

#plt.figure()
#plt.grid()
##plt.yticks([0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30])
##plt.plot(masked_Gauss_reg, Poly_3o(masked_Gauss_reg, p[0], p[1], p[2], p[3])/Continuum, \
##         zorder=4,color = 'red', label = "Poly")
#plt.plot(Gauss_wav_list,Gauss_flux_list,'x')
#plt.plot(masked_Gauss_reg, Continuum)
#plt.show()
##%%
plt.figure()
plt.plot(masked_Gauss_reg, norm_Gauss_spectra, label = f"{filename}") #plot the normalised spectrum
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Normalised Flux", size = "15")
plt.ylim(0.5,2)
plt.legend()
plt.show()
#%%
x1_feature = []
y1_feature = []
err1_feature = []

x2_feature = []
y2_feature = []
err2_feature = []

x3_feature = []
y3_feature = []
err3_feature = []

for a in range(len(masked_Gauss_reg)):
    if masked_Gauss_reg[a] > (start_Gauss_1) and masked_Gauss_reg[a] < (end_Gauss_1):
        x1_feature.append(masked_Gauss_reg[a])
        y1_feature.append(norm_Gauss_spectra[a])
        err1_feature.append(masked_err[a])
    if masked_Gauss_reg[a] > (start_Gauss_2) and masked_Gauss_reg[a] < (end_Gauss_2):
        x2_feature.append(masked_Gauss_reg[a])
        y2_feature.append(norm_Gauss_spectra[a])
        err2_feature.append(masked_err[a])
    if masked_Gauss_reg[a] > (start_Gauss_3) and masked_Gauss_reg[a] < (end_Gauss_3):
        x3_feature.append(masked_Gauss_reg[a])
        y3_feature.append(norm_Gauss_spectra[a])
        err3_feature.append(masked_err[a])
        
x1_feature = np.array(x1_feature)
y1_feature = np.array(y1_feature)
err1_feature = np.array(err1_feature)

x2_feature = np.array(x2_feature)
y2_feature = np.array(y2_feature)
err2_feature = np.array(err2_feature)

x3_feature = np.array(x3_feature)
y3_feature = np.array(y3_feature)
err3_feature = np.array(err3_feature)

def Gaussian(x,mu,sig,A):
    gaus = (-(A)*((1/np.sqrt((2*sp.pi)*sig))*(sp.exp(-(x-mu)**2/(2*sig**2)))))+1
    return gaus

p1 = [5122,50,3]
p2 = [6315,100,3]
p3 = [7090,50,3]


popt_Gauss_1, cov_Gauss_1 = opt.curve_fit(Gaussian, x1_feature, y1_feature, p1, sigma = err1_feature)
popt_Gauss_2, cov_Gauss_2 = opt.curve_fit(Gaussian, x2_feature, y2_feature, p2, sigma = err2_feature)
popt_Gauss_3, cov_Gauss_3 = opt.curve_fit(Gaussian, x3_feature, y3_feature, p3, sigma = err3_feature)

for c in zip(popt_Gauss_1, np.sqrt(np.diag(cov_Gauss_1))):
    print("%.8f pm %.3g" % (c[0], c[1]))
    
for c in zip(popt_Gauss_2, np.sqrt(np.diag(cov_Gauss_2))):
    print("%.8f pm %.3g" % (c[0], c[1]))
    
for c in zip(popt_Gauss_3, np.sqrt(np.diag(cov_Gauss_3))):
    print("%.8f pm %.3g" % (c[0], c[1]))   
    
plt.figure("First Gaussian")
plt.grid()
plt.plot(x1_feature,y1_feature)
plt.plot(x1_feature,Gaussian(x1_feature,*popt_Gauss_1))
plt.savefig('Intermediate Field Gaussian Fit Plots/'+f'{filename}'+'_sigma_minus.png')

plt.figure("Second Gaussian")
plt.grid()
plt.plot(x2_feature,y2_feature)
plt.plot(x2_feature,Gaussian(x2_feature,*popt_Gauss_2))
plt.savefig('Intermediate Field Gaussian Fit Plots/'+f'{filename}'+'_pi.png')

plt.figure("Third Gaussian")
plt.grid()
plt.plot(x3_feature,y3_feature)
plt.plot(x3_feature,Gaussian(x3_feature,*popt_Gauss_3))
plt.savefig('Intermediate Field Gaussian Fit Plots/'+f'{filename}'+'_sigma_plus.png')

#%%
start_B = 20E6 #starting B-field for spline 
end_B = 80E6 #end B-field for spline 

for i in range(len(b_alpha_list)):
    start_list.append(int(np.where(b_alpha_list[i] == min(b_alpha_list[i], key=lambda x:abs(x-start_B)))[0]))
    end_list.append(int(np.where(b_alpha_list[i] == min(b_alpha_list[i], key=lambda x:abs(x-end_B)))[0]))
    
cs_list = [] #list of splines

for i in range(len(b_alpha_list)):
    cs_list.append( spi.CubicSpline( b_alpha_list[i][start_list[i]-3:end_list[i]+3],wavelength_alpha_list[i][start_list[i]-3:end_list[i]+3] ) )
    #perform spline with B-field as x-axis and Wavelength as y-axis (inverted like in high-field fitting)

B_range = np.arange(start_B,end_B,1E4)
wavelength_output=[]

for i in range(len(b_alpha_list)):
    wavelength_output.append(cs_list[i](B_range))
    
pi_list = []

for i in range(6,9):
    pi_list.append(wavelength_output[i])

pi_means = np.array([np.mean(k) for k in zip(*pi_list)])

target_lambda_1 = popt_Gauss_1[0]
target_lambda_2 = popt_Gauss_2[0]
target_lambda_3 = popt_Gauss_3[0]

"""
CS for the 11th Halpha transition (lambda=7090)
"""
start_7090_11 = int(np.where(b_alpha_list[10] == min(b_alpha_list[10], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_7090_11 = int(np.where(b_alpha_list[10] == min(b_alpha_list[10], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_7090_11, end_7090_11)
cs5 = spi.CubicSpline(b_alpha_list[10][start_7090_11-1:end_7090_11+1], wavelength_alpha_list[10][start_7090_11-1:end_7090_11+1])
B_output_7090_11 = cs5(B_range)
B_7090_11_sq = np.square((B_output_7090_11-target_lambda_2))
B_7090_11_chi = (1/target_lambda_2)*(B_7090_11_sq)

"""
CS for the 8th Halpha transition (lambda=6322)
"""
start_6322_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_6322_8 = int(np.where(b_alpha_list[7] == min(b_alpha_list[7], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_6322_8, end_6322_8)
cs2 = spi.CubicSpline(b_alpha_list[7][start_6322_8-1:end_6322_8+1], wavelength_alpha_list[7][start_6322_8-1:end_6322_8+1])
B_output_6322_8 = cs2(B_range)
B_6322_8_sq = np.square((B_output_6322_8-target_lambda_2))
B_6322_8_chi = (1/target_lambda_2)*(B_6322_8_sq)

"""
CS for the 16th Hbeta transition (lambda=5122)
"""
start_5122_16 = int(np.where(b_beta_list[15] == min(b_beta_list[15], key=lambda x:abs(x-start_B)))[0]) #change depending on target_lambda
end_5122_16 = int(np.where(b_beta_list[15] == min(b_beta_list[15], key=lambda x:abs(x-end_B)))[0]) #change depending on target_lambda
print(start_5122_16, end_5122_16)
cs3 = spi.CubicSpline(b_beta_list[15][start_5122_16-1:end_5122_16+1], wavelength_beta_list[15][start_5122_16-1:end_5122_16+1])
B_output_5122_16 = cs3(B_range)
B_5122_16_sq = np.square((B_output_5122_16-target_lambda_3))
B_5122_16_chi = (1/target_lambda_3)*(B_5122_16_sq)

"""
CS for the mean of the 7th, 8th and 9th Halpha transitions (lambda=6322)
"""
pi_sq = np.square(pi_means-target_lambda_2)
pi_chi = (1/(target_lambda_2))*pi_sq

"""
Chi-squared fitting
"""

chi_tot = B_7090_11_chi + B_6322_8_chi + B_5122_16_chi
chi_tot_2 = B_7090_11_chi + pi_chi + B_5122_16_chi

chi_min = chi_tot.min()
chi_min_idx = int((np.where(chi_tot == chi_tot.min()))[0])

chi_min_2 = chi_tot_2.min()
chi_min_idx_2 = int((np.where(chi_tot_2 == chi_tot_2.min()))[0])

chi_B_guess = B_range[chi_min_idx]
chi_B_guess_2 = B_range[chi_min_idx_2]

print(chi_B_guess/1E6)
print(chi_B_guess_2/1E6)
print(chi_min)
print(chi_min_2)
#%%
"""
Chi-squared confidence intervals

Minimum chi-squared + 2.30 gives 68.3% of data
Minimum chi-squared + 6.17 gives 95.45% of data
Minimum chi-squared + 11.8 gives 99.73% of data
"""
confidence_range_idx = np.where((chi_tot <= (chi_min + 2.30)))[0]
B_chi_range = np.array([B_range[i] for i in confidence_range_idx])
B_chi_range = (1E-6)*B_chi_range
B_chi_lower = min(B_chi_range)
B_chi_upper = max(B_chi_range)
print(B_chi_lower, B_chi_upper)
#%%
"""
Chi-squared plotting
"""
x_arr = np.arange(0, len(chi_tot))
chi_con = [chi_tot_2[i] for i in confidence_range_idx]
plt.figure()
plt.grid()
plt.plot(x_arr,chi_tot)
plt.plot(confidence_range_idx, chi_con, linewidth=5)
plt.plot(chi_min_idx, chi_min, 'x', markersize=10, linewidth=20)
plt.show()

