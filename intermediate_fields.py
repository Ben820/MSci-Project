# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 15:00:13 2022

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
""" NEED TO RUN hydrogen_transitions_v2.py BEFORE RUNNING ANY CODE IN THIS FILE"""

#%%
""" Part 1: Load data

Notes: data - MWD spectrum
       wavelength - x-values 
       flux - y-values
"""
#load data and sort into appropriate variables
filename = "DESI_WDJ180813.60+652321.24_bin0p2.dat"
data = np.genfromtxt(f'{filename}', delimiter=' ')

wavelength = data[:,0]
flux = data[:,1]
error = data[:,2]

#wavelength = wavelength[0:28251]
#flux = flux[0:28251]
#error = error[0:28251]

plt.figure("Whole spectrum")
plt.errorbar(wavelength,flux, yerr = error ,label = f"{filename}", fmt ='')
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Flux", size = "15")
plt.xlim(5000, 9000)
plt.ylim(4,20)
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
start_Gauss_1 = 5810
end_Gauss_1 = 6250
start_Gauss_2 = 6465
end_Gauss_2 = 6555
start_Gauss_3 = 6835
end_Gauss_3 = 7150
			

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

p1 = [6040,150,3]
p2 = [6515,50,3]
p3 = [6995,100,3]


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
start_B = 14E6 #starting B-field for spline 
end_B = 60E6 #end B-field for spline 

start_list = [] #list of indices for the start of the spline
end_list = [] #list of indices for the end of the spline

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
#%%
# =============================================================================
# plt.figure()
# plt.grid()
# plt.xlabel('B [G]')
# plt.ylabel('Wavelength $[\AA]$')
# for i in range(len(alpha_array)):
#     plt.plot(b_alpha_list[i],wavelength_alpha_list[i],'x',markersize=15)
# for i in range(len(b_alpha_list)):
#     plt.scatter(B_range,wavelength_output[i],marker='+')
# 
# # =============================================================================
# # plt.scatter(B_range,wavelength_output[0],marker='+')
# # plt.scatter(B_range,wavelength_output[4],marker='+')
# # =============================================================================
# plt.gca().invert_yaxis()
# plt.show()
# =============================================================================
#%%
left_sigma_list = []  #group1 - left sigma
pi_list = [] #group2 - pi
right_sigma_list = [] #group3 - right sigma

for i in range(5):
    left_sigma_list.append(wavelength_output[i])
for i in range(6,9):
    pi_list.append(wavelength_output[i])
for i in range(9,len(b_alpha_list)):
    right_sigma_list.append(wavelength_output[i])
#%%
left_sigma_means = np.array([np.mean(k) for k in zip(*left_sigma_list)])
pi_means = np.array([np.mean(k) for k in zip(*pi_list)])
right_sigma_means = np.array([np.mean(k) for k in zip(*right_sigma_list)])

left_sigma_diff = np.subtract(left_sigma_list[4],left_sigma_list[0])
pi_diff = np.subtract(pi_list[2],pi_list[0])
right_sigma_diff = np.subtract(right_sigma_list[4],right_sigma_list[0])

#%%
left_sigma_peak = int(np.where(left_sigma_means == min(left_sigma_means, key=lambda x:abs(x-popt_Gauss_1[0])))[0])
left_sigma_spread = int(np.where(left_sigma_diff == min(left_sigma_diff, key=lambda x:abs(x-popt_Gauss_1[1])))[0])

pi_peak = int(np.where(pi_means == min(pi_means, key=lambda x:abs(x-popt_Gauss_2[0])))[0])
pi_spread = int(np.where(pi_diff == min(pi_diff, key=lambda x:abs(x-popt_Gauss_2[1])))[0])

right_sigma_peak = int(np.where(right_sigma_means == min(right_sigma_means, key=lambda x:abs(x-popt_Gauss_3[0])))[0])
right_sigma_spread = int(np.where(right_sigma_diff == min(right_sigma_diff, key=lambda x:abs(x-popt_Gauss_3[1])))[0])
#%%
left_B_guess = B_range[left_sigma_peak]
pi_B_guess = B_range[pi_peak]
right_B_guess = B_range[right_sigma_peak]

print(left_B_guess/1E6)
print(pi_B_guess/1E6) 
print(right_B_guess/1E6)

#%%
# =============================================================================
# """
# Chi-squared fitting
# """
# 
# left_sigma_sq = np.square(left_sigma_means-popt_Gauss_1[0])
# pi_sq = np.square(pi_means-popt_Gauss_2[0])
# right_sigma_sq = np.square(right_sigma_means-popt_Gauss_3[0])
# 
# left_sigma_chi = (1/(popt_Gauss_1[1]))*left_sigma_sq
# pi_chi = (1/(popt_Gauss_2[1]))*pi_sq
# right_sigma_chi = (1/(popt_Gauss_3[1]))*right_sigma_sq
# 
# left_sigma_chi_2 = (1/(popt_Gauss_1[0]))*left_sigma_sq
# pi_chi_2 = (1/(popt_Gauss_2[0]))*pi_sq
# right_sigma_chi_2 = (1/(popt_Gauss_3[0]))*right_sigma_sq
# 
# chi_tot = left_sigma_chi+pi_chi+right_sigma_chi
# chi_tot_2 = left_sigma_chi_2+pi_chi_2+right_sigma_chi_2
# 
# chi_min = chi_tot.min()
# chi_min_idx = int((np.where(chi_tot == chi_tot.min()))[0])
# 
# chi_min_2 = chi_tot_2.min()
# chi_min_idx_2 = int((np.where(chi_tot_2 == chi_tot_2.min()))[0])
# 
# chi_B_guess = B_range[chi_min_idx]
# chi_B_guess_2 = B_range[chi_min_idx_2]
# 
# print(chi_B_guess/1E6)
# print(chi_B_guess_2/1E6)
# =============================================================================
#%%
#%%
# =============================================================================
# """
# RUN hydrogen_transitions_v2.py FIRST TO GET ALL THE TRANSITIONS
# """
# start_B = 14E6 #starting B-field for spline 
# end_B = 30E6 #end B-field for spline 
# 
# start_list = [] #list of indices for the start of the spline
# end_list = [] #list of indices for the end of the spline
# 
# for i in range(len(b_alpha_list)):
#     start_list.append(int(np.where(b_alpha_list[i] == min(b_alpha_list[i], key=lambda x:abs(x-start_B)))[0]))
#     end_list.append(int(np.where(b_alpha_list[i] == min(b_alpha_list[i], key=lambda x:abs(x-end_B)))[0]))
# #%%
# sorted_wavelength_list = [] #wavelengths for first 9 transitions
# sorted_B_list = [] #B-fields for first 9 transitions
# remaining_wavelength_list = [] #wavelengths for transitions that are already in a suitable format
# remaining_B_list = [] #B-fields for transitions that are already in a suitable format
# 
# for i in range(0,8):
#     sorted_wavelength_list.append(np.sort(wavelength_alpha_list[i][start_list[i]-3:end_list[i]+3]))
#     sorted_B_list.append(np.sort(b_alpha_list[i][start_list[i]-3:end_list[i]+3]))
# for i in range(9,14):
#     remaining_wavelength_list.append(wavelength_alpha_list[i][start_list[i]-3:end_list[i]+3])
#     remaining_B_list.append(b_alpha_list[i][start_list[i]-3:end_list[i]+3])
# #%%
# sorted_cs_list = [] #splines for the first 9 transition
# remaining_cs_list = [] #splines for the transitions already in a suitable format
# 
# for i in range(len(sorted_wavelength_list)):
#     sorted_cs_list.append( spi.CubicSpline( sorted_wavelength_list[i][start_list[i]-3:end_list[i]+3],sorted_B_list[i][start_list[i]-3:end_list[i]+3] ) )
# for i in range(len(remaining_wavelength_list)):
#     remaining_cs_list.append( spi.CubicSpline( remaining_wavelength_list[i][start_list[i]-3:end_list[i]+3],remaining_B_list[i][start_list[i]-3:end_list[i]+3] ) )
# 
# #%%  
# target_wav_range = np.arange(6400, 5800, -1E1) #find target_lambda from the Gaussian fit to the peak
# B_output = sorted_cs_list[0](target_wav_range )
# plt.figure()
# plt.grid()
# plt.xlabel('B [G]')
# plt.ylabel('Wavelength $[\AA]$')
# plt.plot(wavelength_alpha_list[0],b_alpha_list[0])
# #plt.plot(wavelength_alpha_list[0], sorted_cs_list[0])
# plt.scatter(target_wav_range, B_output, marker='+')
# plt.show()
# 
# #%%
# #cs_list.append(spi.CubicSpline(wavelength_alpha_list[i][start_list[i]-2:end_list[i]+2], b_alpha_list[i][start_list[i]-2:end_list[i]+2]))
# """    
# NEED TO REFLECT THE TRANSITIONS WITH DECREASING X IN THE VERTICAL AXIS THAT GOES THROUGH
# THE MINIMUM VALUE OF THE TRANSITION, DO THE INTERPOLATION AND THEN FLIP BACK (PASSING THE
# WAVELENGTH VALUES THROUGH MAY HAVE TO BE DONE EITHER BEFORE OR AFTER THE FLIP BACK --> NOT SURE)
# """
# =============================================================================
