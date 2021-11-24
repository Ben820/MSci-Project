# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 16:53:42 2021

@author: 44743
"""
"""
Triplet Voigt Testing Script 
"""
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm
#%%
""" Method to clean up noise 
Needs WDspectra py file to run
Problems: Removes details in peaks e.g amplitude and width differences """

import numpy as np
 
kernel_size = 50
kernel = np.ones(kernel_size) / kernel_size
data_convolved = np.convolve(norm_spectra, kernel, mode='same')

plt.figure()
plt.plot(masked_Ha_reg, data_convolved, 'x')
plt.xlim(6400, 6700)
#%%
""" Gaussian simulated data """
x_array = np.linspace(1,100,50)

amp1 = 100
cen1 = 50
sig1 = 10

y_array_gauss = amp1*(1/(sig1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x_array-cen1)/sig1)**2)))

# creating some noise to add the the y-axis data
y_noise_gauss = (np.exp((np.random.ranf(50))))/5
y_array_gauss += y_noise_gauss
##%%
plt.plot(x_array,y_array_gauss, 'x')
plt.show()

#%%
""" Simulated Voigt profile _3voigt (18 params); _3voigt c_fit """
x_array = np.linspace(6000,7200,5000)

ampG1, cenG1, sigG1 = .1, 6200, 50
ampL1, cenL1, widL1 = .1, 6200, 50
ampG2, cenG2, sigG2 = .3, 6500, 50
ampL2, cenL2, widL2 = .3, 6500, 50
ampG3, cenG3, sigG3 = .1, 6800, 50
ampL3, cenL3, widL3 = .1, 6800, 50


y_array_3voigt = -((ampG1*(1/(sigG1*(np.sqrt(2*np.pi))))*(np.exp(-((x_array-cenG1)**2)/((2*sigG1)**2)))) +\
                ((ampL1*widL1**2/((x_array-cenL1)**2+widL1**2)) ) +\
                (ampG2*(1/(sigG2*(np.sqrt(2*np.pi))))*(np.exp(-((x_array-cenG2)**2)/((2*sigG2)**2)))) +\
                ((ampL2*widL2**2/((x_array-cenL2)**2+widL2**2)) ) +\
                (ampG3*(1/(sigG3*(np.sqrt(2*np.pi))))*(np.exp(-((x_array-cenG3)**2)/((2*sigG3)**2)))) +\
                ((ampL3*widL3**2/((x_array-cenL3)**2+widL3**2)) ))+1

# creating some noise to add the the y-axis data
y_noise_3voigt = (((np.random.ranf(5000))))/500
y_array_3voigt += y_noise_3voigt

#plt.figure()
#plt.plot(x_array,y_array_3voigt,'x')
#plt.show()
##%%
""" Curvefitting triplet Voigt; _3voigt (18 params) """

def _3voigt(x, ampG1, cenG1, sigG1, ampL1, cenL1, widL1, ampG2, cenG2, sigG2, ampL2, cenL2, widL2, \
      ampG3, cenG3, sigG3, ampL3, cenL3, widL3):
    return -((ampG1*(1/(sigG1*(np.sqrt(2*np.pi))))*(np.exp(-((x_array-cenG1)**2)/((2*sigG1)**2)))) +\
      ((ampL1*widL1**2/((x_array-cenL1)**2+widL1**2)) ) +\
      (ampG2*(1/(sigG2*(np.sqrt(2*np.pi))))*(np.exp(-((x_array-cenG2)**2)/((2*sigG2)**2)))) +\
      ((ampL2*widL2**2/((x_array-cenL2)**2+widL2**2)) ) +\
      (ampG3*(1/(sigG3*(np.sqrt(2*np.pi))))*(np.exp(-((x_array-cenG3)**2)/((2*sigG3)**2)))) +\
      ((ampL3*widL3**2/((x_array-cenL3)**2+widL3**2)) ))+1

p0 = [ampG1+100, cenG1, sigG1+10, ampL1, cenL1+100, widL1, ampG2+10, cenG2, sigG2+10, ampL2, cenL2, widL2, \
      ampG3+10, cenG3, sigG3+10, ampL3, cenL3, widL3]
# Reasonable variation in p0 still produces a good fit 
popt_3voigt, cov_3voigt = opt.curve_fit(_3voigt, x_array, y_array_3voigt, p0)

for c in zip(popt_3voigt, np.sqrt(np.diag(cov_3voigt))):
    print("%.8f pm %.3g" % (c[0], c[1]))

# Gives both gaussian and lorentzian component parameters (not general voigt parameters!)

plt.figure()
plt.plot(x_array, y_array_3voigt, 'x')
plt.plot(x_array, _3voigt(x_array, *popt_3voigt))
plt.show()

#%%
""" Simulated Lorentzian Triplet; Lorentzian c_fit """
x_array = np.linspace(1,300,250)

amp1 = 50
cen1 = 100
wid1 = 5

amp2 = 100
cen2 = 150
wid2 = 10

amp3 = 50
cen3 = 200
wid3 = 5

y_array_3lorentz = -((amp1*wid1**2/((x_array-cen1)**2+wid1**2)) + \
                    (amp2*wid2**2/((x_array-cen2)**2+wid2**2)) +\
                     (amp3*wid3**2/((x_array-cen3)**2+wid3**2)))+1

# creating some noise to add the the y-axis data
y_noise_3lorentz = (((np.random.ranf(250))))*5
y_array_3lorentz += y_noise_3lorentz

plt.figure()
plt.plot(x_array, y_array_3lorentz, 'x', label = "Simulated Lorentzian profile")
##%%
def _3Lorentzian(x, amp1, cen1, wid1, amp2,cen2,wid2, amp3,cen3,wid3):
    return -((amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
            (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                (amp3*wid3**2/((x-cen3)**2+wid3**2)))+1

popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, x_array, y_array_3lorentz, \
                                            p0=[amp1, cen1, wid1, amp2, cen2, wid2, amp3, cen3, wid3])

for c in zip(popt_3lorentz, np.sqrt(np.diag(cov_3lorentz))):
    print("%.8f pm %.3g" % (c[0], c[1]))

plt.plot(x_array, _3Lorentzian(x_array, *popt_3lorentz), label = "Lorentzian c_fit")
plt.legend()
plt.show()
#%%
""" Simulated Voigt profile; Lorentzian c_fit """
x_array = np.linspace(1,300,250)

ampG1, cenG1, sigG1 = 100, 100, 25
ampL1, cenL1, widL1 = 50, 100, 25
ampG2, cenG2, sigG2 = 100, 150, 25
ampL2, cenL2, widL2 = 50, 150, 25
ampG3, cenG3, sigG3 = 100, 200, 5
ampL3, cenL3, widL3 = 50, 200, 5


y_array_3voigt = -((ampG1*(1/(sigG1*(np.sqrt(2*np.pi))))*(np.exp(-((x_array-cenG1)**2)/((2*sigG1)**2)))) +\
                ((ampL1*widL1**2/((x_array-cenL1)**2+widL1**2)) ) +\
                (ampG2*(1/(sigG2*(np.sqrt(2*np.pi))))*(np.exp(-((x_array-cenG2)**2)/((2*sigG2)**2)))) +\
                ((ampL2*widL2**2/((x_array-cenL2)**2+widL2**2)) ) +\
                (ampG3*(1/(sigG3*(np.sqrt(2*np.pi))))*(np.exp(-((x_array-cenG3)**2)/((2*sigG3)**2)))) +\
                ((ampL3*widL3**2/((x_array-cenL3)**2+widL3**2)) ))+1

# creating some noise to add the the y-axis data
y_noise_3voigt = (((np.random.ranf(250))))/500
y_array_3voigt += y_noise_3voigt

#p0 = [100, 100, 15, 100, 15, 50, 100, 20, 50]
popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, x_array, y_array_3voigt, \
                                            p0=[65, 105, 10, 60, 150, 5, 60, 200, 5])

for c in zip(popt_3lorentz, np.sqrt(np.diag(cov_3lorentz))):
    print("%.8f pm %.3g" % (c[0], c[1]))


plt.figure()
plt.plot(x_array,y_array_3voigt,'x', label = "Voigt simulated data")
plt.plot(x_array, _3Lorentzian(x_array, *popt_3lorentz), label = "Lorentzian c_fit")
plt.legend()
plt.show()
#%%
""" Simulated Lorentzian Triplet; VoigtNew c_fit """
x_array = np.linspace(1,300,250)

amp1 = 50
cen1 = 100
wid1 = 5

amp2 = 100
cen2 = 150
wid2 = 10

amp3 = 50
cen3 = 200
wid3 = 5

y_array_3lorentz = -((amp1*wid1**2/((x_array-cen1)**2+wid1**2)) + \
                    (amp2*wid2**2/((x_array-cen2)**2+wid2**2)) +\
                     (amp3*wid3**2/((x_array-cen3)**2+wid3**2)))+1

# creating some noise to add the the y-axis data
y_noise_3lorentz = (((np.random.ranf(250))))*5
y_array_3lorentz += y_noise_3lorentz

plt.figure()
plt.plot(x_array, y_array_3lorentz, 'x', label = "Simulated Lorentzian profile")
##%%
def VoigtNew(x, A1, centre1, sigma1, gamma1, A2 , centre2, sigma2, gamma2, A3 , centre3, sigma3, gamma3):
    return -((lm.models.voigt(x, A1, centre1, sigma1, gamma1)) \
             + (lm.models.voigt(x, A2, centre2, sigma2, gamma2)) \
             + (lm.models.voigt(x, A3 , centre3, sigma3, gamma3)))+1

p0 = np.array([50, 100, 10, 1.85032, 100, 150, 10, 1.85032, 50, 200, 10, 1.85032])
popt_VoigtNew, cov_VoigtNew = opt.curve_fit(VoigtNew, x_array, y_array_3lorentz, p0)

for c in zip(popt_VoigtNew, np.sqrt(np.diag(cov_VoigtNew))):
    print("%.8f pm %.3g" % (c[0], c[1]))

plt.plot(x_array, VoigtNew(x_array, *popt_VoigtNew), label = "VoigtNew c_fit")
plt.legend()
plt.show()













""" Ran straight after Part 4 cell; Overlays a voigt profile (Not a fit, just a separate function) """
##plt.figure()
#x_test = np.arange(6000,7200,0.01)
#
#p0Vt1 = np.array([1.5, 6530, 2.6, 1.85032])
#p0Vt2= np.array([2.5, 6565, 2.6, 1.85032])
#p0Vt3= np.array([1.5, 6600, 2.6, 1.85032])
#
#y_test = VoigtNew(x_test, p0Vt1[0], p0Vt1[1], p0Vt1[2], p0Vt1[3], p0Vt2[0], p0Vt2[1], p0Vt2[2], p0Vt2[3], p0Vt3[0], p0Vt3[1], p0Vt3[2], p0Vt3[3])
#plt.plot(x_test, y_test,label='data')
#plt.show()
###%% # Curvefits the above function 
##p1 = np.array([1, 6535, 3, 2])
##p2= np.array([2, 6550, 3, 2])
##p3= np.array([1, 6575, 3, 2])
##p, cov = opt.curve_fit(VoigtNew, x_test, y_test, p0 = [p1, p2, p3])
##plt.plot(x_test, VoigtNew(x_test, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], \
##        p[8], p[9], p[10], p[11]), linewidth=2, color = "royalblue", label = "Voigt fit")
##for c in zip(p, np.sqrt(np.diag(cov))):
##    print("%.8f pm %.3g" % (c[0], c[1]))
###plt.legend()
##plt.show()
##%%

