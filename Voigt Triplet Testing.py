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

#%%
""" Simulated Voigt profile _3voigt (18 params) """
x_array = np.linspace(6000,7200,5000)

ampG1, cenG1, sigG1 = .1, 6200, 5
ampL1, cenL1, widL1 = .1, 6200, 5
ampG2, cenG2, sigG2 = .1, 6500, 5
ampL2, cenL2, widL2 = .1, 6500, 5
ampG3, cenG3, sigG3 = .1, 6800, 5
ampL3, cenL3, widL3 = .1, 6800, 5


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

p0 = [ampG1, cenG1, sigG1, ampL1, cenL1, widL1, ampG2, cenG2, sigG2, ampL2, cenL2, widL2, \
      ampG3, cenG3, sigG3, ampL3, cenL3, widL3]
popt_3voigt, cov_3voigt = opt.curve_fit(_3voigt, x_array, y_array_3voigt, p0)

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

ampG1, cenG1, sigG1 = 100, 100, 5
ampL1, cenL1, widL1 = 50, 100, 5
ampG2, cenG2, sigG2 = 100, 150, 5
ampL2, cenL2, widL2 = 50, 150, 5
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














