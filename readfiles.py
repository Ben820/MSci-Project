# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 12:26:07 2021

@author: rr3
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm
import os
import re
#%%
entries = os.listdir('C:\Local\Akshay Laptop Backup 29Dec2019\Imperial\MSci Project\DESI\DESI_DxH\DESI_DxH')
#%%
files = entries.copy()
del files[0:3] #remove DA spectrum that was accidentally included
del files[113:-1] #remove other .py files in this directory
#%%
one = files[0:20]
two = files[20:40]
three = files[40:60]
four = files[60:80]
five = files[80:100]
six = files[100:len(files)]
del six[-1]

#%%
"""
Comment out the 'SourceIdx', 'imgIdx', 'imgfile' and 'plt.savefig' lines 
unless you are saving PNG files of the matplotlib spectra plots
"""
plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]
for i in range(len(six)):
    data = np.genfromtxt(six[i], delimiter=' ') #load data
    wavelength = data[:,0] #unpack wavelength
    flux = data[:,1] #unpack flux
    error = data[:,2] #unpack errors
    SourceIdx = str(six[i]) #gets the source name
    imgIdx = re.sub('\.dat$', '', SourceIdx) #formats the source name
    imgfile = imgIdx+'.png' #frames the source name in a format where an image can be saved
    plt.figure()
    plt.grid()
    plt.plot(wavelength,flux) #plot the spectra
    plt.savefig('DESI Spectra/6/'+imgfile) #saves image to directory 
#%%
"""
test_data = np.genfromtxt('DESI_WDJ111940.70+722847.58_bin0p2.dat', delimiter=' ')
test_wavelength = test_data[:,0]
test_flux = test_data[:,1]
test_error = test_data[:,2]

plt.figure("Whole spectrum")
plt.grid()
plt.plot(test_wavelength,test_flux,'x')
plt.rcParams["figure.figsize"] = (1,1)
plt.show()
"""
    