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

entries = os.listdir('C:\Local\Akshay Laptop Backup 29Dec2019\Imperial\MSci Project\DESI\DESI_DxH\DESI_DxH')

files = entries.copy()
del files[0] #remove DA spectrum that was accidentally included
del files[113:114] #remove other .py files in this directory

one = files[0:20]
two = files[20:40]
three = files[40:60]
four = files[60:80]
five = files[80:100]
six = files[100:len(files)]
#%%

for i in range(len(one)):
    data = np.genfromtxt(one[i], delimiter=' ')
    wavelength = data[:,0]
    flux = data[:,1]
    error = data[:,2]
    plt.figure()
    plt.grid()
    plt.plot(wavelength,flux)
    
    