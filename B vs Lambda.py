# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 18:26:26 2021

@author: 44743
"""
import numpy as np
import matplotlib.pyplot as plt

#%%
# uses wavelength to give B
def delta_lambda(wavelength):
    A = np.square(6562.8/6564)
    B = wavelength/(20.2*A)
    return B

colour = 'royalblue'
label = 'Linear zeeman'
x1 = np.arange(-500, 500, 1)
plt.figure('X')
plt.plot(x1, delta_lambda(abs(x1)), color = colour, label = label)
plt.plot(np.zeros(10000), np.linspace(0,25,10000), color = colour)#, label = label)
plt.xlabel('Wavelength $[\AA]$', fontsize = "13")
plt.ylabel('Magnetic field strength $[MG]$', fontsize = "13")
plt.legend()
#plt.grid()
plt.show()
#%%
# Rudimentary plotting 
# uses B to give wavelength
def delta_lambda(B):
    wav_list = []
    for i in range(len(B)):
        b = B[i]
        
        A = np.square(6562.8/6564)
        ml = np.array([-1,0,1])
        wavelength = [k*b*(20.2*A) for k in ml]
        wav_list.append(np.array(wavelength))
        #wav_list = np.array(wav_list)
    return wav_list

# uses B as input to give wavelength
def delta_lambda_quad(B):
    wav_list = []
    for i in range(len(B)):
        b = B[i]
        A = -(4.98*10**-11)*np.square(6562.8)*(3**4)
        ml = np.array([-2,-1,1*10**-13,1,2])
#        C = [A*(1+(i/abs(i))*np.square(i))*b for i in ml if i != 0]
#        C = [A*(1+np.square(i))*b for i in ml if i == 0]
        C = [A*(1+(i/abs(i))*np.square(i))*np.square(b) for i in ml]

        wav_list.append(np.array(C))
    return wav_list
        
##%%
colour = 'orange'
label = 'Quadratic zeeman'
colour = 'royalblue'
label = 'Linear zeeman'

x1 = np.arange(0, 25, 0.5)

plt.plot(delta_lambda_quad((x1)), np.zeros(50)+x1)
plt.plot(delta_lambda(x1), np.zeros(50)+x1) # This is for delta_lambda linear outputting wavelength
    
A = delta_lambda_quad((x1))
B = delta_lambda(x1)

#plt.plot(np.zeros(10000), np.linspace(0,25,10000), color = red, label = 'y axis')
plt.xlabel('Wavelength $[\AA]$', fontsize = "13")
plt.ylabel('Magnetic field strength $[MG]$', fontsize = "13")
plt.legend()
#plt.grid()
plt.show()

#%%
# Complete plotting Linear and Quadratic effects 
x1 = np.arange(0, 25, 0.5)

# Linear Zeeman effect
def delta_lambda(B):
    wav_list = []
    for i in range(len(B)):
        b = B[i]
        
        A = np.square(6562.8/6564)
        ml = np.array([-1,0,1])
        wavelength = [k*b*(20.2*A) for k in ml]
        wav_list.append(np.array(wavelength))
        #wav_list = np.array(wav_list)
    return wav_list

# Quadratic Zeeman effect
def delta_lambda_quad(B):
    wav_list = []
    for i in range(len(B)):
        b = B[i]
        A = -(4.98*10**-11)*np.square(6562.8)*(3**4)
        ml = np.array([-2,-1,1*10**-13,1,2])
#        C = [A*(1+(i/abs(i))*np.square(i))*b for i in ml if i != 0]
#        C = [A*(1+np.square(i))*b for i in ml if i == 0]
        C = [A*(1+(i/abs(i))*np.square(i))*np.square(b)/20 for i in ml] # scaling factor of 20 applied

        wav_list.append(np.array(C))
    return wav_list

# Total/ combined Zeeman effect
def delta_lambda_total(B):
    wav_list = []
    for i in range(len(B)):
#        b = B[i]

        ml = np.array([0,1,2])
        
#        for j in x1:
        sublist = []
        for k in ml:
            wavelength = delta_lambda(x1)[i][k] + delta_lambda_quad(x1)[i] # an array of 5 elements
            sublist.append(wavelength)
            
        wav_list.append((sublist)) # np.array(sublist)
        #wav_list = np.array(wav_list)
    return wav_list

A = delta_lambda_quad((x1))
B = delta_lambda(x1)
C = delta_lambda_total(x1)

delt_lambda_total_list = []
for i in range(len(delta_lambda_total(x1))):
    compressed_array = np.concatenate(C[i], axis = 0) # compressed array of all 15 values per B value
    delt_lambda_total_list.append(compressed_array+6562.8) # compressed array centered by default to zero
                                                    # needs translation to centre at lambda0 = 6562.8

plt.figure()
plt.plot(delt_lambda_total_list, np.zeros(50)+x1)
#plt.plot(np.zeros(10000), np.linspace(0,25,10000), color = red, label = 'y axis')
plt.xlabel('Wavelength $[\AA]$', fontsize = "13")
plt.ylabel('Magnetic field strength $[MG]$', fontsize = "13")
#plt.legend()
plt.grid()
plt.show()    
#%%
plt.figure()
plt.plot(xp_triplet, abs(popt_3lorentz[1])*yp_triplet,'x', label = "WD Ha data")
plt.plot(xp_triplet, abs(popt_3lorentz[1])*_3Lorentzian(xp_triplet, *popt_3lorentz), linewidth=2, color = "orange", \
                              label = "Lorentzian c_fit")
plt.plot(delt_lambda_total_list, np.zeros(50)+x1)
plt.xlabel('Wavelength $[\AA]$', fontsize = "13")
plt.ylabel('Magnetic field strength $[MG]$', fontsize = "13")
#plt.legend()
plt.grid()
plt.show()    









