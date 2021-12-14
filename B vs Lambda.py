# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 18:26:26 2021

@author: 44743
"""
#%%
import numpy as np

lambda0 = 6562.8

def quadB(ml, B, lambda0):
    A = np.square(lambda0/6564)
    delt_lam_linear = 20.2*A*B
    C = 4.98*10**-23
    delt_lam_quad = C*(3**4)*np.square(lambda0)*(1+np.square(ml))*np.square(B)
    delt_lam = delt_lam_linear+delt_lam_quad
    return delt_lam

ml = np.array([-2,-1,0,1,2])
B = np.arange(0,5,1)
QUad = quadB(ml, 20, lambda0)
#%%
import matplotlib.pyplot as plt

def delta_lambda(y):
    A = np.square(6562.8/6564)
    B = y/(20.2*A)
    return B

x1 = np.arange(-100, 0, 1)
x2 = np.arange(0,100, 1)
plt.plot(x1, -delta_lambda(x1))
plt.plot(x2, delta_lambda(x2))
#%%
plt.figure()
plt.errorbar()