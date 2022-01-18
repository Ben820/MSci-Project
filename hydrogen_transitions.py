# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 12:37:09 2022

@author: rr3
"""

import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm
import pylab
import glob
import os
import pandas as pd
import re
import astropy.io.fits as fits

#%%
column_names = ['NumDataPoints', 'Bfield', 'Wavelength', 'Final State']
transitions_df = pd.read_csv(r'C:\Local\Akshay Laptop Backup 29Dec2019\Imperial\MSci Project\DESI\hydrogencolumns.csv', skiprows=1, names=column_names)
#transitions_array = transitions_df[["Wavelength"]].to_numpy() 
#%%
#transitions_array = float(transitions_array)
#Halphawlengths = str(np.where(transitions_array == min(transitions_array, key=lambda x:abs(x-6562.8)))[0])
#%%
#transitions_df = transitions_df.drop(4299,axis=0)
#transitions_df = transitions_df.drop(4219,axis=0)
#%%
notnan = transitions_df['NumDataPoints'].notnull()
#%%
#some trial code to test out the method before applying to real data
"""
wavelength_trial = transitions_df.loc[1:18,'Wavelength']
b_trial = transitions_df.loc[1:18,'Bfield']
wavelength_trial = wavelength_trial.tolist()
b_trial = b_trial.tolist()
wavelength_trial = np.array(wavelength_trial)
b_trial = np.array(b_trial)
wav_array = [float(i) for i in wavelength_trial]
b_array = [float(i) for i in b_trial]
b_array = np.log10(b_array)
"""
#%%
"""
Indices of the H-alpha transitions (start of the list)
"""
alpha_a = (transitions_df[(transitions_df['Bfield']  == '2p-1') & (transitions_df['Final State'] == '3s0')].index.tolist())
alpha_b = (transitions_df[(transitions_df['Bfield']  == '2s0') & (transitions_df['Final State'] == '3p+1')].index.tolist())
alpha_c = (transitions_df[(transitions_df['Bfield']  == '2p+1') & (transitions_df['Final State'] == '3d+2')].index.tolist())
alpha_d = (transitions_df[(transitions_df['Bfield']  == '2p0') & (transitions_df['Final State'] == '3d+1')].index.tolist())
alpha_e = (transitions_df[(transitions_df['Bfield']  == '2p-1') & (transitions_df['Final State'] == '3d0')].index.tolist())
alpha_f = (transitions_df[(transitions_df['Bfield']  == '2p0') & (transitions_df['Final State'] == '3s0')].index.tolist())
alpha_g = (transitions_df[(transitions_df['Bfield']  == '2p-1') & (transitions_df['Final State'] == '3d-1')].index.tolist())
alpha_h = (transitions_df[(transitions_df['Bfield']  == '2s0') & (transitions_df['Final State'] == '3p0')].index.tolist())
alpha_i = (transitions_df[(transitions_df['Bfield']  == '2p0') & (transitions_df['Final State'] == '3d0')].index.tolist())
alpha_j = (transitions_df[(transitions_df['Bfield']  == '2p+1') & (transitions_df['Final State'] == '3s0')].index.tolist())
alpha_k = (transitions_df[(transitions_df['Bfield']  == '2s0') & (transitions_df['Final State'] == '3p-1')].index.tolist())
alpha_l = (transitions_df[(transitions_df['Bfield']  == '2p-1') & (transitions_df['Final State'] == '3d-2')].index.tolist())
alpha_m = (transitions_df[(transitions_df['Bfield']  == '2p0') & (transitions_df['Final State'] == '3d-1')].index.tolist())
alpha_n = (transitions_df[(transitions_df['Bfield']  == '2p+1') & (transitions_df['Final State'] == '3d0')].index.tolist())
#%%
"""
Number of data points in each transition
"""
num_a = np.int(transitions_df.iloc[alpha_a[0]][0])
num_b = np.int(transitions_df.iloc[alpha_b[0]][0])
num_c = np.int(transitions_df.iloc[alpha_c[0]][0])
num_d = np.int(transitions_df.iloc[alpha_d[0]][0])
num_e = np.int(transitions_df.iloc[alpha_e[0]][0])
num_f = np.int(transitions_df.iloc[alpha_f[0]][0])
num_g = np.int(transitions_df.iloc[alpha_g[0]][0])
num_h = np.int(transitions_df.iloc[alpha_h[0]][0])
num_i = np.int(transitions_df.iloc[alpha_i[0]][0])
num_j = np.int(transitions_df.iloc[alpha_j[0]][0])
num_k = np.int(transitions_df.iloc[alpha_k[0]][0])
num_l = np.int(transitions_df.iloc[alpha_l[0]][0])
num_m = np.int(transitions_df.iloc[alpha_m[0]][0])
num_n = np.int(transitions_df.iloc[alpha_n[0]][0])
#%%
"""
End index of data describing a Halpha transition
"""
end_a = alpha_a[0]+num_a
end_b = alpha_b[0]+num_b
end_c = alpha_c[0]+num_c
end_d = alpha_d[0]+num_d
end_e = alpha_e[0]+num_e
end_f = alpha_f[0]+num_f
end_g = alpha_g[0]+num_g
end_h = alpha_h[0]+num_h
end_i = alpha_i[0]+num_i
end_j = alpha_j[0]+num_j
end_k = alpha_k[0]+num_k
end_l = alpha_l[0]+num_l
end_m = alpha_m[0]+num_m
end_n = alpha_n[0]+num_n
#%%
"""
Wavelength and B-field array for all the Halpha transitions
"""
wavelength_a = transitions_df.loc[alpha_a[0]+2:end_a,'Wavelength']
b_a = transitions_df.loc[alpha_a[0]+2:end_a,'Bfield']
wavelength_b = transitions_df.loc[alpha_b[0]+2:end_b,'Wavelength']
b_b = transitions_df.loc[alpha_b[0]+2:end_b,'Bfield']
wavelength_c = transitions_df.loc[alpha_c[0]+2:end_c,'Wavelength']
b_c = transitions_df.loc[alpha_c[0]+2:end_c,'Bfield']
wavelength_d = transitions_df.loc[alpha_d[0]+2:end_d,'Wavelength']
b_d = transitions_df.loc[alpha_d[0]+2:end_d,'Bfield']
wavelength_e = transitions_df.loc[alpha_e[0]+2:end_e,'Wavelength']
b_e = transitions_df.loc[alpha_e[0]+2:end_e,'Bfield']
wavelength_f = transitions_df.loc[alpha_f[0]+2:end_f,'Wavelength']
b_f = transitions_df.loc[alpha_f[0]+2:end_f,'Bfield']
wavelength_g = transitions_df.loc[alpha_g[0]+2:end_g,'Wavelength']
b_g = transitions_df.loc[alpha_g[0]+2:end_g,'Bfield']
wavelength_h = transitions_df.loc[alpha_h[0]+2:end_h,'Wavelength']
b_h = transitions_df.loc[alpha_h[0]+2:end_h,'Bfield']
wavelength_i = transitions_df.loc[alpha_i[0]+2:end_i,'Wavelength']
b_i = transitions_df.loc[alpha_i[0]+2:end_i,'Bfield']
wavelength_j = transitions_df.loc[alpha_j[0]+2:end_j,'Wavelength']
b_j = transitions_df.loc[alpha_j[0]+2:end_j,'Bfield']
wavelength_k = transitions_df.loc[alpha_k[0]+2:end_k,'Wavelength']
b_k = transitions_df.loc[alpha_k[0]+2:end_k,'Bfield']
wavelength_l = transitions_df.loc[alpha_l[0]+2:end_l,'Wavelength']
b_l = transitions_df.loc[alpha_l[0]+2:end_l,'Bfield']
wavelength_m = transitions_df.loc[alpha_m[0]+2:end_m,'Wavelength']
b_m = transitions_df.loc[alpha_m[0]+2:end_m,'Bfield']
wavelength_n = transitions_df.loc[alpha_n[0]+2:end_n,'Wavelength']
b_n = transitions_df.loc[alpha_n[0]+2:end_n,'Bfield']
#%%
#wav_list = [wavelength_a, wavelength_e, wavelength_f, wavelength_g, wavelength_h, wavelength_i, wavelength_k, wavelength_l, wavelength_m]
#b_list = [b_a, b_e, b_f, b_g, b_h, b_i, b_k, b_l, b_m]
#%%
wavelength_a = wavelength_a.tolist()
wavelength_b = wavelength_b.tolist()
wavelength_c = wavelength_c.tolist()
wavelength_d = wavelength_d.tolist()
wavelength_e = wavelength_e.tolist()
wavelength_f = wavelength_f.tolist()
wavelength_g = wavelength_g.tolist()
wavelength_h = wavelength_h.tolist()
wavelength_i = wavelength_i.tolist()
wavelength_j = wavelength_j.tolist()
wavelength_k = wavelength_k.tolist()
wavelength_l = wavelength_l.tolist()
wavelength_m = wavelength_m.tolist()
wavelength_n = wavelength_n.tolist()
b_a = b_a.tolist()
b_b = b_b.tolist()
b_c = b_c.tolist()
b_d = b_d.tolist()
b_e = b_e.tolist()
b_f = b_f.tolist()
b_g = b_g.tolist()
b_h = b_h.tolist()
b_i = b_i.tolist()
b_j = b_j.tolist()
b_k = b_k.tolist()
b_l = b_l.tolist()
b_m = b_m.tolist()
b_n = b_n.tolist()
#%%
wav_array_a = [float(i) for i in wavelength_a]
b_array_a = [float(i) for i in b_a]
#b_array_a = np.log10(b_array_a)

wav_array_b = [float(i) for i in wavelength_b]
b_array_b = [float(i) for i in b_b]
#b_array_b = np.log10(b_array_b)

wav_array_c = [float(i) for i in wavelength_c]
b_array_c = [float(i) for i in b_c]
#b_array_c = np.log10(b_array_c)

wav_array_d = [float(i) for i in wavelength_d]
b_array_d = [float(i) for i in b_d]
#b_array_d = np.log10(b_array_d)

wav_array_e = [float(i) for i in wavelength_e]
b_array_e = [float(i) for i in b_e]
#b_array_e = np.log10(b_array_e)

wav_array_f = [float(i) for i in wavelength_f]
b_array_f = [float(i) for i in b_f]
#b_array_f = np.log10(b_array_f)

wav_array_g = [float(i) for i in wavelength_g]
b_array_g = [float(i) for i in b_g]
#b_array_g = np.log10(b_array_g)

wav_array_h = [float(i) for i in wavelength_h]
b_array_h = [float(i) for i in b_h]
#b_array_h = np.log10(b_array_h)

wav_array_i = [float(i) for i in wavelength_i]
b_array_i = [float(i) for i in b_i]
#b_array_i = np.log10(b_array_i)

wav_array_j = [float(i) for i in wavelength_j]
b_array_j = [float(i) for i in b_j]
#b_array_j = np.log10(b_array_j)

wav_array_k = [float(i) for i in wavelength_k]
b_array_k = [float(i) for i in b_k]
#b_array_k = np.log10(b_array_k)

wav_array_l = [float(i) for i in wavelength_l]
b_array_l = [float(i) for i in b_l]
#b_array_l = np.log10(b_array_l)

wav_array_m = [float(i) for i in wavelength_m]
b_array_m = [float(i) for i in b_m]
#b_array_m = np.log10(b_array_m)

wav_array_n = [float(i) for i in wavelength_n]
b_array_n = [float(i) for i in b_n]
#b_array_n = np.log10(b_array_n)
#%%
labels=['a','b','c','d','e','f','g','h','i','j','k','l','m','n']
plt.grid()
#plt.xlim(0,100)
#plt.xticks(np.arange(min(wavelength_trial), max(wavelength_trial)+1, 1.0))
plt.xlabel('Wavelength $[\AA]$')
plt.ylabel('logB [G]')
plt.plot(wav_array_a,b_array_a)
plt.plot(wav_array_b,b_array_b)
plt.plot(wav_array_c,b_array_c)
plt.plot(wav_array_d,b_array_d)
plt.plot(wav_array_e,b_array_e)
plt.plot(wav_array_f,b_array_f)
plt.plot(wav_array_g,b_array_g)
plt.plot(wav_array_h,b_array_h)
plt.plot(wav_array_i,b_array_i)
plt.plot(wav_array_j,b_array_j)
plt.plot(wav_array_k,b_array_k)
plt.plot(wav_array_l,b_array_l)
plt.plot(wav_array_m,b_array_m)
plt.plot(wav_array_n,b_array_n)
plt.legend(labels=labels)
plt.show()