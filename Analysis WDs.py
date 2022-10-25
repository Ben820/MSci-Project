# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 19:12:42 2022

@author: 44743
"""
""" Script for analysis """
import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm
import pylab
import glob
import pandas as pd 
from collections import Counter
#%%
column_names = ["Filename", "Class", "B_alpha", "B error alpha", "lambda0", "lambda0 error", "begin", "finish", "start1", "end1", "start2", "end2", "start3", "end3", "Notes"]
datasets = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Detailed class catalogue 2.csv', skiprows = 1, names = column_names)
filename_list = []
#tot_list_ = []

begin_list = []
finish_list = []
start1_list = []
end1_list = []
start2_list = []
end2_list = []
start3_list = []
end3_list = []

tot__list = [begin_list, finish_list, start1_list, end1_list, start2_list, end2_list, start3_list, end3_list]   

for i in range(len(datasets)):
    if datasets.Class[i] == 1:
        filename_list.append(datasets.Filename[i])
#        for j in range(len(cut_regions)):
#            tot__list[j].append(datasets.begin[j])
        begin_list.append(datasets.begin[i])
        finish_list.append(datasets.finish[i])
        start1_list.append(datasets.start1[i])
        end1_list.append(datasets.end1[i])
        start2_list.append(datasets.start2[i])
        end2_list.append(datasets.end2[i])
        start3_list.append(datasets.start3[i])
        end3_list.append(datasets.end3[i])

Bvalues = np.array(datasets.B_alpha)

for i in range(len(Bvalues)):
    if Bvalues[i] == np.isnan:
        print('yes')

nans = np.argwhere(np.isnan(Bvalues))[:,0]

Bvalues_notnan = np.delete(Bvalues, nans)


#%%        
#plt.figure()
#plt.hist(abs(Bvalues_notnan),bins = 50)
#plt.show()

plt.figure()
plt.hist(abs(Bvalues),bins = 50)
#plt.xscale('log', basex = 10) # bins = 200 for log scale (50 otherwise)
plt.xlabel('Magnetic Field B $[MG]$', fontsize = "13")
plt.ylabel('Number of WDs', fontsize = "13")
plt.show()


#%%
data = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Intermediate Fields w out left pi.csv', \
                   skiprows = 0)
##%%
B_lin = []
B_trans = []
[B_lin.append(data['LinearFieldEstimate'][i]) for i in range(0, len(data))]
[B_trans.append(data['B (sigma)'][i]) for i in range(0, len(data))]
B_lin = np.array(B_lin)
B_trans = np.array(B_trans)

B_lin.sort()
B_trans.sort()

a = np.arange(15, 33, 0.05)
b = np.linspace(15, 33, num = len(B_lin))


y_err = 0.08*B_trans
x_err = 0.08*B_lin

plt.figure()
plt.plot(abs(B_lin), abs(B_trans), 'x')
plt.plot(a, a)
plt.plot(b, b, 'o')
plt.show()

residuals = abs(B_trans) - abs(B_lin)
#plt.figure()
#plt.plot()
#%%
fig, axs = plt.subplots(2, gridspec_kw={'height_ratios': [2.5, 1]})
#fig.suptitle(f"Lorentzian fit {filename} \n B = {popt_3lorentz[1]} +/- {np.sqrt(np.diag(cov_3lorentz))[1]}", \
#             fontsize = "13")
axs[0].errorbar(abs(B_lin), abs(B_trans), yerr = y_err, xerr = x_err, fmt = 'x', color = 'darkblue')#, label = "WD Ha data")
axs[0].plot(a,a, linewidth=2, color = "tomato", \
                              label = "y = x")
#axs[0].text(r"B = {popt_3lorentz}")
for ax in axs.flat:
#    ax.set_xlabel('$B_{linear}$ (MG)', fontsize = "15", labelpad = 5)
#    ax.set_ylabel('$B_{inter}$ (MG)', fontsize = "15", labelpad = 13)
    ax.set_xlabel('Linear field (MG)', fontsize = "15", labelpad = 5)
    ax.set_ylabel('Interpolated field (MG)', fontsize = "15", labelpad = 13)
    ax.tick_params(axis='both', labelsize=14)

for ax in axs.flat:
    ax.label_outer()
axs[1].set_ylabel('Deviation', fontsize = "15", labelpad = 23)
axs[0].grid()
#""" Voigt Residuals Plot """
axs[1].plot(abs(B_lin), residuals, 'x', linewidth=2, color = 'darkblue')#, label = "Lorentzian fit Residuals")
axs[1].plot(a, a*0/a, linewidth = 2, color = 'tomato')
axs[0].legend()
#axs[1].legend()
axs[1].grid()
#plt.savefig("Gaussian Residuals", dpi = 1000)
#plt.savefig(f'{filename}.png', dpi = 1000)
plt.show()
#%%
""" Plots alpha beta straight line comparison """
alphabeta_comp = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\H beta alpha comparison csv.csv', \
                   skiprows = 0)

B_alpha = []
B_beta = []
B_alpha_err = []
B_beta_err = []

[B_alpha.append(alphabeta_comp['B alpha'][i]) for i in range(0, len(alphabeta_comp))]
[B_beta.append(alphabeta_comp['B alpha BETA'][i]) for i in range(0, len(alphabeta_comp))]
[B_alpha_err.append(alphabeta_comp['B error alpha '][i]) for i in range(0, len(alphabeta_comp))]
[B_beta_err.append(alphabeta_comp['B error alpha BETA'][i]) for i in range(0, len(alphabeta_comp))]

Balph = abs(np.array(B_alpha))
Bbet = abs(np.array(B_beta))
Balpherr = abs(np.array(B_alpha_err))
Bbeterr = abs(np.array(B_beta_err))

x = Balph
y = Bbet

y = y[np.logical_not(np.isnan(x))]
x = x[np.logical_not(np.isnan(x))]

#plt.figure()
#plt.plot(abs(Balph), abs(Bbet), 'x')
#plt.plot(abs(Bbet), abs(Bbet), color = 'orange')
#plt.show()

residuals = abs(Bbet) - abs(Balph)

res = residuals[np.logical_not(np.isnan(residuals))]



markerline = np.zeros((1, len(Bbet)))+7

fig, axs = plt.subplots(2, gridspec_kw={'height_ratios': [2.5, 1]})
#fig.suptitle(f"Lorentzian fit {filename} \n B = {popt_3lorentz[1]} +/- {np.sqrt(np.diag(cov_3lorentz))[1]}", \
#             fontsize = "13")
axs[0].errorbar(abs(Balph), abs(Bbet), yerr = np.array(B_beta_err), xerr = np.array(B_alpha_err), fmt = 'x', color = 'darkblue')#, label = "WD Ha data")
axs[0].plot(y, y, linewidth=2, color = "tomato", \
                              label = "y = x")
axs[0].plot(markerline[0], np.linspace(0,12, num = len(markerline[0])), 'royalblue', label = 'Field strength\ncut-off')
#axs[0].text(r"B = {popt_3lorentz}")
axs[0].tick_params(axis='both', labelsize=11)
axs[1].tick_params(axis='both', labelsize=11)

for ax in axs.flat:
#    ax.set_xlabel('$B_{linear}$ (MG)', fontsize = "15", labelpad = 5)
#    ax.set_ylabel('$B_{inter}$ (MG)', fontsize = "15", labelpad = 13)
    ax.set_xlabel(r'$\mathrm{H \alpha}$ field strength (MG)', fontsize = "13", labelpad = 5)
    ax.set_ylabel(r'$\mathrm{H \beta}$ field strength (MG)', fontsize = "13", labelpad = 5)
#    ax.tick_params(axis='both', labelsize=12)

for ax in axs.flat:
    ax.label_outer()
axs[1].set_ylabel('Deviation', fontsize = "13", labelpad = 5)
axs[0].grid()
#""" Voigt Residuals Plot """
axs[1].errorbar(abs(Balph), residuals, yerr = np.array(B_beta_err), xerr = np.array(B_alpha_err), \
   fmt = 'x', linewidth=2, color = 'darkblue')#, label = "Lorentzian fit Residuals")
axs[1].plot(abs(Bbet), abs(Bbet)*0/Bbet, linewidth = 2, color = 'tomato')
axs[0].legend(prop={'size': 11})
#axs[1].legend()
axs[1].grid()
#plt.savefig("Gaussian Residuals", dpi = 1000)
#plt.savefig(f'{filename}.pdf'), dpi = 1000)
#plt.savefig('H Alpha Beta comparison errors 2.pdf', bbox_inches = 'tight')
plt.show()
#%%
#\mathrm{\AA}






