# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 23:17:38 2022

@author: 44743
"""
#%%
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
#%%
""" Bfield vs Mass """
x = datalist[2]
y = datalist[0]
xerr = datalist[3]
yerr = datalist[1]
#y = np.array(datalist[0])
#Bfield = y[np.logical_not(np.isnan(x))]
#x = np.array(datalist[2])
#Mass = x[np.logical_not(np.isnan(x))]


fig = plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(4, 4)
ax_main = plt.subplot(gs[1:3, 0:2])
ax_xDist = plt.subplot(gs[0, 0:2], sharex=ax_main)
ax_yDist = plt.subplot(gs[1:3, 2], sharey=ax_main)
    
ax_main.errorbar(x, y, yerr = yerr, xerr = xerr, fmt = 'x', color = 'darkblue', alpha = 0.2)
ax_main.scatter(x, y, marker = 'x', color = 'darkblue')

ax_main.set_xlabel(r'Mass $ (\mathrm{M_{\odot}}) $' , fontsize = "13")
ax_main.set_ylabel('Magnetic Field Strength' + ' (MG)', fontsize = "13")

ax_main.set_xlim(0.23, 1.46)

# tick params 
ax_main.tick_params(axis='both', labelsize=11)
ax_main.tick_params(axis='both', labelsize=11)
ax_xDist.tick_params(axis='both', labelsize=11)
ax_yDist.tick_params(axis='both', labelsize=11)

# Hides ticks on histogram plots 
plt.setp(ax_xDist.get_xticklabels(), visible=False)
plt.setp(ax_yDist.get_yticklabels(), visible=False)

# Plots histogram on the x axis (top of scatter)
ax_xDist.hist(x, bins = 20, align = 'mid', color = 'darkblue', alpha = 1)
ax_xDist.set_ylabel('Number of MWDs', fontsize = "13")
#ax_xCumDist = ax_xDist.twinx()
#ax_xCumDist.hist(x,bins=100,cumulative=True,histtype='step',density=True,color='r',align='mid')
#ax_xCumDist.tick_params('y', colors='r')
#ax_xCumDist.set_ylabel('cumulative',color='r')

# Plots histogram on the y axis (side of scatter)

# LOG SCALE CASE
hist, bins = np.histogram(y, bins=17) # bins = 21
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
#plt.hist(x, bins=logbins, color = 'darkblue', alpha = 1) # , density = True
ax_main.set_yscale('log')
ax_yDist.hist(y, bins = logbins, orientation = 'horizontal', align = 'mid', color = 'darkblue', alpha = 1)
ax_yDist.set_xlabel('Number of MWDs', fontsize = "13")

# LINEAR SCALE CASE
#ax_yDist.hist(y, bins = 20, orientation = 'horizontal', align='mid')
#ax_yDist.set_xlabel('Number of MWDs', fontsize = "13")

#ax_yCumDist = ax_yDist.twiny()
#ax_yCumDist.hist(y,bins=100,cumulative=True,histtype='step',density=True,color='r',align='mid',orientation='horizontal')
#ax_yCumDist.tick_params('x', colors='r')
#ax_yCumDist.set_xlabel('cumulative',color='r')

plt.show()

#plt.savefig(f'Bfieldmass 2.pdf', bbox_inches = 'tight')

#%%
""" B field vs Temperature """
x = np.array(datalist[4])/1000
y = np.array(datalist[0])
xerr = np.array(datalist[5])/1000
yerr = np.array(datalist[1])
#y = np.array(datalist[0])
#Bfield = y[np.logical_not(np.isnan(x))]
#x = np.array(datalist[2])
#Mass = x[np.logical_not(np.isnan(x))]


fig = plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(4, 4)
ax_main = plt.subplot(gs[1:3, 0:2])
ax_xDist = plt.subplot(gs[0, 0:2], sharex=ax_main)
ax_yDist = plt.subplot(gs[1:3, 2], sharey=ax_main)
    
ax_main.errorbar(x, y, yerr = yerr, xerr = xerr, fmt = 'x', color = 'darkblue', alpha = 0.2)
ax_main.scatter(x, y, marker = 'x', color = 'darkblue')

ax_main.set_xlabel(r'Effective Temperature ($10^3$ K) ' , fontsize = "13") # r'$T_{eff} $
ax_main.set_ylabel('Magnetic Field Strength' + ' (MG)', fontsize = "13")

ax_main.set_xlim(3328/1000, 41365/1000)

# tick params 
ax_main.tick_params(axis='both', labelsize=11)
ax_main.tick_params(axis='both', labelsize=11)
ax_xDist.tick_params(axis='both', labelsize=11)
ax_yDist.tick_params(axis='both', labelsize=11)

# Hides ticks on histogram plots 
plt.setp(ax_xDist.get_xticklabels(), visible=False)
plt.setp(ax_yDist.get_yticklabels(), visible=False)

# Plots histogram on the x axis (top of scatter)
ax_xDist.hist(x, bins = 20, align = 'mid', color = 'darkblue', alpha = 1)
ax_xDist.set_ylabel('Number of MWDs', fontsize = "13")
#ax_xCumDist = ax_xDist.twinx()
#ax_xCumDist.hist(x,bins=100,cumulative=True,histtype='step',density=True,color='r',align='mid')
#ax_xCumDist.tick_params('y', colors='r')
#ax_xCumDist.set_ylabel('cumulative',color='r')

# Plots histogram on the y axis (side of scatter)

# LOG SCALE CASE
hist, bins = np.histogram(y, bins=17) # bins = 21
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
#plt.hist(x, bins=logbins, color = 'darkblue', alpha = 1) # , density = True
ax_main.set_yscale('log')
ax_yDist.hist(y, bins = logbins, orientation = 'horizontal', align = 'mid', color = 'darkblue', alpha = 1)
ax_yDist.set_xlabel('Number of MWDs', fontsize = "13")

# LINEAR SCALE CASE
#ax_yDist.hist(y, bins = 20, orientation = 'horizontal', align='mid')
#ax_yDist.set_xlabel('Number of MWDs', fontsize = "13")

#ax_yCumDist = ax_yDist.twiny()
#ax_yCumDist.hist(y,bins=100,cumulative=True,histtype='step',density=True,color='r',align='mid',orientation='horizontal')
#ax_yCumDist.tick_params('x', colors='r')
#ax_yCumDist.set_xlabel('cumulative',color='r')

plt.show()

#plt.savefig(f'BfieldTemp.pdf', bbox_inches = 'tight')

#%%
""" Mass vs Temperature """
x = np.array(datalist[2])
y = np.array(datalist[4])/1000
xerr = np.array(datalist[3])
yerr = np.array(datalist[5])/1000
#y = np.array(datalist[0])
#Bfield = y[np.logical_not(np.isnan(x))]
#x = np.array(datalist[2])
#Mass = x[np.logical_not(np.isnan(x))]


fig = plt.figure(figsize=(8,8))
gs = gridspec.GridSpec(4, 4)
ax_main = plt.subplot(gs[1:3, 0:2])
ax_xDist = plt.subplot(gs[0, 0:2], sharex=ax_main)
ax_yDist = plt.subplot(gs[1:3, 2], sharey=ax_main)
    
ax_main.errorbar(x, y, yerr = yerr, xerr = xerr, fmt = 'x', color = 'darkblue', alpha = 0.2)
ax_main.scatter(x, y, marker = 'x', color = 'darkblue')

ax_main.set_xlabel(r'Mass $ (\mathrm{M_{\odot}}) $' , fontsize = "13")
ax_main.set_ylabel(r'Effective Temperature ($10^3$ K) ' , fontsize = "13") # r'$T_{eff} $

ax_main.set_xlim(0.23, 1.46)
ax_main.set_ylim(3328/1000, 41365/1000)

# tick params 
labelsize = 11
ax_main.tick_params(axis='both', labelsize=labelsize)
ax_main.tick_params(axis='both', labelsize=labelsize)
ax_xDist.tick_params(axis='both', labelsize=labelsize)
ax_yDist.tick_params(axis='both', labelsize=labelsize)

# Hides ticks on histogram plots 
plt.setp(ax_xDist.get_xticklabels(), visible=False)
plt.setp(ax_yDist.get_yticklabels(), visible=False)

# Plots histogram on the x axis (top of scatter)
ax_xDist.hist(x, bins = 20, align = 'mid', color = 'darkblue', alpha = 1)
ax_xDist.set_ylabel('Number of MWDs', fontsize = "13")
#ax_xCumDist = ax_xDist.twinx()
#ax_xCumDist.hist(x,bins=100,cumulative=True,histtype='step',density=True,color='r',align='mid')
#ax_xCumDist.tick_params('y', colors='r')
#ax_xCumDist.set_ylabel('cumulative',color='r')

# Plots histogram on the y axis (side of scatter)

# LOG SCALE CASE
#hist, bins = np.histogram(y, bins=17) # bins = 21
#logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
##plt.hist(x, bins=logbins, color = 'darkblue', alpha = 1) # , density = True
#ax_main.set_yscale('log')
#ax_yDist.hist(y, bins = logbins, orientation = 'horizontal', align = 'mid', color = 'darkblue', alpha = 1)
#ax_yDist.set_xlabel('Number of MWDs', fontsize = "13")

# LINEAR SCALE CASE
ax_yDist.hist(y, bins = 20, orientation = 'horizontal', align='mid', color = 'darkblue', alpha = 1)
ax_yDist.set_xlabel('Number of MWDs', fontsize = "13")

#ax_yCumDist = ax_yDist.twiny()
#ax_yCumDist.hist(y,bins=100,cumulative=True,histtype='step',density=True,color='r',align='mid',orientation='horizontal')
#ax_yCumDist.tick_params('x', colors='r')
#ax_yCumDist.set_xlabel('cumulative',color='r')

plt.show()

#plt.savefig(f'TempMass sub.pdf', bbox_inches = 'tight')








