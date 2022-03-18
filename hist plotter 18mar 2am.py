# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 01:58:36 2022

@author: rr3
"""

import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as spi
from scipy.misc import derivative
import lmfit as lm
import pylab
import glob
import os
import pandas as pd
import re
import astropy.io.fits as fits

#%%
fields_data = pd.read_csv('data for trial hist 18mar 2am v2.csv')
#%%
fields_array = abs(fields_data.to_numpy())
#%%
plt.figure()
plt.xlabel('MWD Field Strength' + ' (MG)')
plt.ylabel('Number of MWDs')
#plt.hist(fields_array,bins=100)

# histogram on linear scale
plt.subplot(211)
hist, bins, _ = plt.hist(fields_array,bins=30)

# histogram on log scale. 
# Use non-equal bin sizes, such that they look equal on log scale.
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
plt.subplot(212)
plt.hist(fields_array, bins=logbins)
plt.xscale('log')
plt.show()