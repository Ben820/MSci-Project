# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 19:03:39 2022

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
"""
Retrieve Gaia data
"""
gaia_data = pd.read_csv('Gaia catalogue.csv')
#%%
"""
Cell that plots a histogram of MWD mass vs MWD population
"""
masses = gaia_data[gaia_data.columns[1]]
masses = masses.dropna()
masses_array = masses.to_numpy()
plt.figure()
plt.xlabel('MWD Mass' + ' (' '$M_{\odot}$' + ')')
plt.ylabel('Number of MWDs')
plt.hist(masses_array, bins=21, color='k')
