# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 13:22:03 2021

@author: rr3
"""

import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm
import pickle
from pathlib import Path
#%%

SourceID = 'DESI_WDJ143122.24+353420.07'
ParamID = 'Params_'+SourceID
CovID = 'Cov_'+SourceID
ErrorID = 'Errors_'+SourceID
file = SourceID+'_bin0p2'
filename = file+'.dat'
imgfile = file+'.png'

#%%

with open (ParamID, 'rb') as fp:
    params = pickle.load(fp)
#%%
with open (ErrorID, 'rb') as lp:
    errors = pickle.load(lp)

