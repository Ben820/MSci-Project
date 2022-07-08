# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 13:31:08 2021

Directories: 
    Users\44743\Documents\Imperial Year 4\MSci Project\DESI_DxH\DESI_DxH
    
Method: 
    Part 1 - Loads data
    Part 2 - Performs cuts on the data to isolate the H-alpha region 
    Part 3 - Removes the absorption region to isolate the continuum
    Part 4 - Fits a polynomial to the continuum and normalises the spectrum
    
"""
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
column_names = ["Filename", 'Identifier', "Class", 'Cut_class', "B alpha",\
                "B error alpha", "lambda0", "lambda0 error", 'sigma1', 'sigma2', \
                "begin", "finish", "start1", "end1", "start2", "end2", "start3", "end3", \
                'Notes', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6']
cut_regions = ["begin", "finish", "start1", "end1", "start2", "end2", "start3", "end3"]
#datasets = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Detailed class catalogue 2_mod.csv', skiprows = 1, names = column_names)
#datasets = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Detailed class catalogue 2.csv', skiprows = 1, names = column_names)
datasets2 = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Detailed class catalogue 3.csv', skiprows = 1, names = column_names)
datasets = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Comprehensive catalogue 1 csv.csv', skiprows = 1, names = column_names)
datasets = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Comprehensive catalogue 2 csv.csv', skiprows = 1, names = column_names)
datasets = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Comprehensive catalogue 3 csv.csv', skiprows = 1, names = column_names)
datasets = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Comprehensive catalogue 5 csv.csv', skiprows = 1, names = column_names)

datasets = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\UROP Year 4\UROP catalogue.csv', skiprows = 1, names = column_names)
datasets = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\UROP Year 4\New DAH MWDs.csv', skiprows = 1, names = column_names)

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
cut_class = []

cut_regions = [begin_list, finish_list, start1_list, end1_list, start2_list, end2_list, start3_list, end3_list,
               cut_class]   

for i in range(len(datasets)):
    if datasets.Class[i] == 1:#'Linear':# or datasets.Class[i] == 2 or datasets.Class[i] == 3 or datasets.Class[i] == 4:
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
        cut_class.append(datasets.Cut_class[i])

#file_name, lin, quad, undec, da, comm = np.loadtxt(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\First categorisation 113 systems.csv',skiprows = 1, delimiter = ',', unpack = True)
""" re engineer this - has to be faster way to write this esp for 200 + systems """
datafolder = [filename_list, cut_regions[0], cut_regions[1], cut_regions[2], cut_regions[3], cut_regions[4], \
                cut_regions[5], cut_regions[6], cut_regions[7], cut_regions[8]]
##%%
# datafolder has system name and 8 different cuts 
# one_ etc. only contains 8 different cuts!!
number_sections = 15

subsections = [[[] for i in range(0, len(datafolder))] for i in range(number_sections)]
#one_  = [[], [], [], [],[], [], [], [], []]

intervals = [0]
[intervals.append(intervals[i]+20) for i in range(0,number_sections)]

for i in range(0,number_sections):
    for j in range(0,len(datafolder)):
        subsections[i][j] = datafolder[j][intervals[i]:intervals[i+1]]
##%% RUN AFTER BIG CELL (just dont want down bottom - cumbersome to edit )
#    
#tot_list = [Bvaluelist, Bvalueerr, lambdalist, lambdaerrlist, beginlist, finishlist, \
#            start1list, end1list, start2list, end2list, start3list, end3list, ]   
#
#aexceldata = np.zeros((len(filenamelist),len(tot_list)))
#for i in range(len(filenamelist)):
#    for j in range(0,len(tot_list)):
#        aexceldata[i][j] = tot_list[j][i]
##%% SPECIFICALLY JUST FOR ROTATING SYSTEM!!!!
#datasetsRotating = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Rotating system analysis csv.csv', skiprows = 1, names = column_names)
#filename_list = []
#for i in range(len(datasetsRotating)):
#    filename_list.append(datasets.Filename[i])
#
#%%
aexceldata = np.zeros((len(filenamelist),len(datalist)))

for i in range(0, len(datalist)):
    aexceldata[:,i] = datalist[i]
#%%
""" Part 1: Load data

Notes: data - MWD spectrum
       wavelength - x-values 
       flux - y-values
"""
filenamelist = []
datalist = [[] for i in range(0, 14)] # lengt(# data params you want (4 ))

subfolder_ = subsections[0] #one_ #[1]# combdata # if want to read in batches of 20 files
subfolder_ = datafolder # if want to read in 216 files at once (all files)
for j in range(119, 124):#len(subfolder_[0][0:])): # 0, len(subfolder_[0][0:]) # 12,20 # 13,14 for Na I DZ
    # range(7,8) subsection[1] for linear Zeeman split viva
    filename = subfolder_[0][j]
#    filename = filename_list[j] # rotating system
    #filename = 'DESI_WDJ160604.06+550924.93_bin0p2.dat'

    #load data and sort into appropriate variables
    #filename = "DESI_WDJ235958.62+301622.51_bin0p2.dat"
    data = np.genfromtxt(f'{filename}', delimiter=' ')
    
    wavelength = data[:,0]
    flux = data[:,1]
    error = data[:,2]
    
#    plt.figure()
#    plt.errorbar(wavelength,flux, yerr = error ,label = f"{filename}", fmt ='')
#    plt.xlabel("Wavelength $\lambda$ $[\mathrm{\AA}]$" , size = "15")
#    plt.ylabel("Flux", size = "15")
#    plt.axvline(6562.8, ls = '--', color = 'green', lw = 1)
#    plt.axvline(4861.35, ls = '--', color = 'green', lw = 1)
#    plt.axvline(4340, ls = '--', color = 'green', lw = 1)
#
##    plt.xlim(3300, 9000)
##    plt.ylim(-20,190)
#    plt.grid()
#    plt.legend()
#    plt.show()
    ##%%
    """ Part 2: Performs cuts on the data to isolate the H-alpha region
    Notes: start/start_Ha - beginning of cut
           last/last_Ha - end of cut 
           masked_Ha_flux - y-values included within the cut
           masked_Ha_reg - x-values included within the cut
    begin/ finish define the whole region including the triplet feature 
    startx/endx define the specific region to be cut out (the absorption feature) """
    
#    filenamelist = []
#    datalist = [[] for i in range(0, 14)] # lengt(# data params you want (4 ))
    
#    p = [[] for k in range(8)]
    p = [[6250, 5800, 5000, 4800, 6300, subfolder_[1][j]], [6900, 7200, 8000, 8200, 6850, subfolder_[2][j]], 
         [6400, 6100, 5500, 5400, 6480, subfolder_[3][j]],
      [6750, 6900, 7500, 7600, 6640, subfolder_[4][j]]]
    
    h = 7
    if subfolder_[9][j] == 1:
        h = 0
    if subfolder_[9][j] == 2:
        h = 1
    if subfolder_[9][j] == 3:
        h = 2
    if subfolder_[9][j] == 4:
        h = 3
    if subfolder_[9][j] == 5:
        h = 4
    if subfolder_[9][j] == 6:
        h = 5
        
    begin = p[0][h]
    finish = p[1][h]
    start1 = p[2][h]
    end1 = p[3][h]
    start2 = 9032
    end2 = 9034
    start3 = 9040
    end3 = 9049

    lambda_guess = 6562.8# 6562.8 #6711#4227 # This is the central wavelength of the transition! (Be it H alpha or Ni, Mg etc)
    # Na I 5893; Mg I 5171; 8545; Fe I: 7445, 6712
    if begin < 4861:
        lambda_guess = 4861.35

#    begin = 5000
#    finish = 8000
#    start1 = 6000
#    end1 = 7000
#    start2 = 9032
#    end2 = 9034
#    start3 = 9040
#    end3 = 9049

#    lambda_guess = 7065
#    begin = subfolder_[1][j]
#    finish = subfolder_[2][j]
#    start1 = subfolder_[3][j]
#    end1 = subfolder_[4][j]
#    start2 = subfolder_[5][j]
#    end2 = subfolder_[6][j]
#    start3 = subfolder_[7][j]
#    end3 = subfolder_[8][j]
    
#    beginlist.append(begin)
#    finishlist.append(finish)
#    start1list.append(start1)
#    end1list.append(end1)
#    start2list.append(start2)
#    end2list.append(end2)
#    start3list.append(start3)
#    end3list.append(end3)
    
    start_Ha = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-begin)))[0])
    end_Ha = int(np.where(wavelength == min(wavelength, key=lambda x:abs(x-finish)))[0])
    masked_Ha_flux = flux[start_Ha:end_Ha]
    masked_Ha_reg = wavelength[start_Ha:end_Ha]
    masked_Ha_err = error[start_Ha:end_Ha]
    
    dip_start1 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-start1)))[0])
    dip_end1 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-end1)))[0])
    dip_start2 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-start2)))[0])
    dip_end2 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-end2)))[0])
    dip_start3 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-start3)))[0])
    dip_end3 = int(np.where(masked_Ha_reg == min(masked_Ha_reg, key=lambda x:abs(x-end3)))[0])
    
    masked_flux = list(masked_Ha_flux)
    masked_wavelength = list(masked_Ha_reg)
    masked_err = list(masked_Ha_err)
    
    ##%%
    flxframe = pd.DataFrame(masked_flux)
    flxframe.loc[dip_start1:dip_end1] = np.nan
    flxframe.loc[dip_start2:dip_end2] = np.nan
    flxframe.loc[dip_start3:dip_end3] = np.nan
    flxframe = (flxframe.dropna()).to_numpy()
    
    wavframe = pd.DataFrame(masked_wavelength)
    wavframe.loc[dip_start1:dip_end1] = np.nan
    wavframe.loc[dip_start2:dip_end2] = np.nan
    wavframe.loc[dip_start3:dip_end3] = np.nan
    wavframe = (wavframe.dropna()).to_numpy()
    
    """ Part 4: Fits a polynomial to the continuum and normalises the spectrum
    Notes: Poly_3o is a third-order polynomial function which fits the continuum using a least-squares method
          
    """
    flxframe = np.array(flxframe)
    wavframe = np.array(wavframe)
    flux_list = flxframe[:,0]
    wav_list = wavframe[:,0]
    ###%%
    #plt.figure()
    #plt.plot(wav_list,flux_list,'x')
    #plt.show()
    ##%%
    def Poly_3o(x, a, b, c, d):
        y = a*x**3 + b*x**2 + c*x + d
        return y
    
    p0 = np.array([1, 1, 1, 1]) #fit parameters
    p, cov = opt.curve_fit(Poly_3o, wav_list,flux_list, p0) # do not change
    Continuum = Poly_3o(masked_Ha_reg, p[0], p[1], p[2], p[3]) # use masked_Ha_reg since more data points in array
    
    #for c in zip(p, np.sqrt(np.diag(cov))):# zips root of diag of cov matrix with related value in curve fit
    #    print('%.15f pm %.4g' % (c[0], c[1]))# prints value and uncertainty, f is decimal places and G is sig figs
    
    norm_spectra = masked_Ha_flux/Continuum
    
    #plt.figure()
    #plt.grid()
    ##plt.yticks([0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30])
    ##plt.plot(masked_Ha_reg, Poly_3o(masked_Ha_reg, p[0], p[1], p[2], p[3])/Continuum, \
    ##         zorder=4,color = 'red', label = "Poly")
    #plt.plot(wav_list,flux_list,'x')
    #plt.plot(masked_Ha_reg, Continuum)
    #plt.show()
    ##%%
    #plt.figure()
    #plt.plot(masked_Ha_reg, norm_spectra, label = f"{filename}") #plot the normalised spectrum
    #plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
    #plt.ylabel("Normalised Flux", size = "15")
    #plt.legend()
    #plt.show()
    ##%%
    """ Triplet Lorentzian Profile 
    
    Notes: 
    """
    xp_triplet = []
    yp_triplet = []
    err_triplet = []
    for a in range(len(masked_Ha_reg)):
        if masked_Ha_reg[a] > (begin) and masked_Ha_reg[a] < (finish):
            xp_triplet.append(masked_Ha_reg[a])
            yp_triplet.append(norm_spectra[a])
            err_triplet.append(masked_err[a])
    xp_triplet = np.array(xp_triplet)
    yp_triplet = np.array(yp_triplet)
    err_triplet = np.array(err_triplet)
    
    """ Ha Data; Lorentzian fit _3Lorentzian; Run after clipping data xp_triplet yp_triplet """
        
    ### SEE WIKI: 'In physics, a three-parameter Lorentzian function is often used'
    # wid = gamma and FWHM = 2*gamma
    
    def _3Lorentzian(x, lambda0, B, amp1, wid1, amp2, wid2):
    #    lambda0 = 6562.8
        A = np.square(lambda0/6564)
        delt_lam = 20.2*A*B
        lambda_minus = lambda0 - delt_lam
        lambda_plus = lambda0 + delt_lam     
        return -((amp1*wid1**2/((x-lambda_minus)**2+wid1**2)) +\
                (amp2*wid2**2/((x-lambda0)**2+wid2**2)) +\
                    (amp1*wid1**2/((x-lambda_plus)**2+wid1**2)))+1
    
    """ SECOND FUNCTION IS FOR DZ!!!!! 
    Make sure you change lambda_guess!! """
#    def _3Lorentzian(x, lambda0, B, amp1, wid1, amp2, wid2):
#    #    lambda0 = 6562.8
#        A = np.square(lambda0)
#        delt_lam = (4.67e-7)*A*B
#        lambda_minus = lambda0 - delt_lam
#        lambda_plus = lambda0 + delt_lam     
#        return -((amp1*wid1**2/((x-lambda_minus)**2+wid1**2)) +\
#                (amp2*wid2**2/((x-lambda0)**2+wid2**2)) +\
#                    (amp1*wid1**2/((x-lambda_plus)**2+wid1**2)))+1

#    lambda_guess = 6562.8# 6562.8 #6711#4227 # This is the central wavelength of the transition! (Be it H alpha or Ni, Mg etc)
#    # Na I 5893; Mg I 5171; 8545; Fe I: 7445, 6712
    Data = {}
    Res_list = []
    B_list = []
    lambda0_list = []
        
    rangeBval = np.arange(0,50.5,0.5)
    for i in rangeBval:
        p0 = np.array([lambda_guess, i, 0.8, 10, 0.8, 10])
        
        try:
            popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, p0, sigma = err_triplet)
        except:
            pass
        
    #    popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, p0)
    #    for c in zip(popt_VoigtNew, np.sqrt(np.diag(cov_VoigtNew))):
    #        print("%.8f pm %.3g" % (c[0], c[1]))
        Residuals= _3Lorentzian(xp_triplet, *popt_3lorentz)-yp_triplet
        Res_list.append(sum(np.square(Residuals)))
        B_list.append(popt_3lorentz[1])
        lambda0_list.append(popt_3lorentz[0])
        
    Data[1] = [Res_list, B_list, lambda0_list]
    
    index = np.where(Res_list == np.amin(Res_list))[0][0]
    
    # Method 
    # The determined B value is BACK INTO curvefit to re-estimate the final parameters of the fit 
    popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, \
                                                p0=[lambda_guess, B_list[index], 0.8, 10, 0.8, 10], \
                                                sigma = err_triplet)#, \
                                                #bounds = ((-np.inf, 0, -np.inf, -np.inf, -np.inf, -np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
    
    ## ALTERNATIVE METHOD: Just take first curvefit estimate - prevents over fitting and finding a 
    # wrong local minima in least squares
    #p0 = np.array([6562.8, rangeBval[2], 0.8, 10, 0.8, 10])
    #popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, p0, sigma = err_triplet)
    
    ## ALTERNATIVE METHOD 2: (Supercedes Method):
    # Puts the list of returned B values from first curvefit into another curvefit, scans through and 
    # looks to see which has the lowest residuals (on the resulting curvefit); 
    # selects the corresponding new index for that fit, and then uses that index to identify the correct
    # B value from B_list (the first! returned set of B values) (NOTE THIS IS THE SECOND TIME IT IS CURVEFITTED)
    Res_list2 = []
    B_list2 = []
    lambda0_list2 = []
    
    for k in range(len(B_list)):
        p0 = np.array([lambda_guess, B_list[k], 0.8, 10, 0.8, 10])
        
        try:
            popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, p0, sigma = err_triplet)
        except:
            pass
        Residuals2= _3Lorentzian(xp_triplet, *popt_3lorentz)-yp_triplet
        Res_list2.append(sum(np.square(Residuals2)))
        B_list2.append(popt_3lorentz[1])
        lambda0_list2.append(popt_3lorentz[0])
    
    index2 = np.where(Res_list2 == np.amin(Res_list2))[0][0]
    
    # Take the minimum residual of both residual lists combined - thus doesnt double iterate 
    # if not requied/ beneficial 
    A = Res_list+Res_list2
    index3 = np.where(A == np.amin(A))[0][0]
    
    # If index3 > len(B_list) then the double iteration is preferable, and so we use index2 
    if index3 >= len(B_list)-1:
        index3 = index2
    ## If any of the initial guesses for B are the best, this will ensure these are selected 
    #if Res_list[index3] >= Res_list[index3]:
    #    B_list[index3] = rangeBval[index3]
    #
    ## If indices do not match up exactly with B values --> consequence of using only B_list below 
    ## Resolves the case if index3 does not correspond to the same B value in B_list and B_list2
    #if Res_list[index3] >= Res_list2[index3]:
    #    B_list[index3] = B_list2[index3]
    #
    #
    ###%%
    ## The determined B value is BACK INTO curvefit to re-estimate the final parameters of the fit 
    ## NOTE THIS IS STILL ONLY THE SECOND CURVEFIT - NOTE USE OF B_list NOT! B_list2
    #popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, \
    #                                            p0=[6562.8, B_list[index3], 0.8, 10, 0.8, 10], \
    #                                            sigma = err_triplet)#, \
    #                                            #bounds = ((-np.inf, 0, -np.inf, -np.inf, -np.inf, -np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
    
    # In the (2 instance so far) case that the lowest residuals outright does not conform to the best fit 
    # This is different to the reason for B_list2 and the secod curvefit; the second curvefit is used when
    # the first curvefit does not find the correct B, so finer precision is required
    # this case is when the correct B value IS identified, but is not selected because another (wrong)
    # B value has lower residuals 
    
    print("B = ", popt_3lorentz[1], 'pm', np.sqrt(np.diag(cov_3lorentz))[1])
    
    #if np.around(np.array([popt_3lorentz[1]]), 2) != Counter(np.around(np.array(B_list),4)).most_common(1)[0][0]:
    #    B_list[index3] = Counter(np.around(np.array(B_list),4)).most_common(1)[0][0]
    
    popt_3lorentz, cov_3lorentz = opt.curve_fit(_3Lorentzian, xp_triplet, yp_triplet, \
                                                p0=[lambda_guess, B_list[index3], 0.8, 10, 0.8, 10], \
                                                sigma = err_triplet)#, \
                                                #bounds = ((-np.inf, 0, -np.inf, -np.inf, -np.inf, -np.inf),(np.inf, np.inf, np.inf, np.inf, np.inf, np.inf)))
    
    for c in zip(popt_3lorentz, np.sqrt(np.diag(cov_3lorentz))):
        print("%.8f pm %.3g" % (c[0], c[1]))
    
    print("B = ", popt_3lorentz[1], 'pm', np.sqrt(np.diag(cov_3lorentz))[1])
    
    Residuals= _3Lorentzian(xp_triplet, *popt_3lorentz)-yp_triplet
    print("Lorentzian Residual sum of squares = ", sum(np.square(Residuals)))
        
    
    filenamelist.append(filename)
    
    datavalues = [popt_3lorentz[1], np.sqrt(np.diag(cov_3lorentz))[1], popt_3lorentz[0], \
                  np.sqrt(np.diag(cov_3lorentz))[0], popt_3lorentz[3], popt_3lorentz[5], \
                  begin, finish, start1, end1, start2, end2, start3, end3]   
    
    for i in range(0,len(datavalues)):
        datalist[i].append(datavalues[i])

   #%% PLotting Cell
    fig, axs = plt.subplots(2, gridspec_kw={'height_ratios': [2.5, 1]})
    fig.suptitle(f"Lorentzian fit {filename} \n B = {popt_3lorentz[1]} +/- {np.sqrt(np.diag(cov_3lorentz))[1]}", \
                 fontsize = "13")
    axs[0].plot(xp_triplet, yp_triplet, color = 'darkblue')#, label = "WD data") #, color = 'midnightblue')#
    axs[0].plot(xp_triplet, _3Lorentzian(xp_triplet, *popt_3lorentz), linewidth=2, \
                                   color = "tomato", label = "Triplet Lorentzian fit")
    axs[0].tick_params(axis='both', labelsize=11)
    axs[1].tick_params(axis='both', labelsize=11)

#    axs[0].text(r"B = {popt_3lorentz}")
    for ax in axs.flat:
        ax.set_xlabel(r'Wavelength $(\mathrm{\AA})$', fontsize = "13")#, labelpad = 5)
        ax.set_ylabel('Normalised flux', fontsize = "13", labelpad = 15)   
#        ax.tick_params(axis = 'both', which = 'both', fontsize = 14)#, pad = 2)
#        ax.set_xticklabels(Fontsize = 14)
#        ax.tick_params(axis='both', labelsize=14)


    for ax in axs.flat:
        ax.label_outer()
    axs[1].set_ylabel(r'$\mathrm{Residuals}$', fontsize = "13", labelpad = 5)
    axs[0].grid()
    #""" Voigt Residuals Plot """
    axs[1].plot(xp_triplet, Residuals, linewidth=2, label = "Lorentzian fit Residuals", color = "darkblue")#
    axs[1].plot(xp_triplet, xp_triplet*0/xp_triplet, linewidth = 2, color = 'tomato')#"orange")
    axs[0].legend(prop={'size': 11})
    #axs[1].legend()
    axs[1].grid()
#    plt.tick_params(axis='both', which='major', labelsize=14)
    
    #plt.savefig("Gaussian Residuals", dpi = 1000)
#    plt.savefig(f'{filename} fit.png', dpi = 1000, bbox_inches = 'tight')
#    plt.savefig('Lineartriplet2.pdf', bbox_inches = 'tight')#,pad_inches = 0)
#    plt.xticks(fontsize=14)
#    plt.yticks(fontsize=14)

    plt.show()

    #%% # Unfit spectra
    plt.figure()
    plt.plot(xp_triplet, yp_triplet, color = 'darkblue')
    plt.tick_params(axis='both', labelsize=14)
    plt.xlabel(r'Wavelength $\lambda$ $(\mathrm{\AA})$', fontsize = "15")
    plt.ylabel('Normalised flux', fontsize = "15", labelpad = 15)   
    plt.grid()
    #plt.savefig("Gaussian Residuals", dpi = 1000)
#    plt.savefig(f'{filename} unfit darkblue.png', dpi = 1000, bbox_inches = 'tight')
#    plt.savefig('VivaFig1.pdf', bbox_inches = 'tight')#,pad_inches = 0)
#    plt.xticks(fontsize=14)
#    plt.yticks(fontsize=14)

    plt.show()

#    #top=0.98,
#bottom=0.135,
#left=0.16,
#right=0.985,
#hspace=0.2,
#wspace=0.2
    #%%
    
    
    
    
    
    
    
    
    
    
    
    
    


























