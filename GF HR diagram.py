"""
This script works with the Gentile Fusilo (GF) catalogue, extracting quantities for our DESI analysis,
and processing the data so we can present values such as effective temperature and WD surface gravity
to accompany our determined B-field values for the DESI MWDs
"""
import pylab as pl
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lmfit as lm
import pickle
from pathlib import Path
import pandas as pd
import astropy.io.fits as fits

#%%
"""
Part 1:
    Open and read the GF FITS file
    Extract quantities from the FITS file
"""
hdulist = fits.open("{:s}final_catalogue_v6.fits".format(r"C:\Users\44743\Documents\Imperial Year 4\MSci Project\\"))

GFdata = hdulist[1].data # List of WD names in GF catalogue 
GFdata0 = hdulist[0].data

header = hdulist[1].header
header0 = hdulist[0].header

#Extract quantities from the overarching GF catalogue
T_eff = GFdata['Teff'] #WD effective temperature
eT_eff = GFdata['eTeff'] #WD effective temperature error
log_g = GFdata['log_g'] #WD surface gravity
elog_g = GFdata['elog_g'] #WD surface gravity error
WD_name = GFdata['WDJ_name'] #GF identifier for the WD
G_abs = GFdata['absG'] #the absolute Gaia magnitude for the WD
BPRP = GFdata['bp_rp'] #the difference between the BP and RP filters used when observing
parallax = GFdata['parallax'] #Gaia parallax of source
parallax_err = GFdata['parallax_error'] #error on the parralax of the source
SN = GFdata['S/N'] #signal to noise
bp_flux = GFdata['phot_bp_mean_flux']
rp_flux = GFdata['phot_rp_mean_flux']
bp_flux_err = GFdata['phot_bp_mean_flux_error']
rp_flux_err = GFdata['phot_rp_mean_flux_error']
Mass = GFdata['mass']
Err_mass = GFdata['emass']

##%%
"""
Part 2:
    Compile all the quantities which will be used to characterise the DESI MWDs into one Pandas dataframe
    Process the data to remove NaN values (wherever there is a NaN for T_eff there is a NaN for log_g)
    Tweak the formatting of the data to make it more suitable for our analysis
"""
names = ["G_abs", "BP-RP", "Parallax", 'Parallax_error', "T_eff", 'eT_eff', "log_g", 'elog_g', "SN", 'Mass', 'Err_mass']
info_df = pd.DataFrame(np.array([G_abs, BPRP, parallax, parallax_err, T_eff, eT_eff, log_g, elog_g, SN, Mass, Err_mass]).transpose(), WD_name, columns=names)
dataset = info_df.reset_index()
dataset = dataset.rename(columns={"index": "WDJ Name"}) #relabel the relic column header from the original dataframe

#%% ESSENTIAL
""" Plotting script for MWD HRD """
# Selects WDs with reciprocal parallaxes/1000 < 100 --> says WDs must be x distance from us (or close to us?)

dist_pc = np.reciprocal(parallax/1000)
idx = np.where(dist_pc < 100)[0]

dist_pc_constrained = [dist_pc[i] for i in idx]
BPRP_sel = [BPRP[i] for i in idx]
G_sel = [G_abs[i] for i in idx]

plt.figure()
plt.gca().invert_yaxis()
plt.xlabel('BP-RP')
plt.ylabel('G_abs')

plt.plot(BPRP_sel,G_sel,'o', markersize=0.25)

#%% DO NOT RUN
nans = info_df[info_df['T_eff'].isna()]

dataset = info_df.dropna().reset_index()
##%%
dataset = dataset.rename(columns={"index": "WDJ Name"}) #relabel the relic column header from the original dataframe
#dataset.columns = dataset.columns.str.rstrip() #remove all white spaces in the strings
#%%
info2 = ['WDJ014202.20-003332.01 ', 
'WDJ075804.43+205957.86 '
]
info2 = filename_list
df2 = dataset[dataset['WDJ Name'].isin(info2)]
print(df2)
#%% ESSENTIAL
""" Plotting script for MWD HRD """
# Selects WDs with reciprocal parallaxes/1000 < 100 --> says WDs must be x distance from us (or close to us?)

dist_pc = np.reciprocal(parallax/1000)
idx = np.where(dist_pc < 100)[0]

dist_pc_constrained = [dist_pc[i] for i in idx]
BPRP_sel = [BPRP[i] for i in idx]
G_sel = [G_abs[i] for i in idx]

plt.figure()
plt.gca().invert_yaxis()
plt.xlabel('BP-RP')
plt.ylabel('G_abs')

plt.plot(BPRP_sel,G_sel,'o', markersize=0.25)
#%%
filename = "DESI_WDJ172329.14+540755.79_bin0p2.dat"
data = np.genfromtxt(f'{filename}', delimiter=' ')
wavelength = data[:,0]
flux = data[:,1]
error = data[:,2]
plt.figure("Whole spectrum")
plt.plot(wavelength,flux, label = f"{filename}")
plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
plt.ylabel("Flux", size = "15")
plt.grid()
plt.legend()
plt.show()
#%%
info_trial = ['WDJ172329.14+540755.79 ']
df_trial = dataset[dataset['WDJ Name'].isin(info_trial)]
print(df_trial)
#%%
T_val = df_trial.iloc[0]['T_eff']
log_g_val = df_trial.iloc[0]['log_g']
G_abs_val = df_trial.iloc[0]['G_abs']
parallax_val = df_trial.iloc[0]['Parallax']
SN_val = df_trial.iloc[0]['SN']
BPRP_val = df_trial.iloc[0]['BP-RP']
name_val = info_trial[0]
print(T_val)
#%%
textstr = "\n".join((
    'Name: %s' % (name_val),
    'G = %2f' % (G_abs_val, ),
    'S/N = %2f' % (SN_val, ),
    'Parallax = %2f' % (parallax_val, ),
    'T_eff = %2f' % (T_val, ),
    'log_g = %2f' % (log_g_val, )))
#%%
f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
a0.set_xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
a0.set_ylabel("Flux", size = "15")
a0.plot(wavelength, flux, label = f"{filename}")
a1.invert_yaxis()
a1.set_xlabel('BP-RP')
a1.set_ylabel('G_abs')
a1.plot(BPRP_sel,G_sel,'o', markersize=0.25)
a1.plot(BPRP_val,G_abs_val,'o',color='red',markersize=5)
props = dict(boxstyle='square', alpha=0.5)
a1.text(0.05, 0.05, textstr, transform=a1.transAxes, fontsize=9,
        verticalalignment='bottom', bbox=props)
#f.tight_layout()


#%%
""" Reading in batches of files """

column_names = ["Filename", "Class", "B alpha", "B error alpha", "lambda0", "lambda0 error", "begin", "finish", "start1", "end1", "start2", "end2", "start3", "end3", "Notes"]
cut_regions = ["begin", "finish", "start1", "end1", "start2", "end2", "start3", "end3"]
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

cut_regions = [begin_list, finish_list, start1_list, end1_list, start2_list, end2_list, start3_list, end3_list]   

for i in range(len(datasets)):
    if datasets.Class[i] == 1 or datasets.Class[i] == 2 or datasets.Class[i] == 3 or datasets.Class[i] == 4:
        filename_list.append((datasets.Filename[i] + ' ')) # NOTE DESI AND GAIA HAVE DIFFERENT NAMES
                                                            # therefore cannot call as same variable!!
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

#file_name, lin, quad, undec, da, comm = np.loadtxt(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\First categorisation 113 systems.csv',skiprows = 1, delimiter = ',', unpack = True)
""" re engineer this - has to be faster way to write this esp for 200 + systems """
datafolder = [filename_list, cut_regions[0], cut_regions[1], cut_regions[2], cut_regions[3], cut_regions[4], \
                cut_regions[5], cut_regions[6], cut_regions[7]]
##%%
number_sections = 6 # for class 1 this is only 6 (was 11 for generality but was empty lists past this point)

subsections = [[[] for i in range(0, len(datafolder))] for i in range(number_sections)]

intervals = [0]
[intervals.append(intervals[i]+20) for i in range(0,number_sections)]

for i in range(0,number_sections):
    for j in range(0,len(datafolder)):
        subsections[i][j] = datafolder[j][intervals[i]:intervals[i+1]]


#%% 
Gaia_filename_list = [filename_list[i].replace('DESI_','') for i in range(len(filename_list))]
Gaia_filename_list = [Gaia_filename_list[i].replace('_bin0p2.dat','') for i in range(len(filename_list))]

sub_Gaia_filename = [[] for i in range(0, number_sections)]
subsections_Gaia_file = [Gaia_filename_list[intervals[i]:intervals[i+1]] for i in range(0, number_sections)]

#tot_list = [Bvaluelist, Bvalueerr, lambdalist, lambdaerrlist, beginlist, finishlist, start1list, end1list, start2list, end2list, start3list, end3list]   
#
#aexceldata = np.zeros((len(filenamelist),len(tot_list)))
#for i in range(len(filenamelist)):
#    for j in range(0,len(tot_list)):
#        aexceldata[i][j] = tot_list[j][i]

#%%
# SPENT AGES getting combdata but combdata is the SAME AS datafolder !!!
combdata = [[] for i in range(0,9)] # range(0,x) x has to match length of Gaiavalues
#completesection = subsections[0] + subsections[1] + subsections[2] + subsections[3] + subsections[4] + \
#                    subsections[5] + subsections[6] + subsections[7] + subsections[8] + subsections[1]

for j in range(0, len(subsections[0])):
    A = []
    for i in range(0, len(subsections)):
        A.append(subsections[i][j])
    conc_list = np.concatenate(A, axis = 0).tolist()
    combdata[j] = conc_list

#A = []
#for i in range(0, len(subsections)):
#    A.append(subsections[i][0])
#conc_list = np.concatenate(A, axis = 0).tolist()
#%%
#aexceldata = np.zeros((len(Exceldata),len(Exceldata[0])))
#for i in range(len(Exceldata)):
#    for j in range(0,len(Exceldata[0])):
#        aexceldata[i][j] = Exceldata[0][j][i]
#
aexceldata = np.zeros((len(Exceldata[0]),len(Exceldata)))

#aexceldata[:,i] = [Exceldata[i] for i in range(1, len(Exceldata))]
#%%
aexceldata = np.zeros((len(Exceldata[0]),len(Exceldata)))

for i in range(1, len(Exceldata)):
    aexceldata[:,i] = Exceldata[i]
    
#%%

Exceldata = [[] for i in range(0,12)] # range(0,x) x has to match length of Gaiavalues


alpha = 5
subfolder_ = subsections[alpha] #[1]# combdata
for j in range(len(subfolder_[0][0:])):
    filename = subfolder_[0][j]

    #filename = "DESI_WDJ153349.03+005916.21_bin0p2.dat"
    data = np.genfromtxt(f'{filename}', delimiter=' ')
    
    wavelength = data[:,0]
    flux = data[:,1]
    error = data[:,2]
    
    #plt.figure("Whole spectrum")
#    plt.figure()
#    plt.plot(wavelength,flux, label = f"{filename}")
#    plt.xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
#    plt.ylabel("Flux", size = "15")
#    plt.grid()
#    plt.legend()
#    plt.show()
            
    #Gaia_data = dataset[dataset['WDJ Name'].isin([subsections_Gaia_file[alpha][j]])] # for batches of 20
    Gaia_data = dataset[dataset['WDJ Name'].isin([Gaia_filename_list[j]])] # for all 216 files
    ##%%
#    Gaia_data = dataset[dataset['WDJ Name'].isin([Gaia_filename_list[0]])]

#    info_trial = ['WDJ172329.14+540755.79 ']
#    df_trial = dataset[dataset['WDJ Name'].isin(info_trial)]
    print(Gaia_data)
    ##%%
    T_val = Gaia_data.iloc[0]['T_eff']
    eT_val = Gaia_data.iloc[0]['eT_eff']
    log_g_val = Gaia_data.iloc[0]['log_g']
    elog_g_val = Gaia_data.iloc[0]['elog_g']
    G_abs_val = Gaia_data.iloc[0]['G_abs']
    parallax_val = Gaia_data.iloc[0]['Parallax']
    eparallax_val = Gaia_data.iloc[0]['Parallax_error']
    SN_val = Gaia_data.iloc[0]['SN']
    BPRP_val = Gaia_data.iloc[0]['BP-RP']
    Mass_val = Gaia_data.iloc[0]['Mass']
    Mass_err = Gaia_data.iloc[0]['Err_mass']
    name_val = Gaia_filename_list[j] #Gaia_data.iloc[0][#Gaia_filename_list[j]
    print(T_val)
    
    Gaiavalues = [name_val, Mass_val, Mass_err, BPRP_val, T_val, eT_val, log_g_val, elog_g_val, \
                  G_abs_val, parallax_val, eparallax_val, SN_val]
    
    for i in range(0,len(Gaiavalues)):
        Exceldata[i].append(Gaiavalues[i])
    
    ##%%
#    textstr = "\n".join((
#        'Name: %s' % (name_val),
#        'G_abs = %2f' % (G_abs_val, ),
#        #'$G_{abs}$ = %2f' % (G_abs_val, ),
#        'S/N = %2f' % (SN_val, ),
#        'Parallax = %2f pm %.3g' % (parallax_val, eparallax_val),
#        'T_eff = %2f pm %.3g' % (T_val, eT_val),
#        'log_g = %2f pm %.3g' % (log_g_val, elog_g_val),
#        'Mass = %2f pm %.3g' % (Mass_val, Mass_err)))
#    
#    textstr = "\n".join((
#        'Name: %s' % (name_val),
#        '$G_{abs}$ = %2f' % (G_abs_val, ),
#        'S/N = %2f' % (SN_val, ),
#        'Parallax = %2f pm %.3g' % (parallax_val, eparallax_val),
#        '$T_{eff}$ = %2f pm %.3g' % (T_val, eT_val),
#        'log(g) = %2f pm %.3g' % (log_g_val, elog_g_val),
#        'Mass = %2f pm %.3g' % (Mass_val, Mass_err)))
#    ##%%
##    print("%.8f pm %.3g" % (c[0], c[1]))
#    
#    f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
#    a0.set_xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
#    a0.set_ylabel("Flux", size = "15")
#    a0.plot(wavelength, flux, label = f"{filename}")
#    a1.invert_yaxis()
#    a1.set_xlabel(r'$G_{BP} - G_{RP}$' , size = "15")
#    a1.set_ylabel(r'$G_{abs}$, $[mag]$' , size = "15")
#    a1.plot(BPRP_sel,G_sel,'o', markersize=0.25) # Plots the whole HR diagram
#    a1.plot(BPRP_val,G_abs_val,'o',color='red',markersize=5) # Plots the BPRP and abs G for this system 
#    props = dict(boxstyle='square', alpha=0.5)
#    a1.text(0.05, 0.05, textstr, transform=a1.transAxes, fontsize=9,
#            verticalalignment='bottom', bbox=props)
#    #f.tight_layout()




















































