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
from scipy.stats import gaussian_kde
from scipy import stats
from matplotlib import colors


#%%
"""
Part 1:
    Open and read the GF FITS file
    Extract quantities from the FITS file
"""
# final_catalogue_v6 is Gaia data - thousands of WD candidates 
hdulist = fits.open("{:s}final_catalogue_v6.fits".format(r"C:\Users\44743\Documents\Imperial Year 4\MSci Project\Datafiles\\"))

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

#%% ESSENTIAL - plots HR diagram (and gets essential data for it )
""" Plotting script for MWD HRD """
# Selects WDs with reciprocal parallaxes/1000 < 100 --> says WDs must be x distance from us (or close to us?)

dist_pc = np.reciprocal(parallax/1000)
idx = np.where(dist_pc < 100)[0]

dist_pc_constrained = [dist_pc[i] for i in idx]
BPRP_sel = [BPRP[i] for i in idx]
G_sel = [G_abs[i] for i in idx]

#plt.figure()
plt.gca().invert_yaxis()
plt.xlabel('BP-RP')
plt.ylabel('G_abs')

plt.plot(BPRP_sel,G_sel,'o', markersize=0.25)


#%%
""" Reading in batches of files """

column_names = ["Filename", "Class", "B alpha", "B error alpha", "lambda0", "lambda0 error", "begin", "finish", "start1", "end1", "start2", "end2", "start3", "end3", "Notes"]
column_names = ["Filename", 'Identifier', "Class", 'Cut_class', "B alpha",\
                "B error alpha", "lambda0", "lambda0 error", 'sigma1', 'sigma2', \
                "begin", "finish", "start1", "end1", "start2", "end2", "start3", "end3", \
                'Notes', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6']

cut_regions = ["begin", "finish", "start1", "end1", "start2", "end2", "start3", "end3"]
datasets = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Detailed class catalogue 2.csv', skiprows = 1, names = column_names)
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

cut_regions = [begin_list, finish_list, start1_list, end1_list, start2_list, end2_list, start3_list, end3_list]   

condition = 1
for i in range(len(datasets)):
    # Hash out beyond .Class[i] == x if you want to read in subsets (otherwise will have all 216 systems)
    """ Bear in mind number_sections needs to be changed to 11 if doing for 216 systems (from 6) """

    if condition == 1:
#    datasets.Class[i] == 1 or datasets.Class[i] == 2 or datasets.Class[i] == 3 or datasets.Class[i] == 4:
        filename_list.append((datasets.Filename[i] + ' ')) # NOTE DESI AND GAIA HAVE DIFFERENT NAMES
                                                            # therefore cannot call as same variable!!
#        for j in range(len(cut_regions)):
#            tot__list[j].append(datasets.begin[j])
        ##%%
        begin_list.append(datasets.begin[i])
        finish_list.append(datasets.finish[i])
        start1_list.append(datasets.start1[i])
        end1_list.append(datasets.end1[i])
        start2_list.append(datasets.start2[i])
        end2_list.append(datasets.end2[i])
        start3_list.append(datasets.start3[i])
        end3_list.append(datasets.end3[i])

#file_name, lin, quad, undec, da, comm = np.loadtxt(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\First categorisation 113 systems.csv',skiprows = 1, delimiter = ',', unpack = True)
datafolder = [filename_list, cut_regions[0], cut_regions[1], cut_regions[2], cut_regions[3], cut_regions[4], \
                cut_regions[5], cut_regions[6], cut_regions[7]]
##%%
""" Bear in mind number_sections needs to be changed to 11 if doing for 216 systems (from 6) """
number_sections = 11 # for class 1 this is only 6 (was 11 for generality but was empty lists past this point)

subsections = [[[] for i in range(0, len(datafolder))] for i in range(number_sections)]

intervals = [0]
[intervals.append(intervals[i]+20) for i in range(0,number_sections)]

for i in range(0,number_sections):
    for j in range(0,len(datafolder)):
        subsections[i][j] = datafolder[j][intervals[i]:intervals[i+1]]


#%% 
number_sections = 11 # for class 1 this is only 6 (was 11 for generality but was empty lists past this point)
intervals = [0]
[intervals.append(intervals[i]+20) for i in range(0,number_sections)]

Gaia_filename_list = [filename_list[i].replace('DESI_','') for i in range(len(filename_list))]
Gaia_filename_list = [Gaia_filename_list[i].replace('_bin0p2.dat','') for i in range(len(filename_list))]

sub_Gaia_filename = [[] for i in range(0, number_sections)]
subsections_Gaia_file = [Gaia_filename_list[intervals[i]:intervals[i+1]] for i in range(0, number_sections)]

#%% Saves data in Excel spreadsheet form - can then simply copy and paste to excel
# RUN AFTER!! CeLL below
aexceldata = np.zeros((len(Exceldata[0]),len(Exceldata)))

for i in range(1, len(Exceldata)):
    aexceldata[:,i] = Exceldata[i]
    
#%%
""" Cell that creates the data in the Gaia catalogue for 216 systems (or batches of 20)"""
Exceldata = [[] for i in range(0,12)] # range(0,x) x has to match length of Gaiavalues

#alpha = 5
#subfolder_ = subsections[alpha] #[1]# combdata # if want to read in batches of 20 files
subfolder_ = [filename_list] # if want to read in 216 files at once (all files)
for j in range(len(subfolder_[0][0:])):
    filename = subfolder_[0][j]

    #filename = "DESI_WDJ153349.03+005916.21_bin0p2.dat"
    data = np.genfromtxt(f'{filename}', delimiter=' ')
    
    wavelength = data[:,0]
    flux = data[:,1]
    error = data[:,2]
                
    #Gaia_data = dataset[dataset['WDJ Name'].isin([subsections_Gaia_file[alpha][j]])] ## for batches of 20
    Gaia_data = dataset[dataset['WDJ Name'].isin([Gaia_filename_list[j]])] ## for all 216 files
    ##%%

    #print(Gaia_data)
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
    #print(T_val)
    
    Gaiavalues = [name_val, Mass_val, Mass_err, BPRP_val, G_abs_val, T_val, eT_val, log_g_val, elog_g_val, \
                  parallax_val, eparallax_val, SN_val]
    
    for i in range(0,len(Gaiavalues)):
        Exceldata[i].append(Gaiavalues[i])
    
    #%% Plotting cell
    textstr = "\n".join((
        'Name: %s' % (name_val),
        'G_abs = %2f' % (G_abs_val, ),
        #'$G_{abs}$ = %2f' % (G_abs_val, ),
        'S/N = %2f' % (SN_val, ),
        'Parallax = %2f pm %.3g' % (parallax_val, eparallax_val),
        'T_eff = %2f pm %.3g' % (T_val, eT_val),
        'log_g = %2f pm %.3g' % (log_g_val, elog_g_val),
        'Mass = %2f pm %.3g' % (Mass_val, Mass_err)))
    
    textstr = "\n".join((
        'Name: %s' % (name_val),
        '$G_{abs}$ = %2f' % (G_abs_val, ),
        'S/N = %2f' % (SN_val, ),
        'Parallax = %2f pm %.3g' % (parallax_val, eparallax_val),
        '$T_{eff}$ = %2f pm %.3g' % (T_val, eT_val),
        'log(g) = %2f pm %.3g' % (log_g_val, elog_g_val),
        'Mass = %2f pm %.3g' % (Mass_val, Mass_err)))
    ##%%
#    print("%.8f pm %.3g" % (c[0], c[1]))
    
    f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
    a0.set_xlabel("Wavelength $\lambda$, $[\AA]$" , size = "15")
    a0.set_ylabel("Flux", size = "15")
    a0.plot(wavelength, flux, label = f"{filename}")
    a1.invert_yaxis()
    a1.set_xlabel(r'$G_{BP} - G_{RP}$' , size = "15")
    a1.set_ylabel(r'$G_{abs}$, $[mag]$' , size = "15")
    a1.plot(BPRP_sel,G_sel,'o', markersize=0.25) # Plots the whole HR diagram
    a1.plot(BPRP_val,G_abs_val,'o',color='red',markersize=5) # Plots the BPRP and abs G for this system 
    props = dict(boxstyle='square', alpha=0.5)
    a1.text(0.05, 0.05, textstr, transform=a1.transAxes, fontsize=9,
            verticalalignment='bottom', bbox=props)
    #f.tight_layout()

#%%
#########################################
    
    # Analysis from now on - read in Gaia data from spreadsheet so do not need to recollect it
    
#########################################
    
    
#columns = ['name_val', 'Mass_val', 'Mass_err', 'BPRP_val', 'G_abs_val', 'T_val', 'eT_val', 'log_g_val', 'elog_g_val', \
#                  'parallax_val', 'eparallax_val', 'SN_val']
Gaia_data_216 = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Gaia catalogue csv.csv', \
                        skiprows = 0)#, unpack = True)#names = columns)
#B, A = np.loadtxt(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Gaia catalogue csv.csv', \
#                        delimiter = ',', skiprows = 1, unpack = True)
Gaia_data_216 = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\Gaia catalogue 2 csv.csv', \
                        skiprows = 0)#, unpack = True)#names = columns)

Gaia_data_216 = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\MSci Project\Catalogues\AlldataDA csv.csv', \
                        skiprows = 0)#, unpack = True)#names = columns)

Gaia_data_216 = pd.read_csv(r'C:\Users\44743\Documents\Imperial Year 4\UROP Year 4\AlldataDA - UROP inc.csv', \
                        skiprows = 0)#, unpack = True)#names = columns)

##%%
##plt.plot(np.array(Gaia_data['BPRP'].tolist()), np.array(Gaia_data['G_abs'].tolist()))
#plt.plot(Gaia_data[3], Gaia_data[4], 'x')
#plt.show()
##%%
#A = np.array([1,2,3])
#B = np.array([4,5,6])
#plt.plot(A, B, 'x')
#A = np.array(Gaia_data[0].tolist())
#
#A = []
##%% 
# !!!
bprp_216 = []
absG_216 = []

[bprp_216.append(Gaia_data_216['BPRP'][i]) for i in range(0, len(Gaia_data_216))]
[absG_216.append(Gaia_data_216['G_abs'][i]) for i in range(0, len(Gaia_data_216))]

mass = []
bfield = []
mass_err = []
bfield_err = []

[mass.append(Gaia_data_216['Mass'][i]) for i in range(0, len(Gaia_data_216))]
[bfield.append(Gaia_data_216['Bfield'][i]) for i in range(0, len(Gaia_data_216))]
[mass_err.append(Gaia_data_216['eMass'][i]) for i in range(0, len(Gaia_data_216))]
[bfield_err.append(Gaia_data_216['eT_eff'][i]) for i in range(0, len(Gaia_data_216))]

#%%
Mass = np.array(mass)
Bfield = np.array(bfield)
eMass = np.array(mass_err)
eBfield = np.array(bfield_err)
#%%
Bfield = Bfield[np.logical_not(np.isnan(Mass))] # remove nans related to mass nans first
Mass = Mass[np.logical_not(np.isnan(Mass))] # remove mass nans (cannot remove B nans after since mass nans now removed)
Mass = Mass[np.logical_not(np.isnan(Bfield))] # remove mass nans (cannot remove B nans after since mass nans now removed)
Bfield = Bfield[np.logical_not(np.isnan(Bfield))] # remove nans related to mass nans first

eBfield = eBfield[np.logical_not(np.isnan(eMass))] # remove nans related to mass nans first
eMass = eMass[np.logical_not(np.isnan(eMass))] # remove mass nans (cannot remove B nans after since mass nans now removed)
eMass = eMass[np.logical_not(np.isnan(eBfield))] # remove mass nans (cannot remove B nans after since mass nans now removed)
eBfield = eBfield[np.logical_not(np.isnan(eBfield))] # remove nans related to mass nans first

#%%
plt.figure()
plt.errorbar(Mass,Bfield, yerr = eBfield, xerr = eMass, color='darkblue', fmt ='x', alpha = 0.2)#, color='darkblue')#, ms = 1)
plt.plot(Mass, Bfield, 'x', color='darkblue')
plt.yscale('log')
plt.show()

arr2d = np.zeros((len(Bfield),2))
arr2d[:,0] = np.array(Mass)
arr2d[:,1] = np.array(Bfield)
#
A = zip(np.corrcoef(arr2d))

sp.stats.pearsonr(Mass, Bfield)
#%%


#plt.figure()
#plt.plot(bprp_216, absG_216, 'x', color = 'orange')


#plt.gca().invert_yaxis()

##%% ESSENTIAL - plots HR diagram (and gets essential data for it )
""" Plotting script for MWD HRD """
# Selects WDs with reciprocal parallaxes/1000 < 100 --> says WDs must be x distance from us (or close to us?)

dist_pc = np.reciprocal(parallax/1000)
idx = np.where(dist_pc < 100)[0]

dist_pc_constrained = [dist_pc[i] for i in idx]
BPRP_sel = [BPRP[i] for i in idx]
G_sel = [G_abs[i] for i in idx]

#plt.figure()
#plt.gca().invert_yaxis()
#plt.xlabel('BP-RP')
#plt.ylabel('G_abs')

#plt.plot(BPRP_sel,G_sel,'o', markersize=0.25)
#plt.plot(bprp, absG, 'x', color = 'orange')


#%%
from scipy.stats import gaussian_kde

bprp = np.array(bprp)
absG = np.array(absG)

xmin = bprp.min()
xmax = bprp.max()
ymin = absG.min()
ymax = absG.max()

X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positionS = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([bprp, absG])
kernel = gaussian_kde(values)
Z = np.reshape(kernel(positionS).T, X.shape)
#Plot the results:

import matplotlib.pyplot as plt

fig, ax = plt.subplots(gridspec_kw={'height_ratios': [200]})
pos = ax.imshow((Z), cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
pos = ax.imshow(np.rot90(Z), cmap=plt.cm.viridis, extent=[xmin, xmax, ymin, ymax])
ax.plot(bprp, absG, 'k.', markersize=2)
fig.colorbar(pos, ax = ax)
#ax.plot(x_p, y_p, 'r.', markersize=2)
#fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax)
#ax.set_xlim([xmin, xmax])
#ax.set_ylim([ymin, ymax])
pl.xlabel("x position", size = "13")
pl.ylabel("y position", size = "13")
plt.gca().invert_yaxis()
#%%
fig, ax = plt.subplots()
pos = plt.imshow((Z), cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
pos = plt.imshow(np.rot90(Z), cmap=plt.cm.viridis, extent=[xmin, xmax, ymin, ymax])
plt.plot(bprp, absG, 'k.', markersize=2)
plt.colorbar(pos, ax = ax)
#ax.plot(x_p, y_p, 'r.', markersize=2)
#fig.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax)
#ax.set_xlim([xmin, xmax])
#ax.set_ylim([ymin, ymax])
plt.xlabel("x position", size = "13")
plt.ylabel("y position", size = "13")
plt.gca().invert_yaxis()
plt.show()

#%% PLotting Cell
    fig, axs = plt.subplots(2, gridspec_kw={'height_ratios': [2.5, 1]})
    fig.suptitle(f"Lorentzian fit {filename} \n B = {popt_3lorentz[1]} +/- {np.sqrt(np.diag(cov_3lorentz))[1]}", \
                 fontsize = "13")
    axs[0].plot(xp_triplet, yp_triplet,'x', label = "WD Ha data")
    axs[0].plot(xp_triplet, _3Lorentzian(xp_triplet, *popt_3lorentz), linewidth=2, color = "orange", \
                                  label = "Lorentzian c_fit")
    #axs[0].text(r"B = {popt_3lorentz}")
    for ax in axs.flat:
        ax.set_xlabel('Wavelength $[\AA]$', fontsize = "13")
        ax.set_ylabel('Normalised flux', fontsize = "13")
    for ax in axs.flat:
        ax.label_outer()
    axs[1].set_ylabel('Flux residuals')
    axs[0].grid()
    #""" Voigt Residuals Plot """
    axs[1].plot(xp_triplet, Residuals, linewidth=2)#, label = "Lorentzian fit Residuals")
    axs[1].plot(xp_triplet, xp_triplet*0/xp_triplet, linewidth = 1)
    axs[0].legend()
    #axs[1].legend()
    axs[1].grid()
    #plt.savefig("Gaussian Residuals", dpi = 1000)
    #plt.savefig(f'{filename}.png', dpi = 1000)
    plt.show()

#%%
from scipy.stats import gaussian_kde

# Generate fake data
x = np.random.normal(size=1000)
y = x * 3 + np.random.normal(size=1000)

# Calculate the point density
xy = np.vstack([bprp,absG])
z = gaussian_kde(xy)(xy)

fig, ax = plt.subplots()
plt.scatter(bprp, bprp, c=z, s=100)
plt.colorbar()
plt.show()



#%%
import matplotlib.pyplot as plt
from matplotlib import colors

#plt.rc('text', usetex=True)

fig, ax = plt.subplots(figsize=(6, 6))
# only show 2D-histogram for bins with more than 10 stars in them
h = ax.hist2d(bprp, absG, bins=300, cmin=10, norm=colors.PowerNorm(0.5), zorder=0.5)
# fill the rest with scatter (set rasterized=True if saving as vector graphics)
ax.scatter(bprp, absG, alpha=0.05, s=1, color='k', zorder=0)

h = ax.hist2d(BPRP_sel, G_sel, bins=150, cmin=10, norm=colors.PowerNorm(0.5), zorder=0.5)
# fill the rest with scatter (set rasterized=True if saving as vector graphics)
ax.scatter(BPRP_sel, G_sel, alpha=0.05, s=1, color='k', zorder=0)

ax.invert_yaxis()
cb = fig.colorbar(h[3], ax=ax, pad=0.02)
ax.set_xlabel(r'$G_{BP} - G_{RP}$')
ax.set_ylabel(r'$M_G$')
cb.set_label(r"$\mathrm{Stellar~density}$")
plt.show()





#%%
from scipy import stats
def measure(n):
    "Measurement model, return two coupled measurements."
    m1 = np.random.normal(size=n)
    m2 = np.random.normal(scale=0.5, size=n)
    return m1+m2, m1-m2

BPRP_sel = np.array(BPRP_sel)
G_sel = np.array(G_sel)

bprp_216 = np.array(bprp_216)
absG_216 = np.array(absG_216)


xmin = BPRP_sel.min()
xmax = BPRP_sel.max()
ymin = G_sel.min()
ymax = G_sel.max()

#m1, m2 = measure(2000)
#xmin = m1.min()
#xmax = m1.max()
#ymin = m2.min()
#ymax = m2.max()

X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([BPRP_sel, G_sel])
kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)

import matplotlib.pyplot as plt
#fig, ax = plt.subplots()
plt.imshow(np.rot90(Z), extent=[xmin, xmax, ymin, ymax])
plt.plot(BPRP_sel, G_sel, 'k.', markersize=2)
#ax.set_xlim([xmin, xmax])
#ax.set_ylim([ymin, ymax])
plt.colorbar()
plt.show()
#%%
from scipy import stats
def measure(n):
    "Measurement model, return two coupled measurements."
    m1 = np.random.normal(size=n)
    m2 = np.random.normal(scale=0.5, size=n)
    return m1+m2, m1-m2
m1, m2 = measure(2000)
xmin = m1.min()
xmax = m1.max()
ymin = m2.min()
ymax = m2.max()
#Perform a kernel density estimate on the data:

X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([X.ravel(), Y.ravel()])
values = np.vstack([m1, m2])
kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)
#Plot the results:

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.imshow(np.rot90(Z), cmap=plt.cm.viridis)#plt.cm.gist_earth_r)#, extent=[xmin, xmax, ymin, ymax])
ax.plot(m1, m2, 'k.', markersize=2)
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
plt.show()
#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
#%%
""" Plots HRD (no density spectrum) with 216 systems (with density) """
# Generate fake data
x = np.random.normal(size=1000)
y = x * 3 + np.random.normal(size=1000)

x = np.array(bprp_216)
y = np.array(absG_216)
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

# HRD positions of 216 Wds
bprp_216 = []
absG_216 = []

[bprp_216.append(Gaia_data_216['BPRP'][i]) for i in range(0, len(Gaia_data_216))]
[absG_216.append(Gaia_data_216['G_abs'][i]) for i in range(0, len(Gaia_data_216))]

fig, ax = plt.subplots()
#ax.scatter(BPRP_sel, G_sel, s = 1, color = 'darkgrey')
ax.scatter(BPRP_sel, G_sel, alpha=1, s=0.05, color='grey', zorder=0) #s = 0.01 or 0.05 or somewhere between
#ax.scatter(BPRP_sel, G_sel, alpha=1, s=0.05, color='grey', zorder=0) #s = 0.01 or 0.05 or somewhere between


ax.scatter(x, y, c=z, s=20, edgecolor='') # MARKER SIZE FOR MWD MARKERS (DENSITY COLOURS)


h = ax.hist2d(bprp_216, absG_216, bins=150, cmin=10, norm=colors.PowerNorm(0.5), zorder=0.5)
# fill the rest with scatter (set rasterized=True if saving as vector graphics)
ax.scatter(bprp_216, absG_216, alpha=0.05, s=0.01, color='k', zorder=0)

cb = fig.colorbar(h[3], ax=ax, pad=0.02)
#cb = fig.colorbar(plt.cm.ScalarMappable(norm=colors.PowerNorm(0.5), cmap=plt.viridis), ax=ax)

#ax.set_xlim([np.amin(BPRP_sel)-0.15, np.amax(BPRP_sel)+0.05])
#ax.set_ylim([np.amin(G_sel)-1, np.amax(G_sel)+1])
ax.set_xlim([np.amin(BPRP_sel)+0.2, np.amax(BPRP_sel)+0.05])
ax.set_ylim([np.amin(G_sel), np.amax(G_sel)])

ax.invert_yaxis()

ax.set_xlabel(r'$G_\mathrm{BP} - G_\mathrm{RP}$ (mag)', fontsize = '13')
ax.set_ylabel(r'$G_\mathrm{abs}$ (mag)', fontsize = '13')
ax.tick_params(axis='both', labelsize=11)
cb.ax.tick_params(axis='both', labelsize=11)

cb.set_label(r"$\mathrm{Relative~density}$", fontsize = '13')
plt.show()

#plt.savefig(f'HRDwithdensityMWDs.pdf', bbox_inches = 'tight')

#%% 
""" HRD diagram with density colourbar """
# Generate fake data
x = np.random.normal(size=1000)
y = x * 3 + np.random.normal(size=1000)

x = np.array(BPRP_sel)
y = np.array(G_sel)
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

### Hash this out for just Gaia HRD density 
# HRD positions of 216 Wds
bprp_216 = []
absG_216 = []

[bprp_216.append(Gaia_data_216['BPRP'][i]) for i in range(0, len(Gaia_data_216))]
[absG_216.append(Gaia_data_216['G_abs'][i]) for i in range(0, len(Gaia_data_216))]
###
##%%
fig, ax = plt.subplots()

ax.scatter(x, y, c=z, s=1, edgecolor='')

ax.scatter(bprp_216, absG_216, s = 5, color = 'tomato') ##


h = ax.hist2d(BPRP_sel, G_sel, bins=150, cmin=10, norm=colors.PowerNorm(0.5), zorder=0.5)
# fill the rest with scatter (set rasterized=True if saving as vector graphics)
ax.scatter(BPRP_sel, G_sel, alpha=0.05, s=1, color='k', zorder=0)

cb = fig.colorbar(h[3], ax=ax, pad=0.02)


#ax.set_xlim([np.amin(BPRP_sel)-0.15, np.amax(BPRP_sel)+0.05])
#ax.set_ylim([np.amin(G_sel)-1, np.amax(G_sel)+1])
ax.set_xlim([np.amin(BPRP_sel)+0.2, np.amax(BPRP_sel)+0.05])
ax.set_ylim([np.amin(G_sel), np.amax(G_sel)])

ax.invert_yaxis()

ax.set_xlabel(r'Colour, $G_{BP} - G_{RP}$', fontsize = '15')
ax.set_ylabel(r'Absolute magnitude, $M_G$', fontsize = '15')
cb.set_label(r"$\mathrm{Stellar~density}$", fontsize = '15')
ax.tick_params(axis='both', labelsize=14)
cb.ax.tick_params(labelsize =14)

#plt.savefig(f'Gaia HRD density.png', dpi = 1000, bbox_inches = 'tight')
plt.show()


#%% 
""" HRD diagram with red markers for 216 systems (no density) """
# Generate fake data
x = np.random.normal(size=1000)
y = x * 3 + np.random.normal(size=1000)

x = np.array(BPRP_sel)
y = np.array(G_sel)
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]
#%%
""" HRD diagram with red markers for 216 systems (no density) """

### Hash this out for just Gaia HRD density 
# HRD positions of 216 Wds
bprp_216 = []
absG_216 = []

[bprp_216.append(Gaia_data_216['BPRP'][i]) for i in range(0, len(Gaia_data_216))]
[absG_216.append(Gaia_data_216['G_abs'][i]) for i in range(0, len(Gaia_data_216))]
###

fig, ax = plt.subplots()
#ax.scatter(x, y, c=z, s=5, edgecolor='')
#ax.scatter(BPRP_sel, G_sel, s = 1, color = 'darkblue')
ax.scatter(BPRP_sel, G_sel, alpha=1, s=0.03, color='darkblue', zorder=0) #s = 0.01 or 0.05 or somewhere between

#ax.scatter(BPRP_sel, G_sel, s=5, edgecolor='') ##


#ax.scatter(bprp_216, absG_216, s = 5, color = 'tomato') ##


#h = ax.hist2d(BPRP_sel, G_sel, bins=150, cmin=10, norm=colors.PowerNorm(0.5), zorder=0.5)
# fill the rest with scatter (set rasterized=True if saving as vector graphics)
#ax.scatter(BPRP_sel, G_sel, alpha=0.05, s=1, color='k', zorder=0)

#cb = fig.colorbar(h[3], ax=ax, pad=0.02)


#ax.set_xlim([np.amin(BPRP_sel)-0.15, np.amax(BPRP_sel)+0.05])
#ax.set_ylim([np.amin(G_sel)-1, np.amax(G_sel)+1])
ax.set_xlim([np.amin(BPRP_sel)+0.2, np.amax(BPRP_sel)+0.05])
ax.set_ylim([np.amin(G_sel), np.amax(G_sel)])

ax.invert_yaxis()

ax.set_xlabel(r'$G_\mathrm{BP} - G_\mathrm{RP}$ (mag)', fontsize = '13')
ax.set_ylabel(r'$G_\mathrm{abs}$ (mag)', fontsize = '13')
ax.tick_params(axis='both', labelsize=11)

#, fontsize = "15", labelpad = 10
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)

cb.set_label(r"$\mathrm{Stellar~density}$") # Stellar since not necessarily WDs (Pwd high but not 1 !!)
#plt.savefig(f'Gaia HRD overplot.png', dpi = 1000, bbox_inches = 'tight')

#plt.savefig(f'HRDnoMWD.pdf', bbox_inches = 'tight')


plt.show()

##%%
#""" Plots HRD (no density spectrum) with 216 systems (with density) """
## Generate fake data
#x = np.random.normal(size=1000)
#y = x * 3 + np.random.normal(size=1000)
#
#x = np.array(bprp_216)
#y = np.array(absG_216)
## Calculate the point density
#xy = np.vstack([x,y])
#z = gaussian_kde(xy)(xy)
#
## Sort the points by density, so that the densest points are plotted last
#idx = z.argsort()
#x, y, z = x[idx], y[idx], z[idx]
#
## HRD positions of 216 Wds
#bprp_216 = []
#absG_216 = []
#
#[bprp_216.append(Gaia_data_216['BPRP'][i]) for i in range(0, len(Gaia_data_216))]
#[absG_216.append(Gaia_data_216['G_abs'][i]) for i in range(0, len(Gaia_data_216))]
#
#fig, ax = plt.subplots()
##ax.scatter(BPRP_sel, G_sel, s = 1, color = 'darkgrey')
#ax.scatter(BPRP_sel, G_sel, alpha=1, s=0.05, color='firebrick', zorder=0) #s = 0.01 or 0.05 or somewhere between
#
#
#ax.scatter(x, y, c=z, s=25, edgecolor='')
#
#
#h = ax.hist2d(bprp_216, absG_216, bins=150, cmin=10, norm=colors.PowerNorm(0.5), zorder=0.5)
## fill the rest with scatter (set rasterized=True if saving as vector graphics)
#ax.scatter(bprp_216, absG_216, alpha=0.05, s=1, color='k', zorder=0)
#
#cb = fig.colorbar(h[3], ax=ax, pad=0.02)
##cb = fig.colorbar(plt.cm.ScalarMappable(norm=colors.PowerNorm(0.5), cmap=plt.viridis), ax=ax)
#
##ax.set_xlim([np.amin(BPRP_sel)-0.15, np.amax(BPRP_sel)+0.05])
##ax.set_ylim([np.amin(G_sel)-1, np.amax(G_sel)+1])
#ax.set_xlim([np.amin(BPRP_sel)+0.2, np.amax(BPRP_sel)+0.05])
#ax.set_ylim([np.amin(G_sel), np.amax(G_sel)])
#
#ax.invert_yaxis()
#
#ax.set_xlabel(r'$G_{BP} - G_{RP}$')
#ax.set_ylabel(r'$M_G$')
#cb.set_label(r"$\mathrm{Stellar~density}$")
#plt.show()
#%%
""" Plots HRD (no density spectrum) with 216 systems (with density) """
# Generate fake data
x = np.random.normal(size=1000)
y = x * 3 + np.random.normal(size=1000)

bprp_216 = np.array(bprp_216)
absG_216 = np.array(absG_216)

bprp_216 = bprp_216[np.logical_not(np.isnan(absG_216))] # remove nans related to mass nans first
absG_216 = absG_216[np.logical_not(np.isnan(absG_216))] # remove mass nans (cannot remove B nans after since mass nans now removed)

x = np.array(bprp_216)
y = np.array(absG_216) # whatever is in x, y is the density plot part!
# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

# HRD positions of 216 Wds
#bprp_216 = []
#absG_216 = []
#
#[bprp_216.append(Gaia_data_216['BPRP'][i]) for i in range(0, len(Gaia_data_216))]
#[absG_216.append(Gaia_data_216['G_abs'][i]) for i in range(0, len(Gaia_data_216))]

fig, ax = plt.subplots()
#ax.scatter(BPRP_sel, G_sel, s = 1, color = 'darkgrey')
ax.scatter(BPRP_sel, G_sel, alpha=1, s=0.05, color='grey', zorder=0) #s = 0.01 or 0.05 or somewhere between
#ax.scatter(BPRP_sel, G_sel, alpha=1, s=0.05, color='grey', zorder=0) #s = 0.01 or 0.05 or somewhere between


ax.scatter(x, y, c=z, s=20, edgecolor='') # MARKER SIZE FOR MWD MARKERS (DENSITY COLOURS)


h = ax.hist2d(bprp_216, absG_216, bins=150, cmin=10, norm=colors.PowerNorm(0.5), zorder=0.5)
# fill the rest with scatter (set rasterized=True if saving as vector graphics)
ax.scatter(bprp_216, absG_216, alpha=0.05, s=0.01, color='k', zorder=0)

cb = fig.colorbar(h[3], ax=ax, pad=0.02)
#cb = fig.colorbar(plt.cm.ScalarMappable(norm=colors.PowerNorm(0.5), cmap=plt.viridis), ax=ax)

#ax.set_xlim([np.amin(BPRP_sel)-0.15, np.amax(BPRP_sel)+0.05])
#ax.set_ylim([np.amin(G_sel)-1, np.amax(G_sel)+1])
ax.set_xlim([np.amin(BPRP_sel)+0.2, np.amax(BPRP_sel)+0.05])
ax.set_ylim([np.amin(G_sel), np.amax(G_sel)])

ax.invert_yaxis()

ax.set_xlabel(r'$G_\mathrm{BP} - G_\mathrm{RP}$ (mag)', fontsize = '13')
ax.set_ylabel(r'$G_\mathrm{abs}$ (mag)', fontsize = '13')
ax.tick_params(axis='both', labelsize=11)
cb.ax.tick_params(axis='both', labelsize=11)

cb.set_label(r"$\mathrm{Relative~density}$", fontsize = '13')
plt.show()

#plt.savefig(f'HRDwithdensityMWDs.pdf', bbox_inches = 'tight')

#%%
""" HRD with colourbar for magnetic field (everything else no density) """

import matplotlib.pyplot as plt
cm = plt.cm.get_cmap('RdYlBu')
#xy = range(20)
#z = xy

bprp_216 = []
absG_216 = []

[bprp_216.append(Gaia_data_216['BPRP'][i]) for i in range(0, len(Gaia_data_216))]
[absG_216.append(Gaia_data_216['G_abs'][i]) for i in range(0, len(Gaia_data_216))]

bprp_216 = np.array(bprp_216)
absG_216 = np.array(absG_216)
Bfieldlist = np.array(datalist[0])
#Bfieldlist = np.array(bfield)
Bfieldlist = Bfieldlist[np.logical_not(np.isnan(absG_216))]
bprp_216 = bprp_216[np.logical_not(np.isnan(absG_216))] # remove nans related to mass nans first
absG_216 = absG_216[np.logical_not(np.isnan(absG_216))] # remove mass nans (cannot remove B nans after since mass nans now removed)

x = bprp_216
y = absG_216
z = Bfieldlist

idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

plt.figure()


plt.xlim([np.amin(BPRP_sel)+0.2, np.amax(BPRP_sel)+0.05])
plt.ylim([np.amin(G_sel), np.amax(G_sel)])

plt.scatter(BPRP_sel, G_sel, alpha=1, s=0.5, color='grey', zorder=0) #s = 0.01 or 0.05 or somewhere between
plt.scatter(x, y, s = 5)
plt.gca().invert_yaxis()
cmap = plt.cm.winter_r
sc = plt.scatter(x, y, c= z, s=10, norm=colors.LogNorm(vmin=z.min(), vmax=z.max()), cmap=cmap) # vmin=0, vmax=20, 
plt.colorbar(sc)

plt.xlabel(r'$G_\mathrm{BP} - G_\mathrm{RP}$ (mag)')#, fontsize = '13')
plt.ylabel(r'$G_\mathrm{abs}$ (mag)')#, fontsize = '13')
#ax.tick_params(axis='both', labelsize=11)
#cb.ax.tick_params(axis='both', labelsize=11)

#plt.colorbar.set_label(r"$\mathrm{B field}$")#, fontsize = '13')

plt.show()















