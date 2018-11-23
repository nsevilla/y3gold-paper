import os,sys
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
import smatch
import warnings
warnings.filterwarnings('ignore')

workdir = '/Users/nsevilla/y3gold-paper/'
datadir = '/Volumes/NO NAME/'
figsdir = '/Users/nsevilla/y3gold-paper/figs/'

def match_cat(tdata_1,tdata_2,radius=0.5/3600):
    maxmatch = 1 
    NSIDE = 4096
    
    ra_1 = tdata_1['ra']
    dec_1 = tdata_1['dec']
    #class_1 = tdata_1['extended_class_mash_sof']
    #mag_1 = tdata_1['mag_auto_i']
    ra_2 = tdata_2['ra']
    dec_2 = tdata_2['dec']
    #class_2 = tdata_2['iclassification_extendedness']
    #mag_2 = tdata_2['imag_kron']

    matches = smatch.match(ra_1, dec_1, radius, ra_2, dec_2, nside=NSIDE, maxmatch=maxmatch)
    
    return matches

def plot_eff_cont(tdata_1,tdata_2,matches,minmax,binning,binvar,field):
    success = np.empty(binning)
    dsuccess = np.empty(binning)
    midbins = np.empty(binning)
    interval = float(minmax[1]-minmax[0])/float(binning)
    ref_class = 'iclassification_extendedness'
    data_class = 'extended_class_mash_sof'
    ref_idx = 'i2'
    truth = 1
    th = 2.5
    for i in range(binning):
        lo = minmax[0]+i*interval
        midbins[i] = minmax[0]+(i*interval)/2.0
        mask = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval) & (tdata_1[data_class][matches['i1']] > th)
        #print(len(tdata_2[ref_class][matches['i2']]))
        good = sum(tdata_2[ref_class][matches['i2']][mask] == truth)
        bad = sum(tdata_2[ref_class][matches['i2']][mask] != truth)
        print(lo,'-',lo+interval,good,bad)
        success[i] = float(good)/float(good+bad)
        dsuccess[i] = 1/float(good+bad)
        dsuccess[i] = dsuccess[i]*np.sqrt(float(good)*(1-dsuccess[i]))
        print(success[i],dsuccess[i])
    print(midbins,success)
    #plt.plot(midbins,success)
    plt.errorbar(midbins,success,yerr=dsuccess,marker='o',label=data_class+' classifier')
    plt.xlabel(binvar,fontsize=14)
    plt.ylabel('Galaxy purity',fontsize=14)
    plt.ylim(0.0,1.0)
    plt.title('Purity vs HSC classification: '+ field, fontsize=16)
    plt.legend(loc='lower right',fontsize=14)
    plt.savefig('purity_vs_HSC_'+field+'.png')

    
radius = 0.5/3600 #in degrees

#fields = ['sxds','deep23','vvds']
fields = ['sxds']

for field in fields:
    print('Matching in field',field)
    hdulist = fits.open(datadir+'field_'+field+'_y3gold.fits',memmap=True)
    tdata_1 = hdulist[1].data
    hdulist = fits.open(datadir+'field_'+field+'_hsc.fits',memmap=True)
    tdata_2 = hdulist[1].data    
    matches = match_cat(tdata_1,tdata_2,radius=radius)
    plot_eff_cont(tdata_1,tdata_2,matches,minmax=[20,25],binning=10,binvar='mag_auto_i',field=field)
