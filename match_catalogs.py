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
    ppv = np.empty(binning)
    dppv = np.empty(binning)
    tpr = np.empty(binning)
    dtpr = np.empty(binning)
    midbins = np.empty(binning)
    interval = float(minmax[1]-minmax[0])/float(binning)
    ref_class = 'iclassification_extendedness'
    data_class = 'extended_class_mash_sof'
    ref_idx = 'i2'
    truth = 1
    th = 2.5
    for i in range(binning):
        lo = minmax[0]+i*interval
        midbins[i] = lo + interval*0.5
        mask = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval) & (tdata_1[data_class][matches['i1']] > th)
        #print(len(tdata_2[ref_class][matches['i2']]))
        tp = sum(tdata_2[ref_class][matches['i2']][mask] == truth)
        fp = sum(tdata_2[ref_class][matches['i2']][mask] != truth)
        print(lo,'-',lo+interval,tp,fp)
        ppv[i] = float(tp)/float(tp+fp)
        dppv[i] = 1/float(tp+fp)
        dppv[i] = dppv[i]*np.sqrt(float(tp)*(1-dppv[i])) #binomial error, temporary
        mask = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval) & (tdata_1[data_class][matches['i1']] < th)
        tn = sum(tdata_2[ref_class][matches['i2']][mask] == truth)        
        tpr[i] = float(tp)/float(tp+tn)
        #mask = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval)
        #gals =  sum(tdata_2[ref_class][matches['i2']][mask] == truth)
        #tpr[i] = float(tp)/float(gals)
        dtpr[i] = 1/float(tp+fp)
        dtpr[i] = dtpr[i]*np.sqrt(float(tp)*(1-dtpr[i])) #binomial error, temporary
        print(midbins[i],1.0-ppv[i],tpr[i])
    #plt.plot(midbins,success)
    #plt.figure()
    plt.errorbar(midbins,1.0-ppv,yerr=dppv,marker='o',label='Contamination '+field)
    plt.errorbar(midbins,tpr,yerr=dtpr,marker='+',label='Efficiency '+field)
    plt.xlabel(binvar,fontsize=14)
    plt.ylabel('Galaxy efficiency/contamination',fontsize=14)
    plt.ylim(0.0,1.0)
    #plt.title('Efficiency/contamination vs HSC classification: '+ field, fontsize=16)
    plt.title('Efficiency/contamination vs HSC classification', fontsize=16)
    plt.legend(loc='lower left',fontsize=14)
    #plt.savefig('effpur_vs_HSC_'+field+'.png')
    plt.savefig('effpur_vs_HSC.png')

    
radius = 0.5/3600 #in degrees

fields = ['sxds','deep23','vvds']
#fields = ['sxds']

for field in fields:
    print('Matching in field',field)
    hdulist = fits.open(datadir+'field_'+field+'_y3gold.fits',memmap=True)
    tdata_1 = hdulist[1].data
    hdulist = fits.open(datadir+'field_'+field+'_hsc.fits',memmap=True)
    tdata_2 = hdulist[1].data    
    matches = match_cat(tdata_1,tdata_2,radius=radius)
    plot_eff_cont(tdata_1,tdata_2,matches,minmax=[20,25],binning=10,binvar='mag_auto_i',field=field)
