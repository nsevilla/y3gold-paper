#!/usr/bin/env python
'''
This script will use a pre-matched catalog to truth to plot detection completeness or efficiency/contamination as a function of magnitude
Author: Nacho Sevilla
Usage: python plot_matched.py
'''
import os,sys
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from descolors import BAND_COLORS
from optparse import OptionParser

matplotlib.style.use('des_dr1')

workdir = '/Users/nsevilla/y3gold-paper/'
datadir = '/Users/nsevilla/y3gold-paper/data/'
figsdir = '/Users/nsevilla/y3gold-paper/figs/'

def plot_detection_completeness_matched(tdata, binvar, minmax, binning):
    compl = np.empty(binning)
    dcompl = np.empty(binning)
    midbins = np.empty(binning)
    interval = float(minmax[1]-minmax[0])/float(binning)
    for i in range(binning):
        lo = minmax[0]+i*interval
        midbins[i] = lo + interval*0.5
        mask_match = (tdata[binvar] > lo) & (tdata[binvar] < lo + interval) & (tdata['detected'] > 0.5)
        mask = (tdata[binvar] > lo) & (tdata[binvar] < lo + interval)
        detgal = len(tdata[mask_match])
        allgal = len(tdata[mask])
        print(lo,lo+interval,float(detgal),float(allgal),float(detgal)/float(allgal))
        compl[i] = float(detgal)/float(allgal)
        dcompl[i] = 1/float(allgal)
        dcompl[i] = dcompl[i]*np.sqrt(float(detgal)*(1-dcompl[i])) #binomial error, temporary
    plt.errorbar(midbins,compl,yerr=dcompl,marker='o',color='red',mfc='red',ecolor='red')
    print(repr(compl))
    print(repr(dcompl))
    plt.xticks(np.arange(minmax[0], minmax[1]+1, 0.5))
    plt.hlines(0.90,19,25)
    plt.xlabel(binvar,fontsize=14)
    plt.ylabel('Completeness',fontsize=14)
    plt.ylim(0.0,1.0)
    plt.grid(True)
    plt.title('Completeness of Balrog objects', fontsize=16)
    plt.savefig(figsdir+'completeness_galaxies_with_balrog_test.png')

def plot_eff_cont_matched(tdata, binvar, minmax, binning):
    ppv = np.empty(binning)
    dppv = np.empty(binning)
    dppv_lo = np.empty(binning)
    dppv_hi = np.empty(binning)
    tpr = np.empty(binning)
    dtpr = np.empty(binning)
    midbins = np.empty(binning)
    interval = float(minmax[1]-minmax[0])/float(binning)
    ref_class = 'extragalactic'
    #data_class = 'EXTENDED_CLASS_MASH_SOF'
    data_class = 'EXTENDED_CLASS_COADD'
    truth_th = 0.5 
    ths = [0.5,1.5,2.5] #for Y3 GOLD
    ### note that the following procedure is constructed for galaxies
    colors = [BAND_COLORS['u'],BAND_COLORS['g'],BAND_COLORS['r']]
    for t,th in enumerate(ths):
        for i in range(binning):
            lo = minmax[0]+i*interval
            midbins[i] = lo + interval*0.5
            mask = (tdata[binvar] > lo) & (tdata[binvar] < lo + interval) & (tdata[data_class] > th) # > th Y3 GOLD 
            #print(len(tdata_2[ref_class][matches['i2']]))
            tp = sum(tdata[ref_class][mask] > truth_th) #> truth_th HSC, == 1 for ACS
            fp = sum(tdata[ref_class][mask] < truth_th) #< truth_th HSC, == 2 for ACS
            #print(lo,'-',lo+interval,tp,fp)
            ppv[i] = float(tp)/float(tp+fp)
            dppv[i] = 1/float(tp+fp) 
            dppv[i] = dppv[i]*np.sqrt(float(tp)*(1-ppv[i])) #binomial error, temporary
            #einterval = efficiencyError.efficiencyError(float(tp+fp),float(tp),0.95).calculate() # See Paterno 2004 Fermilab note
            #print interval[1],'-',interval[2]
            #dppv_lo[i] = einterval[0]-einterval[1]
            #dppv_hi[i] = einterval[2]-einterval[0]
            mask = (tdata[binvar] > lo) & (tdata[binvar] < lo + interval) & (tdata[data_class] < th) # < th Y3 GOLD
            tn = sum(tdata[ref_class][mask] > truth_th) #> truth_th HSC, == 1 for ACS      
            tpr[i] = float(tp)/float(tp+tn)
            #mask = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval)
            #gals =  sum(tdata_2[ref_class][matches['i2']][mask] == truth)
            #tpr[i] = float(tp)/float(gals)
            dtpr[i] = 1/float(tp+fp)
            dtpr[i] = dtpr[i]*np.sqrt(float(tp)*(1-tpr[i])) #binomial error, temporary
            print(midbins[i],float(tp+fp),(1.0-ppv[i])*100,dppv[i]*100,tpr[i]*100,dtpr[i]*100)
        #plt.errorbar(midbins,1.0-ppv,yerr=[dppv_lo,dppv_hi],marker='o',label='Contamination MASH >= '+str(th+0.5))
        plt.errorbar(midbins,1.0-ppv,yerr=dppv,marker='o',label='Contamination MASH $\geq$ '+str(th+0.5),color=colors[t])
        plt.errorbar(midbins,tpr,yerr=dtpr,marker='+',label='Efficiency MASH $\geq$ '+str(th+0.5),color=colors[t],ls='dashed')
    plt.xlabel('i-band magnitude')
    plt.ylabel('Galaxy efficiency/contamination')
    plt.hlines(0.95,minmax[0],minmax[1])
    plt.ylim(0.0,1.0)
    plt.title('Efficiency/contamination vs VHS classification')
    plt.legend(loc='center right')
    plt.savefig(figsdir+'effpur_vs_vhs_coadd_test.png')
    return tpr,dtpr,ppv,dppv_lo,dppv_hi ### this will return the values for the highest threshold th


def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--detection_completeness",action="store_true",dest="measure_completeness",help="Toggle detection completeness",default=False)
    parser.add_option("--efficiency_contamination",action="store_true",dest="measure_effcont",help="Toggle efficiency/contamination",default=False)
    (options, args) = parser.parse_args()

    if options.measure_completeness:
        hdulist = fits.open(datadir+'matched_balrog_x3_v2.fits',memmap=True)
        tdata = hdulist[1].data
        binvar = 'bdf_mag_i'
        minmax = [19,25]
        binning = 10
        plot_detection_completeness_matched(tdata, binvar, minmax, binning)
    if options.measure_effcont:
        hdulist = fits.open(datadir+'goldvhs_y3sgsep_v3.fits',memmap=True)
        tdata = hdulist[1].data
        binvar = 'SOF_CM_MAG_I'
        minmax = [15,21]
        binning = 10
        plot_eff_cont_matched(tdata, binvar, minmax, binning)

if __name__ == "__main__":
    main()
