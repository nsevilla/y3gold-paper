#!/usr/bin/env python
'''
This script will match catalogs for detection completeness and efficiency/contamination measurements
Author: Nacho Sevilla
Usage: python match_catalogs.py
'''
import os,sys
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
import smatch
import warnings
warnings.filterwarnings('ignore')
from optparse import OptionParser

workdir = '/Users/nsevilla/y3gold-paper/'
datadir = '/Users/nsevilla/y3gold-paper/data/'
figsdir = '/Users/nsevilla/y3gold-paper/figs/'

def match_cat(tdata_1,tdata_2,radius):
    print('Matching catalogs')
    maxmatch = 1 
    NSIDE = 4096
    
    ra_1 = tdata_1['alphawin_j2000']
    dec_1 = tdata_1['deltawin_j2000']
    #class_1 = tdata_1['extended_class_mash_sof']
    #mag_1 = tdata_1['mag_auto_i']
    #ra_1 = tdata_1['RA']
    #dec_1 = tdata_1['DEC']
    ra_2 = tdata_2['ra']
    dec_2 = tdata_2['dec']
    #class_2 = tdata_2['iclassification_extendedness']
    #mag_2 = tdata_2['imag_kron']
    print(radius)
    matches = smatch.match(ra_1, dec_1, radius, ra_2, dec_2, nside=NSIDE, maxmatch=maxmatch)
    
    return matches

def plot_detection_completeness(tdata_1,tdata_2,matches,minmax,binning,binvar,field,reference):
    print('Plotting detection completeness')
    compl = np.empty(binning)
    dcompl = np.empty(binning)
    midbins = np.empty(binning)
    interval = float(minmax[1]-minmax[0])/float(binning)
    if reference == 'deep':
        ref_class = 'bdf_T'
        general_mask = (tdata_2['mask_flags'] == 0) & (tdata_2['flags'] == 0) & (tdata_2['bdf_T'] < 30)
        general_mask_matched = (tdata_2['mask_flags'][matches['i2']] == 0) & (tdata_2['flags'][matches['i2']] == 0)  & (tdata_2['bdf_T'][matches['i2']] < 30)
    elif reference == 'HSC':
        ref_class = 'i_extendedness_value'
        general_mask = (tdata_2['i_extendedness_value'] > -1)
        general_mask_matched = (tdata_2['i_extendedness_value'][matches['i2']] > -1)
        #general_mask = (tdata_1['EXTENDED_CLASS_MASH_SOF'] > -1)
        #general_mask_matched = (tdata_1['EXTENDED_CLASS_MASH_SOF'][matches['i1']] > -1)
    else:
        ref_class = None
    thresholds = {"deep":0.01,"Balrog":0,"HSC":0.5}
 
    for typ in ['Stars','Galaxies']:
        for i in range(binning):
            lo = minmax[0]+i*interval
            midbins[i] = lo + interval*0.5
            if field == 'vvds':
                mask_match = (tdata_2[binvar][matches['i2']] > lo) & (tdata_2[binvar][matches['i2']] < lo + interval) \
                & (tdata_2['ra'][matches['i2']] > 337) & general_mask_matched
                mask = (tdata_2[binvar] > lo) & (tdata_2[binvar] < lo + interval) & (tdata_2['ra'] > 337) & general_mask
            else:
                #mask_match = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval)
                #mask = (tdata_1[binvar] > lo) & (tdata_1[binvar] < lo + interval)
                mask_match = (tdata_2[binvar][matches['i2']] > lo) & (tdata_2[binvar][matches['i2']] < lo + interval) & general_mask_matched
                mask = (tdata_2[binvar] > lo) & (tdata_2[binvar] < lo + interval) & general_mask
            if typ == 'Galaxies':
                detgal = sum(tdata_2[ref_class][matches['i2']][mask_match] > thresholds[reference])
                allgal = sum(tdata_2[ref_class][mask] > thresholds[reference])
                #detgal = len(tdata_1[matches['i1']][mask_match])
                #allgal = len(tdata_1[mask])
                #print(lo,lo+interval,detgal,allgal,float(detgal)/float(allgal))
            else:
                detgal = sum(tdata_2[ref_class][matches['i2']][mask_match] < thresholds[reference]) 
                allgal = sum(tdata_2[ref_class][mask] < thresholds[reference])
                #detgal = len(tdata_1[matches['i1']][mask_match])
                #allgal = len(tdata_1[mask])
            compl[i] = float(detgal)/float(allgal)
            dcompl[i] = 1/float(allgal)
            dcompl[i] = dcompl[i]*np.sqrt(float(detgal)*(1-dcompl[i])) #binomial error, temporary
            if (i > 0) and (compl[i] < 0.90) and (compl[i-1] > 0.90):
                mag90 = midbins[i-1]+(midbins[i]-midbins[i-1])*(0.90-compl[i-1])/(compl[i]-compl[i-1])
                print('Completeness at ',mag90,'is ~90% for',typ)
            if (i > 0) and (compl[i] < 0.95) and (compl[i-1] > 0.95):
                mag95 = midbins[i-1]+(midbins[i]-midbins[i-1])*(0.95-compl[i-1])/(compl[i]-compl[i-1])
                print('Completeness at ',mag95,'is ~95% for',typ)
            
        if typ == 'Galaxies':
            plt.errorbar(midbins,compl,yerr=dcompl,color='red',marker='o',label=typ+' '+field)
            print(repr(compl))
            print(repr(dcompl))
        else:
            plt.errorbar(midbins,compl,yerr=dcompl,color='blue',marker='o',label=typ+' '+field,linestyle = '--')            
    plt.xticks(np.arange(minmax[0], minmax[1]+1, 0.5))
    plt.hlines(0.90,19,25)
    plt.xlabel(binvar,fontsize=14)
    plt.ylabel('Completeness',fontsize=14)
    plt.ylim(0.0,1.0)
    plt.title('Completeness vs '+reference+' objects', fontsize=16)
    plt.legend(loc='lower left',fontsize=14)
    plt.grid(True)
    plt.savefig(figsdir+'completeness_galaxies_vs_'+reference+'_test.png')

    # for stars

def plot_eff_cont(tdata_1,tdata_2,matches,minmax,binning,binvar,field,reference):
    ppv = np.empty(binning)
    dppv = np.empty(binning)
    tpr = np.empty(binning)
    dtpr = np.empty(binning)
    midbins = np.empty(binning)
    interval = float(minmax[1]-minmax[0])/float(binning)
    ref_class = 'i_extendedness_value'
    data_class = 'extended_class_mash_sof'
    ref_idx = 'i2'
    truth = 0 #extendedness for HSC
    ths = [0.5,1.5,2.5]
    for th in ths:
        print('Threshold <=',th-0.5)
        for i in range(binning):
            lo = minmax[0]+i*interval
            midbins[i] = lo + interval*0.5
            mask = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval) & (tdata_1[data_class][matches['i1']] < th)
            #print(len(tdata_2[ref_class][matches['i2']]))
            tp = sum(tdata_2[ref_class][matches['i2']][mask] == truth)
            fp = sum(tdata_2[ref_class][matches['i2']][mask] != truth)
            #print(lo,'-',lo+interval,tp,fp)
            ppv[i] = float(tp)/float(tp+fp)
            dppv[i] = 1/float(tp+fp)
            dppv[i] = dppv[i]*np.sqrt(float(tp)*(1-dppv[i])) #binomial error, temporary
            mask = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval) & (tdata_1[data_class][matches['i1']] > th)
            tn = sum(tdata_2[ref_class][matches['i2']][mask] == truth)        
            tpr[i] = float(tp)/float(tp+tn)
            #mask = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval)
            #gals =  sum(tdata_2[ref_class][matches['i2']][mask] == truth)
            #tpr[i] = float(tp)/float(gals)
            dtpr[i] = 1/float(tp+fp)
            dtpr[i] = dtpr[i]*np.sqrt(float(tp)*(1-dtpr[i])) #binomial error, temporary
            print(midbins[i],(1.0-ppv[i])*100,tpr[i]*100)
        plt.errorbar(midbins,1.0-ppv,yerr=dppv,marker='o',label='Contamination '+field+' MASH <= '+str(th+0.5))
        plt.errorbar(midbins,tpr,yerr=dtpr,marker='+',label='Efficiency '+field+' MASH <= '+str(th+0.5))
    plt.xlabel('SOF i-band magnitude',fontsize=14)
    plt.ylabel('Galaxy efficiency/contamination',fontsize=14)
    plt.hlines(0.95,19,25)
    plt.ylim(0.0,1.0)
    plt.title('Efficiency/contamination vs HSC classification: '+ field, fontsize=16)
    plt.legend(loc='lower left',fontsize=10)
    plt.savefig(figsdir+'effpur_vs_'+reference+'_'+field+'_test.png')

def colorterm1(x,A1,B):
    return A1*x[0]+x[1]+B
def colorterm2(x,A1,A2,B):
    return A1*x[0]+A2*x[1]+x[2]+B

def estimate_color_terms(tdata_1,tdata_2,matches):

    g_1 = tdata_1['sof_psf_mag_updated_g'][matches['i1']]
    g_2 = tdata_2['g_psfflux_mag'][matches['i2']]
    r_1 = tdata_1['sof_psf_mag_updated_r'][matches['i1']]
    r_2 = tdata_2['r_psfflux_mag'][matches['i2']]
    i_1 = tdata_1['sof_psf_mag_updated_i'][matches['i1']]
    i_2 = tdata_2['i_psfflux_mag'][matches['i2']]
    z_1 = tdata_1['sof_psf_mag_updated_z'][matches['i1']]
    z_2 = tdata_2['z_psfflux_mag'][matches['i2']]
    gr_2 = tdata_2['g_psfflux_mag'][matches['i2']]-tdata_2['r_psfflux_mag'][matches['i2']]
    ri_2 = tdata_2['r_psfflux_mag'][matches['i2']]-tdata_2['i_psfflux_mag'][matches['i2']]
    iz_2 = tdata_2['i_psfflux_mag'][matches['i2']]-tdata_2['z_psfflux_mag'][matches['i2']]
    zy_2 = tdata_2['z_psfflux_mag'][matches['i2']]-tdata_2['y_psfflux_mag'][matches['i2']]
    pg = curve_fit(colorterm1,[gr_2,g_2],g_1)
    pr = curve_fit(colorterm1,[ri_2,r_2],r_1)
    pi = curve_fit(colorterm1,[iz_2,i_2],i_1)
    #pz = curve_fit(colorterm2,[iz_2,zy_2,z_2],z_1)
    pz = curve_fit(colorterm1,[zy_2,z_2],z_1)
    print(pg[0],np.std(pg[0][0]*gr_2+g_2+pg[0][1]-g_1))
    print(pr[0],np.std(pr[0][0]*ri_2+r_2+pr[0][1]-r_1))
    print(pi[0],np.std(pi[0][0]*iz_2+i_2+pi[0][1]-i_1))
    #print(pz[0],np.std(pz[0][0]*iz_2+pz[0][1]*zy_2+i_2+pz[0][2]-z_1))   
    print(pz[0],np.std(pz[0][0]*zy_2+z_2+pz[0][1]-z_1))
    plt.hist(pg[0][0]*gr_2+g_2+pg[0][1]-g_1,bins=100,histtype='step',range=[-0.5,0.5])
    plt.savefig('testg.png')
    plt.clf()
    plt.hist(pr[0][0]*ri_2+r_2+pr[0][1]-r_1,bins=100,histtype='step',range=[-0.5,0.5])
    plt.savefig('testr.png')
    plt.clf()
    plt.hist(pi[0][0]*iz_2+i_2+pi[0][1]-i_1,bins=100,histtype='step',range=[-0.5,0.5])
    plt.savefig('testi.png')
    plt.clf()
    plt.hist(pz[0][0]*zy_2+z_2+pz[0][1]-z_1,bins=100,histtype='step',range=[-0.5,0.5])
    plt.savefig('testz.png')

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("--detection_completeness",action="store_true",dest="measure_completeness",help="Toggle detection completeness",default=False)
    parser.add_option("--efficiency_contamination",action="store_true",dest="measure_effcont",help="Toggle efficiency/contamination",default=False)
    parser.add_option("--estimate_color_terms",action="store_true",dest="color_terms",help="Toggle color term calibration",default=False)
    parser.add_option("--reference",dest="reference",help="Reference catalog to use (HSC/deep/Balrog)",default='HSC')
    parser.add_option("--band",dest="band",help="Reference band",default='i')
    (options, args) = parser.parse_args()

    #fields = ['sxds','deep23','vvds']
    fields = ['snx3']
    #fields = ['w05']
    
    if options.measure_effcont:
        print('Measuring efficiency/contamination')
        for field in fields:
            hdulist = fits.open(datadir+field+'_y3gold_completeness_Arctmasked_S18grizmasked.fits',memmap=True)
            tdata_1 = hdulist[1].data
            hdulist = fits.open(datadir+field+'_hscw02_completeness_goldmasked_Arctmasked_S18grizmasked_cut.fits',memmap=True)
            tdata_2 = hdulist[1].data    
            print('Matching in field',field)
            matches = match_cat(tdata_1,tdata_2,radius=0.5/3600)
            plot_eff_cont(tdata_1,tdata_2,matches,minmax=[19,22.5],binning=1,binvar='SOF_CM_MAG_CORRECTED_I',field=field,reference=options.reference)
    if options.measure_completeness:
        print('Measuring detection completeness')
        for field in fields:
            if options.reference == 'HSC':
                hdulist1 = fits.open(datadir+field+'_y3gold_y1y3goldmasked_Arctmasked_S18grizmasked.fits',memmap=True)
                hdulist2 = fits.open(datadir+field+'_hsc_y1y3goldmasked_Arctmasked_S18grizmasked_sofconv_cut.fits',memmap=True)
            elif options.reference == 'deep':
                hdulist1 = fits.open(datadir+field+'_y3gold_completeness_dfmasked.fits',memmap=True)
                hdulist2 = fits.open(datadir+field+'_df_masked_goldmasked_selectedcols.fits',memmap=True)
            else:
                print('Data not found',options.reference)
                sys.exit()
            tdata_1 = hdulist1[1].data
            tdata_2 = hdulist2[1].data    
            print('Matching in field',field)
            matches = match_cat(tdata_1,tdata_2,radius = 2.0/3600)
            if options.reference == 'HSC':
                #binvar = options.band+'_cmodel_mag'
                binvar = 'sof_converted_'+options.band
                #binvar='sof_cm_mag_corrected_'+options.band
            elif options.reference == 'deep':
                binvar='bdf_mag_'+options.band
            else:
                binvar='sof_cm_mag_corrected_'+options.band
            print(binvar,field,options.reference)
            plot_detection_completeness(tdata_1,tdata_2,matches,minmax=[19,25],binning=30,binvar=binvar,field=field,reference=options.reference)
    if options.color_terms:
        print('Calibrating color terms')
        hdulist = fits.open(datadir+'y3gold_hscreg_stars_calibration_b.fits',memmap=True)
        tdata_1 = hdulist[1].data
        hdulist = fits.open(datadir+'hsc_stars_calibration.fits',memmap=True)
        tdata_2 = hdulist[1].data    
        matches = match_cat(tdata_1,tdata_2,radius=0.5/3600)
        estimate_color_terms(tdata_1,tdata_2,matches)

if __name__ == "__main__":
    main()
