#!/usr/bin/env python
'''
This script will match catalogs for detection completeness and efficiency/contamination measurements
Author: Nacho Sevilla
Usage: python match_catalogs.py [OPTIONS]
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
import efficiencyError
warnings.filterwarnings('ignore')
from descolors import BAND_COLORS
from optparse import OptionParser

matplotlib.style.use('des_dr1')

def match_cat(tdata_1,tdata_2,radius,data1cat,data2cat):
    print('Matching catalogs')
    maxmatch = 1 
    NSIDE = 4096

    if data1cat == 'y3gold':
        ra_1 = tdata_1['ALPHAWIN_J2000'] 
        dec_1 = tdata_1['DELTAWIN_J2000']
    elif data1cat == 'y1gold':
        ra_1 = tdata_1['RA']
        dec_1 = tdata_1['DEC'] 
    elif data1cat == 'deep':
        ra_1 = tdata_1['ra']  
        dec_1 = tdata_1['dec'] 
    else:
        ra_1 = tdata_1['RA']
        ra_1 = tdata_1['DEC']
    if data2cat == 'deep':
        ra_2 = tdata_2['ra'] 
        dec_2 = tdata_2['dec']
    elif data2cat == 'hsc':
        ra_2 = tdata_2['ra'] 
        dec_2 = tdata_2['dec'] 
    else:
        ra_1 = tdata_1['ra']
        ra_1 = tdata_1['dec']
   
    matches = smatch.match(ra_1, dec_1, radius, ra_2, dec_2, nside=NSIDE, maxmatch=maxmatch)
    
    return matches

def plot_detection_completeness(tdata_1,tdata_2,matches,minmax,binning,binvar,field,reference,objtyp='All'):
    print('Plotting detection completeness')
    compl = np.empty(binning)
    dcompl = np.empty(binning)
    dcompl_lo = np.empty(binning)
    dcompl_hi = np.empty(binning)    
    midbins = np.empty(binning)
    interval = float(minmax[1]-minmax[0])/float(binning)
    mag_bins = np.arange(20.,28.6,0.2)
    if reference == 'deep':
        ref_class = 'bdf_T'
        general_mask = (tdata_2['mask_flags'] == 0) & (tdata_2['flags'] == 0) & (tdata_2['bdf_T'] < 30)
        general_mask_matched = (tdata_2['mask_flags'][matches['i2']] == 0) & (tdata_2['flags'][matches['i2']] == 0)  & (tdata_2['bdf_T'][matches['i2']] < 30) 
    elif reference == 'HSC' or reference == 'deepvsHSC':
        ref_class = 'i_extendedness_value'
        #ref_class = 'iclassification_extendedness'
        general_mask = (tdata_2[ref_class] > -1)
        general_mask_matched = (tdata_2[ref_class][matches['i2']] > -1)
        #general_mask = (tdata_1['EXTENDED_CLASS_MASH_SOF'] > -1)
        #general_mask_matched = (tdata_1['EXTENDED_CLASS_MASH_SOF'][matches['i1']] > -1)
    else:
        ref_class = None
    thresholds = {"deep":0.01,"Balrog":0,"HSC":0.5,"deepvsHSC":0.5}

    for i in range(binning):
        lo = minmax[0]+i*interval
        midbins[i] = lo + interval*0.5
        print(midbins[i])
        if field == 'vvds':
            mask_match = (tdata_2[binvar][matches['i2']] > lo) & (tdata_2[binvar][matches['i2']] < lo + interval) \
                & (tdata_2['ra'][matches['i2']] > 337) & general_mask_matched
            mask = (tdata_2[binvar] > lo) & (tdata_2[binvar] < lo + interval) & (tdata_2['ra'] > 337) & general_mask
        else:
            mask_match = (tdata_2[binvar][matches['i2']] > lo) & (tdata_2[binvar][matches['i2']] < lo + interval) & general_mask_matched
            mask = (tdata_2[binvar] > lo) & (tdata_2[binvar] < lo + interval) & general_mask
        if objtyp == 'Galaxies':
            detgal = sum(tdata_2[ref_class][matches['i2']][mask_match] > thresholds[reference])
            allgal = sum(tdata_2[ref_class][mask] > thresholds[reference])
        elif objtyp == 'Stars':
            detgal = sum(tdata_2[ref_class][matches['i2']][mask_match] < thresholds[reference]) 
            allgal = sum(tdata_2[ref_class][mask] < thresholds[reference])        
        else:
            detgal = sum(tdata_2[ref_class][matches['i2']][mask_match]) 
            allgal = sum(tdata_2[ref_class][mask])
        print(float(detgal),float(allgal))
        compl[i] = float(detgal)/float(allgal)
        dcompl[i] = 1/float(allgal)
        dcompl[i] = dcompl[i]*np.sqrt(float(detgal)*(1-dcompl[i])) #binomial error, temporary
        einterval = efficiencyError.efficiencyError(float(allgal),float(detgal),0.95).calculate() #2-sigma errors
        dcompl_lo[i] = einterval[0]-einterval[1]
        dcompl_hi[i] = einterval[2]-einterval[0]
        if (i > 0) and (compl[i] < 0.90) and (compl[i-1] > 0.90):
            mag90 = midbins[i-1]+(midbins[i]-midbins[i-1])*(0.90-compl[i-1])/(compl[i]-compl[i-1])
            print('Completeness at ',mag90,'is ~90% for',objtyp)
        if (i > 0) and (compl[i] < 0.95) and (compl[i-1] > 0.95):
            mag95 = midbins[i-1]+(midbins[i]-midbins[i-1])*(0.95-compl[i-1])/(compl[i]-compl[i-1])
            print('Completeness at ',mag95,'is ~95% for',objtyp)
    print(midbins)
    print(compl)        
    plt.errorbar(midbins,compl,yerr=[dcompl_lo,dcompl_hi],color='red',marker='o',label=objtyp+' '+field)
    #plt.xticks(np.arange(minmax[0], minmax[1]+1, 0.5))
    plt.xticks(np.arange(minmax[0], minmax[1]+1, 1.0))
    plt.hlines(0.90,19,27)
    plt.xlabel(binvar,fontsize=14)
    #plt.xlabel(binvar,fontsize=14)
    plt.ylabel('Completeness',fontsize=14)
    plt.ylim(0.0,1.0)
    plt.title('Completeness vs '+reference+' objects', fontsize=16)
    #plt.title('Completeness vs deep HSC objects (SN-X3)', fontsize=16)
    #plt.legend(loc='lower left',fontsize=14)
    plt.grid(True)
    #plt.savefig('completeness_galaxies_vs_'+reference+'_test.png')

def plot_eff_cont(tdata_1,tdata_2,matches,minmax,binning,binvar,field,reference,figsdir,do_binning):
    print('Computing extendedness performance in mag range ',minmax)
    ppv = np.empty(binning)
    dppv = np.empty(binning)
    dppv_lo = np.empty(binning)
    dppv_hi = np.empty(binning)
    tpr = np.empty(binning)
    dtpr = np.empty(binning)
    dtpr_lo = np.empty(binning)
    dtpr_hi = np.empty(binning)
    midbins = np.empty(binning)
    interval = float(minmax[1]-minmax[0])/float(binning)
    ref_class = 'i_extendedness_value'
    data_class = 'extended_class_mash_sof'
    #ref_class = 'mu_class_acs'
    #data_class = 'iclassification_extendedness'
    ref_idx = 'i2'
    truth_th = 0.5 #extendedness for HSC
    #truth_th = 1.5 #mu_acs for ACS
    ths = [0.5,1.5,2.5] #for Y3 GOLD
    #ths = [0.5] #for HSC
    ### note that the following procedure is constructed for galaxies
    colors = [BAND_COLORS['u'],BAND_COLORS['g'],BAND_COLORS['r']]
    for t,th in enumerate(ths):
        mask = (tdata_1[binvar][matches['i1']] > minmax[0]) & (tdata_1[binvar][matches['i1']] < minmax[1]) & (tdata_1[data_class][matches['i1']] > th)
        tp = sum(tdata_2[ref_class][matches['i2']][mask] > truth_th)
        fp = sum(tdata_2[ref_class][matches['i2']][mask] < truth_th)
        mask = (tdata_1[binvar][matches['i1']] > minmax[0]) & (tdata_1[binvar][matches['i1']] < minmax[1]) & (tdata_1[data_class][matches['i1']] < th)
        tn = sum(tdata_2[ref_class][matches['i2']][mask] > truth_th)
        print('Purity for extended sample is >',float(tp)/float(tp+fp),' for threshold >= ',th+0.5)
        print('Efficiency for extended sample is >',float(tp)/float(tp+tn),' for threshold >= ',th+0.5)
        mask = (tdata_1[binvar][matches['i1']] > minmax[0]) & (tdata_1[binvar][matches['i1']] < minmax[1]) & (tdata_1[data_class][matches['i1']] < th)
        tp = sum(tdata_2[ref_class][matches['i2']][mask] < truth_th)
        fp = sum(tdata_2[ref_class][matches['i2']][mask] > truth_th)
        mask = (tdata_1[binvar][matches['i1']] > minmax[0]) & (tdata_1[binvar][matches['i1']] < minmax[1]) & (tdata_1[data_class][matches['i1']] > th)
        tn = sum(tdata_2[ref_class][matches['i2']][mask] < truth_th)
        print('Purity for point source sample is >',float(tp)/float(tp+fp),' for threshold <= ',th-0.5)
        print('Efficiency for point source sample is >',float(tp)/float(tp+tn),' for threshold <= ',th-0.5)
    if do_binning:
        for t,th in enumerate(ths):
            for i in range(binning):
                lo = minmax[0]+i*interval
                midbins[i] = lo + interval*0.5
                mask = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval) & (tdata_1[data_class][matches['i1']] > th) # > th Y3 GOLD 
                #print(len(tdata_2[ref_class][matches['i2']]))
                tp = sum(tdata_2[ref_class][matches['i2']][mask] > truth_th) #> truth_th HSC, == 1 for ACS
                fp = sum(tdata_2[ref_class][matches['i2']][mask] < truth_th) #< truth_th HSC, == 2 for ACS
                #print(lo,'-',lo+interval,tp,fp)
                ppv[i] = float(tp)/float(tp+fp)
                dppv[i] = 1/float(tp+fp) 
                dppv[i] = dppv[i]*np.sqrt(float(tp)*(1-ppv[i])) #binomial error, temporary
                einterval = efficiencyError.efficiencyError(float(tp+fp),float(tp),0.95).calculate() # See Paterno 2004 Fermilab note
                #print interval[1],'-',interval[2]
                dppv_lo[i] = einterval[0]-einterval[1]
                dppv_hi[i] = einterval[2]-einterval[0]
                mask = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval) & (tdata_1[data_class][matches['i1']] < th) # < th Y3 GOLD
                tn = sum(tdata_2[ref_class][matches['i2']][mask] > truth_th) #> truth_th HSC, == 1 for ACS      
                tpr[i] = float(tp)/float(tp+tn)
                #mask = (tdata_1[binvar][matches['i1']] > lo) & (tdata_1[binvar][matches['i1']] < lo + interval)
                #gals =  sum(tdata_2[ref_class][matches['i2']][mask] == truth)
                #tpr[i] = float(tp)/float(gals)
                dtpr[i] = 1/float(tp+fp)
                dtpr[i] = dtpr[i]*np.sqrt(float(tp)*(1-tpr[i])) #binomial error, temporary
                einterval = efficiencyError.efficiencyError(float(tp+tn),float(tp),0.95).calculate() # See Paterno 2004 Fermilab note
                #print interval[1],'-',interval[2]
                dtpr_lo[i] = einterval[0]-einterval[1]
                dtpr_hi[i] = einterval[2]-einterval[0]
                print(midbins[i],float(tp+fp),(1.0-ppv[i])*100,dppv[i]*100,tpr[i]*100,dtpr[i]*100)
            #plt.errorbar(midbins,1.0-ppv,yerr=[dppv_lo,dppv_hi],marker='o',label='Contamination '+field.upper()+' MASH >= '+str(th+0.5))
            print(str(int(th+0.5)))
            plt.errorbar(midbins,1.0-ppv,yerr=[dppv_lo,dppv_hi],marker='.',label='Contamination  MASH $\geq$ '+str(int(th+0.5)),color=colors[t])
            plt.errorbar(midbins,tpr,yerr=[dtpr_lo,dtpr_hi],marker='.',label='Efficiency  MASH $\geq$ '+str(int(th+0.5)),color=colors[t],ls='dashed')
            #plt.errorbar(midbins,1.0-ppv,yerr=dppv,marker='o',label='Contamination') # ACS
            #plt.errorbar(midbins,tpr,yerr=dtpr,marker='+',label='Efficiency') # ACS
        #plt.xlabel('SOF i-band magnitude')
        plt.xlabel('i-band magnitude')
        plt.ylabel('Galaxy efficiency/contamination')
        plt.hlines(0.95,19,24)
        plt.ylim(0.0,1.0)
        plt.title('Efficiency/contamination vs HSC galaxy classification: '+ field.upper())
        #plt.title('Efficiency/contamination vs ACS galaxy classification: '+ field.upper(), fontsize=16)
        plt.legend(loc='center left')
        plt.savefig(figsdir+'effpur_vs_'+reference+'_'+field+'_test.png')
        print(repr(ppv))
        print(repr(dppv_lo))
        print(repr(dppv_hi))
        print(repr(dtpr_lo))
        print(repr(dtpr_hi))
    return tpr,dtpr,ppv,dppv ### this will return the values for the highest threshold th

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
    parser.add_option("--reference",dest="reference",help="Reference catalog to use (HSC/deep/Balrog/deepvsHSC)",default='deepvsHSC')
    parser.add_option("--band",dest="band",help="Reference band",default='i')
    parser.add_option("--field",dest="field",help="Field",default='snx3')
    (options, args) = parser.parse_args()

    #fields = ['sxds','deep23','vvds']
    #fields = ['cosmosud']
    #fields = ['w05']
    #fields = ['cosmos']
    #fields = ['snx3']
    fields = [options.field]
    
    workdir = '/Users/nsevilla/y3gold-paper/'
    datadir = '/Users/nsevilla/y3gold-paper/data/'
    figsdir = '/Users/nsevilla/y3gold-paper/figs/'

    print(datadir)
    
    if options.measure_effcont:
        print('Measuring efficiency/contamination')
        for field in fields:
            hdulist = fits.open(datadir+field+'_y3gold_Arctmasked_S18grizmasked.fits',memmap=True)
            #hdulist = fits.open('/Users/nsevilla/data/HSC_14119_cosmos_wide_sg.fits',memmap=True)
            tdata_1 = hdulist[1].data
            hdulist = fits.open(datadir+field+'_hsc_goldmasked_Arctmasked_S18grizmasked_cut.fits',memmap=True)
            #hdulist = fits.open('/Users/nsevilla/data/cosmos_acs_iphot_200709_topcat.fits',memmap=True)
            tdata_2 = hdulist[1].data    
            print('Matching in field',field)
            matches = match_cat(tdata_1,tdata_2,0.5/3600,'y3gold','hsc')
            print(len(matches))
            tpr,dtpr,ppv,dppv = plot_eff_cont(tdata_1,tdata_2,matches,minmax=[19,22.5],binning=10,binvar='SOF_CM_MAG_CORRECTED_I',field=field,reference=options.reference,figsdir=figsdir,do_binning=False)
            #tpr,dtpr,ppv,dppv = plot_eff_cont(tdata_1,tdata_2,matches,minmax=[19,24],binning=10,binvar='imag_kron',field=field,reference=options.reference)
            #print(tpr)
    if options.measure_completeness:
        print('Measuring detection completeness')
        for field in fields:
            if options.reference == 'HSC':
                #hdulist1 = fits.open(datadir+field+'_y3gold_y1y3goldmasked_Arctmasked_S18grizmasked.fits',memmap=True)
                #hdulist2 = fits.open(datadir+field+'_hsc_y1y3goldmasked_Arctmasked_S18grizmasked_sofconv_cut.fits',memmap=True)
                hdulist1 = fits.open(datadir+field+'_y3gold_Arctmasked_S18grizmasked.fits')
                hdulist2 = fits.open(datadir+field+'_hsc_y3goldmasked_Arctmasked_S18grizmasked_sofconv_cut.fits')
                data1cat = 'y3gold'
                data2cat = 'hsc'
            elif options.reference == 'deep':
                hdulist1 = fits.open(datadir+field+'_y3gold_completeness_dfmasked.fits',memmap=True)
                hdulist2 = fits.open(datadir+field+'_dfmasked_goldmasked_selectedcols.fits',memmap=True)
                data1cat = 'y3gold'
                data2cat = 'deep'
            elif options.reference == 'deepvsHSC':
                datadir = '/Users/nsevilla/des/deep_fields/'
                hdulist1 = fits.open(datadir+field+'_df_Arctmasked_S18grizmasked_dfmasked_bounds.fits',memmap=True)
                hdulist2 = fits.open(datadir+field+'_hsc_167138_Arctmasked_S18grizmasked_dfmasked_bounds.fits',memmap=True) #add sofconv
                data1cat = 'deep'
                data2cat = 'hsc'
            else:
                print('Data not found',options.reference)
                sys.exit()
            tdata_1 = hdulist1[1].data
            tdata_2 = hdulist2[1].data    
            print('Matching in field',field)
            matches = match_cat(tdata_1,tdata_2,1.0/3600,data1cat,data2cat)
            print(len(matches))
            if options.reference == 'HSC' or options.reference == 'deepvsHSC':
                #binvar = options.band+'_cmodel_mag' #'cmodel_mag' #'_cmodel_mag'
                binvar = 'sof_converted_'+options.band
                #binvar='sof_cm_mag_corrected_'+options.band
            elif options.reference == 'deep':
                binvar='bdf_mag_'+options.band
            else:
                binvar='sof_cm_mag_corrected_'+options.band
            print(binvar,field,options.reference)
            plot_detection_completeness(tdata_1,tdata_2,matches,minmax=[19,25],binning=10,binvar=binvar,field=field,reference=options.reference,objtyp='All')
    if options.color_terms:
        print('Calibrating color terms')
        hdulist = fits.open(datadir+'y3gold_hscreg_stars_calibration_b.fits',memmap=True)
        tdata_1 = hdulist[1].data
        hdulist = fits.open(datadir+'hsc_stars_calibration.fits',memmap=True)
        tdata_2 = hdulist[1].data    
        matches = match_cat(tdata_1,tdata_2,radius=0.5/3600,data1cat='y3gold',data2cat='hsc')
        estimate_color_terms(tdata_1,tdata_2,matches)

if __name__ == "__main__":
    main()
