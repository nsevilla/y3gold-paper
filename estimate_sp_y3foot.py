#!/usr/bin/env python
'''
This script will create maps and histograms of survey properties
Author: Nacho Sevilla
'''
import os,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from descolors import BAND_COLORS
from astropy.io import fits
from astropy.io.fits.hdu.hdulist import HDUList
import fitsio
import healpy as hp
import skymap
from skymap import Skymap,McBrydeSkymap,OrthoSkymap
from skymap import SurveySkymap,SurveyMcBryde,SurveyOrtho
from skymap import DESSkymap
import warnings
warnings.filterwarnings('ignore')
from optparse import OptionParser

matplotlib.style.use('/Users/nsevilla/.matplotlib/des_dr1.mplstyle')
nside = 4096

def plot_sp_info(datadir,maptype,xlabel,insetlabel,threshold,generate_maps,histmin,histmax):
    footprint = fitsio.read(datadir+'masks/y3a2_footprint_griz_1exp_v2.0.fits.gz',ext=1)['I'].ravel()
    tdataSP_foot = np.empty(hp.nside2npix(nside))
    tdataSP_foot.fill(hp.UNSEEN)
    PIX_TO_ARCSEC = 0.263

    for i,(band,v) in enumerate(BAND_COLORS.items()):
        if band == 'u':
            continue
        if band == 'Y' and (maptype == 'sof' or maptype == 'mof'):
            continue
        print('Processing band',band)
        if maptype == 'maglim':
            filename = 'sp_maps/y3a2_'+band+'_o.4096_t.32768_maglim_EQU.fits'
            tdataSP = fitsio.read(datadir+filename) 
        elif maptype == 'sof':
            filename = 'sp_maps/y3a2_gold_2_2_1_sof_nside4096_nest_'+band+'_depth.fits'
            tdataSP = fitsio.read(datadir+filename)['I'].ravel()
        elif maptype == 'auto':
            filename = 'sp_maps/y3a2_gold_2_2_1_auto_nside4096_nest_'+band+'_depth.fits.gz'
            tdataSP = fitsio.read(datadir+filename)['I'].ravel()
        elif maptype == 'fwhm':
            filename = 'sp_maps/y3a2_'+band+'_o.4096_t.32768_FWHM.WMEAN_EQU.fits.gz'
            tdataSP = fitsio.read(datadir+filename)
        elif maptype == 'skybrite':
            filename = 'sp_maps/y3a2_'+band+'_o.4096_t.32768_SKYBRITE.WMEAN_EQU.fits.gz'
            tdataSP = fitsio.read(datadir+filename)
        elif maptype == 'airmass':
            filename = 'sp_maps/y3a2_'+band+'_o.4096_t.32768_AIRMASS.WMEAN_EQU.fits.gz'
            tdataSP = fitsio.read(datadir+filename)            
        elif maptype == 'exptime':
            filename = 'sp_maps/y3a2_'+band+'_o.4096_t.32768_EXPTIME.SUM_EQU.fits.gz'
            tdataSP = fitsio.read(datadir+filename)
        elif maptype == 'skyvar':
            filename = 'sp_maps/y3a2_'+band+'_o.4096_t.32768_SKYVAR.UNCERTAINTY_EQU.fits'
            tdataSP = fitsio.read(datadir+filename)
        elif maptype == 'sbcontrast':
            filename = 'sp_maps/Y3A1_SURFACE_BRIGHTNESS-v1.fits'
            tdataSP = fitsio.read(datadir+filename)
            if band == 'g':
                sel = (tdataSP['band'] == b'g') & (tdataSP['size'] == 10)
            elif band == 'r':
                sel = (tdataSP['band'] == b'r') & (tdataSP['size'] == 10)
            elif band == 'i':
                sel = (tdataSP['band'] == b'i') & (tdataSP['size'] == 10)
            elif band == 'z':
                sel = (tdataSP['band'] == b'z') & (tdataSP['size'] == 10)
            else:
                sel = (tdataSP['band'] == b'Y') & (tdataSP['size'] == 10)
            tdataSP = tdataSP[sel]
        else:
            print('Map type',maptype,'not found')
            sys.exit(1)

        
        if len(tdataSP) != hp.nside2npix(nside):
            tdataSP_foot[tdataSP['PIXEL']] = tdataSP['SIGNAL']
            tdataSP_foot = tdataSP_foot*footprint
        else:
            tdataSP_foot = tdataSP*footprint

        tdataSP_foot[tdataSP_foot>threshold] = hp.UNSEEN
        tdataSP_foot[tdataSP_foot<0] = hp.UNSEEN
        mask = (tdataSP_foot!=hp.UNSEEN)

        if generate_maps:
            fig = plt.figure(figsize=(12.,4.))
            gs = plt.GridSpec(1,4,wspace=0.001)

            fig.add_subplot(gs[:2])
            smap = DESSkymap()
            smap.draw_hpxmap(tdataSP['SIGNAL'],hp.nest2ring(nside,tdataSP['PIXEL']),nside)
            if maptype == 'skybrite':
                smap.draw_inset_colorbar(fontsize=10, format='%.1f', label=insetlabel)
            else:
                smap.draw_inset_colorbar(fontsize=10, label=insetlabel)
            smap.draw_des()
    
            fig.add_subplot(gs[3])
        
        upper,lower = np.percentile(tdataSP_foot[mask],[84.075,15.825])
        med = np.median(tdataSP_foot[mask])
        print('Median',med,'+',upper-med,'-',med-lower)
        if maptype == 'skyvar':
            medmag = -2.5 * np.log10(med / PIX_TO_ARCSEC**2 ) + 30
            print('Median (mag/arcsec**2)',medmag)
        
        value = tdataSP_foot[mask]
        plt.hist(value,normed=True,histtype='step',bins=100,linewidth=2,
                 color=v,label=band, range=[histmin,histmax])
        plt.xlabel(xlabel)
        plt.xticks(fontsize=12)
        plt.ylabel('PDF')
        plt.tight_layout()
        plt.legend(loc='upper right')
        if generate_maps:
            #plt.savefig('/Users/nsevilla/y3gold-paper/figs/y3gold_'+maptype+'_'+band+'_map.png')
            plt.savefig('/Users/nsevilla/y3gold-paper/figs/y3gold_'+maptype+'_'+band+'_map.pdf')
 
    if generate_maps == False:
        #plt.savefig('/Users/nsevilla/y3gold-paper/figs/y3gold_'+maptype+'_hist.png')
        plt.savefig('/Users/nsevilla/y3gold-paper/figs/y3gold_'+maptype+'_hist.pdf')

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d","--datadir",type="string",dest="datadir",help="Directory containing DES info",default="/Users/nsevilla/des/")
    parser.add_option("-m","--maptype",type="string",dest="maptype",help="Type of SP maps",default="skybrite")
    parser.add_option("-x","--xlabel",type="string",dest="xlabel",help="Label for histograms",default="Sky Brightness (electrons/pixel)")
    parser.add_option("-i","--insetlabel",type="string",dest="insetlabel",help="Label for inset color bar (use '$' to include carriage return)",default="Sky Brightness \n (electrons/pixel)")
    parser.add_option("-t","--threshold",type="int",dest="threshold",help="Threshold",default=12000)
    parser.add_option("--histmin",type="float",dest="histmin",help="Histogram minimum",default=0)
    parser.add_option("--histmax",type="float",dest="histmax",help="Histogram maximum",default=12500)
    parser.add_option("--generate_maps",action="store_true",dest="generate_maps",help="Toggle map generation",default=False)
    (options, args) = parser.parse_args()

    print('Analyzing',options.maptype,'files')
    print('using data in',options.datadir)
    print('generate maps',options.generate_maps)

    plot_sp_info(options.datadir,options.maptype,options.xlabel,options.insetlabel,options.threshold,options.generate_maps,options.histmin,options.histmax)

if __name__ == "__main__":
    main()
