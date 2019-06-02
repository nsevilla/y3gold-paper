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

nside = 4096
datadir = '/Volumes/NO NAME/' 
maptype = 'fwhm'
print('Analyzing',maptype,'files') 
footprint = fitsio.read(datadir+'y3a2_footprint_griz_1exp_v2.0.fits.gz',ext=1)['I'].ravel()
tdataSP_foot = np.empty(hp.nside2npix(nside))
tdataSP_foot.fill(hp.UNSEEN)

for i,(band,v) in enumerate(BAND_COLORS.items()):
    if band == 'u':
        continue
    if band == 'Y' and (maptype == 'sof' or maptype == 'mof'):
        continue
    print('Processing band',band)
    if maptype == 'maglim':
        filename = 'y3a2_'+band+'_o.4096_t.32768_maglim_EQU.fits.gz'
        tdataSP = fitsio.read(datadir+filename)
    elif maptype == 'sof':
        filename = 'y3a2_gold_2_2_1_sof_nside4096_nest_'+band+'_depth.fits.gz'
        tdataSP = fitsio.read(datadir+filename)['I'].ravel()
    elif maptype == 'auto':
        filename = 'y3a2_gold_2_2_1_auto_nside4096_nest_'+band+'_depth.fits.gz'
        tdataSP = fitsio.read(datadir+filename)['I'].ravel()
    elif maptype == 'fwhm':
        filename = 'y3a2_'+band+'_o.4096_t.32768_FWHM.WMEAN_EQU.fits.gz'
        tdataSP = fitsio.read(datadir+filename)
    else:
        print('Map type',maptype,'not found')
        sys.exit(1)
    
    #if len(tdataSP) == hp.nside2npix(nside):
    #    mask_all = (tdataSP != hp.UNSEEN)
        #print('Median',np.median(tdataSP[mask_all]),np.median(tdataSP_foot[mask]))
    if len(tdataSP) != hp.nside2npix(nside):
        tdataSP_foot[tdataSP['PIXEL']] = tdataSP['SIGNAL']
        tdataSP_foot = tdataSP_foot*footprint
    else:
        tdataSP_foot = tdataSP*footprint

    tdataSP_foot[tdataSP_foot>30] = hp.UNSEEN
    tdataSP_foot[tdataSP_foot<0] = hp.UNSEEN
    mask = (tdataSP_foot!=hp.UNSEEN)

    upper,lower = np.percentile(tdataSP_foot[mask],[84.075,15.825])
    med = np.median(tdataSP_foot[mask])
    print('Median',med,'+',upper-med,'-',med-lower)

    #n, bins, patches = plt.hist(tdataSP['SIGNAL'],label=band,bins=150)
    #print('Mean',np.mean(tdataSP['SIGNAL']),np.mean(tdataSP_foot[mask]))
    #print('Median',np.median(tdataSP['SIGNAL']),np.median(tdataSP_foot[mask]))
    #print('Histogram maximum',bins[np.argmax(n)])
