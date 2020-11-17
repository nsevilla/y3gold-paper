#!/usr/bin/env python
'''
This script will create a bad region map file using the Y3 Gold footprint 
Author: Nacho Sevilla (based Eli Rykoff's and Alex Drlica-Wagner's code)
'''
import healpy as hp
import numpy as np
from astropy.io import fits
from astropy import coordinates as coo
import sys
import os
import glob
import matplotlib.pyplot as plt
from optparse import OptionParser
from collections import OrderedDict as odict
import copy
import re
from y3a2_mask_table import make_mask_table
import y3a2_mask_bits
import fitsio

def make_badmask_file(workdir,nside,bits):

    print 'Reading footprint'    
    footprint = fitsio.read(workdir+'y3a2_footprint_griz_1exp_v2.0.fits.gz',ext=1)['I'].ravel()

    objfiles = glob.glob(workdir+'cat4badregions_2_0_0*.fits')

    filelist = []
    badmask = np.zeros(hp.nside2npix(nside),dtype=np.int32)

    for b,d in bits.items():
        if d['name'] == 'Total area':
            continue
        filelist.append(workdir+d['filename'])
        
    if all([os.path.isfile(f) for f in filelist]):
        print "Found previous bad region files",filelist
        for b,d in bits.items():
            if d['name'] == 'Total area':
                continue
            hdu = fits.open(workdir+d['filename'])
            d['dat'] = hdu[1].data
    else:
        pix = None
        for f in objfiles:
            print 'Reading',f
            cat_chunk = fitsio.read(f)
            pix_chunk = hp.ang2pix(nside,cat_chunk['ALPHAWIN_J2000'],cat_chunk['DELTAWIN_J2000'],lonlat=True,nest=True)
            print 'Appending...'
            if pix is None:
                cat = cat_chunk.copy()
                pix = pix_chunk.copy()
            else:
                cat = np.append(cat,cat_chunk)
                pix = np.append(pix,pix_chunk)

        unique,unique_counts = np.unique(pix, return_counts=True)
        print 'Average number of objects in pixel',np.mean(unique_counts)

        #infoot = np.zeros(len(cat['RA']),dtype=bool)
        #pix = hp.ang2pix(nside,cat['RA'],cat['DEC'],lonlat=True,nest=True)    
        #for counter,p in enumerate(pix):
        #    if counter%1000000 == 0:
        #        print 'Analyzed',counter,'objects for footprint masking'
        #    if (footprint[p] > 0): 
        #        infoot[counter] = True
        #cat=cat[infoot]
        #pixgd = hp.ang2pix(nside,cat['RA'],cat['DEC'],nest=True,lonlat=True)

        # useful quantities
        #positions_g = coo.SkyCoord(cat['RA_G'],cat['DEC_G'],unit="deg")
        #positions_i = coo.SkyCoord(cat['RA_I'],cat['DEC_I'],unit="deg")
        #color_separation = positions_g.separation(positions_i)
        #psfres = cat['WAVG_MAG_PSF_I']-cat['MAG_PSF_I'] ###CHANGE TO THIS WHEN MAG_PSF_I BECOMES AVAILABLE
        catlen = len(cat['ALPHAWIN_J2000'])
        print 'Total nb of objects',catlen

        for b,d in bits.items():
            mask_name = d['name']
            file_name = d['filename']
            if mask_name == 'Total area':
                continue
            print 'Analyzing',mask_name
            #if re.search('crazy',mask_name, re.IGNORECASE):
                #mask,=np.where(((cat['GMR'] < -1) | (cat['GMR'] > 4.0) | (cat['IMZ'] < -1) | (cat['IMZ'] > 4.0) | (cat['RMI'] < -1) | (cat['RMI'] > 4.0)) & (abs(cat['GMR'])<20) & (abs(cat['RMI'])<20) &  (abs(cat['IMZ'])<20))                
                #mask,=np.where((((cat['GMR'] < -1) & (cat['MAG_AUTO_I']>20)) | ((cat['GMR'] > 4.0) & (cat['MAG_AUTO_I']>20)) | (cat['IMZ'] < -1) | (cat['IMZ'] > 4.0)) & (abs(cat['GMR'])<20 & (abs(cat['IMZ'])<20))
            if re.search('astrometric',mask_name, re.IGNORECASE):
                mask,=np.where((color_separation > d['tolerance']) & (cat['MAGERR_AUTO_G'] < d['cut']))
            if re.search('psf',mask_name, re.IGNORECASE):
                mask,= np.where((cat['PSFRES'] > d['tolerance']) & (abs(cat['PSFRES']) < 3))
                #mask,= np.where(cat['MAGERR_AUTO_I'] < d['cut'])
            #if re.search('shallow',mask_name, re.IGNORECASE):
                #    mask,= np.where((footprint[pixgd]==1)&(footprint2exp[pixgd]<0))
            if re.search('mof',mask_name,re.IGNORECASE) or re.search('crazy',mask_name,re.IGNORECASE):
                #tilemask = hp.read_map(workdir+'y3v02+a_mof001_badtiles_nside4096_nest.fits',nest=True)
                tilemask = hp.read_map(workdir+file_name,nest=True)
                good = np.zeros(catlen,dtype=bool)
                pix = hp.ang2pix(4096,cat['ALPHAWIN_J2000'],cat['DELTAWIN_J2000'],lonlat=True,nest=True)
                for counter,p in enumerate(pix):
                    if counter%1000000 == 0:
                        print 'Processed',counter,'objects'
                    if (tilemask[p] > 0): # & (data['mag_auto_i'][counter] < 23.5):
                        #print 'FOUND'
                        good[counter] = True
                mask, = np.where(good == True)
                print 'Mask length',len(mask)

            print len(cat[mask]),workdir+d['filename']
            fitsio.write(workdir+d['filename'],cat[mask],clobber=True)
            d['dat'] = cat[mask]

    # make the mask (integer for now; all 0s)
    for b,d in bits.items():
        if d['name'] == 'Total area':
            continue
        badmask = add_mask(workdir,nside,badmask,d['dat'],b,d['threshold'],d['name'])

    print 'Increasing mask resolution...'
    badmask_hires_tmp = np.zeros(hp.nside2npix(4096),dtype=np.int32)
    badmask_hires = np.zeros(hp.nside2npix(4096),dtype=np.int32)
    # only copy over pixels that may be in the mask
    #gpix,=np.where(footprint > 0)
    gpix,=np.where(footprint < 1E40)
    theta,phi = hp.pix2ang(4096,gpix,nest=True)
    ip_mask = hp.ang2pix(nside,theta,phi,nest=True)
    # and set the high res pixels
    badmask_hires_tmp[gpix] = badmask[ip_mask]
    # deal with the MOF avoidance (temporary fix)
    print 'temporary fix for MOF tile mask'
    mofmask = np.zeros(hp.nside2npix(4096),dtype=np.int32)
    tilemask = fitsio.read(workdir+'y3v02+a_mof001_badtiles_nside4096_nest.fits',ext=1)['PIXEL']
    for p in tilemask:
        mofmask[p] = 2.
    badmask_hires[gpix] = badmask_hires_tmp[gpix] + mofmask[gpix]

    print 'Writing...'
    write_badregions(badmask_hires,footprint,workdir)
 
def write_badregions(badmask,footprint,workdir):

    #badmask_float = np.zeros(badmask.size,dtype=np.float32) + hp.UNSEEN
    #maskzero, = np.where(footprint>=0)
    #badmask_float[maskzero] = 0.0 
    #bad, = np.where((footprint>0)&(badmask>0))
    #badmask_float[bad] = badmask[bad].astype(np.float32)
    badmask_float = badmask.astype(np.float32)
    outfilename = '%sy3a2_badregions_mask_v2.1.fits' % (workdir)
    hp.write_map(outfilename,badmask_float,nest=True)

def add_mask(workdir,nside,badmask,cat,code,th,name):
    print "Working on",name,"mask..."
    if re.search('crazy',name, re.IGNORECASE):
        pix = hp.ang2pix(nside,cat['RA'],cat['DEC'],nest=True,lonlat=True)
    else:
        pix = hp.ang2pix(nside,cat['ALPHAWIN_J2000'],cat['DELTAWIN_J2000'],nest=True,lonlat=True)
    pixarea = hp.nside2pixarea(nside,degrees=True)
    if re.search('shallow',name, re.IGNORECASE) or re.search('crazy',name, re.IGNORECASE):
        badpix = pix
    else:
        unique,unique_counts = np.unique(pix, return_counts=True)
        #hist=plt.hist(unique_counts,range=[0,300],bins=300)#bins=np.max(unique_counts),)
        # make the figure
        #fig=plt.figure(1)

        #ax=fig.add_subplot(111)
        #ax.set_xlabel("Number of "+name+" objects")
        #ax.set_ylabel("Number of pixels")
        #ax.set_yscale('log')

        #fig.savefig('%s%s_histogram.png' % (workdir,name))
        #fig.clf()
        mask,=np.where(unique_counts >= th) 
        badpix = pix[mask]

    badmask[badpix] = badmask[badpix] | code        
    print '...done'

    #if re.search('crazy',name, re.IGNORECASE):
    #    for thtest in range(2,40):
    #        mask,=np.where(unique_counts >= thtest)
    #        print thtest,len(pix[mask])*pixarea,float(unique_counts[unique_counts >= thtest].sum())/float(unique_counts.sum())

    return badmask 

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)

    # Parse command line
    (options, args) = parser.parse_args()

    workdir = '/data1/des/y3_validation/'
    nside = 512
    y3a2_mask_bits.init()
    bits = copy.deepcopy(y3a2_mask_bits.BAD_BITS)
   
    print 'Using nside = ',nside
    make_badmask_file(workdir,nside,bits)
    #make_mask_table(workdir,workdir+'y3a2_footprint_griz_1exp.fits',workdir+'y3a2_badregions_mask.fits.gz',workdir+'y3a1_griz_o.4096_t.32768_coverfoot_EQU.fits',nside,'bad')

if __name__ == '__main__':
    sys.exit(main())
