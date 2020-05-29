import numpy as np
import healpy as hp
import fitsio
from astropy.io import fits
from astropy.io import ascii

nside = 4096
datadir = '/Volumes/NO NAME/'

print('Reading mask...')
maskname = 'y3a2_footprint_griz_1exp_v2.0.fits.gz'
#hdu = fits.open(datadir+'y3a2_footprint_griz_1exp_v2.0.fits.gz')
#mask_y3 = hdu[1].data
mask_y3 = fitsio.read(datadir+maskname,ext=1)['I'].ravel()
mask = np.where(mask_y3>0)

print('Reading data...')
filename = 'y3a2_i_o.4096_t.32768_maglim_EQU.fits.gz'
hduSP = fits.open(datadir+filename,memmap=True)
tdataSP = hduSP[1].data
ra = data['ra']
dec = data['dec']
theta = (90.0 - dec)*np.pi/180.
phi = ra*np.pi/180.
pix = hp.ang2pix(nside,theta,phi,nest=False)

mask_y1_ext = np.zeros(12*nside*nside)
print('Filling mask_y1_ext',len(mask_y1))
for i in range(len(mask_y1)):
    mask_y1_ext[mask_y1['pixel_id'][i]] = mask_y1['fraction'][i] 

good = np.zeros(len(data),dtype=bool)
for counter,p in enumerate(pix):
    if counter%1000000 == 0:
        print('Masked',counter,'objects')
    if (mask_y1_ext[p] > 0.8):# & (data['mag_auto_i'][counter] < 23.5):
   
        good[counter] = True

print('Unmasked number of objects = ',len(data))
print('Sum of good objects = ',sum(good))

hdu_masked = fits.BinTableHDU(data=data[good])

print('Writing to file...')

hdu_masked.writeto('/pc/desdsk01/des/sg_challenge/external_catalogs/gaia_dr2_spt_2669_masked.fits',clobber=True)



