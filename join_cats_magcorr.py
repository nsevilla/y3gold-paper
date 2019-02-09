import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

MATCHFLAG_, MAG_AUTO_I_Y1_CORRECTED_, MAG_AUTO_I_Y3_CORRECTED_ = [], [], []

for j in range(1,6):
    hdu = pyfits.open("NSEVILLAY3Y1MATCH_INNERJOIN_0000"+str('%02d' % (j))+".fits", memmap=True)
    MATCHFLAG = hdu[1].data.field('MATCH_FLAG')
    MAG_AUTO_I_Y1_CORRECTED = hdu[1].data.field('MAG_AUTO_I_Y1_CORRECTED')
    MAG_AUTO_I_Y3_CORRECTED = hdu[1].data.field('MAG_AUTO_I_Y3_CORRECTED')
    hdu.close()
    MATCHFLAG_.extend(MATCHFLAG)
    MAG_AUTO_I_Y1_CORRECTED_.extend(MAG_AUTO_I_Y1_CORRECTED)
    MAG_AUTO_I_Y3_CORRECTED_.extend(MAG_AUTO_I_Y3_CORRECTED)

col1 = pyfits.Column(name='MATCH_FLAG', format='K', array=MATCHFLAG_)
col2 = pyfits.Column(name='MAG_AUTO_I_Y1_CORRECTED', format='E', array=MAG_AUTO_I_Y1_CORRECTED_)
col3 = pyfits.Column(name='MAG_AUTO_I_Y3_CORRECTED', format='D', array=MAG_AUTO_I_Y3_CORRECTED_)
cols = pyfits.ColDefs([col1, col2, col3])
tbhdu = pyfits.BinTableHDU.from_columns(cols)
tbhdu.writeto('NSEVILLAY1Y3CORR.fits')
