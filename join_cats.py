import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

cor_min = 0.0  # corte minimo na cor
cor_max = 1.8   # corte maximo na cor
#mag_min = 17.0  # corte minimo na magnitude
#mag_max = 24.  # corte maximo na magnitude

ra_, dec_, match_flag_, mag_auto_i_, diff_mag_ = [],[],[],[],[]

for j in range(1,3):
    print(j)
    hdu = pyfits.open("hpx4096_0000"+str('%02d' % (j))+".fits", memmap=True)
    print(hdu[1].columns)
    ra = hdu[1].data.field('ALPHAWIN_J2000')
    dec = hdu[1].data.field('DELTAWIN_J2000')
    match_flag = hdu[1].data.field('MATCH_FLAG')
    mag_auto_i = hdu[1].data.field('MAG_AUTO_I')
    diff_mag = hdu[1].data.field('Y3Y1_AUTO_RESIDUAL')
    hdu.close()

    ra_.extend(ra); dec_.extend(dec); match_flag_.extend(match_flag); mag_auto_i_.extend(mag_auto_i); diff_mag_.extend(diff_mag)
            
print(len(ra_))

col1 = pyfits.Column(name='ALPHAWIN_J2000', format='E', array=ra_)
col2 = pyfits.Column(name='DELTAWIN_J2000', format='E', array=dec_)
col3 = pyfits.Column(name='MATCH_FLAG', format='E', array=match_flag_)
col4 = pyfits.Column(name='MAG_AUTO_I', format='E', array=mag_auto_i_)
col5 = pyfits.Column(name='DIFF_MAG', format='E', array=diff_mag_)
cols = pyfits.ColDefs([col1, col2, col3, col4, col5])
tbhdu = pyfits.BinTableHDU.from_columns(cols)
tbhdu.writeto('Y3Y1_match_mag.fits')
