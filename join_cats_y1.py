import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

HPX4096_, HPX16384_, counts_ = [], [], []

for j in range(1,3):
    hdu = pyfits.open("hpx4096_y1_0000"+str('%02d' % (j))+".fits", memmap=True)
    print(hdu[1].columns)
    exit()
    HPX4096 = hdu[1].data.field('HPIX_4096')
    HPX16384 = hdu[1].data.field('HPIX_16384')
    counts = hdu[1].data.field('COUNT(COADD_OBJECTS_ID)')
    hdu.close()
    HPX4096_.extend(HPX4096); counts_.extend(counts); HPX16384_.extend(HPX16384)

col1 = pyfits.Column(name='HPIX_4096', format='J', array=HPX4096_)
col2 = pyfits.Column(name='HPIX_16384', format='K', array=HPX16384_)
col3 = pyfits.Column(name='COUNT(COADD_OBJECTS_ID)', format='I', array=counts_)
cols = pyfits.ColDefs([col1, col2, col3])
tbhdu = pyfits.BinTableHDU.from_columns(cols)
tbhdu.writeto('Y1_HPIX_4096_16384.fits')
