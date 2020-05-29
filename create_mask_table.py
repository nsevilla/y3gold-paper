import healpy as hp
import numpy as np
import fitsio
import sys
import os
from optparse import OptionParser
import mask_bits as mbits

def make_mask_table(workdir,footprint_filename,mask_filename,fracdet_filename,nside,identifier):
    mb = mbits.mask_bits() 
    print('Making mask table')
    if identifier == 'bad':
        bt = mb.BAD_BITS
    elif identifier == 'foreground':
        bt = mb.FOREGROUND_BITS
    else:
        print('Identifier',identifier,'not found when making map')
        sys.exit()
    print('Reading footprint')
    footprint = fitsio.read(footprint_filename,ext=1)['I'].ravel()
    print('Reading mask')
    regions_map = fitsio.read(mask_filename,ext=1)['I'].ravel()
    print('Reading detection fraction')
    fracdet = hp.read_map(fracdet_filename,nest=True)

    pixarea = hp.nside2pixarea(nside,degrees=True)
    griz = (footprint >= 1)
    for b in bt.keys():
        mask, = np.where((np.bitwise_and(regions_map.astype(int),int(b))>0) & (footprint>=1))
        print(np.mean(fracdet[mask]))
        selection = np.zeros(hp.nside2npix(4096))
        selection[mask] = 1
        print('| %4s | %8.2f | %8.2f |'%(b,selection[mask].sum()*pixarea,(selection[mask]*fracdet[mask]).sum()*pixarea))      
    
def main():
    '''
    Run code with options
    '''
    nside = 4096
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d","--workdir",dest="syspath",help="Directory to read/store maps",default='/Users/nsevilla/des/masks/')
    parser.add_option("--maptype", dest="maptype", default="foreground", help="Map type")
    parser.add_option("--footprint_filename",dest="footprint_filename",help="Footprint filename (input, fits format)",default='y3a2_footprint_griz_1exp_v2.0.fits.gz')
    parser.add_option("--mask_filename",dest="mask_filename",help="Foreground mask filename (output, fits format)",default='y3a2_foreground_mask_v2.1.fits')
    
    # Parse command line
    (options, args) = parser.parse_args()

    workdir = options.syspath
    footprint_filename = workdir + options.footprint_filename
    mask_filename = workdir + options.mask_filename
    fracdet_filename = workdir + 'y3a2_griz_o.4096_t.32768_coverfoot_EQU.fits'

    make_mask_table(workdir,footprint_filename,mask_filename,fracdet_filename,nside,options.maptype)

    
if __name__ == '__main__':
    sys.exit(main())

