import healpy as hp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import matplotlib.image as mpimg
from matplotlib import rc
plt.rc('font',**{'family':'serif','serif':['Helvetica'], 'size': 12})
mpl.rcParams['legend.numpoints'] = 1
plt.rc('text', usetex=True)
cmap = mpl.cm.get_cmap("inferno_r")
cmap.set_under('darkgray')
cmap.set_bad('darkgray')

'''
hdu = fits.open("Y3Y1_match_mag.fits", memmap=True)
ra = hdu[1].data.field('ALPHAWIN_J2000')
dec = hdu[1].data.field('DELTAWIN_J2000')
angdist = hdu[1].data.field('ANGDIST')
match_flag = hdu[1].data.field('MATCH_FLAG')
mag_auto_i = hdu[1].data.field('MAG_AUTO_I')
diff_mag = hdu[1].data.field('DIFF_MAG')
hdu.close()
'''
hdu = fits.open("NSEVILLAY1Y3CORR.fits", memmap=True)
match_flag = hdu[1].data.field('MATCH_FLAG')
mag_auto_iy1 = hdu[1].data.field('MAG_AUTO_I_Y1_CORRECTED')
mag_auto_i = hdu[1].data.field('MAG_AUTO_I_Y3_CORRECTED')
hdu.close()

diff_mag = mag_auto_i - mag_auto_iy1

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist2d(mag_auto_i, diff_mag, bins=200, range=[[15.5, 25],[-0.8, 1.2]], cmap='inferno_r')
plt.title(r'$\mathrm{Y3-Y1\ mag\_auto(i)}$')
plt.xlabel(r'$\mathrm{mag\_auto\ i\ [Y3]}$')
plt.ylabel(r'$\mathrm{mag\_auto\ i\ [Y3-Y1]}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([15.5,25])
plt.ylim([-0.8,1.2])
cbaxes = fig.add_axes([0.908, 0.11, 0.02, 0.773])
cb = plt.colorbar(cax=cbaxes, cmap='inferno_r', orientation='vertical')
plt.savefig('Y3-Y1corr_diff_mag_hist2d_total_range.png')
plt.close()
plt.clf()
