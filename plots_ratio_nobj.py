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
cmap = mpl.cm.get_cmap("gist_rainbow")
cmap.set_under('darkgray')
cmap.set_bad('darkgray')

def plot_nobj(m, label, filename, labelcb):

    nside = hp.npix2nside(len(m))

    ra_DES, dec_DES = np.loadtxt('data/des_round.dat', usecols=(0, 1), unpack=True)

    size = 800
    hp.visufunc.mollview(map=m, flip='geo', format='%.3g', cmap=cmap, coord='C', hold=True, xsize=size,
                         nest=True, cbar=False, notext=True, margins=None)
    hp.projplot(ra_DES, dec_DES, lonlat=True, coord=['C'], zorder=10, color='b', lw=0.2)

    hp.graticule(dpar=15.,dmer=30., lw=0.1)
    plt.savefig('primary.png', dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()

    d_op = mpimg.imread('primary.png')
    w, h = len(d_op[0,:]), len(d_op[:,0])
    d_op = d_op[int(0.08156*h):int(0.9095*h),int(0.082*w):int(0.982*w)]
    w, h = len(d_op[0,:]), len(d_op[:,0])

    fig = plt.figure(figsize = [9., 6.4])
    ax1 = fig.add_subplot(111)
    ax1.set_facecolor('lightgray')
    # range_ = 0.2
    plt.imshow(d_op[int(0.4293*h):int(0.9497*h), int(0.3611*w):int(0.75*w)], extent=[-45., 90.,-75.,10.],
               aspect='auto', origin='upper', interpolation=None, vmin=np.min(m),
               vmax=np.max(m), cmap=cmap)
    plt.title(r'$\mathrm{%s}$' % label)
    plt.xlabel(r'$\mathrm{\alpha}$')
    plt.ylabel(r'$\mathrm{\delta}$')
    plt.ylim([-75,10])
    plt.xlim([90,-45])


    y_ = [-1.5, -18, -34, -50, -64, -75]
    labels_y=[r'$\mathrm{0^o}$', r'$\mathrm{-15^o}$', r'$\mathrm{-30^o}$', r'$\mathrm{-45^o}$', r'$\mathrm{-60^o}$', r'$\mathrm{-75^o}$']

    x = [67, 54, 41.5, 29, 16, 3, -9, -22, -35]
    labels=[r'$\mathrm{150^o}$', r'$\mathrm{120^o}$', r'$\mathrm{90^o}$', r'$\mathrm{60^o}$', r'$\mathrm{30^o}$', r'$\mathrm{0^o}$', r'$\mathrm{-30^o}$', r'$\mathrm{-60^o}$', r'$\mathrm{-90^o}$']

    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=True)
    plt.tick_params(
        axis='y',
        which='both',
        right=False,
        left=False,
        labelleft=True)

    plt.xticks(x, labels, rotation='horizontal')
    plt.yticks(y_, labels_y)
    cbaxes = fig.add_axes([0.6, 0.62, 0.275, 0.035])
    cb = plt.colorbar(cax=cbaxes, cmap=cmap, orientation='horizontal', label=r'$\mathrm{%s}$' % labelcb)
    cb.ax.tick_params(labelsize=8) 
    plt.savefig('figs/HP_EQU_' + filename + '.png', dpi=150, bbox_inches='tight')
    plt.close()
    plt.clf()


nside1 = 4096
nside2 = 16384

hdu = fits.open("data/Y1_HPIX_4096_16384.fits", memmap=True)
HPIX_4096_Y1 = hdu[1].data.field('HPIX_4096')
HPIX_16384_Y1 = hdu[1].data.field('HPIX_16384')
COUNT_Y1 = hdu[1].data.field('COUNT(COADD_OBJECTS_ID)')
hdu.close()

HPX4096_Y1 = np.bincount(HPIX_4096_Y1, weights=COUNT_Y1, minlength=hp.nside2npix(nside1))
HPX16384_Y1 = np.bincount(HPIX_16384_Y1, weights=COUNT_Y1, minlength=hp.nside2npix(nside2))

hdu2 = fits.open("data/Y3_HPIX_4096_16384.fits", memmap=True)
HPIX_4096_Y3 = hdu2[1].data.field('HPIX_4096')
hdu2.close()

HPX4096_Y3 = np.bincount(HPIX_4096_Y3, minlength=hp.nside2npix(nside1))

plt.hist(HPX4096_Y1, bins=np.arange(1,70,1), color='seagreen', alpha = 0.5, lw=1, density=True, stacked=True,
                                     log=False, label = r'$\mathrm{Y1\ Gold}$')
plt.hist(HPX4096_Y3, bins=np.arange(1,70,1), color='royalblue', alpha = 0.5, lw=1, density=True, stacked=True,
                                     log=False, label = r'$\mathrm{Y3\ Gold}$')
n_y1, bins_y1, patches_y1 = plt.hist(HPX4096_Y1, bins=np.arange(1,70,1), color='seagreen', histtype='step', density=True, stacked=True,
                                     log=False)
n_y3, bins_y3, patches_y3 = plt.hist(HPX4096_Y3, bins=np.arange(1,70,1), color='royalblue', histtype='step', density=True, stacked=True,
                                     log=False)
plt.legend(loc=1, numpoints=1, scatterpoints=1)
plt.xlabel(r'$\mathrm{Objects\ per\ Pixel}$')
plt.title(r'$\mathrm{Objects\ per\ Pixel\ (Nside=4096)}$')
plt.xlim([0.,50.])
plt.ylabel(r'$\mathrm{Fraction\ of\ pixels}$')
plt.savefig('figs/histY1Y3_4096.png')
plt.close()

hdu3 = fits.open("data/Y3_HPIX_16384_bins_not_equal_zero.fits", memmap=True)
HPX16384_Y3 = hdu3[1].data.field('HPIX_16384')
hdu3.close()

plt.hist(HPX16384_Y1, bins=np.arange(1,10,1), color='seagreen', alpha = 0.5, lw=1, density=True, stacked=True,
                                     log=False, label = r'$\mathrm{Y1\ Gold}$')
plt.hist(HPX16384_Y3, bins=np.arange(1,10,1), color='royalblue', alpha = 0.5, lw=1, density=True, stacked=True,
                                     log=False, label = r'$\mathrm{Y3\ Gold}$')
n2_y1, bins2_y1, patches2_y1 = plt.hist(HPX16384_Y1, bins=np.arange(1,10,1), color='royalblue', histtype='step', density=True, stacked=True, log=False)

n2_y3, bins2_y3, patches2_y3 = plt.hist(HPX16384_Y3, bins=np.arange(1,10,1), color='seagreen', histtype='step', density=True, stacked=True, log=False)

plt.legend(loc=1, numpoints=1, scatterpoints=1)
plt.title(r'$\mathrm{Objects\ per\ Pixel\ (Nside=16384)}$')
plt.xlabel(r'$\mathrm{Objects\ per\ Pixel}$')
plt.xlim([0.,10.])
plt.ylim([0.,0.7])
plt.ylabel(r'$\mathrm{Fraction\ of\ Pixels}$')
plt.savefig('figs/histY1Y3_16384.png')
plt.close()
exit()
mask = np.zeros(hp.nside2npix(nside1), dtype=np.bool)

HPX4096 = np.zeros(hp.nside2npix(nside1))

for i in range(len(HPX4096_Y1)):
    if (HPX4096_Y1[i] == 0):
        mask[i] = 1
    else:
        HPX4096[i] = HPX4096_Y3[i] / HPX4096_Y1[i]

# Force limits if you want
# HPX4096[HPX4096 >= 10.] = 10.

m = hp.ma(HPX4096)
m.mask = mask

plot_nobj(m, 'Y3/Y1\ Nobj\ ratio [15<i<25],\ NSIDE=4096', 'y3y1ratio', 'Y3/Y1\ obj\ counts\ [magauto(i)]')
