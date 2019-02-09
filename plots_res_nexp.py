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

def plot_nexp(filey1, filey3, band, lim_inf, lim_sup):

    # Reading Y1 data
    hdu = fits.open(filey1, memmap=True)
    pix_y1 = hdu[1].data.field('PIXEL')
    signal_y1 = hdu[1].data.field('SIGNAL')
    hdu.close()

    # Y1 data is in HP ring format
    pix_y1 = hp.ring2nest(nside, pix_y1)

    # Reading Y3 data
    hdu = fits.open(filey3, memmap=True)
    pix_y3 = hdu[1].data.field('PIXEL')
    signal_y3 = hdu[1].data.field('SIGNAL')
    hdu.close()

    PIX_Y1 = np.bincount(pix_y1, weights=signal_y1, minlength=hp.nside2npix(nside))
    PIX_Y3 = np.bincount(pix_y3, weights=signal_y3, minlength=hp.nside2npix(nside))

    # Y3 fields as same Y1 fields
    PIX_Y3[PIX_Y1 == 0] = 0.

    # Masking empty and bad pixels
    mask = np.zeros(hp.nside2npix(nside), dtype=np.bool)

    mask[PIX_Y1 == 0] = 1

    HPX4096 = PIX_Y3 - PIX_Y1

    # ATTENTION: Forcing limits
    HPX4096[HPX4096 <= lim_inf] = lim_inf
    HPX4096[HPX4096 >= lim_sup] = lim_sup

    m = hp.ma(HPX4096)
    m.mask = mask

    ra_DES, dec_DES = np.loadtxt('data/des_round.dat', usecols=(0, 1), unpack=True)

    size = 800
    hp.visufunc.mollview(map=m, flip='geo', format='%.3g', cmap=cmap, coord='C', hold=True, xsize=size,
                         nest=True, cbar=False, notext=True, margins=None)
    hp.projplot(ra_DES, dec_DES, lonlat=True, coord=['C'], zorder=10, color='b', lw=0.2)

    hp.graticule(dpar=15.,dmer=30., lw=0.1)
    plt.savefig('figs/primary.png', dpi=600, bbox_inches='tight', pad_inches=0)
    plt.close()

    d_op = mpimg.imread('figs/primary.png')
    w, h = len(d_op[0,:]), len(d_op[:,0])
    d_op = d_op[int(0.08156*h):int(0.9095*h),int(0.082*w):int(0.982*w)]
    w, h = len(d_op[0,:]), len(d_op[:,0])

    fig = plt.figure(figsize = [9., 6.4])
    ax1 = fig.add_subplot(111)
    ax1.set_facecolor('lightgray')
    # range_ = 0.2
    plt.imshow(d_op[int(0.4293*h):int(0.9497*h), int(0.3611*w):int(0.75*w)], extent=[-45., 90.,-75.,10.],
               aspect='auto', origin='upper', interpolation=None, vmin=np.min(m), vmax=np.max(m), cmap=cmap)
    plt.title(r'$\mathrm{Y3-Y1\ residual\ exposure\ time\ %s\ band}$' % band)
    plt.xlabel(r'$\mathrm{\alpha}$')
    plt.ylabel(r'$\mathrm{\delta}$')
    plt.ylim([-75,10])
    plt.xlim([90,-45])


    y_ = [-1.5, -18, -34, -50, -64, -75]
    labels_y=[r'$\mathrm{0^o}$', r'$\mathrm{-15^o}$', r'$\mathrm{-30^o}$', r'$\mathrm{-45^o}$', r'$\mathrm{-60^o}$', \
              r'$\mathrm{-75^o}$']

    x = [67, 54, 41.5, 29, 16, 3, -9, -22, -35]
    labels=[r'$\mathrm{150^o}$', r'$\mathrm{120^o}$', r'$\mathrm{90^o}$', r'$\mathrm{60^o}$', r'$\mathrm{30^o}$', \
            r'$\mathrm{0^o}$', r'$\mathrm{-30^o}$', r'$\mathrm{-60^o}$', r'$\mathrm{-90^o}$']

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=True)  # labels along the bottom edge are on
    plt.tick_params(
        axis='y',
        which='both',
        right=False,
        left=False,
        labelleft=True)

    plt.xticks(x, labels, rotation='horizontal')
    plt.yticks(y_, labels_y)
    cbaxes = fig.add_axes([0.6, 0.62, 0.275, 0.035])
    print(m.min(), m.max())
    cb = plt.colorbar(cax=cbaxes, cmap=cmap, orientation='horizontal', label=r'$\mathrm{Y3-Y1\ exposure\ time(s)}$')
    cb.ax.tick_params(labelsize=8) 
    plt.savefig('figs/HP_EQU_Y3-Y1_NEXP_res_' + band + '.png', dpi=150, bbox_inches='tight')
    plt.close()
    plt.clf()

    step = 10.
    plt.hist(m, bins=np.arange(lim_inf,lim_sup+step,step), color='b', alpha = 0.5, lw=1, histtype='step', log=True)
    n, bins, patches = plt.hist(m, bins=np.arange(lim_inf,lim_sup+step,step), color='b', alpha = 0.5, lw=1, log=True)
    # plt.legend(loc=1, numpoints=1, scatterpoints=1)
    plt.title(r'$\mathrm{Y3-Y1\ exposure\ time\ residual\ in\ %s\ band}$' % band)
    plt.xlabel(r'$\mathrm{time\ (s)}$')
    plt.xlim([lim_inf,lim_sup])
    plt.ylim([10.,1.5*np.max(n)])
    plt.grid(color='grey', linestyle='-', linewidth=0.5)
    plt.ylabel(r'$\mathrm{N\ Pixels}$')
    plt.savefig('figs/Y3-Y1_EXPTIME_res_%s.png' % band)
    plt.close()


nside = 4096

plot_nexp("data/Y1A1NEW_COADD_SPT_band_g_nside4096_oversamp4_EXPTIME__total.fits.gz", "y3a2_g_o.4096_t.32768_EXPTIME.SUM_EQU.fits.gz", 'g', -300, 500)
plot_nexp("data/Y1A1NEW_COADD_SPT_band_r_nside4096_oversamp4_EXPTIME__total.fits.gz", "y3a2_r_o.4096_t.32768_EXPTIME.SUM_EQU.fits.gz", 'r', -300, 500)
plot_nexp("data/Y1A1NEW_COADD_SPT_band_i_nside4096_oversamp4_EXPTIME__total.fits.gz", "y3a2_i_o.4096_t.32768_EXPTIME.SUM_EQU.fits.gz", 'i', -300, 500)
plot_nexp("data/Y1A1NEW_COADD_SPT_band_z_nside4096_oversamp4_EXPTIME__total.fits.gz", "y3a2_z_o.4096_t.32768_EXPTIME.SUM_EQU.fits.gz", 'z', -300, 500)
plot_nexp("data/Y1A1NEW_COADD_SPT_band_Y_nside4096_oversamp4_EXPTIME__total.fits.gz", "y3a2_Y_o.4096_t.32768_EXPTIME.SUM_EQU.fits.gz", 'Y', -300, 500)
