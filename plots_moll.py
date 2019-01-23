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
cmap = mpl.cm.get_cmap("seismic")
cmap.set_under('darkgray')
cmap.set_bad('darkgray')


def plot_moll_cut(m, label, filename, labelcb, nside):

    ra_DES, dec_DES = np.loadtxt('des_round.dat', usecols=(0, 1), unpack=True)

    size = 800
    hp.visufunc.mollview(map=m, flip='geo', format='%.3g', cmap=cmap, coord='C', hold=True, xsize=size,
                         nest=True, cbar=False, notext=True, margins=None)
    hp.projplot(ra_DES, dec_DES, lonlat=True, coord=['C'], zorder=10, color='b', lw=0.2)

    hp.graticule(dpar=15.,dmer=30., lw=0.1)
    plt.savefig('primary.png', dpi=600, bbox_inches='tight', pad_inches=0)
    plt.close()

    d_op = mpimg.imread('primary.png')
    w, h = len(d_op[0,:]), len(d_op[:,0])
    d_op = d_op[int(0.08156*h):int(0.9095*h),int(0.082*w):int(0.982*w)]
    w, h = len(d_op[0,:]), len(d_op[:,0])

    fig = plt.figure(figsize = [9., 6.4])
    ax1 = fig.add_subplot(111)
    ax1.set_facecolor('lightgray')
    range_ = 0.2
    plt.imshow(d_op[int(0.4293*h):int(0.9497*h), int(0.3611*w):int(0.75*w)], extent=[-45., 90.,-75.,10.], aspect='auto', origin='upper', interpolation=None, vmin=-1.*range_,
               vmax=range_, cmap=cmap)
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
    print(m.min(), m.max())
    cb = plt.colorbar(cax=cbaxes, cmap=cmap, orientation='horizontal', label=r'$\mathrm{%s}$' % labelcb)
    cb.ax.tick_params(labelsize=8) 
    plt.savefig('HP_EQU_' + filename + '.png', dpi=300, bbox_inches='tight')
    plt.close()
    plt.clf()


hdu = fits.open("Y3Y1_match_mag.fits", memmap=True)
ra = hdu[1].data.field('ALPHAWIN_J2000')
dec = hdu[1].data.field('DELTAWIN_J2000')
angdist = hdu[1].data.field('ANGDIST')
match_flag = hdu[1].data.field('MATCH_FLAG')
mag_auto_i = hdu[1].data.field('MAG_AUTO_I')
diff_mag = hdu[1].data.field('DIFF_MAG')
hdu.close()

ra[ra >=180] -= 360.
fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist2d(ra[(angdist <=0.5)], dec[(angdist <=0.5)], bins=[400,400], cmap='inferno_r')
plt.title(r'$\mathrm{Angular\ separation\ <0.5\ arcsec}$')
plt.xlim([100,-62])
plt.xlabel(r'$\mathrm{RA}$')
plt.ylabel(r'$\mathrm{DEC}$')
cbaxes = fig.add_axes([0.908, 0.11, 0.02, 0.773])
cb = plt.colorbar(cax=cbaxes, cmap='inferno_r', orientation='vertical')
plt.savefig('0.5_arcsec.png')
plt.close()

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist2d(ra[(angdist >=1.5)&(angdist <=3.)], dec[(angdist >=1.5)&(angdist <=3.)], bins=[400,400], cmap='inferno_r')
plt.title(r'$\mathrm{Angular\ separation\ 1.5<r<3\ arcsec}$')
plt.xlabel(r'$\mathrm{RA}$')
plt.ylabel(r'$\mathrm{DEC}$')
plt.xlim([100,-62])
cbaxes = fig.add_axes([0.908, 0.11, 0.02, 0.773])
cb = plt.colorbar(cax=cbaxes, cmap='inferno_r', orientation='vertical')
plt.savefig('1.5_3_arcsec.png')
plt.close()

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist2d(ra[(angdist >=3.5)&(angdist <=5.)], dec[(angdist >=3.5)&(angdist <=5.)], bins=[400,400], cmap='inferno_r')
plt.title(r'$\mathrm{Angular\ separation\ 3.5<r<5\ arcsec}$')
plt.xlabel(r'$\mathrm{RA}$')
plt.ylabel(r'$\mathrm{DEC}$')
plt.xlim([100,-62])
cbaxes = fig.add_axes([0.908, 0.11, 0.02, 0.773])
cb = plt.colorbar(cax=cbaxes, cmap='inferno_r', orientation='vertical')
plt.savefig('3.5_5_arcsec.png')
plt.close()

exit()

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist2d(mag_auto_i, diff_mag, bins=200, range=[[15.5, 25],[-0.8, 1.2]], cmap='inferno_r')
plt.title(r'$\mathrm{Y3-Y1\ mag\_auto(i)}$')
plt.xlabel(r'$\mathrm{mag\_auto\ (i)}$')
plt.ylabel(r'$\mathrm{mag\_auto\ (i)\ [Y3-Y1]}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([15.5,25])
plt.ylim([-0.8,1.2])
cbaxes = fig.add_axes([0.908, 0.11, 0.02, 0.773])
cb = plt.colorbar(cax=cbaxes, cmap='inferno_r', orientation='vertical')
plt.savefig('Y3-Y1_diff_mag_hist2d_total_range.png')
plt.close()
plt.clf()

from scipy.stats import norm
mu, std = norm.fit(diff_mag[(diff_mag > -0.8)&(diff_mag <= 1.2)])

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist(diff_mag, bins=np.arange(-0.8, 1.2, 0.01), color='r', histtype='step', log=False)
plt.title(r'$\mathrm{Y3-Y1\ mag\ auto(i),\ \mu = %s,\ \sigma =%s}$' % (np.round(mu,4), np.round(std,4)))
plt.ylabel(r'$\mathrm{N\_objects}$')
plt.xlabel(r'$\mathrm{mag\_auto\ (i)\ [Y3-Y1]}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([-0.8,1.2])
plt.savefig('Y3-Y1_diff_mag_hist_total_range.png')
plt.close()
plt.clf()

plt.hist(angdist, bins=np.arange(0,10,0.01), color='r', alpha = 0.5, lw=1, histtype='step', log=True)
n, bins, patches = plt.hist(angdist, bins=np.arange(0.,10.,0.01), color='r', alpha = 0.5, lw=1, log=True)
plt.title(r'$\mathrm{Y3-Y1\ separation\ arcsec}$')
plt.xlabel(r'$\mathrm{arcsec}$')
plt.xlim([0,10])
plt.ylim([1e2,1.2*np.max(n)])
plt.grid(color='grey', linestyle='-', linewidth=0.5)
plt.ylabel(r'$\mathrm{N\ objects}$')
plt.savefig('Y3-Y1_angsep_hist.png')
plt.close()
plt.clf()
print(len(match_flag), match_flag)
n, bins, patches = plt.hist(match_flag, bins=np.arange(0.,4.,0.2), color='r', alpha = 0.5, lw=1, log=True)
print(n, patches)
plt.title(r'$\mathrm{Match\ flag\ histogram}$')
plt.xlabel(r'$\mathrm{flag}$')
plt.xlim([0,3.2])
plt.ylim([1e5,1.2e8])
plt.grid(color='grey', linestyle='-', linewidth=0.5)
plt.ylabel(r'$\mathrm{N\ objects}$')
plt.savefig('Y3-Y1_matchflag.png')
plt.close()
plt.clf()

filter_mag = (mag_auto_i >= 17.)&(mag_auto_i <=22.)&(match_flag == 0)
ra, dec, mag_auto_i, diff_mag, angdist = ra[filter_mag], dec[filter_mag], mag_auto_i[filter_mag], diff_mag[filter_mag], angdist[filter_mag]

nside = 512
ipix = hp.ang2pix(nside, (90-dec)/180*np.pi, ra/180*np.pi, nest=True)
HPX_diff_mag_sum = np.bincount(ipix, weights=diff_mag, minlength=hp.nside2npix(nside))

HPX_num = np.bincount(ipix, minlength=hp.nside2npix(nside))

HPX_num = HPX_num.astype(float)

mask = np.zeros(hp.nside2npix(nside), dtype=np.bool)

HPX = np.zeros(hp.nside2npix(nside))

for i in range(len(HPX_num)):
    if (HPX_num[i] == 0):
        mask[i] = 1
    else:
        HPX[i] = HPX_diff_mag_sum[i] / float(HPX_num[i])
# HPX = np.zeros(12*nside**2)
print(HPX.min(), HPX.max())
# for i in range(len(HPX_num)):
#     if HPX_num[i] != 0:
#         HPX[i] = HPX_diff_mag_sum[i] / HPX_num[i]

m = hp.ma(HPX)
m.mask = mask

plot_moll_cut(m, 'Average\ Y3-Y1\ mag\_auto(i)\ [17<i<22],\ 0.5\"\ match', 'diff_magautoi_17_22', 'Y3-Y1\ mag\ auto(i)', nside)
#plot_moll_cut(m, TITLE, FILE, LABEL_CB)
