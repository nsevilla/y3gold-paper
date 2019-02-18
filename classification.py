import healpy as hp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import matplotlib.image as mpimg
from matplotlib import rc
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
plt.rc('font',**{'family':'serif','serif':['DejaVu Sans'], 'size': 12})
mpl.rcParams['legend.numpoints'] = 1
plt.rc('text', usetex=True)
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['lines.linewidth'] = 3
cmap = mpl.cm.get_cmap("inferno")
cmap.set_under('darkgray')
cmap.set_bad('darkgray')

def plot_moll_cut(m, ns, label, filename, obj_k, number, vmin, vmax):
    ra_DES, dec_DES = np.loadtxt('data/des-round17-poly.txt', usecols=(0, 1), unpack=True)
    m_ = np.copy(m)
    m_[(m <= vmin)&(m > 0.)] = vmin
    m_[(m >= vmax)] = vmax
    hp.visufunc.mollview(map=m_, flip='geo', format='%.3g', cmap='inferno', coord='C', hold=True, xsize=800,
                         nest=False, cbar=False, notext=True, margins=None)
    hp.projplot(ra_DES, dec_DES, lonlat=True, coord=['C'], zorder=10, color='b', lw=0.5)
    hp.graticule(dpar=15.,dmer=30., lw=0.1)
    plt.savefig('figs/primary.png', dpi=600, bbox_inches='tight', pad_inches=0)
    plt.close()
    plt.clf()
    
    d_op = mpimg.imread('figs/primary.png')
    w, h = len(d_op[0,:]), len(d_op[:,0])
    # d_op = d_op[int(0.08156*h):int(0.9095*h),int(0.082*w):int(0.982*w)]
    d_op = d_op[int(0.125*h):int(0.88*h),int(0.1*w):int(0.99*w)]
    w, h = len(d_op[0,:]), len(d_op[:,0])
    fig = plt.figure(figsize = [15., 7])
    
    gridspec.GridSpec(12,30)
    plt.subplot2grid((12,30), (0,0), colspan=20, rowspan=12)
    #ax1 = fig.add_subplot(111)
    #plt.set_facecolor('lightgray')
    plt.imshow(d_op[int(0.4293*h):int(0.9497*h), int(0.3611*w):int(0.75*w)], extent=[-45., 90.,-75.,10.],
               aspect='auto', origin='upper', interpolation=None, vmin=0., vmax=0.2, cmap=cmap)
    plt.title(r'$\mathrm{%s}$' % label, fontsize=20)
    plt.xlabel(r'$\mathrm{\alpha}$', fontsize=16)
    plt.ylabel(r'$\mathrm{\delta}$', fontsize=16)
    plt.ylim([-75,10])
    plt.xlim([90,-45])

    y_ = [-0., -15, -34, -50, -66, -75]
    labels_y=[r'$\mathrm{0^o}$', r'$\mathrm{-15^o}$', r'$\mathrm{-30^o}$', r'$\mathrm{-45^o}$', r'$\mathrm{-60^o}$', r'$\mathrm{-75^o}$']
    x = [67, 54, 41.5, 29, 16, 3, -9, -22, -35]
    labels=[r'$\mathrm{150^o}$', r'$\mathrm{120^o}$', r'$\mathrm{90^o}$', r'$\mathrm{60^o}$', r'$\mathrm{30^o}$', r'$\mathrm{0^o}$', r'$\mathrm{-30^o}$', r'$\mathrm{-60^o}$', r'$\mathrm{-90^o}$']

    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
    plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=True)
    plt.xticks(x, labels, rotation='horizontal')
    plt.yticks(y_, labels_y)
    cbaxes = fig.add_axes([0.427, 0.63, 0.2, 0.035])
    cb = plt.colorbar(cax=cbaxes, cmap=cmap, orientation='horizontal', ticks=[np.arange(vmin,vmax,100)]) #, label=r'$\mathrm{%s}$' % labelcb)
    cb.ax.tick_params(labelsize=10)

    plt.subplot2grid((12,30), (0,23), colspan=7, rowspan=12)
    plt.hist(m, bins=np.arange(0,800,10.), color='r', alpha = 0.5, lw=1, histtype='step', log=True, zorder=3)
    n, bins, patches = plt.hist(m, bins=np.arange(0.,800,10.), color='r', alpha = 0.5, lw=1, log=True, zorder=3)
    CA = plt.gca()
    CA.add_patch(Rectangle((vmin, 1.), vmax-vmin, 1.2*n.max(), alpha=0.25, color='b', zorder=1))
    #ar_ = np.linspace(1.,1.3*np.max(n),10)
    #plt.plot(vmin + 0*ar_, ar_, 'k:', lw=0.8)
    #plt.plot(vmax + 0*ar_, ar_, 'k:', lw=0.8)
    plt.xlabel(r'$\mathrm{%s\ ocupation}$' % obj_k)
    plt.xlim([0,800])
    plt.xticks([0,200,400,600,800])
    plt.ylim([1.,1.2*np.max(n)])
    plt.grid(color='grey', linestyle='-', linewidth=0.5, alpha=0.25)
    plt.ylabel(r'$\mathrm{Number\ of\ pixels}$')
    #text1 = str("""Median: \n{0:4.3f} arcsec""".format(np.median(m[np.abs(m) <= 1])))
    #plt.annotate(text1,(0.5,0.95),xycoords='axes fraction',ha='left',va='top',size=12)
    plt.savefig('figs/HP_EQU_sof' + str(number) +'.png', dpi=300, bbox_inches='tight')
    plt.close()
    #plt.clf()
    
def read_map_source(class_number, nside, nside2):
    Y3_star_map = fits.open('data/sources_sof'+str(class_number)+'_16_23.fits',memmap=True)
    Y3_star_data = Y3_star_map[1].data
    #mask = np.zeros(hp.nside2npix(nside), dtype=np.bool)
    #Y3_star_data_fullsky_nest_4096 = np.zeros(hp.nside2npix(nside))
    Y3_star_data_fullsky_nest_4096 = np.full(hp.nside2npix(nside),hp.UNSEEN)
    Y3_star_data_fullsky_nest_4096[Y3_star_data['PIXEL']] = Y3_star_data['SIGNAL']
    Y3_star_data_fullsky_ring_4096 = hp.reorder(Y3_star_data_fullsky_nest_4096,n2r=True)
    Y3_star_data_fullsky = hp.ud_grade(Y3_star_data_fullsky_ring_4096,nside_out=nside2,power=-2)
    pix, = np.where(Y3_star_data_fullsky != hp.UNSEEN)
    mask = Y3_star_data_fullsky != hp.UNSEEN
    #mask[Y3_star_data_fullsky == 0.] = 1
    #m = hp.ma(Y3_star_data_fullsky)
    #m.mask = mask
    #pix = np.where(np.abs(m) < 1)
    nsources = Y3_star_data_fullsky[mask]
    return Y3_star_data_fullsky, nsources

HPX0, nsources0 = read_map_source(0, 4096, 512)
# HPX0[HPX0 >= 300.] = 300

HPX3, nsources3 = read_map_source(3, 4096, 512)
# HPX3[(HPX3 <= 150.)&(HPX3 >= 0.)] = 150

plot_moll_cut(HPX0, nsources0, 'DES\ stellar\ distribution', ' ', 'Star', 0, 70., 500)
plot_moll_cut(HPX3, nsources3, 'DES\ galaxy\ distribution', ' ', 'Galaxy', 3, 100, 400)
