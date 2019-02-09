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

hdu = fits.open("data/NSEVILLAY1Y3CORRgriz.fits", memmap=True)
match_flag = hdu[1].data.field('MATCH_FLAG')
mag_auto_g_y1 = hdu[1].data.field('MAG_AUTO_G_Y1_CORRECTED')
mag_auto_g_y3 = hdu[1].data.field('MAG_AUTO_G_Y3_CORRECTED')
mag_auto_r_y1 = hdu[1].data.field('MAG_AUTO_R_Y1_CORRECTED')
mag_auto_r_y3 = hdu[1].data.field('MAG_AUTO_R_Y3_CORRECTED')
mag_auto_i_y1 = hdu[1].data.field('MAG_AUTO_I_Y1_CORRECTED')
mag_auto_i_y3 = hdu[1].data.field('MAG_AUTO_I_Y3_CORRECTED')
mag_auto_z_y1 = hdu[1].data.field('MAG_AUTO_Z_Y1_CORRECTED')
mag_auto_z_y3 = hdu[1].data.field('MAG_AUTO_Z_Y3_CORRECTED')
hdu.close()

diff_mag_g = mag_auto_g_y3 - mag_auto_g_y1
diff_mag_r = mag_auto_r_y3 - mag_auto_r_y1
diff_mag_i = mag_auto_i_y3 - mag_auto_i_y1
diff_mag_z = mag_auto_z_y3 - mag_auto_z_y1

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist2d(mag_auto_g_y3, diff_mag_g, bins=200, range=[[15.5, 25],[-0.8, 1.2]], cmap='inferno_r')
plt.title(r'$\mathrm{Y3-Y1\ mag\_auto(g)}$')
plt.xlabel(r'$\mathrm{mag\_auto\ g\ [Y3]}$')
plt.ylabel(r'$\mathrm{mag\_auto\ g\ [Y3-Y1]}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([15.5,25])
plt.ylim([-0.8,1.2])
cbaxes = fig.add_axes([0.908, 0.11, 0.02, 0.773])
cb = plt.colorbar(cax=cbaxes, cmap='inferno_r', orientation='vertical')
plt.savefig('figs/Y3-Y1corr_diff_mag_g_hist2d_total_range.png')
plt.close()
plt.clf()

from scipy.stats import norm
mu_g, std_g = norm.fit(diff_mag_g[(diff_mag_g > -0.8)&(diff_mag_g <= 1.2)])

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist(diff_mag_g, bins=np.arange(-0.8, 1.2, 0.01), color='r', histtype='step', log=False)
plt.title(r'$\mathrm{Y3-Y1\ mag\ auto(g),\ \mu = %s,\ \sigma =%s}$' % (np.round(mu_g,4), np.round(std_g,4)))
plt.ylabel(r'$\mathrm{N\_objects}$')
plt.xlabel(r'$\mathrm{mag\_auto\ (g)\ [Y3-Y1]}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([-0.8,1.2])
plt.savefig('figs/Y3-Y1_diff_mag_g_hist_total_range.png')
plt.close()
plt.clf()

###########################################################

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist2d(mag_auto_r_y3, diff_mag_r, bins=200, range=[[15.5, 25],[-0.8, 1.2]], cmap='inferno_r')
plt.title(r'$\mathrm{Y3-Y1\ mag\_auto(r)}$')
plt.xlabel(r'$\mathrm{mag\_auto\ r\ [Y3]}$')
plt.ylabel(r'$\mathrm{mag\_auto\ r\ [Y3-Y1]}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([15.5,25])
plt.ylim([-0.8,1.2])
cbaxes = fig.add_axes([0.908, 0.11, 0.02, 0.773])
cb = plt.colorbar(cax=cbaxes, cmap='inferno_r', orientation='vertical')
plt.savefig('figs/Y3-Y1corr_diff_mag_r_hist2d_total_range.png')
plt.close()
plt.clf()

from scipy.stats import norm
mu_r, std_r = norm.fit(diff_mag_r[(diff_mag_r > -0.8)&(diff_mag_r <= 1.2)])

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist(diff_mag_r, bins=np.arange(-0.8, 1.2, 0.01), color='r', histtype='step', log=False)
plt.title(r'$\mathrm{Y3-Y1\ mag\ auto(r),\ \mu = %s,\ \sigma =%s}$' % (np.round(mu_r,4), np.round(std_r,4)))
plt.ylabel(r'$\mathrm{N\_objects}$')
plt.xlabel(r'$\mathrm{mag\_auto\ (r)\ [Y3-Y1]}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([-0.8,1.2])
plt.savefig('figs/Y3-Y1_diff_mag_r_hist_total_range.png')
plt.close()
plt.clf()

###########################################################


fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist2d(mag_auto_i_y3, diff_mag_i, bins=200, range=[[15.5, 25],[-0.8, 1.2]], cmap='inferno_r')
plt.title(r'$\mathrm{Y3-Y1\ mag\_auto(i)}$')
plt.xlabel(r'$\mathrm{mag\_auto\ i\ [Y3]}$')
plt.ylabel(r'$\mathrm{mag\_auto\ i\ [Y3-Y1]}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([15.5,25])
plt.ylim([-0.8,1.2])
cbaxes = fig.add_axes([0.908, 0.11, 0.02, 0.773])
cb = plt.colorbar(cax=cbaxes, cmap='inferno_r', orientation='vertical')
plt.savefig('figs/Y3-Y1corr_diff_mag_i_hist2d_total_range.png')
plt.close()
plt.clf()

from scipy.stats import norm
mu_i, std_i = norm.fit(diff_mag_i[(diff_mag_i > -0.8)&(diff_mag_i <= 1.2)])

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist(diff_mag_i, bins=np.arange(-0.8, 1.2, 0.01), color='r', histtype='step', log=False)
plt.title(r'$\mathrm{Y3-Y1\ mag\ auto(i),\ \mu = %s,\ \sigma =%s}$' % (np.round(mu_i,4), np.round(std_i,4)))
plt.ylabel(r'$\mathrm{N\_objects}$')
plt.xlabel(r'$\mathrm{mag\_auto\ (i)\ [Y3-Y1]}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([-0.8,1.2])
plt.savefig('figs/Y3-Y1_diff_mag_i_hist_total_range.png')
plt.close()
plt.clf()

###########################################################


fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist2d(mag_auto_z_y3, diff_mag_z, bins=200, range=[[15.5, 25],[-0.8, 1.2]], cmap='inferno_r')
plt.title(r'$\mathrm{Y3-Y1\ mag\_auto(z)}$')
plt.xlabel(r'$\mathrm{mag\_auto\ z\ [Y3]}$')
plt.ylabel(r'$\mathrm{mag\_auto\ z\ [Y3-Y1]}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([15.5,25])
plt.ylim([-0.8,1.2])
cbaxes = fig.add_axes([0.908, 0.11, 0.02, 0.773])
cb = plt.colorbar(cax=cbaxes, cmap='inferno_r', orientation='vertical')
plt.savefig('figs/Y3-Y1corr_diff_mag_z_hist2d_total_range.png')
plt.close()
plt.clf()

from scipy.stats import norm
mu_z, std_z = norm.fit(diff_mag_z[(diff_mag_z > -0.8)&(diff_mag_z <= 1.2)])

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.hist(diff_mag_z, bins=np.arange(-0.8, 1.2, 0.01), color='r', histtype='step', log=False)
plt.title(r'$\mathrm{Y3-Y1\ mag\ auto(z),\ \mu = %s,\ \sigma =%s}$' % (np.round(mu_z,4), np.round(std_z,4)))
plt.ylabel(r'$\mathrm{N\_objects}$')
plt.xlabel(r'$\mathrm{mag\_auto\ (z)\ [Y3-Y1]}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([-0.8,1.2])
plt.savefig('figs/Y3-Y1_diff_mag_z_hist_total_range.png')
plt.close()
plt.clf()

###########################################################

fig = plt.figure(figsize = [9., 6.4])
ax1 = fig.add_subplot(111)
plt.title(r'$\mathrm{Y3-Y1\ mag\ auto}$')
plt.hist(diff_mag_g, bins=np.arange(-0.8, 1.2, 0.01), color='g', histtype='step', log=False, label=r'$\mathrm{g,\ \mu = %s,\ \sigma =%s}$' % (np.round(mu_g,4), np.round(std_g,4)))
plt.hist(diff_mag_r, bins=np.arange(-0.8, 1.2, 0.01), color='r', histtype='step', log=False, label=r'$\mathrm{r,\ \mu = %s,\ \sigma =%s}$' % (np.round(mu_r,4), np.round(std_r,4)))
plt.hist(diff_mag_i, bins=np.arange(-0.8, 1.2, 0.01), color='firebrick', histtype='step', log=False, label=r'$\mathrm{i,\ \mu = %s,\ \sigma =%s}$' % (np.round(mu_i,4), np.round(std_i,4)))
plt.hist(diff_mag_z, bins=np.arange(-0.8, 1.2, 0.01), color='k', histtype='step', log=False, label=r'$\mathrm{z,\ \mu = %s,\ \sigma =%s}$' % (np.round(mu_z,4), np.round(std_z,4)))
plt.ylabel(r'$\mathrm{N\ objects}$')
plt.xlabel(r'$\mathrm{mag\_auto\ (Y3-Y1)}$')
plt.grid(color='grey', linestyle='-', linewidth=0.5, zorder=-1)
plt.xlim([-0.5,0.5])
plt.legend()
plt.savefig('figs/Y3-Y1_diff_mag_griz_hist_05_range.png')
plt.close()
plt.clf()
