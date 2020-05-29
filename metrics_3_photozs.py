import numpy as np
import matplotlib.pylab as plt
import matplotlib
from astropy.io import fits
from astroML.resample import bootstrap
matplotlib.use("Agg")
matplotlib.style.use('des_dr1')

path = '/Users/nsevilla/y3gold-paper/data/'
file = 'photoz_valid_annz2_newdnf_bpz_y1bpz.fits'

def bias(delta_z):
	#delta_z = z_photo-z_spec
	bias = np.mean(delta_z)
	err_bias = np.std(delta_z)/np.sqrt(len(delta_z))
	return bias, err_bias

def sigma68_1z(delta_z,z_spec,axis=None):
	#delta_z = z_photo-z_spec
	upper, lower = np.percentile(delta_z, [84.075, 15.825], axis=axis)
	sigma68_1z = (upper - lower) / 2.0 / (1 + np.mean(z_spec))
	return sigma68_1z
	
	#val_nocorr = sigma_68_1pz(selection_nocorr[refPz]-selection_nocorr[refZ],selection[refZ])
def err_sigma68_1z(delta_z,z_spec,function=sigma68_1z):
	#delta_z = z_photo-z_spec
	n = 50
	err_sigma68_1z = np.std(bootstrap(delta_z, n, function, kwargs=dict(axis=1, z_spec = z_spec)))
	return err_sigma68_1z

data = fits.open(path+file,memmap=True)[1].data

z_photos = ['Z_DNF_NEW','Z_ANNZ2','BPZ_ZMEAN_SOF','BPZY1_MEAN_Z']
z_names = {'Z_DNF_NEW':'DNF','Z_ANNZ2':'ANNz2','BPZ_ZMEAN_SOF':'BPZ','BPZY1_MEAN_Z':'BPZ Y1'}

z_spec_col = 'Z'
z_spec = data[z_spec_col]

zbins = [(0.2,0.3),(0.3,0.4),(0.4,0.5),(0.5,0.6),(0.6,0.7),(0.7,0.8),(0.8,0.9),(0.9,1.0),(1.0,1.1),(1.1,1.2),(1.2,1.3)]

x_offset = 0
plt.figure(figsize=(14,10))
for z_photo_i in z_photos:
	print(z_photo_i)
	z_photo = data[z_photo_i]
	
	zmid_list = []
	bias_i = []
	err_bias_i = []
	for zbin in zbins:
		zmin = zbin[0]
		zmax = zbin[1]
		zmid = (zmin+zmax)/2.
		zmid_list.append(zmid)
		
		z_mask = (z_photo>zmin)*(z_photo<zmax)
		z_photo_sel = z_photo[z_mask]
		z_spec_sel = z_spec[z_mask]
		delta_z = z_photo_sel-z_spec_sel
		
		bias_zbin,err_bias_zbin = bias(delta_z)
		bias_i.append(bias_zbin)
		err_bias_i.append(err_bias_zbin)
		
	plt.errorbar(np.array(zmid_list)+x_offset,bias_i,yerr=err_bias_i,marker='o',ls='',label=z_names[z_photo_i])
	x_offset = x_offset+0.005
plt.grid()
plt.xlabel('photo-z bin')
plt.ylabel(r'$z_{photo}-z_{spec}$',fontsize=30)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.legend(loc="best",fontsize=20)
plt.savefig('/Users/nsevilla/y3gold-paper/figs/photoz_bias_test.png')
#plt.show()
plt.close()

x_offset = 0
plt.figure(figsize=(14,10))
for z_photo_i in z_photos:
	print(z_photo_i)
	z_photo = data[z_photo_i]
	
	zmid_list = []
	sigma68_i = []
	err_sigma68_i = []
	for zbin in zbins:
		zmin = zbin[0]
		zmax = zbin[1]
		zmid = (zmin+zmax)/2.
		zmid_list.append(zmid)
		
		z_mask = (z_photo>zmin)*(z_photo<zmax)
		z_photo_sel = z_photo[z_mask]
		z_spec_sel = z_spec[z_mask]
		delta_z = z_photo_sel-z_spec_sel
		
		sigma68_zbin = sigma68_1z(delta_z,z_spec_sel)
		err_sigma68_zbin = err_sigma68_1z(delta_z,z_spec_sel,function=sigma68_1z)
		sigma68_i.append(sigma68_zbin)
		err_sigma68_i.append(err_sigma68_zbin)
	
	plt.errorbar(np.array(zmid_list)+x_offset,sigma68_i,yerr=err_sigma68_i,marker='o',ls='',label=z_names[z_photo_i])
	x_offset = x_offset+0.005

plt.grid()
plt.xlabel('photo-z bin',fontsize=30)
plt.ylabel(r'$\sigma_{68}/(1+z_{spec})$',fontsize=30)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
plt.legend(loc="best",fontsize=20)
plt.savefig('/Users/nsevilla/y3gold-paper/figs/photoz_sigma68_1z_test.png')
#plt.show()
plt.close()

