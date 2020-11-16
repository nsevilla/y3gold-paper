#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import numpy as np
import pylab as plt
import pandas as pd
import healpy as hp
import fitsio

import skymap

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.axisartist import Subplot
from matplotlib.ticker import MaxNLocator
from descolors import BAND_COLORS


matplotlib.rcParams['xtick.labelsize'] = 12
matplotlib.rcParams['ytick.labelsize'] = 12
matplotlib.rcParams['axes.labelsize'] = 14
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']
matplotlib.rcParams['font.family'] = 'serif'

def set_defaults(kwargs,defaults):
    for k,v in defaults.items():
        kwargs.setdefault(k,v)
    return kwargs

def draw_peak(peak,**kwargs):
    kwargs.setdefault('ls','--')
    kwargs.setdefault('label','%.1f '%(peak))
    ax = plt.gca()
    ax.axvline(peak,**kwargs)

def create_hpxmap_hist_figure():
    #fig = plt.figure(figsize=(10.5,4))
    #gridspec=plt.GridSpec(1, 3)
    fig = plt.figure(figsize=(12.,4.))
    gridspec = plt.GridSpec(1,3,wspace=1.0)
    #gridspec.update(left=0.07,right=0.91,bottom=0.15,top=0.95,wspace=0.08)
    return fig, gridspec

def draw_hist(hpxmap,label,color,**kwargs):
    ax = plt.gca()

    if isinstance(hpxmap,np.ma.MaskedArray):
        pix = np.where(~hpxmap.mask)
    else:
        pix = np.where((np.isfinite(hpxmap)) & (hpxmap != hp.UNSEEN))

    data = hpxmap[pix]

    vmin = kwargs.pop('vmin',np.percentile(data,q=0.1))
    vmax = kwargs.pop('vmax',np.percentile(data,q=99.9))
    nbins = kwargs.pop('nbins',100)
    defaults = dict(bins=np.linspace(vmin,vmax,nbins),
                    histtype='step',density=True,lw=2,
                    peak=False,quantiles=False,color=color,label=label)

    set_defaults(kwargs,defaults)

    do_peak = kwargs.pop('peak')
    do_quantiles = kwargs.pop('quantiles')
    do_overflow = kwargs.pop('overflow',False)
    # Deal with bug: https://github.com/matplotlib/matplotlib/issues/6448/
    if do_overflow:
        data = np.clip(data,kwargs['bins'].min(),kwargs['bins'].max())
    else:
        data = data[(data > kwargs['bins'].min())&(data < kwargs['bins'].max())]

    n,b,p = ax.hist(data,**kwargs)
    ret = dict()

    peak = ((b[1:]+b[:-1])/2.)[np.argmax(n)]
    ret['peak'] = peak
    if do_peak:
        draw_peak(peak,color='k',label='%.1f'%(peak))

    ret['mean'] = np.mean(data)
    ret['std']  = np.std(data)

    quantiles = [5,16,50,84,95]
    percentiles = np.percentile(data,quantiles)
    ret['quantiles']   = quantiles
    ret['percentiles'] = percentiles
    for p,q in zip(percentiles,quantiles):
        ret['q%02d'%q] = p

    if do_quantiles:
        for q,p in zip(quantiles,percentiles):
            draw_peak(p,color='r',label='%.1f (%g%%)'%(p,100-q))

    ax.set_xlim(kwargs['bins'].min(),kwargs['bins'].max())
    return ret

def change_sigma(sb_lim,sigma_out=10.0,sigma_in=1.0):
    """ Convert from one detection significance to another

    Parameters
    ----------
    sb_lim : surface brightness limit
    sigma_out : new significance
    sigma_in  : old significance

    Returns
    -------
    new_sb_lim : surface brightness limit at new significance
    """
    # sb_lim = zeropoint - 2.5 * np.log10(sig_adu / pixel_scale**2)
    delta = -2.5*np.log10(sigma_out/sigma_in)
    return sb_lim + delta

#tiles = pd.read_csv('y3a1_tilenames.csv')
#results = pd.DataFrame(np.load('y3a1_results.npy'))
#merge = results.merge(tiles,on='tilename')
merge = pd.DataFrame(fitsio.read('/Users/nsevilla/des/sp_maps/Y3A1_SURFACE_BRIGHTNESS-v1.fits').byteswap().newbyteorder())

nside=64
bands = ['g','r','i','z','Y']
#bands = ['g']

#sizes = [10,30,60]
sizes = [10]

sigma = 3.0
print(merge['sb'][0])
merge['sb'] = change_sigma(merge['sb'],sigma_out=sigma)
print(merge['sb'][0])
print(merge['size'][0])
print(merge['band'][0])

for i,(band,v) in enumerate(BAND_COLORS.items()):
    print(f"{band}")
    for size in sizes:
        print(f"  {size}")
        if band == 'g':
            sel = (merge['band'] == b'g') & (merge['size'] == size)
        elif band == 'r':
            sel = (merge['band'] == b'r') & (merge['size'] == size)
        elif band == 'i':
            sel = (merge['band'] == b'i') & (merge['size'] == size)
        elif band == 'z':
            sel = (merge['band'] == b'z') & (merge['size'] == size)
        else:
            sel = (merge['band'] == b'Y') & (merge['size'] == size)
        #sel = (merge['band'] == band) & (merge['size'] == size)
        x = merge[sel]
        p16,p50,p84 = np.percentile(x['sb'],[16,50,84])
        print(f"  sb = {p50:.2f} [+{p84 - p50:.2f}, -{p50 - p16:.2f}]")
        print(f"  median(sb_err) = {x['sb_err'].median():.4f}")

        hpxmap = hp.UNSEEN*np.ones(hp.nside2npix(nside))
        pix = hp.ang2pix(64,x['ra_cent'],x['dec_cent'],lonlat=True)
        hpxmap[pix] = x['sb']
        vmin,vmax = np.percentile(x['sb'],[1,99])

        fig,gridspec = create_hpxmap_hist_figure()         
        plt.gca().axis('off')
 
        smap = skymap.DESSkymap(rect=gridspec[0:2])
        #smap.draw_hpxmap(hpxmap,vmin=vmin,vmax=vmax)
        smap.scatter(x['ra_cent'].values,x['dec_cent'].values,c=x['sb'].values,
                     s=4,marker='s',vmin=vmin,vmax=vmax,latlon=True,rasterized=True)
        smap.draw_inset_colorbar(label=r'mag arcsec$^{-2}$')

        ax1 = plt.gca()
        #ax1.annotate(f"{band}-band",(0.05,0.9),xycoords='axes fraction',
        #             fontsize=14)
        ax1.axis['right'].major_ticklabels.set_visible(False)
        #ax1.axis['top'].major_ticklabels.set_visible(False)

        # Plot histogram
        ax2 = Subplot(fig,gridspec[2])
        fig.add_subplot(ax2)
        plt.sca(ax2)
        #fig.add_subplot(gridspec[2])

        bins = np.linspace(x['sb'].min(),x['sb'].max(),35)
        ret = draw_hist(hpxmap,label=band,color=v,peak=False,bins=bins,density=True)
        ax2.yaxis.set_major_locator(MaxNLocator(6,prune='both'))
        ax2.axis['left'].major_ticklabels.set_visible(True)
        ax2.axis['right'].major_ticklabels.set_visible(False)
        ax2.axis['left'].label.set_visible(True)
        #ax2.axis['right'].label.set_text(r'Number of Tiles')
        ax2.axis['left'].label.set_text(r'PDF')
        ax2.axis['bottom'].label.set_visible(True)
        ax2.axis['bottom'].label.set_text(r'Surface Brightness Limit (mag arcsec$^{-2}$)')
        if band == 'z':
            ax2.set_aspect(1.5)
        elif  band == 'Y':
            ax2.set_aspect(1.2)
        else:
            ax2.set_aspect(0.8)

        #forceAspect(ax2,aspect=1)
        #plt.suptitle(f'{band}-band ({size}" x {size}")',y=0.92)

        plt.legend(loc='upper left')

        outfile = f'sbcontrast_{band}_s{size:d}_v2.pdf'
        print(f"  {outfile}")
        plt.savefig(outfile,bbox_inches='tight')

        #outfile = outfile.replace('.pdf','.png')
        #print(f"  {outfile}")
        #plt.savefig(outfile,bbox_inches='tight')
    
