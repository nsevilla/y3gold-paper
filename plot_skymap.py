#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

#import pylab as plt
import matplotlib as mpl
mpl.use("TkAgg") #needed in MacOS
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle
import numpy as np

### https://github.com/kadrlica/skymap
import skymap.survey
from skymap.constants import DECAM
from skymap.utils import gal2cel

plt.rc('text', usetex=True)
plt.rc('font', size=11)
#plt.rc('xtick.major', pad=5)
#plt.rc('xtick.minor', pad=5)
#plt.rc('ytick.major', pad=5)
#plt.rc('ytick.minor', pad=5)

plt.figure(figsize=(10,5))
smap = skymap.survey.SurveyMcBryde()

### Milky Way
smap.draw_milky_way(width=None,ls='-',color='k')
smap.draw_milky_way(width=10,dashes=(5,5),color='k')

### DES Wide footprint
x,y =  smap.draw_des(color='r')
z = np.vstack((x, y)).T
plt.gca().add_artist(PathPatch(Path(z),facecolor='r',alpha=0.2))

### HSC DR2 selected fields
x,y =  smap.draw_survey(hscfile='hsc-pdr2-poly_1.txt',color='b')
z = np.vstack((x, y)).T
plt.gca().add_artist(PathPatch(Path(z),facecolor='b',alpha=0.2))
x,y =  smap.draw_survey(hscfile='hsc-pdr2-poly_2.txt',color='b')
z = np.vstack((x, y)).T
plt.gca().add_artist(PathPatch(Path(z),facecolor='b',alpha=0.2))
x,y =  smap.draw_survey(hscfile='hsc-pdr2-poly_3.txt',color='b')
z = np.vstack((x, y)).T
plt.gca().add_artist(PathPatch(Path(z),facecolor='b',alpha=0.2))
proj = smap.proj(15,10)
plt.text(proj[0],proj[1], 'HSC-SSP DR2', weight='bold',
         fontsize=10, ha='center', va='center', color='b')

### KiDS-S 
#x,y =  smap.draw_survey(hscfile='kids-s-poly.txt',color='g')
#z = np.vstack((x, y)).T
#plt.gca().add_artist(PathPatch(Path(z),facecolor='g',alpha=0.2))
#proj = smap.proj(310,-30)
#plt.text(proj[0],proj[1], 'KiDS-S', weight='bold',
#         fontsize=10, ha='center', va='center', color='g')

### DES SN fields
kwargs = dict(facecolor='none',edgecolor='b',lw=0.75,zorder=10)
for v in skymap.survey.DES_SN.values():
    smap.tissot(v['ra'],v['dec'],DECAM,100,**kwargs)

kwargs['edgecolor'] = 'r'
kwargs['facecolor'] = 'r'
ddec = {'C3':5,'X3':-10,'E2':5,'S1':0}
dra = {'C3':0,'X3':0,'E2':0,'S1':15}
for k in ['C3','X3','E2','S1']:
    v = skymap.survey.DES_SN[k]
    if k != 'S ':
        smap.tissot(v['ra'],v['dec'],DECAM,100,**kwargs)
    ra = v['ra']
    dec = v['dec']
    proj = smap.proj(ra+dra[k],dec+ddec[k])
    plt.text(proj[0],proj[1], 'SN-'+k[0], weight='bold',
         fontsize=10, ha='center', va='center', color='k')
    
### COSMOS
RA_COSMOS = 150.12
DEC_COSMOS = 2.2
proj = smap.proj(RA_COSMOS,DEC_COSMOS+5)
kwargs['edgecolor'] = 'r'
smap.tissot(RA_COSMOS,DEC_COSMOS,DECAM,100,**kwargs)
plt.text(proj[0],proj[1], 'COSMOS', weight='bold',
         fontsize=10, ha='center', va='center', color='k')

### Draw LMC
from skymap.constants import RA_LMC, DEC_LMC, RADIUS_LMC
proj = smap.proj(RA_LMC,DEC_LMC)
smap.tissot(RA_LMC,DEC_LMC,RADIUS_LMC,100,fc='0.7',ec='0.5',alpha=0.5)
plt.text(proj[0],proj[1], 'LMC', weight='bold',
         fontsize=10, ha='center', va='center', color='green')

### Draw SMC
from skymap.constants import RA_SMC, DEC_SMC, RADIUS_SMC
proj = smap.proj(RA_SMC,DEC_SMC)
smap.tissot(RA_SMC,DEC_SMC,RADIUS_SMC,100,fc='0.7',ec='0.5',alpha=0.5)
plt.text(proj[0],proj[1], 'SMC', weight='bold',
         fontsize=10, ha='center', va='center', color='green')

### Galactic
smap.plot(*gal2cel(0, 0),marker='x',color='k',ms=10,latlon=True)
smap.plot(*gal2cel(0, -90),marker='+',color='k',ms=5,latlon=True)

# Inset axis
#llx,lly = smap(-30,-3.5)
#urx,ury = smap(-46,3.5)

#patch = Rectangle((llx,lly),urx-llx,ury-lly,edgecolor='k',facecolor='none',lw=0.5)
#ax = plt.gca()
#ax.add_artist(patch)

# Create the inset axis
#iax = plt.axes([.55, .6, .2, .15], facecolor='w')

# This turns of ticks and labels
#iax.get_xaxis().set_ticks([])
#iax.get_yaxis().set_ticks([])

#iax.set_xlim(llx,urx)
#iax.set_ylim(lly,ury)

# Add the footprint patch
#iax.add_artist(PathPatch(Path(z),facecolor='r',alpha=0.2))

#from mpl_toolkits.axes_grid1.inset_locator import mark_inset
#mark_inset(ax,iax, loc1=3, loc2=4, fc="none",ec='k',lw=0.5)

#plt.savefig('des_footprint.pdf',bbox_inches='tight')
plt.savefig('des_footprint_test.png')
